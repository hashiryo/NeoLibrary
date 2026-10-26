#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <cassert>
#include <chrono>
#include <numeric>
#include <random>
#include <tuple>
#include <vector>
#include "neo/algebra/ZMod.hpp"
#include "neo/linalg/Matrix.hpp"
#include "neo/linalg/det.hpp"
// 疎な n × n 行列の行列式。成分は (行, 列, 値) の列で渡し、同じ位置が複数あれば和を取る。T は ZMod で、法は素数。ほかの型は
// 密な行列にして det を呼ぶ。
// Wiedemann の方法。乱択の対角行列 D を右から掛けた B = A D の列 s_i = u^T B^i v (i < 2n、u と v も乱択) の最小多項式 f を
// Berlekamp–Massey で求める。f は B の最小多項式を割り切るので、次数が n なら f は B の特性多項式で、det(B) = (-1)^n f(0)、
// det(A) = det(B) / det(D)。次数が n 未満でも f(0) = 0 なら B は正則でなく det(A) = 0。どちらでもないとき (A が正則なら
// 1 回あたり n^2 / p 程度の確率で起きる) は乱数を選び直し、3 回続いたら密な行列の掃き出し (neo/linalg/det.hpp) にする。
// 返す値はいつも正しく、乱択で揺れるのは時間だけ。空の行か列があれば 0 を返す。法が偶数 (2) のときは密な掃き出しにする。
// 行と列は成分の数の多い順に同じ並べ替えで付け直し (相似変換なので行列式は変わらない)、行ごとに成分を続けて並べる (CSR)。
// 長さの同じ行が続くので、行の内側のループの終わりの分岐が予測どおりになる。行の和は u64 にため、各行の和は畳んでから
// Montgomery の還元で 2^{-32} 倍して M 未満に直す。毎回 2^{-32} 倍するのは、λ = 2^{-32} として λB の列を作るのと同じなので、
// det(λB) を求めてから λ^{-n} を掛けて戻す。Berlekamp–Massey のずれの内積と、漸化式の更新 (Shoup の方法) も 4 個ずつ。
// procon-judge の yosupo-sparse-matrix-det で書き比べた zmod_wiedemann_csr を写したもの。
namespace sparse_det_internal {
using u32= unsigned;
using u64= unsigned long long;
// x の 4 個の u64 それぞれの 2^32 の位を R = 2^32 mod M に掛けて畳む。結果は (2^32 - 1) M 以下。
inline __m256i fold(__m256i x, __m256i Rv) { return _mm256_add_epi64(_mm256_mul_epu32(_mm256_srli_epi64(x, 32), Rv), _mm256_and_si256(x, _mm256_set1_epi64x(0xffffffffll))); }
template <class Z> u32 hsum(__m256i x, __m256i Rv) {
 alignas(32) u64 t[4];
 _mm256_store_si256((__m256i*)t, fold(x, Rv));
 return (Z(t[0]) + Z(t[1]) + Z(t[2]) + Z(t[3])).val();
}
// sum_{i<n} x[i] y[i] (x, y は M 未満の u32)。
template <class Z> u32 dot(const u32* x, const u32* y, int n, __m256i Rv) {
 __m256i acc= _mm256_setzero_si256();
 int i= 0;
 for(int t= 0; i + 4 <= n; i+= 4) {
  acc= _mm256_add_epi64(acc, _mm256_mul_epu32(_mm256_cvtepu32_epi64(_mm_loadu_si128((const __m128i*)(x + i))), _mm256_cvtepu32_epi64(_mm_loadu_si128((const __m128i*)(y + i)))));
  if(++t == 8) acc= fold(acc, Rv), t= 0;
 }
 Z r= Z::raw(hsum<Z>(acc, Rv));
 for(; i < n; ++i) r+= Z::raw(x[i]) * Z::raw(y[i]);
 return r.val();
}
// CSR (行と列は付け直したあとの番号) の行列 A で 1 回試す。0 なら det に行列式を書き、1 なら正則でないと分かった、2 なら分からない。
template <class Z> int attempt(int N, const std::vector<int>& rp, const std::vector<int>& col, const std::vector<u32>& cv, std::mt19937_64& rng, Z& det) {
 const u32 M= Z::mod();
 const u64 R32= (u64(1) << 32) % M;
 const u32 NIV= [&] {  // -M^{-1} mod 2^32
  u32 x= M;
  for(int i= 0; i < 5; ++i) x*= 2 - M * x;
  return -x;
 }();
 auto rnd= [&]() { return u32(rng() % (M - 1) + 1); };
 const int K= rp[N], N4= (N + 3) & ~3;
 std::vector<u32> d(N), val(K);
 for(auto& x: d) x= rnd();
 for(int e= 0; e < K; ++e) val[e]= (Z::raw(cv[e]) * Z::raw(d[col[e]])).val();
 std::vector<u32> u(N4), w(N4), w2(N4), s(2 * N);
 for(int i= 0; i < N; ++i) u[i]= rnd(), w[i]= rnd();
 std::vector<u64> acc(N4);
 const __m256i Rv= _mm256_set1_epi64x(R32), Mv= _mm256_set1_epi64x(M), Niv= _mm256_set1_epi64x(NIV), M32= _mm256_set1_epi32(M);
 for(int i= 0; i < 2 * N; ++i) {
  s[i]= dot<Z>(u.data(), w.data(), N, Rv);
  if(i + 1 == 2 * N) break;
  for(int q= 0; q < N; ++q) {
   u64 x= 0;
   for(int e= rp[q], k= 0; e < rp[q + 1]; ++e) {
    x+= u64(val[e]) * w[col[e]];
    if(++k == 8) x= (x >> 32) * R32 + (x & 0xffffffffu), k= 0;
   }
   acc[q]= x;
  }
  for(int q= 0; q < N4; q+= 4) {  // 畳んで (2^32 - 1) M 以下にしてから Montgomery の還元 (2M 未満が返る)
   const __m256i x= fold(_mm256_loadu_si256((const __m256i*)&acc[q]), Rv);
   const __m256i t= _mm256_srli_epi64(_mm256_add_epi64(x, _mm256_mul_epu32(_mm256_mul_epu32(x, Niv), Mv)), 32);
   const __m256i r= _mm256_min_epu32(t, _mm256_sub_epi32(t, M32));
   _mm_storeu_si128((__m128i*)&w2[q], _mm256_castsi256_si128(_mm256_permutevar8x32_epi32(r, _mm256_setr_epi32(0, 2, 4, 6, 0, 0, 0, 0))));
  }
  std::swap(w, w2);
 }
 // Berlekamp–Massey。ずれは C と、s を逆順に並べた sr の内積。更新 C[i + m] += (M - coef) B[i] は Shoup の方法で 4 個ずつ。
 std::vector<u32> sr(2 * N), C(N + 4), B(N + 4), T;
 for(int i= 0; i < 2 * N; ++i) sr[i]= s[2 * N - 1 - i];
 C[0]= B[0]= 1;
 int L= 0, LB= 0, m= 1;
 Z bb= Z::raw(1);
 for(int n= 0; n < 2 * N; ++n) {
  const Z dd= Z::raw(dot<Z>(C.data(), &sr[2 * N - 1 - n], L + 1, Rv));
  if(dd == Z()) {
   ++m;
   continue;
  }
  const u32 cf= (-(dd * bb.inv())).val(), cs= u32((u64(cf) << 32) / M);
  const bool grow= 2 * L <= n;
  if(grow) T.assign(C.begin(), C.begin() + L + 1);
  const int top= std::min(LB, N - m);  // C[i + m] (i <= top) を直す
  const __m256i CF= _mm256_set1_epi64x(cf), CS= _mm256_set1_epi64x(cs), M2= _mm256_set1_epi32(2 * M);
  int i= 0;
  for(; i + 4 <= top + 1; i+= 4) {
   const __m256i bv= _mm256_cvtepu32_epi64(_mm_loadu_si128((const __m128i*)&B[i]));
   const __m256i q= _mm256_srli_epi64(_mm256_mul_epu32(bv, CS), 32);
   const __m256i t= _mm256_sub_epi64(_mm256_mul_epu32(bv, CF), _mm256_mul_epu32(q, Mv));               // [0, 2M)
   __m256i x= _mm256_add_epi64(_mm256_cvtepu32_epi64(_mm_loadu_si128((const __m128i*)&C[i + m])), t);  // [0, 3M)
   x= _mm256_min_epu32(x, _mm256_sub_epi32(x, M2));
   x= _mm256_min_epu32(x, _mm256_sub_epi32(x, M32));
   _mm_storeu_si128((__m128i*)&C[i + m], _mm256_castsi256_si128(_mm256_permutevar8x32_epi32(x, _mm256_setr_epi32(0, 2, 4, 6, 0, 0, 0, 0))));
  }
  for(; i <= top; ++i) C[i + m]= (Z::raw(C[i + m]) + Z::raw(cf) * Z::raw(B[i])).val();
  if(grow) {
   std::fill(B.begin(), B.begin() + LB + 1, 0u);
   std::copy(T.begin(), T.end(), B.begin());
   LB= L, L= n + 1 - L, bb= dd, m= 1;
  } else ++m;
 }
 if(L != N) return C[L] ? 2 : 1;  // f(x) = x^L C(1/x) の f(0) = C[L]
 // det(λB) = (-1)^N C[N]、det(A) = det(λB) λ^{-N} / det(D)、λ^{-1} = 2^32
 det= Z::raw(C[N]);
 if(N & 1) det= -det;
 Z dD= Z::raw(1);
 for(u32 x: d) dD*= Z::raw(x);
 det*= Z(R32).pow(N) * dD.inv();
 return 0;
}
template <class T> T dense(int n, const std::vector<std::tuple<int, int, T>>& entries) {
 Matrix<T> A(n, n);
 for(const auto& [i, j, x]: entries) A[i][j]+= x;
 return det(std::move(A));
}
}
template <class T> T sparse_det(int n, const std::vector<std::tuple<int, int, T>>& entries) {
 using namespace sparse_det_internal;
 if constexpr(!linalg_internal::is_zmod_v<T>) return dense(n, entries);
 else {
  const u32 M= T::mod();
  if(n == 0) return T(1);
  if(M == 1) return T();
  {  // 空の行か列があれば正則でない
   std::vector<char> hr(n), hc(n);
   for(const auto& [i, j, x]: entries)
    if(x != T()) hr[i]= hc[j]= 1;
   for(int i= 0; i < n; ++i)
    if(!hr[i] || !hc[i]) return T();
  }
  if(!(M & 1)) return dense(n, entries);
  std::vector<int> len(n), perm(n), rank(n);
  for(const auto& [i, j, x]: entries)
   if(x != T()) ++len[i];
  std::iota(perm.begin(), perm.end(), 0);
  std::stable_sort(perm.begin(), perm.end(), [&](int x, int y) { return len[x] > len[y]; });
  for(int q= 0; q < n; ++q) rank[perm[q]]= q;
  std::vector<int> rp(n + 1);
  for(int q= 0; q < n; ++q) rp[q + 1]= rp[q] + len[perm[q]];
  std::vector<int> col(rp[n]);
  std::vector<u32> cv(rp[n]);
  {
   std::vector<int> pos(rp.begin(), rp.end() - 1);
   for(const auto& [i, j, x]: entries)
    if(x != T()) {
     const int e= pos[rank[i]]++;
     col[e]= rank[j], cv[e]= x.val();
    }
  }
  static std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());
  for(int t= 0; t < 3; ++t) {
   T d;
   const int r= attempt<T>(n, rp, col, cv, rng, d);
   if(r == 0) return d;
   if(r == 1) return T();
  }
  return dense(n, entries);
 }
}
