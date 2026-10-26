#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <cassert>
#include <cmath>
#include <type_traits>
#include <utility>
#include <vector>
#include "neo/algebra/ZMod.hpp"
#include "neo/linalg/Matrix.hpp"
// 正方行列の特性多項式 det(xI - A) の係数 p_0, p_1, ..., p_n (p_n = 1)。T は体で、ZMod なら法は素数。
// 相似変換で上 Hessenberg の形 H に直してから、左上 k × k の特性多項式 p_k を
// p_{k+1} = (x - H[k][k]) p_k - sum_{i<k} H[i][k] (prod_{j=i+1}^{k} H[j][j-1]) p_i で順に求める。
// 列 k の手順では、行 k + 1 以降から列 k の値が 0 でない行を選んで行 k + 1 に置き (行と列を同時に入れ替える)、
// f_i = H[i][k] / H[k + 1][k] (i >= k + 2) として、行 i から行 k + 1 の f_i 倍を引き (左から L^{-1})、列 k + 1 に列 i の f_i 倍を
// 足す (右から L)。
// ZMod は行列を 2M 未満の u32 で持つ。列の操作は、全部の行で「行の列 k + 2 以降と f の内積」を列 k + 1 に足すことになるので、
// 行 k + 2 以降は、行の更新と内積を 1 回の読み書きで済ませる (内積は更新したあとの値で取る)。行の更新は、掛ける数を行ごとに決まった
// 値として Shoup の方法で 2M 未満に直し、内積は u64 にためて 4 回ごとに 2^32 の位を 2^32 mod M に掛けて畳む。漸化式の和は、
// 多項式を u64 で持って 8 本ずつまとめて足す。ほかの体は素朴な形 (浮動小数点数は絶対値が最大の要を選ぶ) にする。
// procon-judge の yosupo-characteristic-polynomial で書き比べた zmod_hess32 を写したもの。
namespace charpoly_internal {
using u32= unsigned;
using u64= unsigned long long;
constexpr int K= 8;
// row[j] += sum_{s<8} f[s] U[s][j] を j in [jb, w) で 4 個ずつ足してから畳む (neo/linalg/det.hpp の写し)。
[[gnu::always_inline]] inline void update_row(u64* __restrict row, const u64* const* U, const u32* f, int jb, int w, u64 R) {
 const __m256i Rv= _mm256_set1_epi64x(R), LO= _mm256_set1_epi64x(0xffffffffll);
 const __m256i f0= _mm256_set1_epi64x(f[0]), f1= _mm256_set1_epi64x(f[1]), f2= _mm256_set1_epi64x(f[2]), f3= _mm256_set1_epi64x(f[3]);
 const __m256i f4= _mm256_set1_epi64x(f[4]), f5= _mm256_set1_epi64x(f[5]), f6= _mm256_set1_epi64x(f[6]), f7= _mm256_set1_epi64x(f[7]);
 const u64 *U0= U[0], *U1= U[1], *U2= U[2], *U3= U[3], *U4= U[4], *U5= U[5], *U6= U[6], *U7= U[7];
 for(int j= jb; j < w; j+= 4) {
  __m256i acc= _mm256_loadu_si256((const __m256i*)(row + j));
  acc= _mm256_add_epi64(acc, _mm256_mul_epu32(_mm256_loadu_si256((const __m256i*)(U0 + j)), f0));
  acc= _mm256_add_epi64(acc, _mm256_mul_epu32(_mm256_loadu_si256((const __m256i*)(U1 + j)), f1));
  acc= _mm256_add_epi64(acc, _mm256_mul_epu32(_mm256_loadu_si256((const __m256i*)(U2 + j)), f2));
  acc= _mm256_add_epi64(acc, _mm256_mul_epu32(_mm256_loadu_si256((const __m256i*)(U3 + j)), f3));
  acc= _mm256_add_epi64(acc, _mm256_mul_epu32(_mm256_loadu_si256((const __m256i*)(U4 + j)), f4));
  acc= _mm256_add_epi64(acc, _mm256_mul_epu32(_mm256_loadu_si256((const __m256i*)(U5 + j)), f5));
  acc= _mm256_add_epi64(acc, _mm256_mul_epu32(_mm256_loadu_si256((const __m256i*)(U6 + j)), f6));
  acc= _mm256_add_epi64(acc, _mm256_mul_epu32(_mm256_loadu_si256((const __m256i*)(U7 + j)), f7));
  acc= _mm256_add_epi64(_mm256_mul_epu32(_mm256_srli_epi64(acc, 32), Rv), _mm256_and_si256(acc, LO));
  _mm256_storeu_si256((__m256i*)(row + j), acc);
 }
}
inline __m256i fold(__m256i x, __m256i Rv) { return _mm256_add_epi64(_mm256_mul_epu32(_mm256_srli_epi64(x, 32), Rv), _mm256_and_si256(x, _mm256_set1_epi64x(0xffffffffll))); }
// 4 個の u64 (どれも 2^62 + 2^32 未満) の和を M 未満で返す。
template <class Z> u32 hsum(__m256i acc) {
 alignas(32) u64 t[4];
 _mm256_store_si256((__m256i*)t, acc);
 return (Z(t[0]) + Z(t[1]) + Z(t[2]) + Z(t[3])).val();
}
// row[j] (2M 未満) に c P[j] mod M を足して 2M 未満に直し (j in [j0, w)、8 個ずつ)、更新したあとの値と f[j] の内積を返す。
// cs = floor(c 2^32 / M)、R = 2^32 mod M。P[j] < 2M、f[j] < M。
template <class Z> u32 fused(u32* __restrict row, const u32* __restrict P, const u32* __restrict f, u32 c, u32 cs, int j0, int w, u64 R) {
 const u32 M= Z::mod();
 const __m256i C= _mm256_set1_epi64x(c), CS= _mm256_set1_epi64x(cs), Mv= _mm256_set1_epi64x(M), MM= _mm256_set1_epi32(2 * M);
 const __m256i Rv= _mm256_set1_epi64x(R);
 __m256i ae= _mm256_setzero_si256(), ao= ae;
 int t= 0;
 for(int j= j0; j < w; j+= 8) {
  __m256i x= _mm256_loadu_si256((const __m256i*)(row + j));
  const __m256i p= _mm256_loadu_si256((const __m256i*)(P + j)), po= _mm256_srli_epi64(p, 32);
  const __m256i qe= _mm256_srli_epi64(_mm256_mul_epu32(p, CS), 32), qo= _mm256_srli_epi64(_mm256_mul_epu32(po, CS), 32);
  const __m256i te= _mm256_sub_epi64(_mm256_mul_epu32(p, C), _mm256_mul_epu32(qe, Mv));
  const __m256i to= _mm256_sub_epi64(_mm256_mul_epu32(po, C), _mm256_mul_epu32(qo, Mv));
  x= _mm256_add_epi32(x, _mm256_blend_epi32(te, _mm256_slli_epi64(to, 32), 0xaa));
  x= _mm256_min_epu32(x, _mm256_sub_epi32(x, MM));
  _mm256_storeu_si256((__m256i*)(row + j), x);
  const __m256i fv= _mm256_loadu_si256((const __m256i*)(f + j));
  ae= _mm256_add_epi64(ae, _mm256_mul_epu32(x, fv));
  ao= _mm256_add_epi64(ao, _mm256_mul_epu32(_mm256_srli_epi64(x, 32), _mm256_srli_epi64(fv, 32)));
  if(++t == 4) ae= fold(ae, Rv), ao= fold(ao, Rv), t= 0;
 }
 return (Z::raw(hsum<Z>(fold(ae, Rv))) + Z::raw(hsum<Z>(fold(ao, Rv)))).val();
}
// sum_j row[j] f[j] (j in [j0, w)、8 個ずつ)。row[j] < 2M、f[j] < M。
template <class Z> u32 dot(const u32* __restrict row, const u32* __restrict f, int j0, int w, u64 R) {
 const __m256i Rv= _mm256_set1_epi64x(R);
 __m256i ae= _mm256_setzero_si256(), ao= ae;
 int t= 0;
 for(int j= j0; j < w; j+= 8) {
  const __m256i x= _mm256_loadu_si256((const __m256i*)(row + j)), fv= _mm256_loadu_si256((const __m256i*)(f + j));
  ae= _mm256_add_epi64(ae, _mm256_mul_epu32(x, fv));
  ao= _mm256_add_epi64(ao, _mm256_mul_epu32(_mm256_srli_epi64(x, 32), _mm256_srli_epi64(fv, 32)));
  if(++t == 4) ae= fold(ae, Rv), ao= fold(ao, Rv), t= 0;
 }
 return (Z::raw(hsum<Z>(fold(ae, Rv))) + Z::raw(hsum<Z>(fold(ao, Rv)))).val();
}
template <class Z> std::vector<Z> zmod(const Matrix<Z>& A) {
 const int n= A.h;
 const u32 M= Z::mod(), M2= 2 * M;
 const u64 R= (u64(1) << 32) % M;
 auto val= [&](u32 x) { return x >= M ? x - M : x; };
 auto addm= [&](u32 x, u32 d) {  // 2M 未満の x に M 未満の d を足して 2M 未満に直す
  x+= d;
  return x >= M2 ? x - M2 : x;
 };
 const int w= (n + 7) & ~7;
 std::vector<u32> H((size_t)n * w), f(w);
 for(int i= 0; i < n; ++i)
  for(int j= 0; j < n; ++j) H[(size_t)i * w + j]= A[i][j].val();
 auto row= [&](int i) { return &H[(size_t)i * w]; };
 for(int k= 0; k + 2 < n; ++k) {
  int p= -1;
  for(int i= k + 1; i < n; ++i)
   if(val(row(i)[k])) {
    p= i;
    break;
   }
  if(p < 0) continue;
  if(p != k + 1) {
   std::swap_ranges(row(p), row(p) + w, row(k + 1));
   for(int i= 0; i < n; ++i) std::swap(row(i)[p], row(i)[k + 1]);
  }
  const Z iv= Z::raw(val(row(k + 1)[k])).inv();
  bool any= false;
  for(int i= k + 2; i < n; ++i) f[i]= (Z::raw(val(row(i)[k])) * iv).val(), any|= f[i] != 0;
  if(!any) continue;
  const u32* P= row(k + 1);
  const int jd= (k + 2) & ~7;
  for(int r= k + 2; r < n; ++r) {
   u32* x= row(r);
   u32 d;
   if(f[r]) {
    const u32 c= M - f[r], cs= u32((u64(c) << 32) / M);
    d= fused<Z>(x, P, f.data(), c, cs, k & ~7, w, R);
   } else d= dot<Z>(x, f.data(), jd, w, R);
   x[k + 1]= addm(x[k + 1], d);
  }
  for(int r= 0; r <= k + 1; ++r) {
   u32* x= row(r);
   x[k + 1]= addm(x[k + 1], dot<Z>(x, f.data(), jd, w, R));
  }
  std::fill(f.begin() + k + 2, f.end(), 0);
 }
 const int wp= (n + 1 + 3) & ~3;
 std::vector<u64> Pl((size_t)(n + 1) * wp), acc(wp), zero(wp);
 Pl[0]= 1;
 for(int k= 0; k < n; ++k) {
  const u64* pk= &Pl[(size_t)k * wp];
  const u32 hk= val(row(k)[k]), c= hk ? M - hk : 0;
  for(int j= 0; j < wp; ++j) acc[j]= (j ? pk[j - 1] : 0) + u64(c) * pk[j];
  const u64* Ub[K];
  u32 cf[K];
  int q= 0;
  Z t= Z::raw(1);
  const int jw= (k + 2 + 3) & ~3;
  for(int i= k - 1; i >= 0; --i) {
   t*= Z::raw(val(row(i + 1)[i]));
   if(t == Z()) break;
   const u32 coef= (t * Z::raw(val(row(i)[k]))).val();
   if(!coef) continue;
   Ub[q]= &Pl[(size_t)i * wp], cf[q++]= M - coef;
   if(q == K) update_row(acc.data(), Ub, cf, 0, jw, R), q= 0;
  }
  if(q) {
   for(int s= q; s < K; ++s) Ub[s]= zero.data(), cf[s]= 0;
   update_row(acc.data(), Ub, cf, 0, jw, R);
  }
  u64* pn= &Pl[(size_t)(k + 1) * wp];
  for(int j= 0; j <= k + 1; ++j) pn[j]= Z(acc[j]).val();
 }
 std::vector<Z> res(n + 1);
 for(int j= 0; j <= n; ++j) res[j]= Z::raw(u32(Pl[(size_t)n * wp + j]));
 return res;
}
template <class T> std::vector<T> field(Matrix<T> A) {
 const int n= A.h;
 const T zero= linalg_internal::zero<T>(), one= linalg_internal::one<T>();
 for(int k= 0; k + 2 < n; ++k) {
  int p= k + 1;
  if constexpr(std::is_floating_point_v<T>) {
   for(int i= k + 2; i < n; ++i)
    if(std::abs(A[i][k]) > std::abs(A[p][k])) p= i;
  } else
   while(p < n - 1 && A[p][k] == zero) ++p;
  if(A[p][k] == zero) continue;
  if(p != k + 1) {
   std::swap_ranges(A[p], A[p] + n, A[k + 1]);
   for(int i= 0; i < n; ++i) std::swap(A[i][p], A[i][k + 1]);
  }
  const T* P= A[k + 1];
  const T iv= one / P[k];
  for(int i= k + 2; i < n; ++i) {
   T* ai= A[i];
   const T fi= ai[k] * iv;
   if(fi == zero) continue;
   for(int j= k; j < n; ++j) ai[j]-= fi * P[j];
   for(int r= 0; r < n; ++r) A[r][k + 1]+= fi * A[r][i];
  }
 }
 std::vector<std::vector<T>> P(n + 1);
 P[0]= {one};
 for(int k= 0; k < n; ++k) {
  std::vector<T>& q= P[k + 1];
  q.assign(k + 2, zero);
  for(int j= 0; j <= k; ++j) q[j + 1]+= P[k][j], q[j]-= A[k][k] * P[k][j];
  T t= one;
  for(int i= k - 1; i >= 0; --i) {
   t*= A[i + 1][i];
   if(t == zero) break;
   const T coef= t * A[i][k];
   for(int j= 0; j <= i; ++j) q[j]-= coef * P[i][j];
  }
 }
 return P[n];
}
}
template <class T> std::vector<T> characteristic_polynomial(Matrix<T> A) {
 assert(A.h == A.w);
 if constexpr(linalg_internal::is_zmod_v<T>) return charpoly_internal::zmod(A);
 else return charpoly_internal::field(std::move(A));
}
