#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <cmath>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>
#include "neo/algebra/ZMod.hpp"
#include "neo/linalg/Matrix.hpp"
// 行列の階数。T は体で、ZMod なら法は素数、浮動小数点数なら絶対値が 1e-8 以下の値を 0 とみなす。
// ZMod は、行の数が min(h, w) になる向きに置き (縦長なら転置する)、行を u64 で持って、要を 8 本ずつ貯めて残りの各行を
// 1 回だけ読み書きする (遅らせた更新)。積は M^2 < 2^60 なので、2^62 + 2^32 未満の値に 8 本足しても u64 に収まり、足したあとは
// 2^32 の位を 2^32 mod M に掛けて畳む。要を探す列の値は、貯めている分をその場で足してから M 未満に直す。
// 法が奇数なら、要を探す列、要の行、行ごとの掛ける数を AVX2 で 4 個ずつ求める (eliminate_v)。組の初めに残りの行の次の 8 列を
// 列ごとの配列 (残りの行の並びの順) に 4 × 4 の転置で写し、要を探す列はそこから読む。行ごとの掛ける数も残りの行の並びの順に
// 列ごとに持ち、要の行を抜くときは末尾の行の分をその位置へ移す。値は 2^32 の位を畳んでから Montgomery の還元で M 未満に直し、
// 要の行を 2^{-32} 倍、掛ける数を 2^32 倍で持つので、積は本当の値の積になり、残りの行は本当の値のまま足される。選ぶ要は
// eliminate と同じ (残りの行の並びで値が 0 でない最初の行)。行列は 2 MB 境界の mmap に置いて huge page を頼む (hvec)。
// max(h, w) <= 16 の ZMod とほかの体は、素朴な掃き出し (浮動小数点数は絶対値が最大の要を選ぶ) にする。
// procon-judge の yosupo-matrix-rank で書き比べた zmod_rank_blk8 を写し、行列を huge page に置いた zmod_rank_hp と、要を探す
// 列をベクトル化した zmod_rank_hpv の形にしたもの。
namespace rank_internal {
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
// n 行 w 列 (w は 4 の倍数) の u64 の行列 a の列 [0, lim) で要を探しながら前進消去し、要の数を返す。要はその列で値が 0 でない
// 最初の行で、まだ要にしていない行の並び rem は、要を抜くたびに末尾の行で埋める。
template <class Z> int eliminate(linalg_internal::hvec<u64>& a, int n, int w, int lim) {
 const u32 M= Z::mod();
 const u64 R= (u64(1) << 32) % M;
 std::vector<int> rem(n);
 std::iota(rem.begin(), rem.end(), 0);
 std::vector<u64> U((size_t)std::min(n, K) * w), zero;  // zero (0 の行) は、組が 8 本に満たないまま足すときに作る
 std::vector<u32> F((size_t)n * K), v(n);
 const u64* Ub[K];
 int rk= 0, col= 0;
 while(col < lim && !rem.empty()) {
  int bt= 0;
  for(; bt < K && col < lim && !rem.empty(); ++col) {
   const int nr= rem.size();
   int pi= -1;
   for(int i= 0; i < nr; ++i) {
    const int r= rem[i];
    const u32* fr= &F[(size_t)r * K];
    u64 x= a[(size_t)r * w + col];
    for(int s= 0; s < bt; ++s) x+= u64(fr[s]) * Ub[s][col];
    v[i]= Z(x).val();
    if(pi < 0 && v[i]) pi= i;
   }
   if(pi < 0) continue;
   const int p= rem[pi];
   u64* up= &U[(size_t)bt * w];
   const u32* fp= &F[(size_t)p * K];
   const u64* ap= &a[(size_t)p * w];
   for(int j= col & ~3; j < w; ++j) {
    u64 x= ap[j];
    for(int s= 0; s < bt; ++s) x+= u64(fp[s]) * Ub[s][j];
    up[j]= Z(x).val();
   }
   const Z iv= Z::raw(v[pi]).inv();
   for(int i= 0; i < nr; ++i) F[(size_t)rem[i] * K + bt]= (-(Z::raw(v[i]) * iv)).val();
   Ub[bt++]= up, ++rk;
   rem[pi]= rem.back(), rem.pop_back();
  }
  if(col >= lim || rem.empty()) break;  // 残りの行をもう読まないので、最後の組は足さない
  if(bt < K && zero.empty()) zero.assign(w, 0);
  for(int s= bt; s < K; ++s) Ub[s]= zero.data();
  const int jb= col & ~3;
  for(int r: rem) {
   u32* fr= &F[(size_t)r * K];
   for(int s= bt; s < K; ++s) fr[s]= 0;
   u32 any= 0;
   for(int s= 0; s < K; ++s) any|= fr[s];
   if(any) update_row(&a[(size_t)r * w], Ub, fr, jb, w, R);  // 掛ける数が全部 0 の行 (疎な行列で多い) は飛ばす
  }
 }
 return rk;
}
// 4 個の u64 (2^63 未満) を、2^32 の位を R = 2^32 mod M に掛けて畳んでから Montgomery の還元で 2^{-32} 倍し、M 未満に直す。
struct Mont {
 __m256i R, LO, NIV, M, M32;
 explicit Mont(u32 m) {
  u32 x= m;
  for(int i= 0; i < 5; ++i) x*= 2 - m * x;
  R= _mm256_set1_epi64x((u64(1) << 32) % m), LO= _mm256_set1_epi64x(0xffffffffll), NIV= _mm256_set1_epi64x(u32(-x));
  M= _mm256_set1_epi64x(m), M32= _mm256_set1_epi32(m);
 }
 [[gnu::always_inline]] __m256i redc(__m256i x) const {
  x= _mm256_add_epi64(_mm256_mul_epu32(_mm256_srli_epi64(x, 32), R), _mm256_and_si256(x, LO));
  const __m256i t= _mm256_srli_epi64(_mm256_add_epi64(x, _mm256_mul_epu32(_mm256_mul_epu32(x, NIV), M)), 32);
  return _mm256_min_epu32(t, _mm256_sub_epi32(t, M32));  // 上の 32 bit は 0 のまま (0 - M は大きいので min で 0)
 }
};
[[gnu::always_inline]] inline __m128i pack_lo(__m256i x) { return _mm256_castsi256_si128(_mm256_permutevar8x32_epi32(x, _mm256_setr_epi32(0, 2, 4, 6, 0, 0, 0, 0))); }
// eliminate と同じ前進消去を、要を探す列、要の行、行ごとの掛ける数を 4 個ずつ求めて行う。法は奇数。a は n * w + 8 個 (行の後ろを
// 8 個まで読む)。
template <class Z> int eliminate_v(linalg_internal::hvec<u64>& a, int n, int w, int lim) {
 const u32 Mod= Z::mod();
 const u64 R= (u64(1) << 32) % Mod;
 const Mont mt(Mod);
 const int n4= (n + 3) & ~3;
 std::vector<int> rem(n);
 std::iota(rem.begin(), rem.end(), 0);
 // C[c][i]: 残りの行の並びで i 番目の行の、列 c0 + c の値。F[s][i]: その行の s 番目の要の掛ける数 (2^32 倍)。V[i]: 今の列の値。
 std::vector<u64> C((size_t)K * n4), U((size_t)std::min(n, K) * w + 8), zero(w + 8);
 std::vector<u32> F((size_t)K * n4), V(n4);
 const u64* Ub[K];
 auto row= [&](int r) { return a.data() + (size_t)r * w; };
 auto fill_cache= [&](int c0) {  // 残りの行の列 [c0, c0 + 8) を C に写す。行の後ろを読んでも、使わない列に入るだけ
  const int nr= rem.size();
  for(int i= 0; i < nr; i+= 4) {
   const u64* r0= row(rem[i]) + c0;
   const u64* r1= i + 1 < nr ? row(rem[i + 1]) + c0 : zero.data();
   const u64* r2= i + 2 < nr ? row(rem[i + 2]) + c0 : zero.data();
   const u64* r3= i + 3 < nr ? row(rem[i + 3]) + c0 : zero.data();
   for(int h= 0; h < 8; h+= 4) {
    const __m256i x0= _mm256_loadu_si256((const __m256i*)(r0 + h)), x1= _mm256_loadu_si256((const __m256i*)(r1 + h));
    const __m256i x2= _mm256_loadu_si256((const __m256i*)(r2 + h)), x3= _mm256_loadu_si256((const __m256i*)(r3 + h));
    const __m256i t0= _mm256_unpacklo_epi64(x0, x1), t1= _mm256_unpackhi_epi64(x0, x1);
    const __m256i t2= _mm256_unpacklo_epi64(x2, x3), t3= _mm256_unpackhi_epi64(x2, x3);
    _mm256_storeu_si256((__m256i*)&C[(size_t)(h + 0) * n4 + i], _mm256_permute2x128_si256(t0, t2, 0x20));
    _mm256_storeu_si256((__m256i*)&C[(size_t)(h + 1) * n4 + i], _mm256_permute2x128_si256(t1, t3, 0x20));
    _mm256_storeu_si256((__m256i*)&C[(size_t)(h + 2) * n4 + i], _mm256_permute2x128_si256(t0, t2, 0x31));
    _mm256_storeu_si256((__m256i*)&C[(size_t)(h + 3) * n4 + i], _mm256_permute2x128_si256(t1, t3, 0x31));
   }
  }
 };
 int rk= 0, col= 0;
 while(col < lim && !rem.empty()) {
  int bt= 0, c0= col;
  fill_cache(c0);
  for(; bt < K && col < lim && !rem.empty(); ++col) {
   if(col == c0 + 8) fill_cache(c0= col);
   const int nr= rem.size(), nr4= (nr + 3) & ~3;
   const u64* cc= &C[(size_t)(col - c0) * n4];
   __m256i us[K];
   for(int s= 0; s < bt; ++s) us[s]= _mm256_set1_epi64x(Ub[s][col]);
   int pi= -1;
   for(int i= 0; i < nr4; i+= 4) {  // V[i] = (C + sum_s F[s] U[s][col]) 2^{-32} mod M と、0 でない最初の i
    __m256i x= _mm256_loadu_si256((const __m256i*)(cc + i));
    for(int s= 0; s < bt; ++s) x= _mm256_add_epi64(x, _mm256_mul_epu32(_mm256_cvtepu32_epi64(_mm_loadu_si128((const __m128i*)&F[(size_t)s * n4 + i])), us[s]));
    const __m256i y= mt.redc(x);
    _mm_storeu_si128((__m128i*)&V[i], pack_lo(y));
    if(pi < 0) {
     const int nz= ~_mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpeq_epi64(y, _mm256_setzero_si256()))) & 15;
     if(nz) pi= i + __builtin_ctz(nz);
    }
   }
   if(pi >= nr) pi= -1;  // 余りの位置 (値は 0) は選ばない
   if(pi < 0) continue;
   const int p= rem[pi];
   u64* up= &U[(size_t)bt * w];
   {  // 要の行 U[j] = (a[p][j] + sum_s F[s][pi] U[s][j]) 2^{-32} mod M
    __m256i fs[K];
    for(int s= 0; s < bt; ++s) fs[s]= _mm256_set1_epi64x(F[(size_t)s * n4 + pi]);
    const u64* ap= row(p);
    for(int j= col & ~3; j < w; j+= 4) {
     __m256i x= _mm256_loadu_si256((const __m256i*)(ap + j));
     for(int s= 0; s < bt; ++s) x= _mm256_add_epi64(x, _mm256_mul_epu32(_mm256_loadu_si256((const __m256i*)(Ub[s] + j)), fs[s]));
     _mm256_storeu_si256((__m256i*)(up + j), mt.redc(x));
    }
   }
   {  // F[bt][i] = -V[i] c mod M (c = V[pi]^{-1} 2^32)。Shoup の方法で 4 個ずつ。
    const u32 c= (Z::raw(V[pi]).inv() * Z::raw(u32(R))).val(), cs= u32((u64(c) << 32) / Mod);
    const __m256i Cv= _mm256_set1_epi64x(c), CS= _mm256_set1_epi64x(cs);
    u32* fb= &F[(size_t)bt * n4];
    for(int i= 0; i < nr4; i+= 4) {
     const __m256i vv= _mm256_cvtepu32_epi64(_mm_loadu_si128((const __m128i*)&V[i]));
     const __m256i q= _mm256_srli_epi64(_mm256_mul_epu32(vv, CS), 32);
     __m256i t= _mm256_sub_epi64(_mm256_mul_epu32(vv, Cv), _mm256_mul_epu32(q, mt.M));  // [0, 2M)
     t= _mm256_min_epu32(t, _mm256_sub_epi32(t, mt.M32));
     t= _mm256_sub_epi32(mt.M32, t);
     t= _mm256_min_epu32(t, _mm256_sub_epi32(t, mt.M32));
     _mm_storeu_si128((__m128i*)(fb + i), pack_lo(t));
    }
   }
   Ub[bt++]= up, ++rk;
   const int last= nr - 1;  // 要の行を抜き、末尾の行の分を位置 pi へ移して、末尾を 0 にする
   rem[pi]= rem[last], rem.pop_back();
   for(int s= 0; s < bt; ++s) F[(size_t)s * n4 + pi]= F[(size_t)s * n4 + last], F[(size_t)s * n4 + last]= 0;
   for(int q= 0; q < K; ++q) C[(size_t)q * n4 + pi]= C[(size_t)q * n4 + last], C[(size_t)q * n4 + last]= 0;
  }
  if(!bt || col >= lim || rem.empty()) break;  // 残りの行をもう読まないので、最後の組は足さない
  for(int s= bt; s < K; ++s) Ub[s]= zero.data();
  const int jb= col & ~3, nr= rem.size();
  for(int i= 0; i < nr; ++i) {
   u32 f[K], any= 0;
   for(int s= 0; s < K; ++s) any|= (f[s]= s < bt ? F[(size_t)s * n4 + i] : 0);
   if(any) update_row(row(rem[i]), Ub, f, jb, w, R);
  }
  for(int s= 0; s < bt; ++s) std::fill(&F[(size_t)s * n4], &F[(size_t)s * n4] + n4, 0u);
 }
 return rk;
}
template <class Z> int large(const Matrix<Z>& A) {
 const bool tr= A.h > A.w;
 const int n= tr ? A.w : A.h, m= tr ? A.h : A.w, w= (m + 3) & ~3;
 linalg_internal::hvec<u64> a((size_t)n * w + 8);
 for(int i= 0; i < A.h; ++i)
  for(int j= 0; j < A.w; ++j) a[tr ? (size_t)j * w + i : (size_t)i * w + j]= A[i][j].val();
 return Z::mod() & 1 ? eliminate_v<Z>(a, n, w, m) : eliminate<Z>(a, n, w, m);
}
template <class T> int field(Matrix<T> A) {
 const int h= A.h, w= A.w;
 const T zero= linalg_internal::zero<T>(), one= linalg_internal::one<T>();
 int r= 0;
 for(int c= 0; c < w && r < h; ++c) {
  int p= r;
  if constexpr(std::is_floating_point_v<T>) {
   for(int i= r + 1; i < h; ++i)
    if(std::abs(A[i][c]) > std::abs(A[p][c])) p= i;
   if(std::abs(A[p][c]) <= T(1e-8)) continue;
  } else {
   while(p < h && A[p][c] == zero) ++p;
   if(p == h) continue;
  }
  if(p != r) std::swap_ranges(A[p] + c, A[p] + w, A[r] + c);
  const T* pr= A[r];
  const T iv= one / pr[c];
  for(int i= r + 1; i < h; ++i) {
   T* ai= A[i];
   const T f= ai[c] * iv;
   if(f == zero) continue;
   for(int j= c + 1; j < w; ++j) ai[j]-= f * pr[j];
  }
  ++r;
 }
 return r;
}
}
// 名前を rank にすると、using namespace std; のもとで std::rank (型の特性) とぶつかって rank(A) が曖昧になるので、matrix_rank にする。
template <class T> int matrix_rank(Matrix<T> A) {
 if constexpr(linalg_internal::is_zmod_v<T>) {
  if(std::max(A.h, A.w) > 16) return rank_internal::large(A);
 }
 return rank_internal::field(std::move(A));
}
