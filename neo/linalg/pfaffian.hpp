#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>
#include "neo/algebra/ZMod.hpp"
#include "neo/linalg/Matrix.hpp"
// 交代行列 (A^T = -A) のパフィアン Pf(A) (Pf(A)^2 = det(A))。辺が奇数なら 0。T は体で、ZMod なら法は素数。
// 列 c = 0, 1, ..., n - 2 の順に、要を行 c + 1 以降から選んで行 c + 1 に置き (行と、同じ番号の列を入れ替えると符号が変わる)、
// 行 c + 2 以降から行 c + 1 の倍数を引いて列 c の値を消す。列の操作はしない。消し終えた行列 H は上 Hessenberg で、引いた操作を
// 下三角の L (対角は 1) とすると H = L A なので、L A L^T = H L^T は交代行列かつ上 Hessenberg、つまり 3 重対角で、その下の対角は
// H と同じになる。Pf(L A L^T) = det(L) Pf(A) = Pf(A) なので、Pf(A) は c が偶数の列の -H[c + 1][c] の積になる。c が奇数の列は
// 積に入らないので、要がなければ (列 c の行 c + 1 以降がすべて 0 なら) 何もせずに次の列へ進む (0 を返さない)。
// ZMod は、行を u64 で持ち、要を 8 本ずつ貯めて残りの各行を 1 回だけ読み書きする (遅らせた更新)。列を入れ替えるときは、残りの行と、
// 貯めている要の行の両方で入れ替える。n <= 16 の ZMod とほかの体は素朴な形 (浮動小数点数は絶対値が最大の要を選ぶ) にする。
// 法が奇数なら、要を探す列、要の行、行ごとの掛ける数を AVX2 で 4 個ずつ求める (zmod_v)。組の初めに、残りの行の次の 8 列を位置で
// 引く列ごとの配列 C に 4 × 4 の転置で写し、掛ける数 F も位置で引く列ごとの配列に持つ。位置 c + 1 と sel の行を入れ替えるときは
// C と F と今の列の値 V も入れ替え、列 c + 1 と sel が両方とも写した 8 列の中なら C の 2 列を入れ替え、c + 1 だけが中なら
// 入れ替えたあとの行から読み直す。値は Montgomery の還元で M 未満に直し、要の行を 2^{-32} 倍、掛ける数を 2^32 倍で持つので、
// 偶数の列で掛ける要の値も 2^{-32} 倍になり、最後に掛けた数だけ 2^32 を掛けて戻す。選ぶ要は zmod と同じ。
// procon-judge の yosupo-pfaffian-of-matrix で書き比べた zmod_pf_blk8 と zmod_pf_hpv を写したもの。
namespace pfaffian_internal {
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
template <class Z> Z zmod(const Matrix<Z>& A) {
 const int n= A.h, w= (n + 3) & ~3;
 const u32 M= Z::mod();
 const u64 R= (u64(1) << 32) % M;
 linalg_internal::hvec<u64> a((size_t)n * w);
 std::vector<u64> U((size_t)K * w), zero(w);
 for(int i= 0; i < n; ++i)
  for(int j= 0; j < n; ++j) a[(size_t)i * w + j]= A[i][j].val();
 std::vector<int> pos(n);  // pos[t] は位置 t にある行
 std::iota(pos.begin(), pos.end(), 0);
 std::vector<u32> F((size_t)n * K), v(n);
 const u64* Ub[K];
 Z res= Z::raw(1);
 int bt= 0, first= 1;  // first は残りの行の最初の位置
 auto flush= [&](int col) {
  if(!bt) return;
  for(int s= bt; s < K; ++s) Ub[s]= zero.data();
  for(int t= first; t < n; ++t) {
   u32* fr= &F[(size_t)pos[t] * K];
   u32 any= 0;
   for(int s= 0; s < K; ++s) {
    if(s >= bt) fr[s]= 0;
    any|= fr[s];
   }
   if(any) update_row(&a[(size_t)pos[t] * w], Ub, fr, col & ~3, w, R);
  }
  bt= 0;
 };
 for(int c= 0; c + 1 < n; ++c) {
  int sel= -1;
  for(int t= first; t < n; ++t) {
   const int r= pos[t];
   const u32* fr= &F[(size_t)r * K];
   u64 x= a[(size_t)r * w + c];
   for(int s= 0; s < bt; ++s) x+= u64(fr[s]) * Ub[s][c];
   v[t]= Z(x).val();
   if(sel < 0 && v[t]) sel= t;
  }
  if(sel < 0) {  // 列 c が消えている。c が奇数なら H[c + 1][c] = 0 は積に入らないので、そのまま次の列へ進む
   if(!(c & 1)) return Z();
   first= c + 2;
   continue;
  }
  if(sel != c + 1) {
   std::swap(pos[c + 1], pos[sel]), std::swap(v[c + 1], v[sel]);
   for(int t= c + 1; t < n; ++t) {
    u64* row= &a[(size_t)pos[t] * w];
    std::swap(row[c + 1], row[sel]);
   }
   for(int s= 0; s < bt; ++s) std::swap(U[(size_t)s * w + c + 1], U[(size_t)s * w + sel]);
   res= -res;
  }
  const int p= pos[c + 1];
  if(!(c & 1)) res*= -Z::raw(v[c + 1]);
  u64* up= &U[(size_t)bt * w];
  const u32* fp= &F[(size_t)p * K];
  const u64* ap= &a[(size_t)p * w];
  for(int j= c & ~3; j < w; ++j) {
   u64 x= ap[j];
   for(int s= 0; s < bt; ++s) x+= u64(fp[s]) * Ub[s][j];
   up[j]= Z(x).val();
  }
  const Z iv= Z::raw(v[c + 1]).inv();
  for(int t= c + 2; t < n; ++t) F[(size_t)pos[t] * K + bt]= (-(Z::raw(v[t]) * iv)).val();
  Ub[bt++]= up;
  first= c + 2;
  if(bt == K) flush(c + 1);
 }
 return res;
}
// 4 個の u64 を、2^32 の位を R = 2^32 mod M に掛けて畳んでから Montgomery の還元で 2^{-32} 倍し、M 未満に直す (M は奇数)。
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
// u64 の 4 個の下の 32 bit を u32 の 4 個に詰める。
[[gnu::always_inline]] inline __m128i pack_lo(__m256i x) { return _mm256_castsi256_si128(_mm256_permutevar8x32_epi32(x, _mm256_setr_epi32(0, 2, 4, 6, 0, 0, 0, 0))); }
template <class Z> Z zmod_v(const Matrix<Z>& A) {
 const int n= A.h, w= (n + 3) & ~3, n4= w;
 const u32 M= Z::mod();
 const u64 R= (u64(1) << 32) % M;
 const Mont mt(M);
 linalg_internal::hvec<u64> a((size_t)n * w + 8);  // 写した 8 列を読むとき、最後の行の後ろを 8 個まで読む
 for(int i= 0; i < n; ++i) {
  u64* r= &a[(size_t)i * w];
  for(int j= 0; j < n; ++j) r[j]= A[i][j].val();
  for(int j= n; j < w; ++j) r[j]= 0;
 }
 std::fill(a.begin() + (size_t)n * w, a.end(), 0);
 auto row= [&](int r) { return &a[(size_t)r * w]; };
 std::vector<int> pos(n);  // pos[t] は位置 t にある行
 std::iota(pos.begin(), pos.end(), 0);
 std::vector<u64> C((size_t)K * n4), U((size_t)K * w + 8), zero(w + 8);
 std::vector<u32> F((size_t)K * n4), V(n4);
 const u64* Ub[K];
 Z res= Z::raw(1);
 int bt= 0, first= 1, c0= -1000, ne= 0;  // first は残りの行の最初の位置、c0 は C に写した最初の列、ne は掛けた要の値の数
 auto fill_cache= [&](int cs) {
  c0= cs;
  for(int t= first & ~3; t < n; t+= 4) {
   const u64* r0= row(pos[t]) + cs;
   const u64* r1= t + 1 < n ? row(pos[t + 1]) + cs : zero.data();
   const u64* r2= t + 2 < n ? row(pos[t + 2]) + cs : zero.data();
   const u64* r3= t + 3 < n ? row(pos[t + 3]) + cs : zero.data();
   for(int h= 0; h < 8; h+= 4) {
    const __m256i x0= _mm256_loadu_si256((const __m256i*)(r0 + h)), x1= _mm256_loadu_si256((const __m256i*)(r1 + h));
    const __m256i x2= _mm256_loadu_si256((const __m256i*)(r2 + h)), x3= _mm256_loadu_si256((const __m256i*)(r3 + h));
    const __m256i t0= _mm256_unpacklo_epi64(x0, x1), t1= _mm256_unpackhi_epi64(x0, x1);
    const __m256i t2= _mm256_unpacklo_epi64(x2, x3), t3= _mm256_unpackhi_epi64(x2, x3);
    _mm256_storeu_si256((__m256i*)&C[(size_t)(h + 0) * n4 + t], _mm256_permute2x128_si256(t0, t2, 0x20));
    _mm256_storeu_si256((__m256i*)&C[(size_t)(h + 1) * n4 + t], _mm256_permute2x128_si256(t1, t3, 0x20));
    _mm256_storeu_si256((__m256i*)&C[(size_t)(h + 2) * n4 + t], _mm256_permute2x128_si256(t0, t2, 0x31));
    _mm256_storeu_si256((__m256i*)&C[(size_t)(h + 3) * n4 + t], _mm256_permute2x128_si256(t1, t3, 0x31));
   }
  }
 };
 auto flush= [&](int col) {
  if(!bt) return;
  for(int s= bt; s < K; ++s) Ub[s]= zero.data();
  for(int t= first; t < n; ++t) {
   u32 f[K], any= 0;
   for(int s= 0; s < K; ++s) any|= (f[s]= s < bt ? F[(size_t)s * n4 + t] : 0);
   if(any) update_row(row(pos[t]), Ub, f, col & ~3, w, R);
  }
  for(int s= 0; s < bt; ++s) std::fill(&F[(size_t)s * n4], &F[(size_t)s * n4] + n4, 0u);
  bt= 0, c0= -1000;  // 行が変わったので写し直す
 };
 for(int c= 0; c + 1 < n; ++c) {
  if(c < c0 || c >= c0 + 8) fill_cache(c);
  const u64* cc= &C[(size_t)(c - c0) * n4];
  __m256i us[K];
  for(int s= 0; s < bt; ++s) us[s]= _mm256_set1_epi64x(Ub[s][c]);
  int sel= -1;
  for(int t= first & ~3; t < n4; t+= 4) {
   __m256i x= _mm256_loadu_si256((const __m256i*)(cc + t));
   for(int s= 0; s < bt; ++s) x= _mm256_add_epi64(x, _mm256_mul_epu32(_mm256_cvtepu32_epi64(_mm_loadu_si128((const __m128i*)&F[(size_t)s * n4 + t])), us[s]));
   const __m256i y= mt.redc(x);
   _mm_storeu_si128((__m128i*)&V[t], pack_lo(y));
   if(sel < 0) {
    int nz= ~_mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpeq_epi64(y, _mm256_setzero_si256()))) & 15;
    if(t < first) nz&= ~((1 << (first - t)) - 1);  // 位置 first より前は終わった行
    if(nz) sel= t + __builtin_ctz(nz);
   }
  }
  if(sel >= n) sel= -1;
  if(sel < 0) {  // 列 c が消えている。c が奇数なら H[c + 1][c] = 0 は積に入らないので、そのまま次の列へ進む
   if(!(c & 1)) return Z();
   first= c + 2;
   continue;
  }
  if(sel != c + 1) {
   std::swap(pos[c + 1], pos[sel]), std::swap(V[c + 1], V[sel]);
   for(int q= 0; q < K; ++q) std::swap(C[(size_t)q * n4 + c + 1], C[(size_t)q * n4 + sel]);
   for(int s= 0; s < bt; ++s) std::swap(F[(size_t)s * n4 + c + 1], F[(size_t)s * n4 + sel]);
   for(int t= c + 1; t < n; ++t) {
    u64* r= row(pos[t]);
    std::swap(r[c + 1], r[sel]);
   }
   for(int s= 0; s < bt; ++s) std::swap(U[(size_t)s * w + c + 1], U[(size_t)s * w + sel]);
   const int q1= c + 1 - c0, q2= sel - c0;
   if(q1 < 8 && q2 < 8) {
    for(int t= first & ~3; t < n4; ++t) std::swap(C[(size_t)q1 * n4 + t], C[(size_t)q2 * n4 + t]);
   } else if(q1 < 8) {
    for(int t= first; t < n; ++t) C[(size_t)q1 * n4 + t]= row(pos[t])[c + 1];
   }
   res= -res;
  }
  const int p= pos[c + 1];
  if(!(c & 1)) res*= -Z::raw(V[c + 1]), ++ne;
  u64* up= &U[(size_t)bt * w];
  {  // 要の行 U[j] = (a[p][j] + sum_s F[s][c + 1] U[s][j]) 2^{-32} mod M
   __m256i fs[K];
   for(int s= 0; s < bt; ++s) fs[s]= _mm256_set1_epi64x(F[(size_t)s * n4 + c + 1]);
   const u64* ap= row(p);
   for(int j= c & ~3; j < w; j+= 4) {
    __m256i x= _mm256_loadu_si256((const __m256i*)(ap + j));
    for(int s= 0; s < bt; ++s) x= _mm256_add_epi64(x, _mm256_mul_epu32(_mm256_loadu_si256((const __m256i*)(Ub[s] + j)), fs[s]));
    _mm256_storeu_si256((__m256i*)(up + j), mt.redc(x));
   }
  }
  {  // F[bt][t] = -V[t] cinv mod M (t >= c + 2、cinv = V[c + 1]^{-1} 2^32)。Shoup の方法で 4 個ずつ
   const u32 ci= (Z::raw(V[c + 1]).inv() * Z::raw(u32(R))).val(), cs= u32((u64(ci) << 32) / M);
   const __m256i Cv= _mm256_set1_epi64x(ci), CS= _mm256_set1_epi64x(cs);
   u32* fb= &F[(size_t)bt * n4];
   for(int t= (c + 2) & ~3; t < n4; t+= 4) {
    const __m256i vv= _mm256_cvtepu32_epi64(_mm_loadu_si128((const __m128i*)&V[t]));
    const __m256i qq= _mm256_srli_epi64(_mm256_mul_epu32(vv, CS), 32);
    __m256i x= _mm256_sub_epi64(_mm256_mul_epu32(vv, Cv), _mm256_mul_epu32(qq, mt.M));
    x= _mm256_min_epu32(x, _mm256_sub_epi32(x, mt.M32));
    x= _mm256_sub_epi32(mt.M32, x);
    x= _mm256_min_epu32(x, _mm256_sub_epi32(x, mt.M32));
    _mm_storeu_si128((__m128i*)(fb + t), pack_lo(x));
   }
  }
  Ub[bt++]= up;
  first= c + 2;
  if(bt == K) flush(c + 1);
 }
 return res * Z(R).pow(ne);  // 掛けた要の値は 2^{-32} 倍
}
template <class T> T field(Matrix<T> A) {
 const int n= A.h;
 const T zero= linalg_internal::zero<T>(), one= linalg_internal::one<T>();
 T res= one;
 for(int i= 1; i < n; ++i) {
  int sel= i;
  if constexpr(std::is_floating_point_v<T>) {
   for(int j= i + 1; j < n; ++j)
    if(std::abs(A[j][i - 1]) > std::abs(A[sel][i - 1])) sel= j;
  } else
   while(sel < n - 1 && A[sel][i - 1] == zero) ++sel;
  if(A[sel][i - 1] == zero) {
   if(i & 1) return zero;
   continue;
  }
  if(sel != i) {
   std::swap_ranges(A[i] + i - 1, A[i] + n, A[sel] + i - 1);
   for(int k= i; k < n; ++k) std::swap(A[k][i], A[k][sel]);
   res= zero - res;
  }
  const T* ai= A[i];
  if(i & 1) res*= zero - ai[i - 1];
  const T iv= one / ai[i - 1];
  for(int j= i + 1; j < n; ++j) {
   T* aj= A[j];
   const T f= aj[i - 1] * iv;
   if(f == zero) continue;
   for(int k= i - 1; k < n; ++k) aj[k]-= f * ai[k];
  }
 }
 return res;
}
}
template <class T> T pfaffian(Matrix<T> A) {
 assert(A.h == A.w);
 if(A.h & 1) return linalg_internal::zero<T>();
 if constexpr(linalg_internal::is_zmod_v<T>) {
  if(A.h > 16) return T::mod() & 1 ? pfaffian_internal::zmod_v(A) : pfaffian_internal::zmod(A);
 }
 return pfaffian_internal::field(std::move(A));
}
