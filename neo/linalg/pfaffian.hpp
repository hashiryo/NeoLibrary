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
// procon-judge の yosupo-pfaffian-of-matrix で書き比べた zmod_pf_blk8 を写したもの。
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
 std::vector<u64> a((size_t)n * w), U((size_t)K * w), zero(w);
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
  if(A.h > 16) return pfaffian_internal::zmod(A);
 }
 return pfaffian_internal::field(std::move(A));
}
