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
// max(h, w) <= 16 の ZMod とほかの体は、素朴な掃き出し (浮動小数点数は絶対値が最大の要を選ぶ) にする。
// procon-judge の yosupo-matrix-rank で書き比べた zmod_rank_blk8 を写したもの。
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
template <class Z> int eliminate(std::vector<u64>& a, int n, int w, int lim) {
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
template <class Z> int large(const Matrix<Z>& A) {
 const bool tr= A.h > A.w;
 const int n= tr ? A.w : A.h, m= tr ? A.h : A.w, w= (m + 3) & ~3;
 std::vector<u64> a((size_t)n * w);
 for(int i= 0; i < A.h; ++i)
  for(int j= 0; j < A.w; ++j) a[tr ? (size_t)j * w + i : (size_t)i * w + j]= A[i][j].val();
 return eliminate<Z>(a, n, w, m);
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
// using namespace std; のもとでは std::rank (型の特性) と名前がぶつかるので、::rank(A) と書く。
template <class T> int rank(Matrix<T> A) {
 if constexpr(linalg_internal::is_zmod_v<T>) {
  if(std::max(A.h, A.w) > 16) return rank_internal::large(A);
 }
 return rank_internal::field(std::move(A));
}
