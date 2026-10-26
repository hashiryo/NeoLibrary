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
#include "neo/number_theory/inv_gcd.hpp"
// 正方行列の行列式。T は ZMod (法は合成数でもよい)、浮動小数点数、割り算を持つ体のどれか。
// ZMod は、列に単元 (法と互いに素な値) があればそれを要にして体と同じように掃き出し、単元のない列だけ、その列の値が 0 でない
// 行どうしに互除法の 2×2 の変換 (行列式 1) を掛けて値を 1 行に集め、その行を要にする (要の値は単元でなくてよい)。法が素数なら
// 0 でない値はどれも単元なので互除法の道は通らず、素数かどうかを調べる必要もない。
// n > 16 では、行を u64 で持ち、要を 8 本ずつ貯めて、8 本たまったら残りの各行を 1 回だけ読み書きして 8 本分の積をまとめて足す
// (遅らせた更新)。積は M^2 < 2^60 なので、2^62 + 2^32 未満の値に 8 本足しても u64 に収まり、足したあとは 2^32 の位を
// 2^32 mod M に掛けて畳む。要を探す列の値は、貯めている分をその場で足してから M 未満に直す。
// n <= 16 では準備をせず、行を u64 で持って要の行を足すだけにする。1 行に足すのは 15 回以下なので、畳まなくても u64 に収まる。
// 浮動小数点数は絶対値が最大の要を選ぶ。ほかの体は 0 でない最初の要で割る。
// procon-judge の yosupo-matrix-det と yosupo-matrix-det-arbitrary-mod で書き比べた zmod_det_blk8 を写したもの。
namespace det_internal {
using u32= unsigned;
using u64= unsigned long long;
using i64= long long;
constexpr int K= 8;
// row[j] += sum_{s<8} f[s] U[s][j] を j in [jb, w) で 4 個ずつ足してから畳む。U[s][j] < M、f[s] < M、row[j] < 2^62 + 2^32。
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
// x (0 < x < M) が単元なら逆元を、そうでなければ 0 を返す。
inline u32 unit_inv(u32 x, u32 M) {
 const auto [g, y]= inv_gcd(x, M);
 return g == 1 ? y : 0;
}
// 置換 p (t 番目の要が行 p[t]) が奇置換か。
inline bool odd_perm(const std::vector<int>& p) {
 const int n= p.size();
 std::vector<char> seen(n);
 int cyc= 0;
 for(int i= 0; i < n; ++i)
  if(!seen[i]) {
   ++cyc;
   for(int j= i; !seen[j]; j= p[j]) seen[j]= 1;
  }
 return (n - cyc) & 1;
}
// 行 P と行 J (列 [c, e) の値は M 未満、P[c] と J[c] は 0 でない) に、x P[c] + y J[c] = g = gcd(P[c], J[c]) となる整数 x, y で
// [[x, y], [-J[c] / g, P[c] / g]] (行列式 1) を掛ける。P[c] は g に、J[c] は 0 になる。
template <class Z> void euclid(u64* P, u64* J, int c, int e) {
 const i64 M= Z::mod();
 i64 g= P[c], b= J[c], x0= 1, x1= 0, y0= 0, y1= 1;
 while(b) {
  const i64 q= g / b;
  g-= q * b, std::swap(g, b);
  x0-= q * x1, std::swap(x0, x1);
  y0-= q * y1, std::swap(y0, y1);
 }
 auto md= [&](i64 t) {
  const i64 r= t % M;
  return u64(r < 0 ? r + M : r);
 };
 const u64 cx= md(x0), cy= md(y0), cu= md(-i64(J[c]) / g), cv= md(i64(P[c]) / g);
 for(int k= c; k < e; ++k) {
  const u64 pk= P[k], jk= J[k];
  P[k]= Z(cx * pk + cy * jk).val(), J[k]= Z(cu * pk + cv * jk).val();
 }
}
// n <= 16。行を u64 で持ち、要の行 (M 未満に直した u32) の f 倍を足す。足す回数が 15 以下なので畳まない。
template <class Z> Z small(const Matrix<Z>& A) {
 const int n= A.h, w= (n + 3) & ~3;
 const u32 M= Z::mod();
 u64 a[16][16];
 u32 pv[16];
 for(int i= 0; i < n; ++i)
  for(int j= 0; j < w; ++j) a[i][j]= j < n ? A[i][j].val() : 0;
 Z d= Z::raw(1);
 bool neg= false;
 for(int c= 0; c < n; ++c) {
  int p= -1, nz= -1;
  u32 iv= 0;
  for(int r= c; r < n; ++r) {
   const u32 x= Z(a[r][c]).val();
   a[r][c]= x;
   if(x && p < 0) {
    if(nz < 0) nz= r;
    if((iv= unit_inv(x, M))) p= r;
   }
  }
  if(nz < 0) return Z();
  if(p < 0) {  // 単元のない列。値が 0 でない行を、互除法で行 nz に集める
   p= nz;
   for(int j= c + 1; j < n; ++j) a[p][j]= Z(a[p][j]).val();
   for(int r= p + 1; r < n; ++r)
    if(a[r][c]) {
     for(int j= c + 1; j < n; ++j) a[r][j]= Z(a[r][j]).val();
     euclid<Z>(a[p], a[r], c, n);
    }
  }
  if(p != c) std::swap(a[p], a[c]), neg= !neg;
  for(int j= 0; j < w; ++j) pv[j]= j < c ? 0 : Z(a[c][j]).val();
  d*= Z::raw(pv[c]);
  if(!iv) {
   if(d == Z()) return Z();
   continue;
  }
  const Z ivz= Z::raw(iv);
  for(int r= c + 1; r < n; ++r) {
   const u32 x= u32(a[r][c]);
   if(!x) continue;
   const u32 f= (-(Z::raw(x) * ivz)).val();
   u64* row= a[r];
   for(int j= c & ~3; j < w; j+= 4)
    for(int l= 0; l < 4; ++l) row[j + l]+= u64(f) * pv[j + l];
  }
 }
 return neg ? -d : d;
}
// n > 16。要を 8 本ずつ貯める。残りの行の並び rem は、要を抜くたびに末尾の行で埋める。
template <class Z> Z large(const Matrix<Z>& A) {
 const int n= A.h, w= (n + 3) & ~3;
 const u32 M= Z::mod();
 const u64 R= (u64(1) << 32) % M;
 std::vector<u64> a((size_t)n * w), U((size_t)K * w), zero(w);
 for(int i= 0; i < n; ++i)
  for(int j= 0; j < n; ++j) a[(size_t)i * w + j]= A[i][j].val();
 std::vector<int> rem(n), prow;
 std::iota(rem.begin(), rem.end(), 0);
 prow.reserve(n);
 std::vector<u32> F((size_t)n * K), v(n);
 const u64* Ub[K];
 Z d= Z::raw(1);
 int bt= 0;
 // 貯めた bt 本の更新を、残りの行の列 [col を 4 の倍数に切り下げた位置, w) に足す。
 auto flush= [&](int col) {
  if(!bt) return;
  for(int s= bt; s < K; ++s) Ub[s]= zero.data();
  for(int r: rem) {
   u32* fr= &F[(size_t)r * K];
   u32 any= 0;
   for(int s= 0; s < K; ++s) {
    if(s >= bt) fr[s]= 0;
    any|= fr[s];
   }
   if(any) update_row(&a[(size_t)r * w], Ub, fr, col & ~3, w, R);  // 掛ける数が全部 0 の行 (疎な行列で多い) は飛ばす
  }
  bt= 0;
 };
 for(int col= 0; col < n; ++col) {
  const int nr= rem.size();
  int pi= -1, nz= -1;
  u32 iv= 0;
  for(int i= 0; i < nr; ++i) {
   const int r= rem[i];
   const u32* fr= &F[(size_t)r * K];
   u64 x= a[(size_t)r * w + col];
   for(int s= 0; s < bt; ++s) x+= u64(fr[s]) * Ub[s][col];
   v[i]= Z(x).val();
   if(v[i] && pi < 0) {
    if(nz < 0) nz= i;
    if((iv= unit_inv(v[i], M))) pi= i;
   }
  }
  if(nz < 0) return Z();
  if(pi < 0) {  // 単元のない列。貯めた更新を足してから、値が 0 でない行を互除法で行 rem[nz] に集める
   flush(col);
   const int p= rem[nz];
   u64* P= &a[(size_t)p * w];
   for(int k= col; k < w; ++k) P[k]= Z(P[k]).val();
   for(int i= nz + 1; i < nr; ++i) {
    if(!v[i]) continue;
    u64* J= &a[(size_t)rem[i] * w];
    for(int k= col; k < w; ++k) J[k]= Z(J[k]).val();
    euclid<Z>(P, J, col, w);
   }
   d*= Z::raw(u32(P[col]));
   if(d == Z()) return Z();
   prow.push_back(p);
   rem[nz]= rem.back(), rem.pop_back();
   continue;
  }
  const int p= rem[pi];
  u64* up= &U[(size_t)bt * w];
  const u32* fp= &F[(size_t)p * K];
  const u64* ap= &a[(size_t)p * w];
  for(int j= col & ~3; j < w; ++j) {
   u64 x= ap[j];
   for(int s= 0; s < bt; ++s) x+= u64(fp[s]) * Ub[s][j];
   up[j]= Z(x).val();
  }
  const Z ivz= Z::raw(iv);
  for(int i= 0; i < nr; ++i) F[(size_t)rem[i] * K + bt]= (-(Z::raw(v[i]) * ivz)).val();
  d*= Z::raw(v[pi]);
  prow.push_back(p);
  Ub[bt++]= up;
  rem[pi]= rem.back(), rem.pop_back();
  if(bt == K) flush(col + 1);
 }
 return odd_perm(prow) ? -d : d;
}
template <class T> T field(Matrix<T> A) {
 const int n= A.h;
 const T zero= linalg_internal::zero<T>(), one= linalg_internal::one<T>();
 T d= one;
 for(int c= 0; c < n; ++c) {
  int p= c;
  if constexpr(std::is_floating_point_v<T>) {
   for(int r= c + 1; r < n; ++r)
    if(std::abs(A[r][c]) > std::abs(A[p][c])) p= r;
  } else
   while(p < n && A[p][c] == zero) ++p;
  if(p == n || A[p][c] == zero) return zero;
  if(p != c) std::swap_ranges(A[p], A[p] + n, A[c]), d= zero - d;
  const T* pc= A[c];
  d*= pc[c];
  const T iv= one / pc[c];
  for(int r= c + 1; r < n; ++r) {
   T* pr= A[r];
   const T f= pr[c] * iv;
   if(f == zero) continue;
   for(int j= c + 1; j < n; ++j) pr[j]-= f * pc[j];
  }
 }
 return d;
}
}
template <class T> T det(Matrix<T> A) {
 assert(A.h == A.w);
 if constexpr(linalg_internal::is_zmod_v<T>) {
  if(A.h <= 16) return det_internal::small(A);
  return det_internal::large(A);
 } else return det_internal::field(std::move(A));
}
