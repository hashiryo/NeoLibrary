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
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>
#include "neo/algebra/ZMod.hpp"
#include "neo/linalg/Matrix.hpp"
// 連立一次方程式 A x = b と、A x = 0 の解の空間 (核)。T は体で、ZMod なら法は素数、浮動小数点数なら絶対値が 1e-8 以下の値を
// 0 とみなす。解は、自由な変数 (要のない列) を 0 にした 1 つ。核の基底は、自由な変数 1 つを 1、ほかの自由な変数を 0 にした
// ものを行に並べた行列で返す。
// ZMod は、[A | b] を前進消去し、要の列が作る上三角の行列 T で、自由な列と b の列をまとめて後退代入する (T X = Y)。
// 前進消去は、行を u64 で持ち、要を 8 本ずつ貯めて残りの各行を 1 回だけ読み書きする (遅らせた更新)。積は M^2 < 2^60 なので、
// 2^62 + 2^32 未満の値に 8 本足しても u64 に収まり、足したあとは 2^32 の位を 2^32 mod M に掛けて畳む。後退代入は下から 8 行ずつの
// 組にして、組より下の行の分を 8 行まとめて足す。max(h, w) <= 16 の ZMod とほかの体は、素朴な掃き出し (浮動小数点数は絶対値が
// 最大の要を選ぶ) で簡約した階段の形にする。
// procon-judge の yosupo-system-of-linear-equations で書き比べた zmod_sle_blk8t を写したもの。
namespace solve_internal {
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
// acc の 8 行 (行の間隔 wx) に X の ns 行 (行の間隔 wx) を足す: acc[q][j] += sum_s c[s * 8 + q] X[s][j]。4 列ずつ 8 行分の和を
// レジスタに持ち、X の各行の 4 列を 1 回読んで 8 行に足す。8 本足すごとに畳む。c[...] < M、X[s][j] < M。
inline void block_update(u64* __restrict acc, const u64* __restrict X, const u32* c, int ns, int wx, u64 R) {
 const __m256i Rv= _mm256_set1_epi64x(R), LO= _mm256_set1_epi64x(0xffffffffll);
 auto fold= [&](__m256i x) { return _mm256_add_epi64(_mm256_mul_epu32(_mm256_srli_epi64(x, 32), Rv), _mm256_and_si256(x, LO)); };
 for(int j= 0; j < wx; j+= 4) {
  u64* ap= acc + j;
  __m256i a0= _mm256_loadu_si256((const __m256i*)ap), a1= _mm256_loadu_si256((const __m256i*)(ap + wx));
  __m256i a2= _mm256_loadu_si256((const __m256i*)(ap + 2 * wx)), a3= _mm256_loadu_si256((const __m256i*)(ap + 3 * wx));
  __m256i a4= _mm256_loadu_si256((const __m256i*)(ap + 4 * wx)), a5= _mm256_loadu_si256((const __m256i*)(ap + 5 * wx));
  __m256i a6= _mm256_loadu_si256((const __m256i*)(ap + 6 * wx)), a7= _mm256_loadu_si256((const __m256i*)(ap + 7 * wx));
  const u64* xp= X + j;
  const u32* cp= c;
  for(int s= 0; s < ns;) {
   for(const int e= std::min(ns, s + 8); s < e; ++s, xp+= wx, cp+= 8) {
    const __m256i x= _mm256_loadu_si256((const __m256i*)xp);
    a0= _mm256_add_epi64(a0, _mm256_mul_epu32(x, _mm256_set1_epi32(cp[0])));
    a1= _mm256_add_epi64(a1, _mm256_mul_epu32(x, _mm256_set1_epi32(cp[1])));
    a2= _mm256_add_epi64(a2, _mm256_mul_epu32(x, _mm256_set1_epi32(cp[2])));
    a3= _mm256_add_epi64(a3, _mm256_mul_epu32(x, _mm256_set1_epi32(cp[3])));
    a4= _mm256_add_epi64(a4, _mm256_mul_epu32(x, _mm256_set1_epi32(cp[4])));
    a5= _mm256_add_epi64(a5, _mm256_mul_epu32(x, _mm256_set1_epi32(cp[5])));
    a6= _mm256_add_epi64(a6, _mm256_mul_epu32(x, _mm256_set1_epi32(cp[6])));
    a7= _mm256_add_epi64(a7, _mm256_mul_epu32(x, _mm256_set1_epi32(cp[7])));
   }
   a0= fold(a0), a1= fold(a1), a2= fold(a2), a3= fold(a3), a4= fold(a4), a5= fold(a5), a6= fold(a6), a7= fold(a7);
  }
  _mm256_storeu_si256((__m256i*)ap, a0), _mm256_storeu_si256((__m256i*)(ap + wx), a1);
  _mm256_storeu_si256((__m256i*)(ap + 2 * wx), a2), _mm256_storeu_si256((__m256i*)(ap + 3 * wx), a3);
  _mm256_storeu_si256((__m256i*)(ap + 4 * wx), a4), _mm256_storeu_si256((__m256i*)(ap + 5 * wx), a5);
  _mm256_storeu_si256((__m256i*)(ap + 6 * wx), a6), _mm256_storeu_si256((__m256i*)(ap + 7 * wx), a7);
 }
}
// n 行 w 列 (w は 4 の倍数) の行列の前進消去。要は列の順に探し、その列で値が 0 でない最初の行を選ぶ (まだ要にしていない行の
// 並び rem は、要を抜くたびに末尾の行で埋める)。
template <class Z> struct Elim {
 int n, w;
 std::vector<u64> a;                // n × w。2^62 + 2^32 未満で、その成分と mod M で等しい
 std::vector<int> rem, prow, pcol;  // まだ要にしていない行、要にした行と列 (要にした順)
 std::vector<u32> pval;             // 要の値
 std::vector<u64> U;                // 要の行を M 未満に直したもの (要にした順に w 個ずつ)。要の列を 4 の倍数に切り下げた位置より左は使わない
 Elim(int n, int w): n(n), w(w), a((size_t)n * w), rem(n) { std::iota(rem.begin(), rem.end(), 0); }
 u64* row(int r) { return a.data() + (size_t)r * w; }
 const u64* urow(int t) const { return U.data() + (size_t)t * w; }
 u32 value(int r, int j) const { return Z(a[(size_t)r * w + j]).val(); }
 // 列 [0, lim) で要を探しながら前進消去する。要のない列は飛ばす。
 void run(int lim) {
  const u32 M= Z::mod();
  const u64 R= (u64(1) << 32) % M;
  U.assign((size_t)std::min(n, lim) * w, 0);
  std::vector<u32> F((size_t)n * K), v(n);
  std::vector<u64> zero(w);
  const u64* Ub[K];
  int col= 0;
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
    u64* up= &U[prow.size() * w];
    const u32* fp= &F[(size_t)p * K];
    const u64* ap= row(p);
    for(int j= col & ~3; j < w; ++j) {
     u64 x= ap[j];
     for(int s= 0; s < bt; ++s) x+= u64(fp[s]) * Ub[s][j];
     up[j]= Z(x).val();
    }
    const Z iv= Z::raw(v[pi]).inv();
    for(int i= 0; i < nr; ++i) F[(size_t)rem[i] * K + bt]= (-(Z::raw(v[i]) * iv)).val();
    prow.push_back(p), pcol.push_back(col), pval.push_back(v[pi]);
    Ub[bt++]= up;
    rem[pi]= rem.back(), rem.pop_back();
   }
   if(!bt) break;
   for(int s= bt; s < K; ++s) Ub[s]= zero.data();
   const int jb= col & ~3;
   for(int r: rem) {
    u32* fr= &F[(size_t)r * K];
    for(int s= bt; s < K; ++s) fr[s]= 0;
    u32 any= 0;
    for(int s= 0; s < K; ++s) any|= fr[s];
    if(any) update_row(row(r), Ub, fr, jb, w, R);  // 掛ける数が全部 0 の行 (疎な行列で多い) は飛ばす
   }
  }
 }
};
// 前進消去のあと、要の列が作る上三角の行列 T (T[t][s] = U_t[pcol[s]]) で T X = Y を解く。Y の行 t は U_t の列 cols。
// X は要の数 × wx (wx は cols の数を 4 の倍数に切り上げたもの) で、M 未満の値を u64 で持つ。下から 8 行ずつの組にして、
// 組より下の行の分は block_update で 8 行まとめて足し、組の中は下の行から順に update_row で足す。
template <class Z> std::vector<u64> tri_solve(const Elim<Z>& e, const std::vector<int>& cols, int& wx) {
 const u32 M= Z::mod();
 const u64 R= (u64(1) << 32) % M;
 const int r= e.prow.size(), nc= cols.size();
 wx= (nc + 3) & ~3;
 std::vector<u64> X((size_t)r * wx), acc((size_t)K * wx), zero(wx);
 std::vector<u32> c((size_t)r * K);
 for(int t1= r; t1 > 0; t1-= K) {
  const int t0= std::max(0, t1 - K), bs= t1 - t0, ns= r - t1;
  for(int q= 0; q < K; ++q) {
   u64* aq= &acc[(size_t)q * wx];
   std::fill(aq, aq + wx, 0);
   if(q < bs) {
    const u64* ut= e.urow(t0 + q);
    for(int k= 0; k < nc; ++k) aq[k]= ut[cols[k]];
   }
  }
  for(int s= t1; s < r; ++s)
   for(int q= 0; q < K; ++q) {
    const u32 x= q < bs ? u32(e.urow(t0 + q)[e.pcol[s]]) : 0;
    c[(size_t)(s - t1) * K + q]= x ? M - x : 0;
   }
  if(ns) block_update(acc.data(), &X[(size_t)t1 * wx], c.data(), ns, wx, R);
  for(int q= bs - 1; q >= 0; --q) {
   const int t= t0 + q;
   u64* aq= &acc[(size_t)q * wx];
   const u64* Ub[K];
   u32 f[K];
   for(int i= 0; i < K; ++i) {
    const int s= t + 1 + i;
    if(s < t1) {
     const u32 x= u32(e.urow(t)[e.pcol[s]]);
     Ub[i]= &X[(size_t)s * wx], f[i]= x ? M - x : 0;
    } else Ub[i]= zero.data(), f[i]= 0;
   }
   update_row(aq, Ub, f, 0, wx, R);
   const Z iv= Z::raw(e.pval[t]).inv();
   u64* xt= &X[(size_t)t * wx];
   for(int k= 0; k < wx; ++k) xt[k]= (Z(aq[k]) * iv).val();
  }
 }
 return X;
}
// A x = b (b が nullptr なら核だけ) を解く。解がなければ false。sol と ker が nullptr でなければ、特解と核の基底を書く。
template <class Z> bool zmod(const Matrix<Z>& A, const std::vector<Z>* b, std::vector<Z>* sol, Matrix<Z>* ker) {
 const int n= A.h, m= A.w, hb= b != nullptr;
 Elim<Z> e(n, (m + hb + 3) & ~3);
 for(int i= 0; i < n; ++i) {
  u64* row= e.row(i);
  for(int j= 0; j < m; ++j) row[j]= A[i][j].val();
  if(hb) row[m]= (*b)[i].val();
 }
 e.run(m);
 if(hb)  // 残りの行は A の部分が 0 なので、b の値が 0 でなければ解がない
  for(int i: e.rem)
   if(e.value(i, m)) return false;
 const int r= e.prow.size();
 std::vector<int> cols;
 if(ker) {
  std::vector<char> isp(m);
  for(int c: e.pcol) isp[c]= 1;
  for(int j= 0; j < m; ++j)
   if(!isp[j]) cols.push_back(j);
 }
 const int nf= cols.size();
 if(sol) cols.push_back(m);
 if(ker) *ker= Matrix<Z>(nf, m);
 if(sol) sol->assign(m, Z());
 if(cols.empty()) return true;
 int wx;
 const auto X= tri_solve(e, cols, wx);
 if(sol)
  for(int t= 0; t < r; ++t) (*sol)[e.pcol[t]]= Z::raw(u32(X[(size_t)t * wx + nf]));
 for(int k= 0; k < nf; ++k) {
  Z* v= (*ker)[k];
  v[cols[k]]= Z::raw(1);
  for(int t= 0; t < r; ++t) {
   const u32 x= u32(X[(size_t)t * wx + k]);
   v[e.pcol[t]]= Z::raw(x ? Z::mod() - x : 0);
  }
 }
 return true;
}
template <class T> bool is_zero(const T& x) {
 if constexpr(std::is_floating_point_v<T>) return std::abs(x) <= T(1e-8);
 else return x == linalg_internal::zero<T>();
}
// 体の上の素朴な形。[A | b] を簡約した階段の形にしてから読む。
template <class T> bool field(const Matrix<T>& A, const std::vector<T>* b, std::vector<T>* sol, Matrix<T>* ker) {
 const int n= A.h, m= A.w, W= m + (b != nullptr);
 const T zero= linalg_internal::zero<T>(), one= linalg_internal::one<T>();
 Matrix<T> a(n, W);
 for(int i= 0; i < n; ++i) {
  std::copy(A[i], A[i] + m, a[i]);
  if(b) a[i][m]= (*b)[i];
 }
 std::vector<int> pcol;
 int r= 0;
 for(int c= 0; c < m && r < n; ++c) {
  int p= r;
  if constexpr(std::is_floating_point_v<T>) {
   for(int i= r + 1; i < n; ++i)
    if(std::abs(a[i][c]) > std::abs(a[p][c])) p= i;
  } else
   while(p < n - 1 && a[p][c] == zero) ++p;
  if(is_zero(a[p][c])) continue;
  if(p != r) std::swap_ranges(a[p] + c, a[p] + W, a[r] + c);
  T* pr= a[r];
  const T iv= one / pr[c];
  for(int j= c; j < W; ++j) pr[j]*= iv;
  for(int i= 0; i < n; ++i) {
   if(i == r) continue;
   T* ai= a[i];
   const T f= ai[c];
   if(f == zero) continue;
   for(int j= c; j < W; ++j) ai[j]-= f * pr[j];
  }
  pcol.push_back(c), ++r;
 }
 if(b)
  for(int i= r; i < n; ++i)
   if(!is_zero(a[i][m])) return false;
 if(sol) {
  sol->assign(m, zero);
  for(int t= 0; t < r; ++t) (*sol)[pcol[t]]= a[t][m];
 }
 if(ker) {
  std::vector<char> isp(m);
  for(int c: pcol) isp[c]= 1;
  std::vector<int> fc;
  for(int j= 0; j < m; ++j)
   if(!isp[j]) fc.push_back(j);
  *ker= Matrix<T>(fc.size(), m);
  for(int k= 0; k < (int)fc.size(); ++k) {
   T* v= (*ker)[k];
   v[fc[k]]= one;
   for(int t= 0; t < r; ++t) v[pcol[t]]= zero - a[t][fc[k]];
  }
 }
 return true;
}
template <class T> bool dispatch(const Matrix<T>& A, const std::vector<T>* b, std::vector<T>* sol, Matrix<T>* ker) {
 if constexpr(linalg_internal::is_zmod_v<T>) {
  if(std::max(A.h, A.w) > 16) return zmod(A, b, sol, ker);
 }
 return field(A, b, sol, ker);
}
}
// A x = b の解を 1 つ返す。解がなければ std::nullopt。
template <class T> std::optional<std::vector<T>> solve(const Matrix<T>& A, const std::vector<T>& b) {
 assert((int)b.size() == A.h);
 std::vector<T> x;
 if(!solve_internal::dispatch(A, &b, &x, (Matrix<T>*)nullptr)) return std::nullopt;
 return x;
}
// A x = b の解 1 つと、核の基底 (行に並べた行列) の組を返す。解がなければ std::nullopt。
template <class T> std::optional<std::pair<std::vector<T>, Matrix<T>>> solve_with_kernel(const Matrix<T>& A, const std::vector<T>& b) {
 assert((int)b.size() == A.h);
 std::pair<std::vector<T>, Matrix<T>> r;
 if(!solve_internal::dispatch(A, &b, &r.first, &r.second)) return std::nullopt;
 return r;
}
// A x = 0 の解の空間の基底を行に並べた行列 ((w - rank) × w)。
template <class T> Matrix<T> kernel(const Matrix<T>& A) {
 Matrix<T> k;
 solve_internal::dispatch(A, (const std::vector<T>*)nullptr, (std::vector<T>*)nullptr, &k);
 return k;
}
