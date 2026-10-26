#pragma once
#include <cassert>
#include <numeric>
#include <vector>
#include "neo/algebra/ZMod.hpp"
#include "neo/linalg/Matrix.hpp"
#include "neo/linalg/inverse.hpp"
// 正方行列の余因子行列 adj(A) (adj(A)[i][j] = (-1)^{i+j} det(A から行 j と列 i を除いたもの)、A adj(A) = det(A) I)。T は体で、
// ZMod なら法は素数、浮動小数点数なら絶対値が 1e-8 以下の要を 0 とみなす。
// 階数が n なら adj(A) = det(A) A^{-1}、n - 2 以下なら 0 になる。n - 1 なら、A x = 0 と y^T A = 0 の 0 でない解で
// adj(A) = c x y^T と書ける。[A | I] を前進消去すると (要のない列は 1 本まで飛ばす)、要にならなかった行 r0 の右半分が
// y (y_{r0} = 1) になり、要のない列 f から後退代入で x (x_f = 1) が求まる。c = adj(A)[f][r0] は、A から行 r0 と列 f を除いた
// 行列の行列式に (-1)^{f + r0} を掛けたもので、その行列式は要の値の積に要の行の並びの符号を掛けたものになる。乱数は使わない。
// 前進消去と後退代入は neo/linalg/inverse.hpp の核 (右半分の付け替えつき) を使う。
namespace adjugate_internal {
// 置換 p (0 から n - 1 の並べ替え) が奇置換か。
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
// 階数 n - 1 のときの c の符号。要の行の並び prow (行 r0 を除いた n - 1 個) と、要のない列 f から。
inline bool minor_sign(std::vector<int> prow, int r0, int f) {
 for(int& r: prow) r-= r > r0;
 return odd_perm(prow) ^ ((f + r0) & 1);
}
inline int free_col(const std::vector<int>& pcol, int n) {
 std::vector<char> isp(n);
 for(int c: pcol) isp[c]= 1;
 int f= 0;
 while(isp[f]) ++f;
 return f;
}
template <class Z> Matrix<Z> zmod(const Matrix<Z>& A) {
 using namespace inverse_internal;
 const int n= A.h;
 Matrix<Z> res(n, n);
 Elim<Z> e(A);
 if(!e.run(1)) return res;
 Z c= Z::raw(1);
 for(u32 p: e.pval) c*= Z::raw(p);
 int wx;
 if((int)e.prow.size() == n) {  // 正則。adj(A) = det(A) A^{-1}
  if(e.odd()) c= -c;
  std::vector<int> cols(n);
  std::iota(cols.begin(), cols.end(), n);
  const auto X= tri_solve(e, cols, wx);
  const auto C= c.fixed();
  for(int t= 0; t < n; ++t) {
   Z* rt= res[t];
   const u64* xt= &X[(size_t)t * wx];
   for(int q= 0; q < n; ++q) rt[e.pos2row[q]]= Z::raw(u32(xt[q])) * C;
  }
  return res;
 }
 const int r0= e.rem[0], f= free_col(e.pcol, n);
 if(minor_sign(e.prow, r0, f)) c= -c;
 const auto X= tri_solve(e, std::vector<int>{f}, wx);
 std::vector<Z> x(n), y(n);
 x[f]= Z::raw(1);
 for(int t= 0; t < n - 1; ++t) x[e.pcol[t]]= -Z::raw(u32(X[(size_t)t * wx]));
 const u64* ar= e.row(r0) + n;
 for(int j= 0; j < n; ++j) y[j]= Z(ar[e.row2pos[j]]);
 for(int i= 0; i < n; ++i) {
  const auto xi= (c * x[i]).fixed();
  Z* ri= res[i];
  for(int j= 0; j < n; ++j) ri[j]= y[j] * xi;
 }
 return res;
}
template <class T> Matrix<T> field(const Matrix<T>& A) {
 using namespace inverse_internal;
 const int n= A.h;
 const T zero= linalg_internal::zero<T>();
 Matrix<T> res(n, n);
 Field<T> e(A);
 if(!e.run(1)) return res;
 T c= linalg_internal::one<T>();
 for(int t= 0; t < (int)e.prow.size(); ++t) c*= e.a[e.prow[t]][e.pcol[t]];
 if((int)e.prow.size() == n) {
  if(odd_perm(e.prow)) c= zero - c;
  std::vector<int> cols(n);
  std::iota(cols.begin(), cols.end(), n);
  const Matrix<T> X= e.back(cols);
  for(int t= 0; t < n; ++t)
   for(int j= 0; j < n; ++j) res[t][j]= c * X[t][j];
  return res;
 }
 const int r0= e.rem[0], f= free_col(e.pcol, n);
 if(minor_sign(e.prow, r0, f)) c= zero - c;
 const Matrix<T> X= e.back(std::vector<int>{f});
 std::vector<T> x(n, zero);
 x[f]= linalg_internal::one<T>();
 for(int t= 0; t < n - 1; ++t) x[e.pcol[t]]= zero - X[t][0];
 for(int i= 0; i < n; ++i) {
  const T xi= c * x[i];
  for(int j= 0; j < n; ++j) res[i][j]= xi * e.a[r0][n + j];
 }
 return res;
}
}
template <class T> Matrix<T> adjugate(const Matrix<T>& A) {
 assert(A.h == A.w);
 if constexpr(linalg_internal::is_zmod_v<T>) {
  if(A.h > 16) return adjugate_internal::zmod(A);
 }
 return adjugate_internal::field(A);
}
