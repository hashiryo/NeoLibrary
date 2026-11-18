#pragma once
#include <algorithm>
#include <vector>
#include "neo/polynomial/convolve_naive.hpp"
// [x^k] P(x)/Q(x) (Q(0) は可逆) と、定数係数の線形漸化式の第 k 項を、Bostan–Mori の方法で O(d^2 log k) で求める
// (d は Q の次数、FFT を使わない)。T は可換環 (ZMod なら法は合成数でもよい) で、割り算は最初に Q(0) の逆元を掛けるときだけ。
// 1 段で Q(x) Q(-x) = V(x^2) と P(x) Q(-x) の偶数部か奇数部を作り、k を半分にする。Q = E(x^2) + x O(x^2)、
// P = Pe(x^2) + x Po(x^2) と分けると、V(y) = E(y)^2 - y O(y)^2、偶数部 = Pe E - y Po O、奇数部 = Po E - Pe O で、
// 半分の大きさの積 4 回で 1 段進む。ZMod の積は poly_internal::conv (短い側が 16 以下ならスカラー、それより大きければ AVX2)。
// procon-judge の self-poly-kth-term で書き比べた bmo を写したもの。
namespace linrec_internal {
using u32= unsigned;
using u64= unsigned long long;
// P と Q を半分の大きさの積で 1 段ずつ進め、k = 0 になったときの P(0) を返す (Q(0) で割る前)。
template <class Z> u32 zmod(std::vector<u32> P, std::vector<u32> Q, u64 k) {
 const u32 M= Z::mod();
 auto sub= [&](u32 x, u32 y) { return x >= y ? x - y : x + M - y; };
 static thread_local std::vector<u32> E, O, Pe, Po, T1, T2;
 for(; k; k>>= 1) {
  while(!P.empty() && P.back() == 0) P.pop_back();
  if(P.empty()) return 0;
  const int lq= int(Q.size()), lp= int(P.size());
  const int ne= (lq + 1) / 2, no= lq / 2, npe= (lp + 1) / 2, npo= lp / 2;
  E.resize(ne), O.resize(no), Pe.resize(npe), Po.resize(npo);
  for(int i= 0; i < ne; ++i) E[i]= Q[2 * i];
  for(int i= 0; i < no; ++i) O[i]= Q[2 * i + 1];
  for(int i= 0; i < npe; ++i) Pe[i]= P[2 * i];
  for(int i= 0; i < npo; ++i) Po[i]= P[2 * i + 1];
  // V(y) = E^2 - y O^2 (長さ lq)
  T1.assign(2 * ne, 0), T2.assign(std::max(2 * no, 1), 0);
  poly_internal::conv<Z>(E.data(), ne, E.data(), ne, T1.data(), 2 * ne - 1);
  if(no) poly_internal::conv<Z>(O.data(), no, O.data(), no, T2.data(), 2 * no - 1);
  for(int i= 0; i < lq; ++i) Q[i]= sub(i < 2 * ne - 1 ? T1[i] : 0, i >= 1 && i - 1 < 2 * no - 1 ? T2[i - 1] : 0);
  if(k & 1) {  // Po E - Pe O
   const int l1= npo ? npo + ne - 1 : 0, l2= no ? npe + no - 1 : 0, len= std::max(l1, l2);
   T1.assign(len, 0), T2.assign(len, 0);
   if(l1) poly_internal::conv<Z>(Po.data(), npo, E.data(), ne, T1.data(), l1);
   if(l2) poly_internal::conv<Z>(Pe.data(), npe, O.data(), no, T2.data(), l2);
   P.resize(len);
   for(int i= 0; i < len; ++i) P[i]= sub(T1[i], T2[i]);
  } else {  // Pe E - y Po O
   const int l1= npe + ne - 1, l2= npo && no ? npo + no - 1 : 0, len= std::max(l1, l2 + 1);
   T1.assign(len, 0), T2.assign(len, 0);
   poly_internal::conv<Z>(Pe.data(), npe, E.data(), ne, T1.data(), l1);
   if(l2) poly_internal::conv<Z>(Po.data(), npo, O.data(), no, T2.data() + 1, l2);
   P.resize(len);
   for(int i= 0; i < len; ++i) P[i]= sub(T1[i], T2[i]);
  }
 }
 return P.empty() ? 0 : P[0];
}
template <class T> std::vector<T> mul(const std::vector<T>& a, const std::vector<T>& b) {
 if(a.empty() || b.empty()) return {};
 std::vector<T> c(a.size() + b.size() - 1, T());
 for(std::size_t i= 0; i < a.size(); ++i)
  for(std::size_t j= 0; j < b.size(); ++j) c[i + j]+= a[i] * b[j];
 return c;
}
template <class T> T generic(std::vector<T> P, std::vector<T> Q, u64 k) {
 for(; k; k>>= 1) {
  if(P.empty()) return T();
  std::vector<T> Qm= Q;
  for(std::size_t i= 1; i < Qm.size(); i+= 2) Qm[i]= -Qm[i];
  const std::vector<T> U= mul(P, Qm), V= mul(Q, Qm);
  P.clear(), Q.clear();
  for(std::size_t i= k & 1; i < U.size(); i+= 2) P.push_back(U[i]);
  for(std::size_t i= 0; i < V.size(); i+= 2) Q.push_back(V[i]);
 }
 return P.empty() ? T() : P[0];
}
}  // namespace linrec_internal
// [x^k] P(x)/Q(x)。Q(0) は可逆であること。
template <class T> T div_at(std::vector<T> P, std::vector<T> Q, unsigned long long k) {
 // Q(0) = 1 にそろえる (段を進めても Q(0) は 1 のまま)。
 if(!(Q[0] == T(1))) {
  const T iv= T(1) / Q[0];
  for(auto& x: P) x*= iv;
  for(auto& x: Q) x*= iv;
 }
 if constexpr(poly_internal::is_zmod_v<T>) return T::raw(linrec_internal::zmod<T>(poly_internal::to_u32(P), poly_internal::to_u32(Q), k));
 else return linrec_internal::generic(P, Q, k);
}
// a_n = Σ_{j=1}^{d} c_j a_{n-j} (n ≥ d) の第 k 項。c = (c_1, ..., c_d)、a は初項 a_0, ..., a_{d-1} (d 個以上あれば先頭の
// d 個を使い、k がその範囲なら a_k を返す)。
template <class T> T linear_recurrence(const std::vector<T>& c, const std::vector<T>& a, unsigned long long k) {
 if(k < a.size()) return a[k];
 const int d= int(c.size());
 if(d == 0) return T();
 if(d == 1) {
  T r= a[0], b= c[0];
  for(; k; k>>= 1, b*= b)
   if(k & 1) r*= b;
  return r;
 }
 std::vector<T> Q(d + 1);
 Q[0]= T(1);
 for(int j= 1; j <= d; ++j) Q[j]= -c[j - 1];
 std::vector<T> A(a.begin(), a.begin() + d), P= convolve_naive(A, Q);
 P.resize(d);
 if constexpr(poly_internal::is_zmod_v<T>) return T::raw(linrec_internal::zmod<T>(poly_internal::to_u32(P), poly_internal::to_u32(Q), k));
 else return linrec_internal::generic(P, Q, k);
}
