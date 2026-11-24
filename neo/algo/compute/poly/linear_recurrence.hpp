#pragma once
#include <algorithm>
#include <vector>
#include "neo/base/core/convolve/naive.hpp"
// [x^k] P(x)/Q(x) (Q(0) は可逆) と、定数係数の線形漸化式の第 k 項を、Bostan–Mori の方法で O(d^2 log k) で求める
// (d は Q の次数、FFT を使わない)。T は可換環 (ZMod なら法は合成数でもよい) で、割り算は最初に Q(0) の逆元を掛けるときだけ。
// 1 段で Q(x) Q(-x) = V(x^2) と P(x) Q(-x) の偶数部か奇数部を作り、k を半分にする。Q = E(x^2) + x O(x^2)、
// P = Pe(x^2) + x Po(x^2) と分けると、V(y) = E(y)^2 - y O(y)^2、偶数部 = Pe E - y Po O、奇数部 = Po E - Pe O で、
// 半分の大きさの積 4 回で 1 段進む。ZMod の積は convolve_internal::conv (短い側が 16 以下ならスカラー、それより大きければ AVX2)。
// procon-judge の self-poly-kth-term で書き比べた bmo を写したもの。
namespace linrec_internal {
using u32= unsigned;
using u64= unsigned long long;
// P と Q を半分の大きさの積で 1 段ずつ進め、k = 0 になったときの P(0) を返す (Q(0) で割る前)。
template <class Z> u32 zmod_small(const u32* P0, int lp, const u32* Q0, int lq, u64 k);
template <class Z> u32 zmod(std::vector<u32> P, std::vector<u32> Q, u64 k) {
 if(Q.size() <= 16 && P.size() <= 15) return zmod_small<Z>(P.data(), int(P.size()), Q.data(), int(Q.size()), k);
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
  convolve_internal::conv<Z>(E.data(), ne, E.data(), ne, T1.data(), 2 * ne - 1);
  if(no) convolve_internal::conv<Z>(O.data(), no, O.data(), no, T2.data(), 2 * no - 1);
  for(int i= 0; i < lq; ++i) Q[i]= sub(i < 2 * ne - 1 ? T1[i] : 0, i >= 1 && i - 1 < 2 * no - 1 ? T2[i - 1] : 0);
  if(k & 1) {  // Po E - Pe O
   const int l1= npo ? npo + ne - 1 : 0, l2= no ? npe + no - 1 : 0, len= std::max(l1, l2);
   T1.assign(len, 0), T2.assign(len, 0);
   if(l1) convolve_internal::conv<Z>(Po.data(), npo, E.data(), ne, T1.data(), l1);
   if(l2) convolve_internal::conv<Z>(Pe.data(), npe, O.data(), no, T2.data(), l2);
   P.resize(len);
   for(int i= 0; i < len; ++i) P[i]= sub(T1[i], T2[i]);
  } else {  // Pe E - y Po O
   const int l1= npe + ne - 1, l2= npo && no ? npo + no - 1 : 0, len= std::max(l1, l2 + 1);
   T1.assign(len, 0), T2.assign(len, 0);
   convolve_internal::conv<Z>(Pe.data(), npe, E.data(), ne, T1.data(), l1);
   if(l2) convolve_internal::conv<Z>(Po.data(), npo, O.data(), no, T2.data() + 1, l2);
   P.resize(len);
   for(int i= 0; i < len; ++i) P[i]= sub(T1[i], T2[i]);
  }
 }
 return P.empty() ? 0 : P[0];
}
// Q の長さが 16 以下で P の長さが 15 以下のとき。スタックの配列で、P(x) Q(-x) と Q(x) Q(-x) の要る偶奇の項だけを
// 出力ごとに u64 にためる (16 項以下なので M^2 < 2^60 の積が溢れない)。次数の小さい問いをたくさん聞くとき、作業の配列と
// 核を呼ぶ手間を省く。
template <class Z> u32 zmod_small(const u32* P0, int lp, const u32* Q0, int lq, u64 k) {
 const u32 M= Z::mod();
 u32 P[32]= {}, Q[32]= {}, Qm[32]= {};
 std::copy(P0, P0 + lp, P), std::copy(Q0, Q0 + lq, Q);
 for(; k; k>>= 1) {
  while(lp && P[lp - 1] == 0) --lp;
  if(!lp) return 0;
  for(int i= 0; i < lq; ++i) Qm[i]= (i & 1) && Q[i] ? M - Q[i] : Q[i];
  u32 nP[32], nQ[32];
  const int par= int(k & 1), np= (lp + lq - 1 - par + 1) / 2;
  for(int j= 0; j < np; ++j) {
   const int t= 2 * j + par;
   u64 s= 0;
   for(int i= std::max(0, t - lq + 1); i <= std::min(lp - 1, t); ++i) s+= u64(P[i]) * Qm[t - i];
   nP[j]= Z(s).val();
  }
  for(int j= 0; j < lq; ++j) {
   const int t= 2 * j;
   u64 s= 0;
   for(int i= std::max(0, t - lq + 1); i <= std::min(lq - 1, t); ++i) s+= u64(Q[i]) * Qm[t - i];
   nQ[j]= Z(s).val();
  }
  lp= std::min(np, 31);
  std::copy(nP, nP + lp, P), std::copy(nQ, nQ + lq, Q);
 }
 return lp ? P[0] : 0;
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
  for(std::size_t i= 1; i < Qm.size(); i+= 2) Qm[i]= T() - Qm[i];  // 単項の - を持たない型 (GF2p64 など) のため
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
 if constexpr(convolve_internal::is_zmod_v<T>) return T::raw(linrec_internal::zmod<T>(convolve_internal::to_u32(P), convolve_internal::to_u32(Q), k));
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
 if constexpr(convolve_internal::is_zmod_v<T>) {
  if(d <= 15) {  // 配列をヒープに取らない: Q = 1 - Σ c_j x^j、P = (A Q) mod x^d
   using linrec_internal::u32, linrec_internal::u64;
   const u32 M= T::mod();
   u32 Q[16], P[16];
   Q[0]= T(1).val();
   for(int j= 1; j <= d; ++j) Q[j]= c[j - 1].val() ? M - c[j - 1].val() : 0;
   for(int i= 0; i < d; ++i) {
    u64 s= 0;
    for(int j= 0; j <= i; ++j) s+= u64(a[i - j].val()) * Q[j];
    P[i]= T(s).val();
   }
   return T::raw(linrec_internal::zmod_small<T>(P, d, Q, d + 1, k));
  }
 }
 std::vector<T> Q(d + 1);
 Q[0]= T(1);
 for(int j= 1; j <= d; ++j) Q[j]= T() - c[j - 1];
 std::vector<T> A(a.begin(), a.begin() + d), P= convolve_naive(A, Q);
 P.resize(d);
 if constexpr(convolve_internal::is_zmod_v<T>) return T::raw(linrec_internal::zmod<T>(convolve_internal::to_u32(P), convolve_internal::to_u32(Q), k));
 else return linrec_internal::generic(P, Q, k);
}
