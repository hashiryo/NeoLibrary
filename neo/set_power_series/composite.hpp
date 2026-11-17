#pragma once
#include "neo/set_power_series/generic.hpp"
#include "neo/set_power_series/ranked_xor.hpp"
// 集合冪級数への合成。f の長さは 2 冪 (2^n)。
//   egf_comp(F, f) = sum_i F_i f^i / i!。f[∅] = 0 が要る。F の n より先は使わない。
//   poly_comp(P, f) = sum_i P_i f^i (f[∅] は何でもよい)。f[∅] のまわりの P の微分 P^{(k)}(f[∅]) を EGF の係数にして egf_comp する。
//   pow(f, k) = f^k (k は 64 bit)。
//   法が奇数の ZMod (5 <= n <= 22) は rank 付きの XOR 変換の領域で 1 bit ずつ伸ばす速い核、ほかの環はランク付きのゼータ変換。
//   どちらも割り算を使わず O(n^2 2^n) (poly_comp は P の長さ M に O(n M) が足される)。
namespace sps {
template <class T> std::vector<T> egf_comp(const std::vector<T>& F, const std::vector<T>& f) {
 const int n= sps_xor::log2_size<T>(f.size());
 if constexpr(sps_xor::is_zmod<T>::value) {
  if(sps_xor::fast<T>(n, 5)) {
   assert(f[0] == T());
   std::vector<sps_xor::u32> e(n + 1, 0);
   for(int i= 0; i <= n && i < int(F.size()); ++i) e[i]= F[i].val();
   std::vector<T> h(f.size());
   const sps_xor::MK k(T::mod());
   sps_xor::dispatch<5>(n, [&]<int N>() { sps_xor::egf_comp_n<N>(sps_xor::raw(f), e.data(), sps_xor::raw(h), k); });
   return h;
  }
 }
 std::vector<T> G(F.begin(), F.begin() + std::min<std::size_t>(F.size(), n + 1));
 if(G.empty()) G.push_back(T());
 return sps_generic::egf_comp(G, f);
}
template <class T> std::vector<T> poly_comp(const std::vector<T>& P, const std::vector<T>& f) {
 const int n= sps_xor::log2_size<T>(f.size());
 if constexpr(sps_xor::is_zmod<T>::value) {
  if(sps_xor::fast<T>(n, 5)) {
   std::vector<T> h(f.size());
   if(P.empty()) return h;
   const sps_xor::MK k(T::mod());
   const std::vector<sps_xor::u32> e= sps_xor::taylor(sps_xor::raw(P), int(P.size()), f[0].val(), n, k);
   sps_xor::dispatch<5>(n, [&]<int N>() { sps_xor::egf_comp_n<N>(sps_xor::raw(f), e.data(), sps_xor::raw(h), k); });  // f[∅] は読まない
   return h;
  }
 }
 if(P.empty()) return std::vector<T>(f.size());
 return sps_generic::poly_comp(P, f);
}
template <class T> std::vector<T> pow(std::vector<T> f, std::uint64_t k) {
 const unsigned N= f.size(), n= sps_xor::log2_size<T>(N);
 std::vector<T> F;
 unsigned i, m;
 if(n < k) {
  F.resize(n + 1), F[m= n]= T(1);
  T x= f[0];
  for(std::uint64_t l= k - n;; x*= x)
   if(l & 1 ? F[n]*= x : T(), !(l>>= 1)) break;
 } else F.resize(k + 1), F[m= unsigned(k)]= T(1);
 for(i= m; i--;) F[i]= F[i + 1] * f[0];
 for(T t= T(1); ++i < m;) F[i + 1]*= (t*= T(k - i));
 f[0]= T();
 return egf_comp(F, f);
}
}  // namespace sps
