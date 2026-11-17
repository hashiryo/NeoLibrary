#pragma once
#include "neo/set_power_series/generic.hpp"
#include "neo/set_power_series/ranked_xor.hpp"
// 集合冪級数の log(f) (exp の逆)。f[∅] = 1 が要る。f の長さは 2 冪。
//   法が奇数の ZMod (5 <= n <= 22): rank 付きの XOR 変換の領域で、上の要素から 1 つずつ割る速い核。
//   ほかの環: ランク付きのゼータ変換。どちらも割り算を使わず O(n^2 2^n)。
namespace sps {
template <class T> std::vector<T> log(const std::vector<T>& f) {
 const int n= sps_xor::log2_size<T>(f.size());
 if constexpr(sps_xor::is_zmod<T>::value) {
  if(sps_xor::fast<T>(n, 5)) {
   assert(f[0] == T(1));
   std::vector<T> h(f.size());
   const sps_xor::MK k(T::mod());
   sps_xor::dispatch<5>(n, [&]<int N>() { sps_xor::log_n<N>(sps_xor::raw(f), sps_xor::raw(h), k); });
   return h;
  }
 }
 return sps_generic::log(f);
}
}  // namespace sps
