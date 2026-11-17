#pragma once
#include "neo/set_power_series/generic.hpp"
#include "neo/set_power_series/ranked_xor.hpp"
// 部分集合畳み込み h[S] = sum_{T ⊆ S} f[T] g[S \ T]。f と g は同じ 2 冪の長さ 2^n。
//   法が奇数の ZMod (4 <= n <= 22): rank 付きの XOR 変換の速い核 (ranked_xor.hpp)。O(n^2 2^n)。
//   ほかの環 (引き算を持つ型): ランク付きのゼータ変換 (割り算を使わない)。O(n^2 2^n)。
//   引き算の無い半環 (MinPlus、MaxPlus など): 3^n で素直に足す。
namespace sps {
template <class T> std::vector<T> convolve(const std::vector<T>& f, const std::vector<T>& g) {
 const int n= sps_xor::log2_size<T>(f.size());
 assert(f.size() == g.size());
 if constexpr(sps_xor::is_zmod<T>::value) {
  if(sps_xor::fast<T>(n, 4)) {
   std::vector<T> h(f.size());
   const sps_xor::MK k(T::mod());
   sps_xor::dispatch<4>(n, [&]<int N>() { sps_xor::convolve_n<N>(sps_xor::raw(f), sps_xor::raw(g), sps_xor::raw(h), k); });
   return h;
  }
 }
 if constexpr(requires(T a, T b) { a - b; }) return sps_generic::convolve(f, g);
 else return sps_generic::semiring_convolve(f, g);
}
}  // namespace sps
