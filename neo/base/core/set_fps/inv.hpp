#pragma once
#include "neo/base/core/set_fps/generic.hpp"
#include "neo/base/core/set_fps/ranked_xor.hpp"
// 集合冪級数の逆元と割り算 (部分集合畳み込みについて)。f の長さは 2 冪。
//   inv(f): 1 / f。f[∅] = 1 が要る。
//   div(f, g): f / g。g[∅] = 1 が要る。
//   法が奇数の ZMod (4 <= n <= 22) は、rank 付きの XOR 変換の領域の各点で割って 1 回だけ逆変換する速い核。ほかの環はランク付きの
//   ゼータ変換。どちらも割り算を使わず O(n^2 2^n)。
namespace sps {
template <class T> std::vector<T> div(const std::vector<T>& f, const std::vector<T>& g) {
 const int n= sps_xor::log2_size<T>(f.size());
 assert(f.size() == g.size());
 if constexpr(sps_xor::is_zmod<T>::value) {
  if(sps_xor::fast<T>(n, 4)) {
   assert(g[0] == T(1));
   std::vector<T> h(f.size());
   const sps_xor::MK k(T::mod());
   sps_xor::dispatch<4>(n, [&]<int N>() { sps_xor::div_n<N>(sps_xor::raw(f), sps_xor::raw(g), sps_xor::raw(h), k); });
   return h;
  }
 }
 return sps_generic::div(f, g);
}
template <class T> std::vector<T> inv(const std::vector<T>& f) {
 const int n= sps_xor::log2_size<T>(f.size());
 if constexpr(sps_xor::is_zmod<T>::value) {
  if(sps_xor::fast<T>(n, 4)) {
   std::vector<T> e(f.size());
   e[0]= T(1);
   return div(e, f);
  }
 }
 return sps_generic::inv(f);
}
}  // namespace sps
