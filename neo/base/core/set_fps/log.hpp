#pragma once
#include "neo/base/core/set_fps/generic.hpp"
#include "neo/base/core/set_fps/ranked_or64.hpp"
#include "neo/base/core/set_fps/ranked_xor.hpp"
// 集合冪級数の log(f) (exp の逆)。f[∅] = 1 が要る。f の長さは 2 冪。
//   法が奇数の ZMod (5 <= n <= 22): rank 付きの XOR 変換の領域で、上の要素から 1 つずつ割る速い核。
//   32 bit と 64 bit の整数 (法 2^64、32 bit は切り詰める): OR のランク付きのゼータ変換の領域で、段ごとに各点で割る速い核。
//   ほかの環: ランク付きのゼータ変換。どれも割り算を使わず O(n^2 2^n)。
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
 if constexpr(sps_or64::fits<T>) {
  assert(f[0] == T(1));
  return sps_or64::run<true>(f, n);
 }
 return sps_generic::log(f);
}
}  // namespace sps
