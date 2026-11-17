#pragma once
#include "neo/set_power_series/generic.hpp"
#include "neo/set_power_series/ranked_xor.hpp"
// 集合冪級数の exp(f) = sum_k f^k / k! (部分集合畳み込みの冪)。f の長さは 2 冪。
//   環: f[∅] = 0 が要る。法が奇数の ZMod (5 <= n <= 22) は rank 付きの XOR 変換の領域で 1 bit ずつ伸ばす速い核、ほかの環は
//   ランク付きのゼータ変換。どちらも割り算を使わず O(n^2 2^n)。
//   引き算の無い半環 (MinPlus など、T::one() が要る): h[∅] = 1、h[S] = sum_{S の最小の要素を含む T ⊆ S} f[T] h[S \ T]
//   (集合を塊に分ける和。MinPlus なら塊の費用の和の最小)。f[∅] は読まない。3^n。
namespace sps {
template <class T> std::vector<T> exp(const std::vector<T>& f) {
 const int n= sps_xor::log2_size<T>(f.size());
 if constexpr(requires(T a, T b) { a - b; }) {
  if constexpr(sps_xor::is_zmod<T>::value) {
   if(sps_xor::fast<T>(n, 5)) {
    assert(f[0] == T());
    std::vector<T> h(f.size());
    const sps_xor::MK k(T::mod());
    sps_xor::dispatch<5>(n, [&]<int N>() { sps_xor::exp_n<N>(sps_xor::raw(f), sps_xor::raw(h), k); });
    return h;
   }
  }
  return sps_generic::exp(f);
 } else return sps_generic::semiring_exp(f, T::one());
}
}  // namespace sps
