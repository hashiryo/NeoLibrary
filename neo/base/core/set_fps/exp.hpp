#pragma once
#include "neo/base/core/set_fps/generic.hpp"
#include "neo/base/core/set_fps/ranked_or64.hpp"
#include "neo/base/core/set_fps/ranked_xor.hpp"
#include "neo/base/core/set_fps/semiring_simd.hpp"
// 集合冪級数の exp(f) = sum_k f^k / k! (部分集合畳み込みの冪)。f の長さは 2 冪。
//   環: f[∅] = 0 が要る。法が奇数の ZMod (5 <= n <= 22) は rank 付きの XOR 変換の領域で 1 bit ずつ伸ばす速い核、32 bit と 64 bit の
//   整数 (法 2^64、32 bit は切り詰める) は OR のランク付きのゼータ変換の領域で伸ばす速い核、ほかの環はランク付きのゼータ変換。
//   どれも割り算を使わず O(n^2 2^n)。
//   引き算の無い半環 (MinPlus など、T::one() が要る): h[∅] = 1、h[S] = sum_{S の最小の要素を含む T ⊆ S} f[T] h[S \ T]
//   (集合を塊に分ける和。MinPlus なら塊の費用の和の最小)。f[∅] は読まない。3^n。MinPlus、MaxPlus、MaxMin で値が 32 bit か
//   64 bit の符号付き整数なら、AVX2 の核 (semiring_simd.hpp) で行う。
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
  if constexpr(sps_or64::fits<T>) {
   assert(f[0] == T());
   return sps_or64::run<false>(f, n);
  }
  return sps_generic::exp(f);
 } else return sps_semi::exp(f);
}
}  // namespace sps
