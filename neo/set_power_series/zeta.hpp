#pragma once
#include "neo/set_power_series/generic.hpp"
// 部分集合と上位集合のゼータ変換・メビウス変換、アダマール変換、OR と AND と XOR の畳み込み。長さは 2 冪で、その場で書き換える。
// 値の型 T の演算子だけを使う。ゼータは T の += だけを使うので、和を min にした MinPlus や、+= を書いた自前の型 (上位 2 個の
// max など) もそのまま使える。メビウスとアダマールは -= も使う。xor_convolve は 2^n で割る (T(1) / 2^n が 0 でなければそれを掛け、
// 0 なら T の割り算で 2^n で割る。法が偶数の ZMod では割れない)。
//   subset_zeta(f):   f[S] <- sum_{T ⊆ S} f[T]
//   subset_mobius(f): subset_zeta の逆
//   supset_zeta(f):   f[S] <- sum_{S ⊆ T} f[T]
//   supset_mobius(f): supset_zeta の逆
//   hadamard(f):      f[S] <- sum_T (-1)^{|S ∩ T|} f[T]
//   or_convolve(f, g)、and_convolve(f, g)、xor_convolve(f, g): h[S] = sum_{U op T = S} f[U] g[T]
namespace sps {
using sps_generic::and_convolve;
using sps_generic::hadamard;
using sps_generic::or_convolve;
using sps_generic::subset_mobius;
using sps_generic::subset_zeta;
using sps_generic::supset_mobius;
using sps_generic::supset_zeta;
using sps_generic::xor_convolve;
}  // namespace sps
