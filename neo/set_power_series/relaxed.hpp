#pragma once
#include "neo/set_power_series/generic.hpp"
// オンライン (relaxed) の部分集合畳み込み。f の長さは 2 冪。今の Library と同じ汎用の道だけで、速い核はまだない。
//   semi_relaxed_convolve(f, init, phi): h[∅] = init、h[S] = phi(S, sum_{T ⊊ S} h[T] f[S \ T])。phi は [](int S, T& x) で x を書き換える。
//   self_relaxed_convolve<T>(phi, n): h[S] = phi(S, (1/2) sum_{∅ ≠ T ⊊ S} h[T] h[S \ T]) (h[∅] = 0)。
//   どちらも O(n^2 2^n)。
namespace sps {
using sps_generic::self_relaxed_convolve;
using sps_generic::semi_relaxed_convolve;
}  // namespace sps
