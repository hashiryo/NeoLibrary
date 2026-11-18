#pragma once
#include <cstddef>
#include <vector>
#include "neo/set_power_series/generic.hpp"
// オンライン (relaxed) の部分集合畳み込みとゼータ変換。f の長さは 2 冪。
//   semi_relaxed_convolve(f, init, phi): h[∅] = init、h[S] = phi(S, sum_{T ⊊ S} h[T] f[S \ T])。phi は [](int S, T& x) で x を書き換える。
//   self_relaxed_convolve<T>(phi, n): h[S] = phi(S, (1/2) sum_{∅ ≠ T ⊊ S} h[T] h[S \ T]) (h[∅] = 0)。
//   どちらも今の Library と同じ汎用の道で O(n^2 2^n)。
//   online_subset_zeta<T>(n, phi): h[S] = phi(S, sum_{T ⊊ S} h[T]) を S の小さい順に求める。phi は [](int S, T& x) で、x は真の部分集合の
//   h の和で入り、h[S] に書き換える。T は T() が 0 で += があればよい。O(n 2^n)。上位集合の和なら添字を補集合にして呼ぶ。
//   区画 [lo, lo + 2^b) の再帰で、作業の配列 Y は入るときに Y[t] = sum_{j < lo, j ⊆ lo + t} h[j]、出るときに区画の中の h のゼータ変換を
//   持つ。下の半分を済ませると下の半分は下の半分のゼータ変換になるので、それを上の半分に足してから上の半分を済ませ、もう一度足すと全体の
//   ゼータ変換になる。足し算はどれも連続した配列どうし。大きさ 8 以下の区画は素朴に求める。
namespace sps_relaxed_internal {
template <class T, class P> void online_zeta_rec(T* Y, T* h, int lo, int b, const P& phi) {
 if(b <= 3) {
  const int w= 1 << b;
  for(int t= 0; t < w; ++t) {
   T x= Y[t];
   for(int u= (t - 1) & t; u != t; u= (u - 1) & t) {  // t の真の部分集合 (t = 0 なら無い)
    x+= h[lo + u];
    if(!u) break;
   }
   phi(lo + t, x);
   h[lo + t]= x;
  }
  for(int t= 0; t < w; ++t) Y[t]= h[lo + t];
  for(int k= 0; k < b; ++k)
   for(int t= 0; t < w; ++t)
    if(t >> k & 1) Y[t]+= Y[t ^ (1 << k)];
  return;
 }
 const std::size_t hw= std::size_t(1) << (b - 1);
 online_zeta_rec(Y, h, lo, b - 1, phi);
 for(std::size_t t= 0; t < hw; ++t) Y[hw + t]+= Y[t];
 online_zeta_rec(Y + hw, h, lo + int(hw), b - 1, phi);
 for(std::size_t t= 0; t < hw; ++t) Y[hw + t]+= Y[t];
}
}  // namespace sps_relaxed_internal
namespace sps {
using sps_generic::self_relaxed_convolve;
using sps_generic::semi_relaxed_convolve;
template <class T, class P> std::vector<T> online_subset_zeta(int n, const P& phi) {
 std::vector<T> h(std::size_t(1) << n), Y(std::size_t(1) << n);
 sps_relaxed_internal::online_zeta_rec(Y.data(), h.data(), 0, n, phi);
 return h;
}
}  // namespace sps
