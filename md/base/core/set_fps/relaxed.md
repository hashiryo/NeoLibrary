---
title: オンラインの部分集合畳み込みとゼータ変換
documentation_of: ../../neo/base/core/set_fps/relaxed.hpp
---

値が真の部分集合の値で決まる部分集合畳み込みとゼータ変換を求めます。

```cpp
template <class T, class P> std::vector<T> sps::semi_relaxed_convolve(const std::vector<T>& f, T init, const P& phi);
template <class T, class P> std::vector<T> sps::self_relaxed_convolve(const P& phi, int n);
template <class T, class P> std::vector<T> sps::online_subset_zeta(int n, const P& phi);
```

<!-- textlint-disable ja-technical-writing/sentence-length -->
- `semi_relaxed_convolve` は $h[\emptyset] = \mathrm{init}$、$h[S] = \phi(S, \sum_{T \subsetneq S} h[T]\, f[S \setminus T])$ を返します。`phi` は `[](int S, T& x)` で、`x` を書き換えます。
- `self_relaxed_convolve<T>` は $h[\emptyset] = 0$、$h[S] = \phi(S, \frac{1}{2} \sum_{\emptyset \ne T \subsetneq S} h[T]\, h[S \setminus T])$ を返します。
- `online_subset_zeta<T>` は $h[S] = \phi(S, \sum_{T \subsetneq S} h[T])$ を返します。`phi` は `[](int S, T& x)` で、真の部分集合の $h$ の和で入る `x` を $h[S]$ に書き換えます。
<!-- textlint-enable ja-technical-writing/sentence-length -->

畳み込みの 2 つは、今の Library と同じ汎用の道 (ランク付きのゼータ変換) で $O(n^2 2^n)$ です。速い核はまだありません。

`online_subset_zeta` は $O(n 2^n)$ です。`phi` は $S$ の小さい順に呼ぶので、前の $S$ の値を `phi` の外に持って使えます。`T` は `T()` が 0 で `+=` があればよく、値の組 (2 つの列を同時に決めるとき) でも使えます。上位集合の和が要るときは、添字を補集合にして呼びます (`phi` の中で `S` を補集合に戻す)。添字の揃った区画 (l から l + 2^b - 1 まで) の再帰で求めます。作業の配列は、入るときに区画より前の値のうち区画の各位置の部分集合のものの和を持ち、出るときに区画の中の値のゼータ変換を持ちます。下の半分を済ませたあと、そのゼータ変換を上の半分に足してから上の半分を済ませ、もう一度足すと区画全体のゼータ変換になります。足し算はどれも連続した配列どうしです。

```cpp
// AtCoder ARC198 E: dp[k] = dp2[k - 1]、dp2[k] = A(k) sum_{j ⊆ k} dp[j] - sum_{j ⊊ k} dp2[j]
struct Pair {
  Z a, b;
  Pair& operator+=(const Pair& o) { return a += o.a, b += o.b, *this; }
};
Z prev = 1;  // dp[0]
sps::online_subset_zeta<Pair>(n, [&](int k, Pair& x) {
  const Z dp = prev, dp2 = A[k] * (x.a + dp) - x.b;
  x = Pair{dp, dp2}, prev = dp2;
});
// 答えは prev (= dp2[2^n - 1])
```
