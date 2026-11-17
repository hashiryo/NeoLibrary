---
title: オンラインの部分集合畳み込み
documentation_of: ../../neo/set_power_series/relaxed.hpp
---

値が真の部分集合の値で決まる部分集合畳み込みを、$O(n^2 2^n)$ で求めます。

```cpp
template <class T, class P> std::vector<T> sps::semi_relaxed_convolve(const std::vector<T>& f, T init, const P& phi);
template <class T, class P> std::vector<T> sps::self_relaxed_convolve(const P& phi, int n);
```

<!-- textlint-disable ja-technical-writing/sentence-length -->
- `semi_relaxed_convolve` は $h[\emptyset] = \mathrm{init}$、$h[S] = \phi(S, \sum_{T \subsetneq S} h[T]\, f[S \setminus T])$ を返します。`phi` は `[](int S, T& x)` で、`x` を書き換えます。
- `self_relaxed_convolve<T>` は $h[\emptyset] = 0$、$h[S] = \phi(S, \frac{1}{2} \sum_{\emptyset \ne T \subsetneq S} h[T]\, h[S \setminus T])$ を返します。
<!-- textlint-enable ja-technical-writing/sentence-length -->

今の Library と同じ汎用の道 (ランク付きのゼータ変換) で、速い核はまだありません。
