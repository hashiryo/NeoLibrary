---
title: 部分集合畳み込み
documentation_of: ../../neo/set_power_series/convolve.hpp
---

部分集合畳み込み $h[S] = \sum_{T \subseteq S} f[T]\, g[S \setminus T]$ を求めます。

```cpp
template <class T> std::vector<T> sps::convolve(const std::vector<T>& f, const std::vector<T>& g);
```

`f` と `g` は同じ長さ $2^n$ です。値の型 `T` で道が変わります。

- 法が奇数の `ZMod` で $4 \le n \le 22$ なら、rank 付きの XOR 変換の速い核 ([ranked_xor.hpp](ranked_xor.md)) を使う。$O(n^2 2^n)$。
- 引き算を持つほかの型 (法が偶数の `ZMod`、整数、mod $2^{64}$、実数など) は、ランク付きのゼータ変換を使う。割り算を使わず $O(n^2 2^n)$。
- 引き算の無い半環 (`MinPlus`、`MaxPlus` など) は、$3^n$ で素直に足す。`MinPlus` なら $h[S] = \min_{T \subseteq S} (f[T] + g[S \setminus T])$ になる。
