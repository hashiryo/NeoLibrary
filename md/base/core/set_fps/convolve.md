---
title: 部分集合畳み込み
documentation_of: ../../neo/base/core/set_fps/convolve.hpp
---

部分集合畳み込み $h[S] = \sum_{T \subseteq S} f[T]\, g[S \setminus T]$ を求めます。

```cpp
template <class M> std::vector<typename M::T> sps::convolve(const std::vector<typename M::T>& f, const std::vector<typename M::T>& g);
template <class T> std::vector<T> sps::convolve(const std::vector<T>& f, const std::vector<T>& g);
```

構造 `M` をテンプレート引数で受けます。値の型は `M::T` で、`f` と `g` は同じ長さ $2^n$ です。`sps::convolve<MinPlus<int>>(f, g)` の `f` と `g` は `std::vector<int>` です。構造を省くと、`PlusTimes<T>` (普通の + と ×) として扱います。

道は構造で変わります。

- 値の型が `ZMod` で法が奇数、$4 \le n \le 22$ なら、rank 付きの XOR 変換の速い核 ([ranked_xor.hpp](ranked_xor.md)) を使う。$O(n^2 2^n)$。
- 構造が引き算 `sub` を持つ (環) なら、ランク付きのゼータ変換を使う。割り算を使わず $O(n^2 2^n)$。算法は値の型の `+`、`-`、`*` で書いてある。
- 引き算の無い半環 (`MinPlus`、`MaxPlus` など) は、$3^n$ で足す。`MinPlus` なら $h[S] = \min_{T \subseteq S} (f[T] + g[S \setminus T])$ になる。`MinPlus`、`MaxPlus`、`MaxMin` で値が 32 bit か 64 bit の符号付き整数なら、AVX2 の核 ([semiring_simd.hpp](semiring_simd.md)) で行う。足し算と掛け算は構造の `add` と `mul` で行う。
