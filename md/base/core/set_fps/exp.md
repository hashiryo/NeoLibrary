---
title: 集合冪級数の exp
documentation_of: ../../neo/base/core/set_fps/exp.hpp
---

集合冪級数の $\exp(f) = \sum_k f^k / k!$ (部分集合畳み込みの冪) を求めます。

```cpp
template <class M> std::vector<typename M::T> sps::exp(const std::vector<typename M::T>& f);
template <class T> std::vector<T> sps::exp(const std::vector<T>& f);
```

構造 `M` をテンプレート引数で受けます。値の型は `M::T` で、`f` の長さは $2^n$ です。構造を省くと、`PlusTimes<T>` (普通の + と ×) として扱います。

- 環 (構造が引き算 `sub` を持つ) では $f[\emptyset] = 0$ が要ります。法が奇数の `ZMod` で $5 \le n \le 22$ なら、rank 付きの XOR 変換の領域で 1 bit ずつ伸ばす速い核 ([ranked_xor.hpp](ranked_xor.md)) を使います。32 bit と 64 bit の整数の型 (符号付きも) は、法 $2^{64}$ の値として OR のランク付きのゼータ変換の領域で伸ばす速い核 ([ranked_or64.hpp](ranked_or64.md)) を使います。32 bit の型は 64 bit で求めてから切り詰めます。ほかの環はランク付きのゼータ変換です。どれも割り算を使わず $O(n^2 2^n)$ です。算法は値の型の `+`、`-`、`*` で書いてあります。
- 引き算の無い半環 (`MinPlus` など、`one` が要る) では、$h[\emptyset] = 1$、$h[S] = \sum f[T]\, h[S \setminus T]$ を $3^n$ で求めます。$T$ は、$S$ の最小の要素を含む $S$ の部分集合です。集合を塊に分ける和で、`MinPlus` なら塊の費用の和の最小になります。$f[\emptyset]$ は読みません。`MinPlus`、`MaxPlus`、`MaxMin` で値が 32 bit か 64 bit の符号付き整数なら、AVX2 の核 ([semiring_simd.hpp](semiring_simd.md)) で行います。足し算と掛け算は構造の `add` と `mul` で行います。
