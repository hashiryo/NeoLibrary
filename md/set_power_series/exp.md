---
title: 集合冪級数の exp
documentation_of: ../../neo/set_power_series/exp.hpp
---

集合冪級数の $\exp(f) = \sum_k f^k / k!$ (部分集合畳み込みの冪) を求めます。

```cpp
template <class T> std::vector<T> sps::exp(const std::vector<T>& f);
```

`f` の長さは $2^n$ です。

- 環では $f[\emptyset] = 0$ が要ります。法が奇数の `ZMod` で $5 \le n \le 22$ なら、rank 付きの XOR 変換の領域で 1 bit ずつ伸ばす速い核 ([ranked_xor.hpp](ranked_xor.md)) を使います。ほかの環はランク付きのゼータ変換です。どちらも割り算を使わず $O(n^2 2^n)$ です。
- 引き算の無い半環 (`MinPlus` など、`T::one()` が要る) では、$h[\emptyset] = 1$、$h[S] = \sum f[T]\, h[S \setminus T]$ を $3^n$ で求めます。$T$ は、$S$ の最小の要素を含む $S$ の部分集合です。集合を塊に分ける和で、`MinPlus` なら塊の費用の和の最小になります。$f[\emptyset]$ は読みません。`MinPlus`、`MaxPlus`、`MaxMin` で値が 32 bit か 64 bit の符号付き整数なら、AVX2 の核 ([semiring_simd.hpp](semiring_simd.md)) で行います。
