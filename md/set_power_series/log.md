---
title: 集合冪級数の log
documentation_of: ../../neo/set_power_series/log.hpp
---

集合冪級数の $\log(f)$ (exp の逆) を求めます。

```cpp
template <class T> std::vector<T> sps::log(const std::vector<T>& f);
```

`f` の長さは $2^n$ で、$f[\emptyset] = 1$ が要ります。法が奇数の `ZMod` で $5 \le n \le 22$ なら、rank 付きの XOR 変換の領域で上の要素から 1 つずつ割る速い核 ([ranked_xor.hpp](ranked_xor.md)) を使います。32 bit と 64 bit の整数の型 (符号付きも) は、法 $2^{64}$ の値として OR のランク付きのゼータ変換の領域で段ごとに各点で割る速い核 ([ranked_or64.hpp](ranked_or64.md)) を使います。32 bit の型は 64 bit で求めてから切り詰めます。ほかの環はランク付きのゼータ変換です。どれも割り算を使わず $O(n^2 2^n)$ です。浮動小数点数の型でも使えます。
