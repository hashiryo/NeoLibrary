---
title: 集合冪級数の log
documentation_of: ../../neo/set_power_series/log.hpp
---

集合冪級数の $\log(f)$ (exp の逆) を求めます。

```cpp
template <class T> std::vector<T> sps::log(const std::vector<T>& f);
```

`f` の長さは $2^n$ で、$f[\emptyset] = 1$ が要ります。法が奇数の `ZMod` で $5 \le n \le 22$ なら、rank 付きの XOR 変換の領域で上の要素から 1 つずつ割る速い核 ([ranked_xor.hpp](ranked_xor.md)) を使います。ほかの環はランク付きのゼータ変換です。どちらも割り算を使わず $O(n^2 2^n)$ です。浮動小数点数の型でも使えます。
