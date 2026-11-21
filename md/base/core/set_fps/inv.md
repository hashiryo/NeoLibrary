---
title: 集合冪級数の逆元と割り算
documentation_of: ../../neo/base/core/set_fps/inv.hpp
---

部分集合畳み込みについての逆元と割り算を求めます。

```cpp
template <class T> std::vector<T> sps::inv(const std::vector<T>& f);                            // 1 / f
template <class T> std::vector<T> sps::div(const std::vector<T>& f, const std::vector<T>& g);   // f / g
```

長さは $2^n$ で、割る側の空集合の値は 1 である必要があります。法が奇数の `ZMod` で $4 \le n \le 22$ なら、rank 付きの XOR 変換の領域の各点で割って 1 回だけ逆変換する速い核 ([ranked_xor.hpp](ranked_xor.md)) を使います。ほかの環はランク付きのゼータ変換です。どちらも割り算を使わず $O(n^2 2^n)$ です。
