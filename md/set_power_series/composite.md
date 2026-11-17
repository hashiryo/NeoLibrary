---
title: 集合冪級数への合成
documentation_of: ../../neo/set_power_series/composite.hpp
---

多項式や EGF を集合冪級数に合成します。

```cpp
template <class T> std::vector<T> sps::egf_comp(const std::vector<T>& F, const std::vector<T>& f);  // sum_i F_i f^i / i!
template <class T> std::vector<T> sps::poly_comp(const std::vector<T>& P, const std::vector<T>& f); // sum_i P_i f^i
template <class T> std::vector<T> sps::pow(std::vector<T> f, std::uint64_t k);                       // f^k
```

`f` の長さは $2^n$ です。`egf_comp` は $f[\emptyset] = 0$ が要り、`F` の $n$ より先は使いません。`poly_comp` と `pow` の $f[\emptyset]$ は何でもよく、$c_0 = f[\emptyset]$ のまわりの微分 $P^{(k)}(c_0)$ を EGF の係数にして `egf_comp` します。

法が奇数の `ZMod` で $5 \le n \le 22$ なら、rank 付きの XOR 変換の領域で 1 bit ずつ伸ばす速い核 ([ranked_xor.hpp](ranked_xor.md)) を使います。ほかの環はランク付きのゼータ変換です。どちらも割り算を使わず $O(n^2 2^n)$ で、`poly_comp` は `P` の長さ $M$ に $O(nM)$ の微分と評価が足されます (速い核では 8 個ずつ)。
