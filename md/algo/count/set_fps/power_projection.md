---
title: power projection と合成の転置
documentation_of: ../../neo/algo/count/set_fps/power_projection.hpp
---

集合冪級数の冪の重み付きの和と、合成の転置を求めます。

```cpp
template <class T> std::vector<T> sps::power_projection(const std::vector<T>& f, const std::vector<T>& w, int M);
template <class T> std::vector<T> sps::egf_T(const std::vector<T>& f);
template <class T> std::vector<T> sps::egf_T(const std::vector<T>& f, const std::vector<T>& g);
```

- `power_projection` は $\mathrm{ans}_m = \sum_S w[S]\, [S] f^m$ ($m = 0, \ldots, M-1$) を返します。$f[\emptyset]$ は何でもかまいません。
- `egf_T(f)` は全体集合の係数 $[X^{\text{全体}}] f^k / k!$ ($k = 0, \ldots, n$) を、`egf_T(f, g)` は $[X^{\text{全体}}] (f^k / k!)\, g$ を返します。どちらも $f[\emptyset]$ を 0 とみなします。

次の値を使います。

- $c_0 = f[\emptyset]$、$g = f - c_0$、$V_k = \sum_S w[S]\, (g^k / k!)[S]$ とおく。
- $\mathrm{ans}_m = \sum_k (m)_k\, c_0^{m-k} V_k$ ($(m)_k$ は下降階乗)。
- $V$ は、合成 ([composite.hpp](composite.md)) の転置を $w$ に当てて求める。
- 答えは $u_k(m+1) = c_0 u_k(m) + k\, u_{k-1}(m)$ で $m$ を 1 つずつ進めるので、$c_0$ に逆元が無くてもよい。

法が奇数の `ZMod` で $5 \le n \le 22$ なら速い核 ([ranked_xor.hpp](ranked_xor.md)) を使い、ほかの環はランク付きのゼータ変換です。時間は $V$ が $O(n^2 2^n)$、答えが $O(nM)$ です。
