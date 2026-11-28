---
title: 多項式のテイラーシフト
documentation_of: ../../neo/base/core/poly/taylor_shift.hpp
---

多項式を $f(x) \mapsto f(x + c)$ とずらします。

```cpp
template <class Conv, class K>
std::vector<K> taylor_shift(const std::vector<K>& f, K c);
```

$f$ の全体 (次数 + 1 項) を返します。手間は $O(M(n))$ で、$M$ は畳み込みの手間です。`K` は体で、$|f|$ は法 (標数) より小さいことが要ります ($1/j!$ を使うため)。

核はテンプレート引数 `Conv` で受け、畳み込みの共通の口の `convolve` に渡します。省くと、include した核のうち `K` を受け持つ既定の核になります。

## 中身

$$f(x+c) = \sum_j \frac{x^j}{j!} \sum_{i \ge j} a_i\, i!\, \frac{c^{i-j}}{(i-j)!}$$

で、後ろの和は $A_i = a_{n-1-i}(n-1-i)!$ と $B_k = c^k/k!$ の 1 回の畳み込みで、$n-1-j$ 次の係数がその和になります。最後に $1/j!$ を掛けます。

$c = 0$ のときは $f$ をそのまま返します。
