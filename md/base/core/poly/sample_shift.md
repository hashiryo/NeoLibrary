---
title: 標本点のずらし
documentation_of: ../../neo/base/core/poly/sample_shift.hpp
---

連続する点での値から、別の連続する点での値を求めます。

```cpp
template <class Conv, class K>
std::vector<K> sample_points_shift(const std::vector<K>& y, K c, int m);
```

$y_i = f(i)$ ($0 \le i < n$、$f$ は $n - 1$ 次以下) から $f(c), f(c+1), \ldots, f(c+m-1)$ の $m$ 個を返します。手間は $O(M(n + m))$ で、$M$ は畳み込みの手間です。`K` は体で、`ZMod` なら法は素数です。

核はテンプレート引数 `Conv` で受け、畳み込みの共通の口の `convolve_range` に渡します。省くと、include した核のうち `K` を受け持つ既定の核になります。

## 中身

ラグランジュの形 $f(x) = P(x) \sum_i a_i / (x - i)$ ($P(x) = \prod_j (x - j)$、$a_i = y_i (-1)^{n-1-i} / (i! (n-1-i)!)$) を使います。$a$ と $1/(c-n+1+k)$ の畳み込みで $\sum_i a_i/(c+j-i)$ を 1 回にまとめ、$P(c+j)$ は $v_k = c-(n-1)+k$ の窓の積で作ります。

$c+j$ が標本点に重なると $P(c+j) = 0$ になり、割り算が壊れます。その窓は $f(c+j) = y_i$ をそのまま返します。$v_k$ の逆数と窓の積は、0 を 1 に置き換えた前積と 1 回の逆数でまとめて作り、0 の位置だけ別に扱います。
