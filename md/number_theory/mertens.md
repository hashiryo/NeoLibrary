---
title: Mertens 関数
documentation_of: ../../neo/number_theory/mertens.hpp
---

Mertens 関数 $M(N) = \sum_{d \le N} \mu(d)$ を求めます。

```cpp
long long mertens(unsigned long long N);
```

$N < 2^{50}$ を受け、範囲の外は assert で止めます。$N = 0$ なら 0 を返します。$\Theta(N^{2/3})$ 時間と $\Theta(N^{1/2})$ の場所で求めます。

$N$ の床の値 $\lfloor N/q \rfloor$ ($q \le \sqrt{N}$) すべてで $M$ を求めます。$1.6 N^{2/3}$ ほどまでは、30 と互いに素な数だけを詰めた区画ごとに $\mu$ を篩い、和を書き留めます。それより大きい値には、$M(v) = 1 - \sum_{j \ge 2} M(\lfloor v/j \rfloor)$ を双曲線の方法で分けた再帰を使います。再帰の $\sum \mu(t) \lfloor v/t \rfloor$ は AVX2 の double で足します。$\lfloor v/t \rfloor$ は、$1/t$ の double との積を切り捨てて割り算と一致させます。手元の M2 で $M(10^{10})$ は 13 ms、$M(10^{12})$ は 0.3 s ほどでした。
