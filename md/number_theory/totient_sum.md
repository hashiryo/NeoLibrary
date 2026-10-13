---
title: トーシェントの和
documentation_of: ../../neo/number_theory/totient_sum.hpp
---

$\Phi(N) = \sum_{i \le N} \varphi(i)$ を求めます。

```cpp
unsigned __int128 totient_sum(unsigned long long N);
```

$N < 2^{50}$ を受け、範囲の外は assert で止めます。$\Phi(N)$ は $N = 10^{10}$ で $2^{64}$ を超えるので、unsigned __int128 で丸めずに返します。法での余りが欲しいときは、返り値を割ってください。$\Theta(N^{2/3})$ 時間と $\Theta(N^{1/2})$ の場所で求めます。

$\Phi(N) = (1 + \sum_{d \le N} \mu(d) \lfloor N/d \rfloor^2) / 2$ を、mertens と同じ手順で作る Mertens 関数の床の値から求めます。$d \le \sqrt{N}$ は $\mu(d)$ を直に足し、それより大きい $d$ は $\lfloor N/d \rfloor$ ごとにまとめます。procon-judge の yosupo-sum-of-totient-function ($N \le 10^{10}$) では、x64 の CI のランナーで最大 12 ms ほどでした。今の Library の multiplicative_sum の 11 倍ほどの速さです。
