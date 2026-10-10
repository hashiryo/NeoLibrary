---
title: 素数判定
documentation_of: ../../neo/number_theory/is_prime.hpp
---

64 bit の符号なし整数が素数かどうかを判定します。

```cpp
constexpr bool is_prime(unsigned long long n);
```

`n` が素数なら true、そうでなければ false を返します。$0 \le n < 2^{64}$ のすべての `n` で正しい答えを返し、constexpr なのでコンパイル時にも使えます。

Baillie-PSW 判定で計算し、表は持ちません。37 以下の素数で割り切れない $n \ge 41^2$ について、底 2 の強擬素数判定と extra strong Lucas 判定の両方に通れば素数とします。Lucas 判定は $Q = 1$ で、$P$ は $3, 4, 5, \dots$ のうち $\left(\frac{P^2 - 4}{n}\right) = -1$ となる最初のものです。この組み合わせに $2^{64}$ 未満の反例が無いことは、Feitsma の底 2 の擬素数の一覧で確かめられています (Math::Prime::Util 0.31 の Changes)。

2 つの判定は 1 つのループで同時に回すので、時間は $n$ の桁数と同じ長さの Montgomery 乗算の鎖 1 本でほぼ決まります。$n \approx 10^{18}$ の素数で、1 回あたり 200〜300 ns ほどです (procon-judge の yosupo-primality-test、x64 の CI のランナー)。
