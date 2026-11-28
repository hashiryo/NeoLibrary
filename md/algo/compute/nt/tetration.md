---
title: テトレーション
documentation_of: ../../neo/algo/compute/nt/tetration.hpp
---

$a$ の $b$ 回のテトレーション $a \uparrow\uparrow b$ を $m$ で割った余りを返します。

```cpp
constexpr unsigned long long mod_tetration(unsigned long long a, unsigned long long b, unsigned long long m);
```

$a \uparrow\uparrow 1 = a$、$a \uparrow\uparrow b = a^{a \uparrow\uparrow (b-1)}$ で、$0^0 = 1$、$b = 0$ は 1、$m = 1$ は 0 です。

指数をそのまま作れないので、指数を $L(m)$ で割った余りに落とします。$L$ は、$2^e$ では $e \le 3$ のとき $\varphi(2^e)$、$e \ge 4$ のとき $\lambda(2^e) = 2^{e-2}$、奇数の素数べきでは $\varphi(p^e)$ とし、素数べきごとの値を lcm で束ねたものです。$L(m)$ は $\log_2 m$ 以上なので、$a$ と $m$ が互いに素でなくても、指数が $L$ 以上なら $a^E \bmod m$ は $E \bmod L$ で決まります。再帰の結果が 0 のときは、真の指数が $L$ の正の倍数なので $L$ を使います。

実行時は、$m < 2^{32}$ の因数分解を 65536 未満の素数の表で確保なしに試し割りします。テトレーションの再帰では法が小さくなるので、この道が大半です。素数は先に `is_prime` で分けます。表は constexpr の段数に収まらないので実行時に作り、定数評価のときは `factorize_pairs` の道を通します。`mod_tetration` は constexpr です。

`yosupo-tetration-mod` ($T = 10^3$、$m \le 10^9$) の手元の Apple M2 Max で、最大ケースは 22 ms、今の Library は 35 ms でした。`self-oeis-mod-tetration` の 3 組の OEIS の値と、`self-constexpr-mod-tetration` の 2 つの static_assert も通しています。
