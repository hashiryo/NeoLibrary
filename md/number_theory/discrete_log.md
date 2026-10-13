---
title: 離散対数
documentation_of: ../../neo/number_theory/discrete_log.hpp
---

$x^k \equiv y \pmod m$ となる最小の $k$ を求めます。

```cpp
long long discrete_log(unsigned long long x, unsigned long long y, unsigned long long m);
```

最小の $k \ge 0$ を返し、無ければ -1 を返します。$0^0 = 1$ とし、$m = 1$ なら 0 を返します。$1 \le m < 2^{32}$ を受け、範囲の外は assert で止めます。$m$ は合成数でもかまいません。`x` と `y` は $m$ 以上でもかまいません ($m$ で割った余りに直します)。

$x$ と $m$ が互いに素でない部分は、小さい $k$ を直に調べてから、$\gcd(x, m)$ で割って外します。残りは Pohlig-Hellman で、$x$ の位数を $\varphi(m)$ の素因数分解から求め、位数の素因数 $r$ ごとに $r$ 進の桁を求めて CRT でまとめます。桁は位数 $r$ の元の離散対数で、baby-step giant-step で引きます。表はほとんど空きにしておき、表に無い鍵を引くたびの分岐が外れないようにしています。$r > 2^{22}$ で群が巡回群 ($m$ が素数冪かその 2 倍) のときは、Pollard の rho 法 (区別点を使い、16 本の歩みを並べる) に切り替えます。

手間は位数のいちばん大きい素因数 $r$ について $\sqrt{r}$ 回ほどの掛け算です。$m \le 10^9$ の 100 問は、procon-judge の yosupo-discrete-logarithm-mod の x64 の CI のランナーで 11〜12 ms ほどでした。今の Library の DiscreteLogarithm の 10 倍ほどの速さです。
