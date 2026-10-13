---
title: k 乗根 mod p
documentation_of: ../../neo/number_theory/kth_root_mod.hpp
---

素数 $p$ を法として、$X^k \equiv y$ となる $X$ を 1 つ求めます。

```cpp
long long kth_root_mod(unsigned long long k, unsigned long long y, unsigned long long p);
```

$0 \le X < p$ の $X$ を 1 つ返し、無ければ -1 を返します。$0^0 = 1$ とします。$2 \le p < 2^{63}$ の素数 `p` を受け、範囲の外は assert で止めます。$p$ が素数であることは確かめません。`y` は $p$ 以上でもかまいません ($p$ で割った余りに直します)。

$m = p - 1$、$g = \gcd(k, m)$ とすると、解があるのは $y^{m/g} = 1$ のときです。$k = gu$ として $v = u^{-1} \bmod (m/g)$ とすると、$y^v$ の任意の $g$ 乗根が答えになります。$g$ の素因数 $r$ ごとに、群を $r$ の Sylow 部分群とそれ以外に分けて $r^e$ 乗根を取ります。Sylow 部分群の中は Pohlig-Hellman で $r$ 進の桁ごとに離散対数を求め、桁は位数 $r$ の元の離散対数を baby-step giant-step で引きます ($r \le 64$ なら表を順に見ます)。

剰余の掛け算は Montgomery で、$p < 2^{30}$ なら 32 bit、それより大きければ 64 bit の道を通ります。手間は $\log p$ 回ほどの累乗と、$r^2 \mid p - 1$ となる素因数 $r$ についての $\sqrt{r}$ 回ほどの掛け算です。$p \le 10^9$ の 5000 問は、procon-judge の yosupo-kth-root-mod の x64 の CI のランナーで 2〜5 ms ほどでした。今の Library の mod_kth_root の 2.4 倍ほどの速さです。
