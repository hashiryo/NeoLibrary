---
title: 平方根 mod p
documentation_of: ../../neo/number_theory/sqrt_mod.hpp
---

素数 $p$ を法として、$X^2 \equiv y$ となる $X$ を 1 つ求めます。

```cpp
long long sqrt_mod(unsigned long long y, unsigned long long p);
```

$0 \le X < p$ の $X$ を 1 つ返し、無ければ -1 を返します。2 つの根のどちらを返すかは決めていません。$2 \le p < 2^{63}$ の素数 `p` を受け、範囲の外は assert で止めます。$p$ が素数であることは確かめません。`y` は $p$ 以上でもかまいません ($p$ で割った余りに直します)。

$p \equiv 3 \pmod 4$ なら $y^{(p+1)/4}$、$p \equiv 5 \pmod 8$ なら Atkin の方法で求め、2 乗して確かめます。$p \equiv 1 \pmod 8$ のときは、奇数 $Q$ で $p - 1 = 2^s Q$ と書き、$s \le 8$ なら Tonelli-Shanks を使います。平方非剰余は、相互法則を使って小さい素数の表から探します。$s > 8$ なら、Jacobi 記号で平方非剰余を除いてから $a^2 - y$ が平方非剰余になる $a = 2^j$ を探し、和が 2、積が $y / a^2$ の Lucas 列の $V_{(p+1)/2}$ を梯子で求めます。Cipolla の方法を Lucas 列に直した形で、梯子の 1 段は掛け算 1 回ぶんの深さです。梯子の選択は分岐にしているので、同じ $p$ が続くと分岐が当たります。

剰余の掛け算は Montgomery で、$p < 2^{30}$ なら 32 bit、それより大きければ 64 bit の道を通ります。手間は $\log p$ 回ほどの掛け算の累乗か梯子と、平方非剰余を探す Jacobi 記号です。$p \le 10^9$ の $10^5$ 問は、procon-judge の yosupo-sqrt-mod の x64 の CI のランナーで最大 16 ms ほどでした ($p = 998244353$ のケース)。今の Library の mod_sqrt の 3 倍ほどの速さです。
