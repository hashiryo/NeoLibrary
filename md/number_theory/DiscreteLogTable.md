---
title: 法が固定の離散対数
documentation_of: ../../neo/number_theory/DiscreteLogTable.hpp
---

素数 $p$ と原始根 $g$ を固定して、$g^x \equiv a \pmod p$ となる $x$ を表で引きます。

```cpp
DiscreteLogTable D(p, g);  // 表を作る
long long x = D(a);        // 0 ≤ x < p - 1。a ≡ 0 なら -1
```

$2 \le p < 2^{30}$ の素数 `p` と、その原始根 `g` を受けます。範囲の外と、`g` が原始根でないときは assert で止めます。$p$ が素数であることは確かめません。`a` は $p$ 以上でもかまいません ($p$ で割った余りに直します)。1 つの $p$ について多くの $a$ を引くときに使い、数問だけなら discrete_log のほうが速いです。

maspy の 1 問 O(1) の方法です。$a \equiv num / den \pmod p$ で、$|num|$ が $1.3 \times 10^6$ ほど以下、$den$ が $2^{16}$ 未満となる分数を表で引きます。分数の表は、$a / p$ を $10^6$ 等分した区画ごとに 1 つずつ持ちます。log の表を $|num|$ の上限まで前もって作っておき、1 問は表を 2 回引くだけで済ませます。log の表は、最初の 10 個の素数を Pohlig-Hellman で、$\sqrt{p}$ までの残りの素数を分数の表を使った歩みで求めます。$\sqrt{p}$ より大きい整数は、$p = (p / i) i + (p \bmod i)$ の関係で埋めます。$p < 2 \times 10^6$ なら $g$ の冪をたどって全部の log を表にします。

表は合わせて 9 MB ほどで、Linux では huge page を頼みます。procon-judge の yosupo-discrete-logarithm-fixed-mod は、$p \le 10^9$ で $10^6$ 問のケースです。表を作るのも含めて、x64 の CI のランナーで最大 20 ms ほどでした。
