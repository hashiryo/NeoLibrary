---
title: 素数の個数
documentation_of: ../../neo/number_theory/prime_pi.hpp
---

$n$ 以下の素数の個数 $\pi(n)$ を数えます。

```cpp
unsigned long long prime_pi(unsigned long long n);
```

$0 \le n \le 10^{14}$ の `n` を受け、範囲の外は assert で止めます。constexpr ではありません。

Deléglise–Rivat の方法で数えます。$y = \alpha n^{1/3}$ 以下の素数で $\phi(n, \pi(y))$ の展開を打ち切り、ordinary leaves、special leaves、$P_2$ の 3 つの和に分けます。special leaves のうち $\pi$ の表で済まないもの (hard leaves) は、$[1, n / y]$ を区間ごとに篩いながら、残った数の個数で数えます。$\pi$ の表は同じ篩の走査で作ります。表で済む easy leaves のうち値が階段状に並ぶ部分は、Gourdon の反転で短い和に直します。

AVX2 が使え、gather の速い CPU では、$\pi$ の表を奇数だけで持ち、easy leaves の表引きと ordinary leaves と $P_2$ を 4 つずつ回します。商は double の掛け算で求め、表は gather で引きます。Intel の Haswell から Rocket Lake の世代は gather が遅いので (Downfall の対策の microcode で遅くなる世代を含む)、`__builtin_cpu_is` で見分けます。その世代では 30 の車輪の表を 1 つずつ引く道を通します。AVX2 が使えないときも同じ道です。

時間はおよそ $n^{2/3}$ の割合で伸び、$10^{11}$ で 2 ms ほどです (procon-judge の yosupo-counting-primes、x64 の CI のランナーで 1.4〜2.2 ms)。手元の Apple M2 (1 つずつ引く道) では、$10^{12}$ が 7 ms、$10^{13}$ が 30 ms、$10^{14}$ が 0.14 s でした。メモリは $n / y$ に比例し、$10^{14}$ で 21 MB (gather の道では 29 MB) ほど使います。
