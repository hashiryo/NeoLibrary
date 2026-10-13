---
title: 線形関数の床の和
documentation_of: ../../neo/number_theory/floor_sum.hpp
---

$\sum_{i=0}^{n-1} \lfloor (a i + b) / m \rfloor$ を求めます。

```cpp
long long floor_sum(long long n, long long m, long long a, long long b);
```

AtCoder Library の floor_sum と同じ約束で、$0 \le n < 2^{32}$、$1 \le m < 2^{32}$ を受け、範囲の外は assert で止めます。$a$ と $b$ は負の値も $m$ 以上の値も受けます。和は $2^{64}$ で割った余りを long long にして返します。

先に $a$ と $b$ を $m$ の剰余に直し、残りは互除法で求めます。1 段では、$a n + b$ を $m$ で割った商と余り、$\lfloor a/m \rfloor$、$\lfloor b/m \rfloor$ の 3 つの割り算をします。商を $a$ と $b$ を割る前の値から求めるので、3 つは互いに待ちません。AtCoder Library の形は $a \bmod m$ と $b \bmod m$ を作ってから割るので、割り算 2 回ぶんが 1 段の待ちになります。手間は $(a \bmod m, m)$ の互除法の段の数で、$O(\log m)$ です。procon-judge の yosupo-sum-of-floor-of-linear ($10^5$ 問) では、x64 の CI のランナーで最大 7 ms ほどでした。AtCoder Library と同じ形の 1.4 倍、今の Library の FloorSumTable の 29 倍ほどの速さです。
