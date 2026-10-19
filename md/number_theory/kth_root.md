---
title: 整数の k 乗根
documentation_of: ../../neo/number_theory/kth_root.hpp
---

$\lfloor a^{1/k} \rfloor$ を求めます。

```cpp
unsigned long long kth_root(unsigned long long a, unsigned k);
```

$0 \le a < 2^{64}$ と $1 \le k \le 64$ を受け、$k$ が範囲の外なら assert で止めます。

浮動小数点で見積もった $r$ を、整数の $r^k$ と $(r+1)^k$ で 1 だけ直します。直しは分岐を使わない算術です。$a$ が $r^k$ のすぐそばにあると、double に丸めた $a$ からの見積もりが 1 大きく出るかどうかが半々に変わり、分岐にすると外れるためです。

- $k = 2$ は sqrt で見積もります。
- $k = 3$ は、$a = 2^e m$ の指数を 3 で割った余りと仮数の上位 7 bit で表を引き、3 次の多項式で立方根を出します。
- $k \ge 4$ は、表と 2 次の多項式で log2 と exp2 を出します。
- $k \ge 5$ は答えが 7132 以下なので、$r^k$ の表を引いて確かめます。

表はコンパイル時に作ります ($r^k$ の表は 1 万項ほど)。procon-judge の yosupo-kth-root-integer ($5 \times 10^5$ 問) では、x64 の CI のランナーで最大 4 ms ほどでした。std::pow で見積もる yosupo の最速提出の写しの 2.5 倍ほどの速さです。
