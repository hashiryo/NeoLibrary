---
title: 冪の表
documentation_of: ../../neo/base/core/nt/pow_table.hpp
---

$i^d$ ($0 \le i \le n$) の表を作ります。

```cpp
template <class T> std::vector<T> pow_table(int n, unsigned long long d);
```

$0^0 = 1$ とします。`T` は掛け算を持つ型 (`ZMod` など) で、整数から作れることが要ります。補間 (`neo/algo/compute/poly/lagrange.hpp`) と合わせて、べき乗和や $\sum r^i i^d$ を求めるときに使います。

## 中身

$i^d$ は完全乗法的なので、素数 $p$ だけ $p^d$ を冪で求めます。合成数 $m$ は、最小の素因数 $p$ について $m^d = (m/p)^d p^d$ とします。$\sqrt{n}$ 以下の素数で $2^{15}$ 個ずつの区画を篩い、合成数 $m$ には $m/p$ と $p$ の番号を記録します。素数は大きい順に書くので、最後に書いた最小の素数が残ります。$m/p$ は区画の前か、区画の中で $m$ より前にあるので、前から埋めればそろっています。

`ZMod` (法は奇数) では、2 番目以降の区画の $m/p$ はどれも前の区画にあります。そこで区画の素数を集め、32 個ずつ AVX2 の Montgomery の積で冪を求めます (8 個ずつの 4 本を並べて待ちを隠す)。値は答えの `vector` に直接書きます。
