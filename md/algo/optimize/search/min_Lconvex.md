---
title: L 凸な関数の最小 (座標ごとの山登り)
documentation_of: ../../neo/algo/optimize/search/min_Lconvex.hpp
---

L 凸な f の最小を、座標を alpha ずつ動かす山登りで求めます。alpha は改善しなくなるまで半分にしていきます。

```cpp
template <class TD, class TR, class F> std::pair<TR, std::vector<TD>> min_Lconvex(const F& f, std::vector<TD> x, TD alpha);
```

返すのは (f の最小値, その x) です。各段で 2^n 通りの部分集合を試すので、座標の数は 30 程度までを想定します。

## 使い方

```cpp
std::vector<long long> x(n, 0);
auto [value, best] = min_Lconvex<long long, long long>(f, x, 1LL << 30);
```
