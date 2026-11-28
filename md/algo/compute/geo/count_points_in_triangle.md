---
title: 三角形の内部の点の数
documentation_of: ../../neo/algo/compute/geo/count_points_in_triangle.hpp
---

点の集合 `a` の 3 点が作る三角形の内部 (境界を除く) にある、点の集合 `b` の点の数を、クエリごとに返します。

```cpp
template <class K>
std::vector<int> count_points_in_triangle(
    const std::vector<Point<K>>& a,
    const std::vector<Point<K>>& b,
    const std::vector<std::array<int, 3>>& qs);
```

各クエリは `a` の 3 点の添字です。3 点が一直線のときは 0 を返します。

各 `(a[x], a[y])` について、`a[x] -> a[y]` の左側にある `b` の点をビットセットで持ちます。クエリでは、三角形を反時計回りに直し、3 辺のビットセットの積を数えます。`b` の点が境界にあるときは、どの辺についても左側でないので数えません。

計算量は $O(N^2 M + Q M / 64)$ 時間、$O(N^2 M / 64)$ メモリです ($N = |a|$、$M = |b|$)。$N$、$M$ が数百のときの形で、もっと大きいときは角度で並べて Fenwick で数える形にします。
