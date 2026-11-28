---
title: ユークリッド距離の MST
documentation_of: ../../neo/algo/compute/geo/euclidean_mst.hpp
---

ユークリッド距離による最小全域木の辺を返します。

```cpp
template <class K> std::vector<std::pair<int, int>> euclidean_mst(const std::vector<Point<K>>& ps);
```

最小全域木の辺はドロネー三角形分割の辺に入るので、ドロネーの辺だけを距離の 2 乗で並べて Kruskal します。距離の 2 乗の順と距離の順は同じなので、2 乗のまま比べてよいです。計算量は $O(n \log n)$ です。

同じ距離の辺が複数あるときは、そのうちの 1 つの MST を返します。同じ座標の点は、ドロネーの長さ 0 の辺で繋がります。
