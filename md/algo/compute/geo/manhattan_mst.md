---
title: マンハッタン距離の MST
documentation_of: ../../neo/algo/compute/geo/manhattan_mst.hpp
---

平面の点の集合の、マンハッタン距離による最小全域木を返します。

```cpp
template <class K> struct ManhattanMst {
  wide_t<K> weight;                    // 辺の重みの和
  std::vector<std::pair<int, int>> edges;  // 辺の両端の点の添字
};
template <class K> ManhattanMst<K> manhattan_mst(std::vector<Point<K>> ps);
```

点を 45 度回す操作を 4 通り試し、各点について $x$ と $y$ の両方が大きくない点のうち $x - y$ が最大のものを Fenwick で探して候補の辺にします。候補の辺はたかだか $4n$ 本で、これを距離で並べて Kruskal で MST を作ります。計算量は $O(n \log n)$ です。同じ距離の辺が複数あるときは、そのうちの 1 つの MST を返します。
