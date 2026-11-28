---
title: 各頂点の最遠点
documentation_of: ../../neo/algo/compute/geo/all_furthest_neighbors.hpp
---

狭義凸多角形の各頂点について、最も遠い頂点を返します。

```cpp
template <class K> std::vector<int> all_furthest_neighbors(const std::vector<Point<K>>& ps);
```

`ps` は狭義凸多角形の頂点を反時計回りに並べたもので、`n >= 3` を前提にします。返り値の `i` 番目は、`ps[i]` から最も遠い頂点の添字です。同じ距離の頂点が複数あるときは、そのうちの 1 つを返します。

距離の行列 $D(i, j) = |p_i - p_{j \bmod n}|^2$ ($i < j < i + n$) は単調なので、行ごとの最大を分割統治で求めます。計算量は $O(n \log n)$ です。

単純な回転カリパスは、頂点からの距離が境界に沿って単峰でない多角形で誤ります。Library Checker の `all_furthest_neighbors` には、その形のテストケース (`non_unimodal`) があります。
