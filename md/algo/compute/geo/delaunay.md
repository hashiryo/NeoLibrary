---
title: ドロネー三角形分割
documentation_of: ../../neo/algo/compute/geo/delaunay.hpp
---

点の集合のドロネー三角形分割の辺を返します。

```cpp
template <class K> std::vector<std::pair<int, int>> delaunay_edges(const std::vector<Point<K>>& ps);
```

返すのは辺の両端の点の添字の組です。同じ座標の点があるときは 1 つにまとめ、まとめた点と元の点を結ぶ長さ 0 の辺も返します。

点を $x$ の順に並べ、左右の三角形分割を共通接線で繋ぎます (Guibas–Stolfi の分割統治)。辺は quad-edge 構造で両方向に持ち、辺を消すときは番号を使い回します。計算量は $O(n \log n)$ です。

座標が整数なら、外接円の内側かの判定は `wide_t<K>` の行列式で exact になります。`wide_t<K>` は、座標の差の 4 乗の 12 倍が入る大きさが要ります。
