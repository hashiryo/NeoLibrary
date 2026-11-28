---
title: 最近点対
documentation_of: ../../neo/base/core/geo/closest_pair.hpp
---

点の集合の最近点対を返します。

```cpp
template <class K> std::pair<Point<K>, Point<K>> closest_pair(std::vector<Point<K>> ps);
```

点を $y$ の順に並べ、分割統治します。統合のときは、$x$ の差が今の最短距離未満の点だけを見ます。計算量は $O(n \log n)$ です。同じ距離の点の組が複数あるときは、そのうちの 1 組を返します。
