---
title: 凸多角形
documentation_of: ../../neo/base/core/geo/Convex.hpp
---

凸多角形を作り、直径・最遠点対・半平面での切断を置きます。

```cpp
template <class K> struct Convex : Polygon<K> {
  Convex();
  Convex(std::vector<Point<K>> ps, bool strict = true);
  std::pair<Point<K>, Point<K>> farthest_pair() const;
  long double diameter() const;
  Convex cut(const Line<K>& l, int side = 1) const;
  Convex cut(K a, K b, K c, int side = 1) const;
};
template <class K> std::pair<Point<K>, Point<K>> farthest_pair(const std::vector<Point<K>>& ps);
```

`Convex(ps, strict)` は点の集合の凸包を作ります。`strict` が `true` なら内角 180 度の点を頂点から外し、`false` なら辺上の点も頂点として残します。頂点は反時計回りです。

`farthest_pair` は最も遠い 2 頂点を回転カリパスで返し、`diameter` はその距離です。`cut(l, side)` は直線 `l` で切った凸多角形を返します。`side` が正なら左側、負なら右側との共通部分です。
