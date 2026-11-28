---
title: 直線
documentation_of: ../../neo/base/core/geo/Line.hpp
---

直線 $p + t d$ と、その上の射影・交点・距離を置きます。

```cpp
template <class K> struct Line {
  Point<K> p, d;
  Line();
  Line(const Point<K>& p, const Point<K>& d);
  Line(K a, K b, K c);          // ax + by + c = 0
  int where(const Point<K>& q) const;   // +1: 左, 0: 上, -1: 右
  Point<K> project(const Point<K>& q) const;
  std::tuple<K, K, K> coef() const;
};
template <class K> Line<K> line_through(const Point<K>& p, const Point<K>& q);
template <class K> bool is_parallel(const Line<K>& l, const Line<K>& m);
template <class K> bool is_orthogonal(const Line<K>& l, const Line<K>& m);
template <class K> std::vector<Point<K>> cross_points(const Line<K>& l, const Line<K>& m);
template <class K> Line<K> bisector(const Point<K>& p, const Point<K>& q);
template <class K> wide_t<K> dist2(const Line<K>& l, const Point<K>& p);
template <class K> long double dist(const Line<K>& l, const Point<K>& p);
template <class K> Affine<K> reflect(const Line<K>& l);
```

`where` は直線の左側を $+1$、直線上を $0$、右側を $-1$ とします。`project` は点の直線への射影です。`cross_points` は 1 点で交わるならその点、平行で交わらないなら空、同じ直線なら 2 点を返します。

`reflect(l)` は直線 `l` に対する折り返しのアフィン変換を返します。点 `p` の折り返しは `reflect(l)(p)` です。

`dist2` は距離の 2 乗で、整数なら切り捨ての割り算を含みます。距離そのものが要るときは `dist` を使ってください。
