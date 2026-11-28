---
title: 線分
documentation_of: ../../neo/base/core/geo/Segment.hpp
---

線分 $p q$ と、その上の点の判定・最近点・交点・距離を置きます。

```cpp
template <class K> struct Segment {
  Point<K> p, q;
  bool on(const Point<K>& r) const;
  Point<K> closest_point(const Point<K>& r) const;
  long double length() const;
};
template <class K> std::vector<Point<K>> cross_points(const Segment<K>& s, const Line<K>& l);
template <class K> std::vector<Point<K>> cross_points(const Segment<K>& s, const Segment<K>& t);
enum class Intersection { crossing, touching, disjoint, overlap };
template <class K> Intersection intersection(const Segment<K>& s, const Segment<K>& t);
template <class K> wide_t<K> dist2(const Segment<K>& s, const Point<K>& p);
template <class K> long double dist(const Segment<K>& s, const Point<K>& p);
template <class K> long double dist(const Segment<K>& s, const Segment<K>& t);
```

`on` は点が線分の上 (両端を含む) かどうかです。`closest_point` は線分上の最も近い点を返します。`cross_points` は 1 点で交わるならその点、交わらないなら空、同じ直線上で重なるなら重なりの端点を返します。

`intersection` は 2 つの線分の位置関係を `crossing`、`touching`、`disjoint`、`overlap` の 4 つに分けます。`dist2` と `dist` は、点・直線・線分との距離です。
