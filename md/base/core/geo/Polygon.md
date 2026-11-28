---
title: 多角形
documentation_of: ../../neo/base/core/geo/Polygon.hpp
---

反時計回りの多角形の面積・内外判定・凸かどうか・距離を置きます。

```cpp
template <class K> class Polygon {
  Polygon();
  Polygon(const std::vector<Point<K>>& ps);
  int prev(int i) const;
  int next(int i) const;
  const Point<K>& operator[](int i) const;
  std::size_t size() const;
  std::vector<Segment<K>> edges() const;
  bool is_convex() const;
  wide_t<K> area() const;
  wide_t<K> area2() const;
  int where(const Point<K>& p) const;   // +1: 内, 0: 周上, -1: 外
  bool in(const Segment<K>& s) const;
  bool out(const Segment<K>& s) const;
};
template <class K> wide_t<K> dist2(const Polygon<K>& g, const Point<K>& p);
template <class K> long double dist(const Polygon<K>& g, const Point<K>& p);
```

頂点は反時計回りに直して持ちます。`area` は 2 倍の面積を 2 で割ったもので、整数なら切り捨てです。`area2` は 2 倍の面積そのものです。`where` は点が内部なら $+1$、周上なら $0$、外部なら $-1$ を返します。`in` と `out` は線分が多角形の内側か外側かです。

`dist2` と `dist` は、点・直線・線分・多角形との距離です。点が多角形の内側か周上なら 0 です。
