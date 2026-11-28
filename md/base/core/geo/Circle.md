---
title: 円
documentation_of: ../../neo/base/core/geo/Circle.hpp
---

円と、その交点・接線・共通接線・外接円・内接円を置きます。浮動小数を受けます。

```cpp
template <class R> struct Circle {
  Point<R> o;  // 中心
  R r;         // 半径
  Circle();
  Circle(const Point<R>& o, R r);
  R area() const;
  int where(const Point<R>& p) const;  // +1: 内, 0: 上, -1: 外
  int where(const Line<R>& l) const;   // +1: 2 点で交わる, 0: 接する, -1: 交わらない
  bool in(const Circle& c) const;
  std::vector<Line<R>> tangent(const Point<R>& p) const;
};
template <class R> std::vector<Point<R>> cross_points(const Circle<R>& c, const Circle<R>& d);
template <class R> std::vector<Point<R>> cross_points(const Circle<R>& c, const Line<R>& l);
template <class R> std::vector<Point<R>> cross_points(const Circle<R>& c, const Segment<R>& s);
template <class R> Circle<R> circumscribed_circle(const Point<R>& A, const Point<R>& B, const Point<R>& C);
template <class R> Circle<R> inscribed_circle(const Point<R>& A, const Point<R>& B, const Point<R>& C);
template <class R> std::vector<Line<R>> common_tangent(const Circle<R>& c, const Circle<R>& d);
```

`where` は点や直線との位置関係を返します。`in` は円 `c` がこの円の中に収まるかです。`tangent` は点 `p` を通る接線を返します。

`cross_points` は、円と円なら 2 点・接点 1 点・空・同じ円なら円周上の 3 点、円と直線や線分なら交点を返します。`common_tangent` は 2 円の共通接線を返します。

`area` は円の面積です。円と円、円と多角形の共通部分の面積は [intersection_area](intersection_area.md) にあります。
