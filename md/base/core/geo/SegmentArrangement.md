---
title: 線分アレンジメント
documentation_of: ../../neo/base/core/geo/SegmentArrangement.hpp
---

線分の集合を交点で分けて平面グラフを作ります。

```cpp
template <class K> class SegmentArrangement {
  explicit SegmentArrangement(const std::vector<Segment<K>>& ss);
  std::size_t vertex_size() const;
  std::size_t face_size() const;
  std::span<const int> out_edges(int v) const;
  Point<K> point(int v) const;
  int vertex(const Point<K>& p) const;
  int to_v(int e) const;
  int from_v(int e) const;
  int twin(int e) const;
  int next(int e) const;
  int prev(int e) const;
  int incident_face(int e) const;
  int component_e(int f) const;
  wide_t<K> area2(int f) const;
  wide_t<K> area(int f) const;
  bool is_outside(int f) const;
};
```

線分の全ペアで交点を調べ、交点で線分を分けて平面グラフを作ります。交差の列挙は $O(n^2)$ です。辺は両方向に持ち、`e ^ 1` が逆方向です。`next(e)` は辺 `e` の行き先で反時計回りに次に来る辺で、面は `next` をたどる閉路です。`area2(f)` は面の面積の 2 倍で、向きが反時計回りなら正、時計回りなら負です。`is_outside(f)` は面積が負かどうかで、外側の面を表します。
