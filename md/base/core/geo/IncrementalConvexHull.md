---
title: 足すだけの動的凸包
documentation_of: ../../neo/base/core/geo/IncrementalConvexHull.hpp
---

点を足すたびに凸包を作り直す、削除の無い動的凸包です。

```cpp
template <class K> class IncrementalConvexHull {
  IncrementalConvexHull();
  std::size_t edge_size() const;
  wide_t<K> area() const;
  wide_t<K> area2() const;
  int where(const Point<K>& p) const;  // +1: 内, 0: 周上, -1: 外
  void insert(const Point<K>& p);
};
```

下側の鎖を `std::set` で持ち、点を足したときに要らなくなった点を左右から消します。上側の鎖は、点を反転して同じ構造で持ちます。`where` は今の凸包に対する点の位置です。`area2` は 2 倍の面積で、整数なら exact です。
