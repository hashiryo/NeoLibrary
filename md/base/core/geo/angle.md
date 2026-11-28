---
title: 角度と偏角の比較
documentation_of: ../../neo/base/core/geo/angle.hpp
---

角度の変換、回転、偏角の比較を置きます。

```cpp
long double radian_to_degree(long double r);
long double degree_to_radian(long double d);
long double normalize_radian(long double r);
template <class K> long double angle(const Point<K>& p);
template <class K> long double angle(const Point<K>& p, const Point<K>& q);
template <class K> Affine<K> rotate(long double theta);
template <class K> Affine<K> rotate(const Point<K>& p, long double theta);
template <class K> Affine<K> rotate90(const Point<K>& p);
template <class K> class AngleComp;
```

`angle(p)` は点の偏角、`angle(p, q)` は 2 点のなす角です。`rotate` は原点または点 `p` のまわりの回転のアフィン変換を返します。`rotate90` は点 `p` のまわりの 90 度回転です。

`AngleComp` は偏角の比較で、$-\pi < \theta \le \pi$ の順に並べます。座標が整数なら外積の符号で比べるので exact です。`SegmentArrangement` の辺の並べ替えに使います。
