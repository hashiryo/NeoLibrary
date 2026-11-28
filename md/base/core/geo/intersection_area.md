---
title: 共通部分の面積
documentation_of: ../../neo/base/core/geo/intersection_area.hpp
---

円と円、円と多角形の共通部分の面積を返します。

```cpp
template <class R> R intersection_area(const Circle<R>& c, const Circle<R>& d);
template <class R> R intersection_area(const Circle<R>& c, const Polygon<R>& g);
```

円と円の面積は、2 つの扇形の和から凧形の面積を引いて作ります。接する近くでは `acos` の引数が 1 に近づいて精度を落とすので、角度は Heron の式で作った三角形の高さと射影から `atan2` で作ります。半径が $10^4$、距離が $10^4$ の接する近くのケースでも、`long double` (arm では 64 bit) で 1e-12 ほどの精度が出ます。

円と多角形の面積は、多角形の各辺について円の中心から見た符号付きの面積を足します。
