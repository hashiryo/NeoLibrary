---
title: 2 次元の点
documentation_of: ../../neo/base/core/geo/Point.hpp
---

2 次元の点 `Point<K>` と、その和・差・スカラー倍、内積・外積、距離の 2 乗を置きます。`K` は整数でも浮動小数でもよく、整数なら判定が exact になります。

```cpp
template <class K> struct Point {
  K x, y;
  Point(K x = K(), K y = K());
  // +, -, 単項 -, スカラー倍, +=, -=, *=, /=, !(90 度回す)
  bool operator<(const Point& p) const;   // 座標の辞書式
  bool operator==(const Point& p) const;  // 座標が同じ
};
template <class K> using wide_t = typename geo_internal::wide<K>::type;
template <class K> constexpr int sgn(K x);
template <class K> constexpr wide_t<K> dot(const Point<K>& p, const Point<K>& q);
template <class K> constexpr wide_t<K> cross(const Point<K>& p, const Point<K>& q);
template <class K> constexpr wide_t<K> norm2(const Point<K>& p);
template <class K> constexpr wide_t<K> dist2(const Point<K>& p, const Point<K>& q);
```

`dot` と `cross` は、`K` より 1 段広い型 `wide_t<K>` で計算します。`int` なら `long long`、`long long` なら `__int128`、`float` なら `double`、`double` なら `long double` です。座標が 32 bit に収まる整数なら、外積は 64 bit で足ります。64 bit の座標でも、外積と三角形の面積の符号は `__int128` で exact に決まります。

`sgn(x)` は、整数なら `x` の符号をそのまま返し、浮動小数なら絶対値が `1e-9` 以下を 0 とみなします。位置の判定は `sgn(cross(b - a, c - a))` のように `sgn` を通して書いてください。

`Point` の比較 `operator<` と `operator==` は座標の比較そのもので、EPS を入れません。並べ替えの順序が推移的であることと、同じ座標の点を重複として扱えることを優先しています。EPS を入れた比較は推移的でないことがあり、`std::sort` の前提を壊します。
