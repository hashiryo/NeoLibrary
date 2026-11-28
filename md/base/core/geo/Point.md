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
template <class K> long double dist(const Point<K>& p, const Point<K>& q);
template <class K> struct Affine { K a00, a01, a10, a11; Point<K> b; Point<K> operator()(const Point<K>& p) const; Affine operator*(const Affine& r) const; };
enum class Ccw { counter_clockwise, clockwise, online_back, online_front, on_segment };
const char* ccw_name(Ccw c);
template <class K> Ccw ccw(const Point<K>& p0, const Point<K>& p1, const Point<K>& p2);
```

`dot` と `cross` は、`K` より 1 段広い型 `wide_t<K>` で計算します。`int` なら `long long`、`long long` なら `__int128`、`float` なら `double`、`double` なら `long double` です。座標が 32 bit に収まる整数なら、外積は 64 bit で足ります。64 bit の座標でも、外積と三角形の面積の符号は `__int128` で exact に決まります。

`sgn(x)` は、整数なら `x` の符号をそのまま返し、浮動小数なら絶対値が `1e-9` 以下を 0 とみなします。位置の判定は `sgn(cross(b - a, c - a))` のように `sgn` を通して書いてください。

`Point` の比較 `operator<` と `operator==` は座標の比較そのもので、EPS を入れません。並べ替えの順序が推移的であることと、同じ座標の点を重複として扱えることを優先しています。EPS を入れた比較は推移的でないことがあり、`std::sort` の前提を壊します。

`ccw(p0, p1, p2)` は 3 点の向きを返します。反時計回り、時計回り、`p2` が線分 `p0 p1` の後ろ側、前側、線分の上、の 5 つです。`operator<<` は `COUNTER_CLOCKWISE` などの名前を出します。

`Affine` は $p \mapsto A p + b$ のアフィン変換です。`operator*` は変換の合成です。`dist` は距離で、距離の 2 乗の平方根を `long double` で取ります。
