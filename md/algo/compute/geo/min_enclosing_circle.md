---
title: 最小包含円
documentation_of: ../../neo/algo/compute/geo/min_enclosing_circle.hpp
---

点の集合の最小包含円を、円を決める 1〜3 点 (支援集合) として返します。

```cpp
template <class K> struct MecSupport {
  int m;          // 点の数 (0〜3)
  Point<K> a, b, c;
  int side(const Point<K>& p) const;  // +1: 内, 0: 上, -1: 外
};
template <class K> MecSupport<K> min_enclosing_circle(std::vector<Point<K>> ps);
```

`m = 1` なら点 `a` の円 (半径 0)、`m = 2` なら `a` と `b` を直径とする円、`m = 3` なら `a`、`b`、`c` の外接円です。`side(p)` は、その円に対する `p` の位置を返します。

点をランダムな順に足し、いまの円の外に出た点があれば、その点を円周に載せた形で作り直します (Welzl の逐次追加)。期待計算量は $O(n)$ です。整数の座標では、直径の内外と外接円の incircle 行列式を `wide_t<K>` で計算するので、円周上の点の判定が exact になります。浮動小数では `sgn` の EPS で判定します。乱数は splitmix64 の固定の種を使うので、同じ入力なら同じ結果になります。

円の中心と半径は、支援集合から `sqrt` で作ります。中心と半径が要るときは、`m` と `a`、`b`、`c` から計算してください。
