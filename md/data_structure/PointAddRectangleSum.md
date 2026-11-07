---
title: 点に足せる長方形の和
documentation_of: ../../neo/data_structure/PointAddRectangleSum.hpp
---

座標が先に決まった 2 次元の点の集合について、点の重みに足すことと、長方形 `[xl, xr) × [yl, yr)` に入る点の重みの和を引くことを、交互に行えます。あとから足す点の座標も、作るときに渡しておきます。

```cpp
template <class T = int, class W = long long> class PointAddRectangleSum;
PointAddRectangleSum(const std::vector<T>& x, const std::vector<T>& y, const std::vector<W>& w);  // 重み w[i] の点 (x[i], y[i]) の集合から作る
int size() const;
void add(int i, W v);                 // 点 i (作るときの添字) の重みに v を足す
W sum(T xl, T xr, T yl, T yr) const;  // xl <= x < xr かつ yl <= y < yr の点の重みの和。xl >= xr や yl >= yr なら 0
```

T は整数の型で、負の値も 64 bit の値も使えます。W は足し引きのできる型です。同じ座標の点が何個あってもかまいません。あとから足す点は、作るときに重み 0 で入れておき、その添字に `add` で足します。計算量は、構築が $O(N \log N)$、`add` と `sum` が $O(\log^2 N)$ です。重みを変えないなら [RectangleSum](RectangleSum.md) の方が速く動きます。

## 中身

点を x の順に並べ、並べた順の y の順位の列から、1 段で 1 bit ずつ分ける 2 分の wavelet matrix を作ります。段ごとに、分けたあとの並びの重みを BIT で持ち、元の並びの重みの BIT も 1 本持ちます。`add` は、点の位置を段ごとに rank で追い、各段の BIT に足します。`sum` は、y の 2 つの境の bit が同じ段を 1 回で下り、違う段から下は 2 本を同じループで下りて、境の bit が 1 の段だけ 0 の側へ行く範囲を BIT で引きます。bit が同じ段では、2 つの境の和が打ち消し合うので BIT を引きません。x の範囲と y の境の順位は、[RectangleCount](RectangleCount.md) と同じ表で引きます。

方式は procon-judge の yosupo-point-add-rectangle-sum ($N, Q \le 10^5$) で比べて選びました。この型を使う提出の最大は、EPYC 7763 の x64-clang で 113.4 ms ほど、x64-gcc で 130.5 ms ほどでした。今の Library の 2 次元セグメント木は、同じ回で 500 ms ほどです。x64-clang では写し元の提出 (wm_bit) と並び、x64-gcc では 1 割ほど遅く出ています (wm_bit は 116.9 ms)。遅いのは点の数の多い、座標のばらばらなケースで、座標の集まったケースではこの型の方が速く出ました。点の並べ替えと y の詰め方 (写し元は 10 bit ずつの基数ソート) の違いと見ています。段ごとの BIT を、塊の中の累積と塊ごとの和の BIT の 2 段に替えた形も試しましたが、x64 では BIT のままの方が速く出ました。記録は algo-notes の notes/wavelet-matrix.md にあります。
