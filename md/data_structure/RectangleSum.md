---
title: 長方形に入る点の重みの和
documentation_of: ../../neo/data_structure/RectangleSum.hpp
---

重み付きの 2 次元の点の集合について、長方形 `[xl, xr) × [yl, yr)` に入る点の重みの和を答えます。点は作ったあと変えられません。

```cpp
template <class T = int, class W = long long> class RectangleSum;
RectangleSum(const std::vector<T>& x, const std::vector<T>& y, const std::vector<W>& w);  // 重み w[i] の点 (x[i], y[i]) の集合から作る
int size() const;
W sum(T xl, T xr, T yl, T yr) const;  // xl <= x < xr かつ yl <= y < yr の点の重みの和。xl >= xr や yl >= yr なら 0
```

T は整数の型で、負の値も 64 bit の値も使えます。W は足し引きのできる型です。同じ座標の点が何個あってもかまいません。計算量は、構築が $O(N \log N)$、1 回の問い合わせが $O(\log N)$ です。点の数だけが要るなら [RectangleCount](RectangleCount.md) の方が速く動きます。

## 中身

点を x の順に並べ、並べた順の y の列と重みの列から [WaveletMatrixSum](WaveletMatrixSum.md) を作ります。x の範囲は、[RectangleCount](RectangleCount.md) と同じ表で位置の範囲に直します。並べ方も同じで、添字の並べ替えの `argsort` (neo/sort/argsort.hpp) を使います。

方式は procon-judge の yosupo-rectangle-sum ($N, Q \le 2 \times 10^5$) で比べて選びました。写し元の提出 (wm_sep_brtbl) の x64-gcc の最大は、EPYC 9V45 で 119.1 ms ほどでした。同じ回で、x の範囲を `std::lower_bound` で引く形は 139.0 ms ほどかかりました。記録は algo-notes の notes/wavelet-matrix.md にあります。
