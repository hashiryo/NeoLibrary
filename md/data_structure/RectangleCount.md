---
title: 長方形に入る点の数
documentation_of: ../../neo/data_structure/RectangleCount.hpp
---

2 次元の点の集合について、長方形 `[xl, xr) × [yl, yr)` に入る点の数を答えます。点は作ったあと変えられません。

```cpp
template <class T = int> class RectangleCount;
RectangleCount(const std::vector<T>& x, const std::vector<T>& y);  // 点 (x[i], y[i]) の集合から作る
int size() const;
int count(T xl, T xr, T yl, T yr) const;  // xl <= x < xr かつ yl <= y < yr の点の数。xl >= xr や yl >= yr なら 0
```

T は整数の型で、負の値も 64 bit の値も使えます。同じ座標の点が何個あってもかまいません。計算量は、構築が $O(N \log N)$、1 回の問い合わせが $O(\log N)$ です。

## 中身

点を x の順に並べ、並べた順の y の列から [WaveletMatrix](WaveletMatrix.md) を作ります。問い合わせでは、x の昇順の列のうち xl 未満のものの数と xr 未満のものの数を引いて、x の範囲を位置の範囲に直します。その範囲で y が `[yl, yr)` に入るものを、WaveletMatrix の `count` で数えます (2 つの境を一緒に下ります)。x の昇順の列から数を引く表は、WaveletMatrix が値から順位を引く表と同じ形です。値の範囲を最大 65536 個の区間に分けた表で区間を決め、その中を cmov の二分探索で探します。

並べ方は添字の並べ替えの `argsort` (neo/sort/argsort.hpp) です。4 byte の整数は基数ソートの形で並べ、ほかの型は添字の列を `std::stable_sort` で並べます。

方式は procon-judge の aoj-2426 ($N \le 5000$、$M \le 5 \times 10^5$) で比べて選びました。この型を使う提出の x64-gcc の最大は、EPYC 9V74 で 41.0 ms ほどで、同じ回の 4 分の形の提出 (wm4) の 60.7 ms より 3 割ほど速く出ました。WaveletMatrix を 8 分の形にした分と、y の 2 つの境を一緒に下りる分です。x の範囲を `std::lower_bound` で引いていたときは、EPYC 9V45 で 92.6 ms ほどかかり、wm4 (43.2 ms) の 2 倍でした。記録は algo-notes の notes/wavelet-matrix.md にあります。
