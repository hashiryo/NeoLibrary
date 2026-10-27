---
title: Wavelet matrix (区間の重みの和)
documentation_of: ../../neo/data_structure/WaveletMatrixSum.hpp
---

静的な整数の列について、区間で値の範囲に入るものの個数と重みの和、小さい方 (大きい方) から k 個の重みの和を答えます。列は作ったあと変えられません。

```cpp
template <class T = int, class W = long long> class WaveletMatrixSum;
explicit WaveletMatrixSum(const std::vector<T>& a);                    // 値そのものを重みにする
WaveletMatrixSum(const std::vector<T>& a, const std::vector<W>& w);    // a[i] の重みを w[i] にする
int size() const;
W sum(int l, int r) const;                                    // a[l, r) の重みの和
std::pair<int, W> count_sum_lt(int l, int r, T x) const;      // 値が x 未満のものの個数と重みの和
std::pair<int, W> count_sum_le(int l, int r, T x) const;      // x 以下
std::pair<int, W> count_sum(int l, int r, T lo, T hi) const;  // 値が [lo, hi) に入るもの。lo >= hi なら {0, 0}
W sum_smallest(int l, int r, int k) const;                    // 値の小さい方から k 個の重みの和。0 <= k <= r - l
W sum_largest(int l, int r, int k) const;                     // 大きい方から k 個
T kth_smallest(int l, int r, int k) const;                    // 小さい方から k 番目 (0 始まり)
T kth_largest(int l, int r, int k) const;
```

T は整数の型で、負の値も 64 bit の値も使えます。W は足し引きのできる型です。同じ値が複数あるとき、小さい方から k 個に入るのは位置の前のものからです。計算量は、構築が $O(N \log \sigma)$、1 回の操作が $O(\log \sigma)$ です。

## 使い分け

和の要らない k 番目、個数、前後の値は、[WaveletMatrix](WaveletMatrix.md) の方が速く動きます。この型は段の数が 2 倍で、累積の配列 (段ごとに $N + 1$ 個の W) を持つので、メモリも多く使います。

2 次元の点の重みの和は、点を x の順に並べ、y の列と重みの列からこの型を作って、x の範囲を二分探索で位置の範囲に直してから `count_sum` で引きます。

## 中身

値の詰め方と値から順位への変換は WaveletMatrix と同じです。1 段で 1 bit ずつ分ける 2 分の形で、bit 列は 64 bit の語の列と語ごとの累積の列で持ちます。段ごとに、分けたあとの並びの重みの累積を持ち、値の上限の bit が 1 の段で、0 の側へ行くものの和を足します。累積は、bit が 1 の段だけ読みます (分岐で飛ばす)。全部の段で読んでマスクで捨てるより、主記憶への読み出しが半分で済み、3 割ほど速く出ました。小さい方から k 個の和は、k 番目と同じく下り、1 の側へ行く段で 0 の側の和を足し、最後の段の並びの累積で残りを足します。2 MB 以上の配列は、Linux では 2 MB 境界に取って huge page を頼みます。

方式は procon-judge の yosupo-static-range-sum-with-upper-bound などで比べて選びました。4 分の形に記号ごとの累積を持たせる形も試しましたが、累積の配列が 1.5 倍になり、x64 では和を引く 4 問ともこちらが速く出ました。atcoder-abc339-g ($N, Q \le 2 \times 10^5$、オンライン) の最大は、Intel Xeon の x64-gcc で 76.1 ms ほどです。今の Library の SegmentTree_2D の 11 分の 1 にあたります。atcoder-abc281-e (小さい方から K 個の和) は 63.0 ms ほどで、2 つの multiset で窓をずらす想定解より 2 割ほど速く出ました。記録は algo-notes の notes/wavelet-matrix.md にあります。
