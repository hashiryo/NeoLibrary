---
title: Wavelet matrix (区間の k 番目と個数)
documentation_of: ../../neo/data_structure/WaveletMatrix.hpp
---

静的な整数の列について、区間の k 番目の値、値の範囲に入る個数、前後の値を答えます。列は作ったあと変えられません。

```cpp
template <class T = int> class WaveletMatrix;
explicit WaveletMatrix(const std::vector<T>& a);
int size() const;
T kth_smallest(int l, int r, int k) const;        // a[l, r) の小さい方から k 番目 (0 始まり)。0 <= k < r - l
T kth_largest(int l, int r, int k) const;         // 大きい方から k 番目
int count_lt(int l, int r, T x) const;            // a[l, r) のうち x 未満の個数
int count_le(int l, int r, T x) const;            // x 以下の個数
int count(int l, int r, T lo, T hi) const;        // 値が [lo, hi) に入る個数。lo >= hi なら 0
int freq(int l, int r, T x) const;                // ちょうど x の個数
std::optional<T> max_lt(int l, int r, T x) const; // x 未満の最大。無ければ nullopt
std::optional<T> max_le(int l, int r, T x) const; // x 以下の最大
std::optional<T> min_ge(int l, int r, T x) const; // x 以上の最小
std::optional<T> min_gt(int l, int r, T x) const; // x より大きい最小
```

T は整数の型で、負の値も 64 bit の値も使えます。区間は $0 \le l \le r \le N$ で使います。空の区間を渡すと、個数が 0、前後の値が nullopt になります。計算量は、構築が $O(N \log \sigma)$、1 回の操作が $O(\log \sigma)$ です ($\sigma$ は異なる値の数)。前後の値は、個数と k 番目を 1 回ずつ引いて求めます。

## 使い分け

重みの和 (値が x 未満のものの和、小さい方から k 個の和) が要るなら [WaveletMatrixSum](WaveletMatrixSum.md) を使います。この型より段の数が倍ですが、和の配列の読み出しが少なく、和を引く問題ではそちらが速く動きます。

2 次元の点の長方形の個数は、点を x の順に並べ、y の列からこの型を作って、x の範囲を二分探索で位置の範囲に直してから `count` で数えます。

ある値の個数 (`freq`) だけを引く問題もあります。procon-judge の yosupo-static-range-frequency では、値ごとに出てくる位置を並べておいて二分探索する形が 1.7 倍ほど速く動きました。

## 中身

値は内部で異なる値の順位に詰めます。32 bit の整数は (値, 位置) の基数ソート、ほかの型は std::sort で並べます。値から順位への変換は、値の範囲を最大 65536 個の区間に分けた表で区間を決め、その中を cmov の二分探索で探します。

1 段で 2 bit ずつ分ける 4 分の形です。順位の上の 2 bit から段ごとに、その 2 bit (0 から 3 の記号) の順に安定に分けるので、段の数は 2 分の形の半分です。段ごとに 64 個ずつのブロックを、記号の上の bit の列と下の bit の列と、ブロックより前の記号ごとの数の 32 byte に詰めます。区間の記号ごとの数は、両端のブロックを 1 つずつ読み、popcount で求めます。k 番目は、記号ごとの累計と k を比べたマスクで記号を選び、分岐を使いません。2 MB 以上の配列は、Linux では 2 MB 境界に取って huge page を頼みます。

方式は procon-judge の yosupo-range-kth-smallest などで比べて選びました。yosupo-range-kth-smallest ($N, Q \le 2 \times 10^5$) の最大は、EPYC 9V45 の x64-gcc で 25.9 ms ほどです。今の Library の WaveletMatrix の 3 分の 1 にあたります。1 段で 1 bit ずつ分ける形より 2 割速く出ました。前後の値の aoj-1549、範囲の個数の aoj-2674、yosupo-static-range-count-distinct でも、今の Library の 1.8 倍から 4.3 倍の速さでした。記録は algo-notes の notes/wavelet-matrix.md にあります。
