---
title: Wavelet matrix (区間の k 番目と個数)
documentation_of: ../../neo/data_structure/WaveletMatrix.hpp
---

静的な整数の列について、区間の k 番目の値、値の範囲に入る個数、前後の値を答えます。列は作ったあと変えられません。

```cpp
template <class T = int, bool Compress = true> class WaveletMatrix;
explicit WaveletMatrix(std::vector<T> a);
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

`Compress = false` のときは、値を詰めず、そのまま順位として使います。値はどれも 0 以上 $2^{31}$ 未満にしてください。$\sigma$ は値の最大 + 1 になり、段の数もこれで決まります。前に同じ値が出た位置の列のように、値がもともと 0 以上 $N$ 以下の列で使うと、作るときに値を詰める手間と、クエリのたびに値から順位を引く手間が省けます。procon-judge の yosupo-static-range-count-distinct では、前に出た位置 + 1 の列をこの形で渡す提出が、EPYC 9V45 の x64-gcc で 52.7 ms ほどでした。同じ回の 4 分の形の提出 (wm4) は 65.4 ms です。値を詰め直す形では、wm4 の 1.3 倍ほどかかっていました。

## 使い分け

重みの和 (値が x 未満のものの和、小さい方から k 個の和) が要るなら [WaveletMatrixSum](WaveletMatrixSum.md) を使います。この型より段の数が倍ですが、和の配列の読み出しが少なく、和を引く問題ではそちらが速く動きます。

2 次元の点の長方形に入る個数は [RectangleCount](RectangleCount.md) で引けます。点を x の順に並べた y の列からこの型を作り、x の範囲を位置の範囲に直す表を添えたものです。

x 未満の最大と x 以上の最小を両方引くなら、`count_lt` で x 未満の個数 c を 1 回数え、`kth_smallest` で c - 1 番目と c 番目を引くと、個数を数え直さずに済みます。`max_lt` と `min_ge` を続けて呼ぶと、個数を 2 回数えます。

ある値の個数 (`freq`) だけを引く問題もあります。procon-judge の yosupo-static-range-frequency では、値ごとに出てくる位置を並べておいて二分探索する形が 1.7 倍ほど速く動きました。

## 中身

値は内部で異なる値の順位に詰めます。詰め方は座標圧縮の `compress` (neo/sort/compress.hpp) です。値から順位への変換は、値の範囲を最大 65536 個の区間に分けた表で区間を決め、その中を cmov の二分探索で探します。

1 段で 3 bit ずつ分ける 8 分の形です。順位の上の 3 bit から段ごとに、その 3 bit (0 から 7 の記号) の順に安定に分けるので、段の数は 2 分の形の 3 分の 1 です。段ごとに 64 個ずつのブロックを、記号の 3 本の bit の列と、ブロックより前の記号 t 未満の数 (t = 0..8) の 64 byte に詰めます。1 つのブロックがキャッシュの 1 行に収まります。区間の記号 t 未満の数は、両端のブロックを 1 つずつ読み、t 未満の印の列の popcount で求めます。k 番目は、両端で記号 1 から 7 未満の数を 7 個ずつ求め、区間の中の数と k を比べた 7 つの比較の和で記号を選び、分岐を使いません。個数を数えるときは、x の記号が段ごとに決まっているので、印の列の作り方を先に決めておけます。`count` の `[lo, hi)` の個数は、2 つの境の記号が同じ段を 1 回で下り、違う段から下は 2 本を同じループで下ります。記号が同じ段では両方の数が打ち消し合うので、数えずに済みます。2 MB 以上の配列は、Linux では 2 MB 境界に取って huge page を頼みます。

方式は procon-judge の yosupo-range-kth-smallest などで比べて選びました。x64-gcc で 1 段で 2 bit ずつ分ける 4 分の形 (wm4) と比べると、個数を数える問題で 2 割から 3 割速く出ました。aoj-2674 で 26.4 ms と 33.5 ms、aoj-2426 で 58.4 ms と 73.7 ms、yosupo-static-range-count-distinct で 94.8 ms と 120.2 ms です。`[lo, hi)` を 1 回で下りる形にすると、aoj-2674 は 24.7 ms、aoj-2426 は 48.4 ms になりました。k 番目は、yosupo-range-kth-smallest ($N, Q \le 2 \times 10^5$) の Xeon 8573C で 40.4 ms ほどで、4 分の形の 49.1 ms より 2 割弱速く出ました。値の種類が少なく段の数が減らない列では、4 分の形より遅くなります (値が 16 種類以下なら段の数は同じで、dense_small_a で 5.8 ms と 4.4 ms)。記録は algo-notes の notes/wavelet-matrix.md にあります。
