---
title: 和 (monoid::Sum)
documentation_of: ../../../neo/data_structure/segtree/Sum.hpp
---

和のモノイド `monoid::Sum<T>` と、それを渡したときの [SegmentTree](../SegmentTree) の中身です。口は SegmentTree と同じで、`mul(i, x)` が 1 点の加算になります。

```cpp
SegmentTree<monoid::Sum<long long>> seg(a);
seg.mul(i, x);                 // a[i] += x
long long s= seg.prod(l, r);   // a[l] + ... + a[r - 1]
int r= seg.max_right(l, [&](long long v) { return v <= k; });
```

`T` が 64 bit 以下の整数なら、節を 8 個に広げた接頭辞和の木になります。段 h の表には、塊 (8 個ずつ) の中でその子より前にある子の和を持ちます。先頭から r 個の和は、段ごとに表を 1 つずつ読んで足すだけで済みます。1 点の加算は、段ごとに塊の後ろ側へ AVX2 でまとめて足します。塊の 8 個はキャッシュの 1 行に収まります。値は u64 の桁あふれのまま足して差を `T` に戻すので、`T` の範囲に収まる和は正しく出ます。`get(i)` は、i 番目だけの区間の和で求めます。元の値の配列を別に持つ形は、1 点の加算のたびにその配列も読み書きするので、x64 で 1 割遅くなりました。`max_right` と `min_left` は、段ごとに塊の 8 個の境目の接頭辞和を見て、条件が切り替わる子へ下ります。それ以外の `T` (浮動小数点数など) は二分木の汎用の形になります。

procon-judge の yosupo-point-add-range-sum ($N = Q = 5 \times 10^5$) で、同じ中身の wide8 は EPYC 9V74 の x64-gcc で最大 7.05 ms でした。今の Library の BinaryIndexedTree の 2.5 分の 1、SegmentTree の 5.3 分の 1 です。
