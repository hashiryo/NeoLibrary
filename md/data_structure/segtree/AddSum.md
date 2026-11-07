---
title: 区間加算と区間和 (monoid::AddSum)
documentation_of: ../../../neo/data_structure/segtree/AddSum.hpp
---

区間加算と区間和の `monoid::AddSum<T>` と、それを渡したときの [SegmentTree](../SegmentTree) の中身です。口は SegmentTree と同じで、`apply(l, r, x)` が区間の加算です。

```cpp
SegmentTree<monoid::AddSum<long long>> seg(a);
seg.apply(l, r, x);              // a[l..r) += x
long long s= seg.prod(l, r);     // a[l] + ... + a[r - 1]
```

`T` が 64 bit 以下の整数なら、区間加算は、差分の配列 $d$ の 2 点の加算 ($d_l \mathrel{+}= x$、$d_r \mathrel{-}= x$) と見ます。中身は、[monoid::Sum](Sum) と同じ 8 分木の接頭辞和を 2 本持つ形です。$a_0 + \cdots + a_{p-1}$ は $S_p + p B_1 - B_2$ で出します。$S_p$ は初めの値の接頭辞和、$B_1 = \sum_{j<p} d_j$、$B_2 = \sum_{j<p} j d_j$ で、BIT を 2 本使う定番の形と同じです。区間の和は、2 つの接頭辞の 4 本の和を同じループで読みます。区間の加算は 2 本の木に 2 点ずつ足します。値は u64 の桁あふれのまま計算し、差を `T` に戻すので、`T` の範囲に収まる和は正しく出ます。`max_right` と `min_left` は区間の和を使った二分探索です。それ以外の `T` は二分木の汎用の形になります。

区間加算と区間和は、遅延セグ木で書くより BIT を 2 本使う形がずっと速い問題です。procon-judge の self-segtree-range-add-range-sum ($N = Q = 5 \times 10^5$) では、手元の M2 で 5 倍の差がありました。
