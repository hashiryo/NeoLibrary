---
title: 区間代入 (monoid::AssignSum、monoid::AssignMin、monoid::AssignMax)
documentation_of: ../../../neo/data_structure/segtree/Assign.hpp
---

区間代入の 3 つの構造体と、それを渡したときの [SegmentTree](../SegmentTree) の中身です。`monoid::AssignSum<T>` は区間和、`monoid::AssignMin<T>` は区間最小、`monoid::AssignMax<T>` は区間最大を取ります。口は SegmentTree と同じで、`apply(l, r, x)` が区間の代入です。

```cpp
SegmentTree<monoid::AssignMin<long long>> seg(a);
seg.apply(l, r, x);              // a[l..r) = x
long long m= seg.prod(l, r);     // min(a[l], ..., a[r - 1])
```

`T` が 64 bit 以下の符号付き整数なら、子を 16 個持つ遅延セグ木になります。段 h ≥ 1 の塊は、子 16 個ぶんの値 (和か最小か最大) と、まだ配っていない代入を並べて持ちます。代入が無いことは i64 の最小で表すので、代入する値に i64 の最小は使えません。区間の代入は、根から両端の道を配り、下から全部掛かるレーンの値と代入を AVX2 でまとめて書き、両端の道の値を子の 16 個から取り直します。代入は上の段ほど新しく、子の部分木全体を同じ値にします。そのため区間の取得と 1 点の取得は、道の途中で代入の残った子に当たれば、その値 (和なら値と個数の積) を返して下りずに済みます。`max_right` と `min_left` は区間の積を使った二分探索です。それ以外の `T` は二分木の汎用の形になります。

procon-judge に自作の 2 問を置いて測りました。self-segtree-range-assign-range-min と self-segtree-range-assign-range-sum で、どちらも $N = Q = 5 \times 10^5$ です。手元の M2 では、どちらも今の Library の SegmentTree の 2.3 倍の速さでした。
