---
title: 区間加算と区間最小・最大 (monoid::AddMin、monoid::AddMax)
documentation_of: ../../../neo/data_structure/segtree/AddMin.hpp
---

区間加算と区間最小の `monoid::AddMin<T>`、区間加算と区間最大の `monoid::AddMax<T>` と、それを渡したときの [SegmentTree](../SegmentTree) の中身です。口は SegmentTree と同じで、`apply(l, r, x)` が区間の加算です。

```cpp
SegmentTree<monoid::AddMin<long long>> seg(a);
seg.apply(l, r, x);              // a[l..r) += x
long long m= seg.prod(l, r);     // min(a[l], ..., a[r - 1])
int r= seg.max_right(l, [&](long long v) { return v >= k; });
```

`T` が 64 bit 以下の符号付き整数なら、子を 16 個持つ、作用を子へ配らない木になります。段 h ≥ 1 の塊は、子の部分木全体に足した量 t と、子の最小に t を足した s を子 16 個ぶん並べて持ちます。加算は最小と可換なので配らずに済みます。区間の加算は、下から全部掛かるレーンの s と t に AVX2 でまとめて足し、両端の道の s を取り直します。区間の最小は、l と r - 1 が同じ子にある間は t を足しながら下ります。分かれたら左右の道を下りて、全部掛かるレーンの最小にそこまでの t を足して比べます。値は i64 で持ちます。`max_right` と `min_left` は、端から塊を広げて最小を集め、条件が崩れる子へ下ります。それ以外の `T` は二分木の汎用の形になります。

procon-judge の yosupo-range-add-range-min ($N = Q = 5 \times 10^5$) で測ると、同じ中身の wide3 は最大 21.97 ms でした。今の Library の SegmentTree の 10.4 分の 1 です (EPYC 9V45 の x64-gcc)。
