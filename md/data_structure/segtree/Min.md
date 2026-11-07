---
title: 最小・最大 (monoid::Min、monoid::Max)
documentation_of: ../../../neo/data_structure/segtree/Min.hpp
---

最小のモノイド `monoid::Min<T>`、最大のモノイド `monoid::Max<T>` と、それを渡したときの [SegmentTree](../SegmentTree) の中身です。口は SegmentTree と同じです。

```cpp
SegmentTree<monoid::Min<long long>> seg(a);
seg.set(i, x);                  // a[i] = x
seg.mul(i, x);                  // a[i] = min(a[i], x)
long long m= seg.prod(l, r);    // min(a[l], ..., a[r - 1])
int r= seg.max_right(l, [&](long long v) { return v >= k; });
```

`T` が 64 bit 以下の符号付き整数なら、子を 16 個持つ木になります。段 h ≥ 1 の塊は、子 16 個ぶんの最小を並べて持ちます ([monoid::AddMin](AddMin) の木から加算の量を除いた形)。1 点の代入は葉を書き換え、段ごとに子の塊の 16 個の最小を AVX2 で取り直し、親の値が変わらなくなった段で止めます。`mul` は値がよくなる向きにしか動かないので、子を読み直さずに上へ min を取るだけで済みます。区間の最小は、l と r - 1 が同じ子にある間は下り、分かれたら左右の道を下りて、全部掛かるレーンの最小を比べます。それ以外の `T` は二分木の汎用の形になります。
