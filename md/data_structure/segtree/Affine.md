---
title: 区間の 1 次関数 (monoid::AffineSum、monoid::Affine)
documentation_of: ../../../neo/data_structure/segtree/Affine.hpp
---

区間の 1 次関数と区間和の `monoid::AffineSum<Z>` (遅延)、区間の 1 次関数と 1 点取得の `monoid::Affine<Z>` (双対) と、それを渡したときの [SegmentTree](../SegmentTree) の中身です。作用 `{b, c}` は $x \mapsto bx + c$ で、`b = 0` なら代入、`b = 1` なら加算になります。

```cpp
using Z= ZMod<998244353>;
SegmentTree<monoid::AffineSum<Z>> seg(a);
seg.apply(l, r, {b, c});   // a[l..r) = b a[i] + c
Z s= seg.prod(l, r);       // a[l] + ... + a[r - 1]
SegmentTree<monoid::Affine<Z>> dual(a);
dual.apply(l, r, {b, c});
Z v= dual.get(i);
```

`Z` が法を定数にした `ZMod<MOD>` なら、どちらも子を 16 個持つ木になります。段 h ≥ 1 の塊は、子 16 個ぶんの和 (遅延のときだけ) と、まだ配っていない作用を並べて持ちます。区間の作用は、根から両端の道を配り、下から全部掛かるレーンへ AVX2 でまとめて当て、両端の道の和を取り直します。操作の初めに、その操作が下りる 2 本の道の塊を先読みします。区間の和は書き込みをせず、l と r - 1 が同じ子にある間は作用を合成しながら下り、分かれたら左右の道を下りて、全部掛かるレーンの和にそこまでに合成した作用を当てます。値は剰余をそのままの値で $2p$ 未満に持ち、掛け算は Shoup の方法で 8 レーンずつします。1 点の取得は、葉の値に下の段から順に作用を当てるだけです。`max_right` と `min_left` は `prod` を使った二分探索です。それ以外の `Z` (実行時の法の `ZMod<0>` など) は二分木の汎用の形になります。

procon-judge の yosupo-range-affine-range-sum ($N = Q = 5 \times 10^5$) で、同じ中身の wide_lazy3 は最大 67.21 ms でした。今の Library の SegmentTree の 7.4 分の 1 です。yosupo-range-affine-point-get では、wide_dual2_pf が最大 35.63 ms (5.0 分の 1) でした。どちらも EPYC 7763 の x64-gcc の値です。
