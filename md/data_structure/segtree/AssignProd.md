---
title: 区間代入と区間の積 (monoid::AssignProd)
documentation_of: ../../../neo/data_structure/segtree/AssignProd.hpp
---

どのモノイドにも使える区間代入の構造体 `monoid::AssignProd<M>` と、それを渡したときの [SegmentTree](../SegmentTree) の中身です。`M` は SegmentTree と同じモノイド (`using T`、`static T ti()`、`static T op(T a, T b)`) で、可換でなくてかまいません。口は SegmentTree と同じで、`apply(l, r, x)` が `T` の値 `x` の区間代入です。

```cpp
using Z= ZMod<998244353>;
using CM= monoid::Composite<Z>;                 // 1 次関数の合成 (Composite.hpp)
SegmentTree<monoid::AssignProd<CM>> seg(f);
seg.apply(l, r, {c, d});                         // f[l..r) = (x -> c x + d)
Z y= CM::eval(seg.prod(l, r), x);                // f[r - 1](...f[l](x)...)
```

整数の和、最小、最大なら、[Assign](Assign) の `monoid::AssignSum` などの方が速く動きます。

中身は二分木の遅延セグ木です。長さ $2^j$ の節に `x` を代入した値は $x^{2^j}$ (`x` を $2^j$ 個並べた積) なので、代入のたびに `x` から 2 乗を重ねた表を作り、節の作用には表の位置を持たせます。表は代入が覆う節のいちばん高い段までだけ作ります。配るときは子の段の冪を表から引くので、作用を当てるたびに冪を求め直しません。1 点の取得は、根から下りて初めに当たった代入の値を返します (根に近い代入ほど新しいため)。表は 1 回の代入ごとに段の数だけ伸びます。32 MB (節の数の 4 倍の個数より小さければそちら) を超えたら、残っている作用が引く分だけを新しい表へ詰め直します。残る分は節の数の 2 倍以下です。`max_right` と `min_left` は ACL と同じで、`check(ti())` が真、条件が区間の伸びに対して単調であることを前提にします。
