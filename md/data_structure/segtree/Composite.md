---
title: 1 次関数の合成 (monoid::Composite)
documentation_of: ../../../neo/data_structure/segtree/Composite.hpp
---

1 次関数の合成のモノイド `monoid::Composite<Z>` です。値 `{a, b}` は $x \mapsto ax + b$ で、`op(f, g)` は `f` を先に、`g` をあとに当てた関数です。`eval(f, x)` で $f(x)$ を返します。

```cpp
using Z= ZMod<998244353>;
using CM= monoid::Composite<Z>;
SegmentTree<CM> seg(f);
seg.set(i, {c, d});                  // f[i] = (x -> c x + d)
Z y= CM::eval(seg.prod(l, r), x);    // f[r - 1](...f[l](x)...)
```

[SegmentTree](../SegmentTree) に渡しても中身は二分木の汎用の形のままで、区間の積は覆う節を分岐なしで集めてから 2 本の鎖で畳みます。procon-judge の yosupo-point-set-range-composite では、この形が Montgomery で書いた書き比べの首位と並びました。区間代入と合わせるときは [AssignProd](AssignProd) に渡します。
