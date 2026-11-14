---
title: 区間 chmin、chmax、加算と区間和 (monoid::ChminChmaxAddSum)
documentation_of: ../../../neo/data_structure/segtree/ChminChmaxAddSum.hpp
---

区間 chmin、区間 chmax、区間加算と、区間の和、最大、最小を取る `monoid::ChminChmaxAddSum<V>` と、それを渡したときの [SegmentTree](../SegmentTree) の中身です (Segment Tree Beats)。`V` は 64 bit 以下の符号付き整数です。作用 `{lb, ub, add}` は $x \mapsto \min(\max(x, lb), ub) + add$ で、`chmin(x)`、`chmax(x)`、`add(x)` で作ります。値は `{sum, max, min}` で、`sum` は long long です。

```cpp
using M= monoid::ChminChmaxAddSum<long long>;
SegmentTree<M> seg(n, [&](int i) { return a[i]; });
seg.apply(l, r, M::chmin(x));       // a[i] = min(a[i], x)
seg.apply(l, r, M::chmax(x));       // a[i] = max(a[i], x)
seg.apply(l, r, M::add(x));         // a[i] += x
long long s= seg.prod(l, r).sum;    // .max と .min もある
```

口は SegmentTree と同じですが、`mul` はありません。`set(i, x)` と `get(i)` は 1 要素の値で、`V` から暗に作れます。`SegmentTree(n)` は 0 で埋めます。`max_right` と `min_left` は `{sum, max, min}` についての条件で使います。値は絶対値 $2^{61}$ 未満で使います ($2^{61}$ を番兵にするため)。

中身はこの 3 つの操作に特化した Beats です。節は最大と 2 番目の最大、最小と 2 番目の最小、和、まだ配っていない加算、最大と最小の個数を 64 byte に並べます。根から再帰で下り、chmin は最大だけを、chmax は最小だけを動かせば済む節で止めます。chmin と chmax の作用は節に持たず、親の最大と最小で表します。配るときは加算を先に当て、子の最大が親の最大を超えていれば子にも chmin を当てます (chmax も同じ)。
