---
title: 永続遅延セグメント木
documentation_of: ../../neo/base/ds/range/PersistentLazySegmentTree.hpp
---

モノイドと作用を渡す永続遅延セグメント木です。`apply` と `copy_range` が新しい版を返し、版の写しは根 1 個です。

```cpp
template <class M> class PersistentLazySegmentTree;
PersistentLazySegmentTree(int n);
PersistentLazySegmentTree(const std::vector<T>& a);
int size() const;
T get(int i) const;
T prod(int l, int r) const;
T all_prod() const;
PersistentLazySegmentTree apply(int l, int r, E x) const;
PersistentLazySegmentTree copy_range(int l, int r, const PersistentLazySegmentTree& src) const;
static void reserve(std::size_t n);
static bool pool_empty();
static void gc(std::vector<PersistentLazySegmentTree>& vs);
```

`M` の決まりは [SegmentTree](SegmentTree.md) と同じです。`using T`、`using E`、`ti`、`op`、`mp`、`cp` を渡します。`mp` が bool を返す Beats の作用は受けません。

```cpp
using PLST = PersistentLazySegmentTree<monoid::AffineSum<ZMod<998244353>>>;
PLST v0(a);
PLST v1 = v0.apply(l, r, {b, c});
PLST v2 = v1.copy_range(l, r, v0);  // v1 の [l, r) を v0 から写す
auto s = v2.prod(l, r);
```

`copy_range` は、写す側と写される側の両方の道の節を写してから遅延を配ります。元の版を変えず、元の祖先の遅延を写した範囲に残さないためです。`apply` も `copy_range` も O(log n) 個の節が増えます。`prod` と `get` は上から遅延を合成して当てるので、版を写しません。

手元の値は algo-notes の `notes/persistent-data-structures.md` にあります。Library Checker の persistent_range_affine_range_sum を通しました。
