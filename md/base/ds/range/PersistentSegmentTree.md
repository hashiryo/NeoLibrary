---
title: 永続セグメント木
documentation_of: ../../neo/base/ds/range/PersistentSegmentTree.hpp
---

モノイドを渡す永続セグメント木です。`set` と `mul` が新しい版を返し、版の写しは根 1 個です。1 回の更新で根から葉への道の節を写すので、O(log n) 個の節が増えます。

```cpp
template <class M> class PersistentSegmentTree;
PersistentSegmentTree(int n);                 // ti() で埋めた 1 版
PersistentSegmentTree(const std::vector<T>& a);
int size() const;
T get(int i) const;
T prod(int l, int r) const;
T all_prod() const;
PersistentSegmentTree set(int i, T x) const;  // 新しい版
PersistentSegmentTree mul(int i, T x) const;  // a[i] = op(a[i], x)
static void reserve(std::size_t n);
static bool pool_empty();
static void gc(std::vector<PersistentSegmentTree>& vs);
```

`M` の決まりは [SegmentTree](SegmentTree.md) と同じです。`using T`、`static T ti()`、`static T op(T a, T b)` を渡します。

```cpp
using PST = PersistentSegmentTree<monoid::Sum<i64>>;
PST v0(n);
PST v1 = v0.set(3, 5);
PST v2 = v1.mul(4, 2);
PST v3 = v1;              // 根を写すだけ (O(1))
i64 s = v2.prod(0, n);
```

## 2 版の差

`T` が整数 (数え上げ) のとき、次の 2 つが使えます。

- `kth(a, b, k)`: 版 b と版 a の差の k 番目 (0 始まり) の位置。
- `max_xor(a, b, x)`: 版 b と版 a の差の値のうち、x との xor が最大になるものの xor。

`prod_diff(a, b, l, r, sub)` は、版 b と版 a の `[l, r)` の差を、覆う節ごとに `sub(版 b の値, 版 a の値)` で求めて `op` でまとめます。`sub` は呼ぶ側が渡します。和なら `[](T b, T a) { return b - a; }`、 (個数, 和) なら両方の引き算です。

## 置き場

節は型ごとに 1 本の置き場に取り、0 番を空の節 (子が自分、値は `ti()`) にします。`reserve` で先に取ると、配列を作り直す分を省けます。`pool_empty()` が真を返したら、生きている版を渡して `gc` を呼んでください。`gc` は渡した版から届く節だけを前へ詰め、版どうしで共有する節は共有したまま残します。渡さなかった版は使えなくなります。

```cpp
vector<PST> ver(q + 1);
ver[i] = ver[j].set(k, x);
if (PST::pool_empty()) PST::gc(ver);
```

## 測った値

静的な range-kth-smallest、static-range-sum-with-upper-bound、rectangle-sum は wavelet matrix の方が速いです。版が要る問い合わせと、`max_xor` のように値の bit を使う問い合わせで使います。手元の値と自作問題の結果は algo-notes の `notes/persistent-data-structures.md` にあります。
