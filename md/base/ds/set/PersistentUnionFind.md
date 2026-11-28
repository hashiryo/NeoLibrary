---
title: 永続 union-find
documentation_of: ../../neo/base/ds/set/PersistentUnionFind.hpp
---

永続 union-find です。親の配列を `PersistentArray<int, 64>` に持ち、`unite` が新しい版を返します。版の写しは根 1 個です。経路圧縮はせず、合併は大きさの大きいほうに寄せるので、`find` は O(log n) 段、各段の `get` は O(log_64 n) ですみます。

```cpp
struct PersistentUnionFind;
PersistentUnionFind();
PersistentUnionFind(int n);
int find(int x) const;
bool same(int a, int b) const;
int size(int x) const;
PersistentUnionFind unite(int a, int b) const;  // 新しい版
static void reserve(std::size_t n);
```

根は `p[x] < 0` で、その絶対値が大きさです。`unite` は同じ成分なら `*this` を返します。

```cpp
PersistentUnionFind uf0(n);
PersistentUnionFind uf1 = uf0.unite(a, b);
bool c = uf1.same(a, b);
```

手元の値は algo-notes の `notes/persistent-data-structures.md` にあります。Library Checker の persistent_unionfind を通しました。
