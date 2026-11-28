---
title: 永続 Li Chao 木
documentation_of: ../../neo/algo/optimize/line/PersistentLiChaoTree.hpp
---

[LiChaoTree](LiChaoTree.md) の永続版です。`add` / `add_line` / `add_segment` が新しい版を返し、`query` はその版の最小 (MAX なら最大) を引きます。版の写しは根 1 個です。

```cpp
template <class T = long long, bool MAX = false, class F = LiChaoLine<T, MAX>> class PersistentLiChaoTree;
PersistentLiChaoTree(long long lo, long long hi);
PersistentLiChaoTree add(const F& f) const;
PersistentLiChaoTree add(const F& f, long long l, long long r) const;
PersistentLiChaoTree add_line(T a, T b) const;
PersistentLiChaoTree add_segment(long long l, long long r, T a, T b) const;
T query(long long x) const;
static void reserve(std::size_t n);
static bool pool_empty();
static void gc(std::vector<PersistentLiChaoTree>& vs);
```

`F` の決まりは [LiChaoTree](LiChaoTree.md) と同じです。x の範囲 `[lo, hi)` を渡します。全体に足す関数は O(log)、線分は O(log^2) 節を写します。`query` は版を写しません。

```cpp
PersistentLiChaoTree<long long> t0(lo, hi);
auto t1 = t0.add_line(a, b);
auto t2 = t1.add_segment(l, r, a, b);
long long y = t2.query(x);
```

版を残さないなら [LiChaoTree](LiChaoTree.md) の方が速いです。手元では時間が 1.4 倍、メモリが 7 倍でした。版を残す使い方でだけ意味があります。手元の値は algo-notes の `notes/persistent-data-structures.md` にあります。
