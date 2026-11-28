---
title: 永続配列
documentation_of: ../../neo/base/ds/range/PersistentArray.hpp
---

M 分木の永続配列です。`set` が新しい版を返し、版の写しは根 1 個です。1 つの節に M 個の子を持たせ、添字 i を M 進法で下ります。M を大きくすると木が浅くなり、節の写しが M 個の子の分だけ重くなります。

```cpp
template <class T, int M = 8> class PersistentArray;
PersistentArray();
PersistentArray(int n, const T& v = T());
PersistentArray(const std::vector<T>& a);
int size() const;
T get(int i) const;
PersistentArray set(int i, const T& x) const;
static void reserve(std::size_t n);
static bool pool_empty();
static void gc(std::vector<PersistentArray>& vs);
```

値の無い節は `T()` を返します。`T` は既定構築できることが必要です。union-find の親の配列には `M = 64` を使います。

```cpp
PersistentArray<int, 64> p(n, -1);
auto q = p.set(3, 7);
```

手元の値は algo-notes の `notes/persistent-data-structures.md` にあります。
