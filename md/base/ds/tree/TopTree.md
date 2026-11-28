---
title: 自己調整 top tree
documentation_of: ../../neo/base/ds/tree/TopTree.hpp
---

森の辺を切り貼りしながら、部分木と道の積・作用を O(log N) で扱う自己調整 top tree です。道と部分木の両方を同時に動的に扱えるので、dynamic tree DP の土台になります。

```cpp
template <class M>
class TopTree {
 public:
  explicit TopTree(const std::vector<T>& a);
  void link(int u, int v);
  void cut(int u, int v);
  void evert(int v);
  void apply_subtree(int v, int p, const E& e);
  T prod_subtree(int v, int p);
  void apply_path(int u, int v, const E& e);
  T prod_path(int u, int v);
  void set(int v, T x);
  T get(int v);
};
```

`M` は可換モノイド (`using T`、`ti()`、`op(T, T)`、`using commute`) で、作用は長さを受ける `mp(T&, E, int)` と `cp(E&, E)` を持ちます。作用の逆元は要りません (構造を変える前に遅延を降ろして子から再計算するため)。値の逆元も要りません。

節は道の splay (`l`, `r`) と rake の splay (`h`) を 1 つの型で持ち、道の頂点数と積 (`len`, `sum`)、部分木の頂点数と積 (`cap`, `sub`)、道への遅延 `lzy`、部分木への遅延 `hzy` を持ちます。expose で x を根のクラスタにし、`x->h` が x の子の部分木 (rake) になります。節は `std::deque` に置くので、増えても番地が変わりません。

部分木の 2 引数形は、辺 (v, p) を親子として見たときの v の側を扱います。道は u から v の順で、値が可換なので向きは問いません。
