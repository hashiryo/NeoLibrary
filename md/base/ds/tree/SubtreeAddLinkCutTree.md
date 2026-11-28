---
title: 部分木への作用付き Link-Cut Tree
documentation_of: ../../neo/base/ds/tree/SubtreeAddLinkCutTree.hpp
---

森の辺を切り貼りしながら、辺 (v, p) について p を親とした v の部分木の積を答え、その部分木に作用を当てる Link-Cut Tree です。

```cpp
template <class M>
class SubtreeAddLinkCutTree {
 public:
  explicit SubtreeAddLinkCutTree(int n);
  void set(int v, T x);
  T get(int v);
  void evert(int v);
  void link(int u, int v);
  void cut(int u, int v);
  bool connected(int u, int v);
  int root(int v);
  int parent(int v);
  int lca(int u, int v);
  T prod_tree(int v);
  T prod_subtree(int v, int p);
  void apply_subtree(int v, int p, const E& e);
};
```

`M` は可換群 (`using T`、`ti()`、`op(T, T)`、`using commute`、`inv(T)`) で、作用は長さを受ける `mp(T&, E, int)` と `cp(E&, E)`、それに作用の逆元 `inve(E)` を持ちます。部分木に作用を当てるときに、仮想の子へ移す遅延を打ち消すために逆元を使うためです。

節は道の積 `sum` と頂点数 `sz` に加えて、仮想の子の積 `vir` と頂点数 `vir_sz`、道の遅延 `lz`、仮想の子への遅延 `vir_lz` を持ちます。access で右の子を仮想へ移すときは `vir` に足して `-vir_lz` を当て、仮想から道へ取るときは `+vir_lz` を当ててから `vir` から引きます。すべての操作は O(log N) です。

辺 (v, p) が存在することを前提にした 2 引数の形です。`prod_subtree(v, p)` と `apply_subtree(v, p, e)` は、その辺を親子として見たときの v の側の部分木を扱います。
