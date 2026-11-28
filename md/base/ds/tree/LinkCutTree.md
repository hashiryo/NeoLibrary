---
title: Link-Cut Tree
documentation_of: ../../neo/base/ds/tree/LinkCutTree.hpp
---

森の辺を切り貼りしながら、道の積と道への作用を答える Link-Cut Tree です。道の積は u から v への順で、可換でないモノイドでもよいです。辺 (v, p) について p を親とした v の部分木の積を答える `SubtreeLinkCutTree` も同じヘッダにあります。

```cpp
template <class M> class LinkCutTree {
public:
 explicit LinkCutTree(int n);
 int size() const;
 T get(int v);                              // 1 点の値
 void set(int v, T x);                      // 1 点の代入
 void mul(int v, const T& d);               // 1 点の値を op で動かす
 void init_value(int v, T x);               // 初期値 (splay しない)
 void init_parent(int v, int p);            // 初期木の親を張る (splay しない)
 void evert(int v);                         // v を根にする
 void link(int u, int v), cut(int u, int v);
 bool connected(int u, int v);
 int root(int v), parent(int v), lca(int u, int v);
 T prod(int u, int v);                      // u から v への道の積 (両端を含む)
 void apply(int u, int v, const E& e);      // 道への作用
};

template <class M> class SubtreeLinkCutTree {
public:
 explicit SubtreeLinkCutTree(int n);
 T get(int v); void set(int v, T x); void add(int v, T d);
 void evert(int v); void link(int u, int v), cut(int u, int v);
 bool connected(int u, int v);
 int root(int v), parent(int v), lca(int u, int v);
 T prod_tree(int v);                        // v の木の全体の積
 T prod_subtree(int v, int p);              // p を親とした v の部分木の積
};
```

`M` は SegmentTree と同じモノイドです (`using T`、`static T ti()`、`static T op(T, T)`)。`using commute= void` があれば可換として扱い、逆向きの積を持ちません。作用は `using E`、`mp(T&, E)`、`cp(E&, E)` で受けます。作用はモノイドの `op` に分配する必要があります。和のモノイドに「各頂点に足す」作用を載せるには、要素数を持つモノイドにするか、長さを受ける作用にしてください。

`SubtreeLinkCutTree` は可換群 (`commute` と `static T inv(T)`) を受けます。部分木は子の順が決まらないので、可換でない `op` では使えません。`prod_subtree(v, p)` は v と p が辺で結ばれていることを前提にします。

## 使い方

初期の木は、値を `init_value` で置き、親を `init_parent` で張ると O(N) で組めます。すべての辺を仮想辺として張るので、`link` を N-1 回呼ぶより速いです。

```cpp
LinkCutTree<monoid::Sum<long long>> lct(n);
for (int v = 0; v < n; ++v) lct.init_value(v, a[v]);
for (int v = 1; v < n; ++v) lct.init_parent(v, par[v]);
...
lct.cut(u, v), lct.link(w, x);
lct.mul(p, x);
long long s = lct.prod(u, v);
```

`prod` の u と v は同じ木にあることが前提です。`lca` は今の根から見た lca で、`root` と `parent` も今の根に依存します。`evert` で根を変えられます。

## 中身

節は配列に置き、親を符号つきの int で持ちます。正は splay の親、負は道の親、0 はなしです。頂点の番号は 1 始まりで、0 番を番兵にします。splay は下から回転し、回転の直前に祖父母・親・自分の 3 節だけ push します。

回転では、x の部分木の区間が p の前の部分木と同じになるので、x の集約を p から写してから p だけ求め直します。集約の計算を回転ごとに 2 回しない分だけ速いです。`access` は標準の形で、最後に x を splay します。道の積は `evert(u)`、`access(v)` のあとの `sum(v)` です。v の右の子は空なので、和は u から v への道の順になります。

`SubtreeLinkCutTree` は、splay の部分木の和 (道の上の頂点と、仮想の子の部分木) と、仮想の子の部分木の和 `vir` を持ちます。`access` で右の子を仮想へ移すとき `vir` に足し、仮の道の子を仮想から外すとき `vir` から引きます。部分木の積は `evert(v)`、`access(p)` のあとの `op(val[v], vir[v])` です。

procon-judge の提出 (`problems/_shared/dyntree/lct.hpp` と同じ実装) の値です。run 46102092510 の x64-gcc (Xeon 8573C) で、入出力込みの最大と合計です。

| 問題 | 提出 | 最大 | 合計 | Library Checker の首位 |
| --- | --- | --- | --- | --- |
| dynamic-tree-vertex-add-path-sum | `lct` | 257 ms | 1248 ms | 174 ms |
| dynamic-tree-vertex-set-path-composite | `lct` | 315 ms | 1543 ms | 141 ms |
| dynamic-tree-vertex-add-subtree-sum | `lct_sub` | 285 ms | 2553 ms | 115 ms (オフラインの森) |

部分木の和は、オフラインでよいなら期限付きの最大全域森の方が速いです。記録は algo-notes の notes/dynamic-tree.md にあります。
