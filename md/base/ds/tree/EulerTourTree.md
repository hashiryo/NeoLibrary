---
title: Euler Tour Tree
documentation_of: ../../neo/base/ds/tree/EulerTourTree.hpp
---

森の辺を切り貼りしながら、連結成分の積と、辺 (v, p) について p を親とした v の部分木の積を答える Euler Tour Tree です。部分木への作用も遅延で当てられます。

```cpp
template <class M> class EulerTourTree {
public:
 explicit EulerTourTree(int n);
 EulerTourTree(int n, const std::vector<std::array<int, 2>>& edges, const std::vector<T>& a);  // 辺の列から O(N)
 int size() const;                          // 頂点の数
 bool edge_exists(int u, int v) const;
 T get(int v); void set(int v, T x);
 void link(int u, int v), cut(int u, int v);
 bool connected(int u, int v);
 T prod_tree(int v);                        // v の連結成分の積
 int tree_size(int v);                      // v の連結成分の頂点数
 void apply_tree(int v, const E& e);        // 連結成分への作用
 T prod_subtree(int v, int p);              // p を親とした v の部分木の積
 void apply_subtree(int v, int p, const E& e);
};
```

`M` はモノイド (`using T`、`static T ti()`、`static T op(T, T)`) です。作用は `using E`、`cp(E&, E)` と、区間の頂点の数を受ける `mp(T&, E, int)` で受けます。辺の節は頂点として数えないので、区間の長さを受ける形が要ります。

`prod_subtree(v, p)` と `apply_subtree(v, p, e)` は、v と p が辺で結ばれていることを前提にします。p を親とした v の部分木は、p から v へ向かう節の直後から、v から p へ向かう節の直前までの区間です。

辺の列が先に分かっているときは、`EulerTourTree(n, edges, a)` で O(N) に組めます。各成分の根から反復 DFS で巡回列を作り、区間の中央を根にした平衡な splay 木を組みます。

## 使い方

```cpp
EulerTourTree<monoid::AffineSum<ZMod<998244353>>> ett(n);
for (int v = 0; v < n; ++v) ett.set(v, a[v]);
for (auto [u, v] : edges) ett.link(u, v);
...
ett.cut(u, v), ett.link(w, x);
ett.apply_subtree(v, p, {b, c});
auto s = ett.prod_subtree(v, p);
```

辺は (u, v) の値で覚えるので、`cut` は `link` したときと逆順でもかまいません。同じ 2 頂点の辺を何度も切って張り直せます。辺の節は free list で使い回すので、節の数は生きた辺の数で決まります。

## 中身

木を、頂点の節と辺ごとの 2 つの向きの節を並べた巡回列で持ち、splay 木で分けます。`link(u, v)` は u の列、u から v へ向かう節、v の列、v から u へ向かう節の順につなぎます。向きは、p から v へ向かう節から前方に進むと v の側、v から p へ向かう節の手前で終わるように決めています。

`cut` は、まず u から v へ向かう節の位置まで列を回してから、その節と v から u へ向かう節で列を 3 つに分けます。巡回列なので、切る前に列の線形化でどちらの節が先に来ているかを気にしなくてよいです。

`prod_subtree` は、p から v へ向かう節の位置まで列を回し、v から p へ向かう節を splay して、その左の部分木の積を返します。辺の節は頂点の数 0 として扱うので、左の部分木には v の側の頂点だけが入ります。`apply_subtree` は同じ区間に遅延を置きます。

splay の回転では、x の部分木の区間が p の前の部分木と同じになるので、x の集約を p から写してから p だけ求め直します。push は回転に参加する 3 節だけにします。

procon-judge の提出 (`problems/_shared/dyntree/ett.hpp` と同じ実装) の値です。run 46102092510 の x64-gcc (Xeon 8573C) で、dynamic-tree-subtree-add-subtree-sum の最大 817 ms、合計 7370 ms、メモリ 65 MB でした。同問題の Library Checker の首位は 192 ms (オフラインのクラスタ木) です。部分木の加算と和は、クエリを先読みできるならオフラインの構造の方が速いです。記録は algo-notes の notes/dynamic-tree.md にあります。
