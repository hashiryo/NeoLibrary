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

辺の列が先に分かっているときは、`EulerTourTree(n, edges, a)` で O(N) に組めます。各成分の根から反復 DFS で巡回列を作り、区間の中央を根にした釣り合った木を組みます。

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

`link(u, v)` は u と v が別の連結成分にあるときだけ呼んでください。同じ成分の 2 頂点をつなぐと巡回列が二重になり、以後の積と `cut` が壊れます (assert で確かめています)。辺の列からの構築も森 (閉路の無いグラフ) を前提にします。

## 中身

木を、頂点の節と辺ごとの 2 つの向きの節を並べた巡回列で持ち、大きさで釣り合いを取る木 (WBT) で分けます。優先度は持たず、併合は大きい側の木へ降ります。`link(u, v)` は u の列、u から v へ向かう節、v の列、v から u へ向かう節の順につなぎます。向きは、p から v へ向かう節から前方に進むと v の側、v から p へ向かう節の手前で終わるように決めています。

大きい側へ降りるだけでは釣り合いが崩れます。降りなかった側の子と継いだ側の大きさを比べ、継いだ側が 4 倍を超えたら、継いだ側の根を持ち上げる回転 1 回で直します。直しを入れないと、鎖を 0-1、1-2、… と末尾に継ぎ足していくだけで木が一本の道になり、高さが N になります (yuki-828 の e_hand が TLE しました)。

`prod_subtree` は、v から p へ向かう節と p から v へ向かう節の位置を根まで登って出し、その間を区間として、木を組み替えずに降りて積を求めます。v から p へ向かう節が先に来るときは、巡回列の折り返しで区間が 2 つに分かれます。辺の節は頂点の数 0 として扱うので、区間には v の側の頂点だけが入ります。`apply_subtree` は同じ区間に遅延を置きます。

`cut(u, v)` は 2 つの節の位置を出し、4 回の分割で u の側と v の側に分けて、u の側だけ継ぎ直します。巡回列なので、切る前にどちらの節が先に来るかを気にしなくてよいです。分割は木を切るだけで、高さは元の木より大きくなりません。

区間の積と作用は、splay 版のように木を組み替えて取らないので、そのぶんクエリが軽いです。逆に分割と併合は splay 版より少し重いです。

値です。procon-judge の提出 (`submissions/lib-neo-build.cpp`) は dynamic-tree-subtree-add-subtree-sum の 18 ケースがすべて AC で、x64-gcc の最大 660 ms、合計 5709 ms、メモリ 84.5 MB です (AMD EPYC 7763。同じ CPU モデルの splay 版の提出は 774 ms、6600 ms、84.3 MB)。Xeon 6975P-C では WBT が 500 ms、splay が 716 ms、素朴な WBT (回転の直しなし) が 582 ms でした。yuki-828 は 46 ケースがすべて AC、x64-gcc で最大 361 ms です。手元の M2 Max の almost_line_00 は 600〜620 ms です (splay 版の下書き `_shared/dyntree/ett.hpp` が 645〜688 ms)。同問題の Library Checker の首位は 192 ms (オフラインのクラスタ木) です。部分木の加算と和は、クエリを先読みできるならオフラインの構造の方が速いです。記録は algo-notes の notes/dynamic-tree.md にあります。
