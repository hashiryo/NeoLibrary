---
title: 部分木加算と部分木和 (オフライン)
documentation_of: ../../neo/base/ds/tree/AntiMonopolyTree.hpp
---

頂点 1 個 1 節のクラスタ木で、成分全体への加算と成分の和を答える、オフラインの構造です。2/3 の大きさ平衡 (anti-monopoly) を保ち、splay を使いません。部分木は辺 (v, p) を切った成分で、イベント列の再生で切る・足す・聞く・戻します。

```cpp
template <class T> class FoldableAntiMonopolyTree {
public:
 explicit FoldableAntiMonopolyTree(const std::vector<T>& a);
 bool link(int u, int v, int w);           // 辺を重み w で入れる。閉路なら最大重みの辺を外す
 bool erase(int u, int v, int w);          // 重み w の辺を外す
 void cut(int v);                          // v の親辺を切る
 void add_vertex(int v, T x);              // 頂点 v の値を x だけ動かす
 void add_all(int v, T x);                 // v の成分全体に x を足す
 T fold(int v);                            // v の成分の和
 std::pair<int, int> path_max(int u, int v);
};

template <class T> class OfflineSubtree {
public:
 explicit OfflineSubtree(int n);
 void link(int u, int v);                  // 辺を張る (木)
 void cut(int u, int v);                   // 生きている辺を消す
 void add_all(int v, T x);                 // 成分全体への加算を記録する
 void fold(int v);                         // 成分の和を記録する
 std::vector<T> solve(const std::vector<T>& a);  // 再生し、fold の答えを順に返す
};
```

`T` は整数のような、`+`、`-`、`int` 倍、`+=`、`-=` が使える型です。`long long` と `ZMod` が載ります。

`OfflineSubtree` は、辺 (v, p) を切って p を親とした v の部分木を、切った成分として扱います。型 1 のクエリは「切る・成分全体への加算・つなぐ」、型 2 のクエリは「切る・成分の和・つなぐ」の 3 イベントに分けて積みます。辺の重みは削除時刻の符号反転で、初期の辺は `link` として積みます。`solve` の中で、辺の寿命の前計算と再生をします。

Library Checker の dynamic_tree_subtree_add_subtree_sum の首位 (anon123 の 401546、192 ms) の核を写したものです。procon-judge の提出の CI (run 46103525567、x64-gcc、Xeon 6975P-C) は、最大 241 ms、合計 1646 ms、メモリ 39 MB でした。記録は algo-notes の notes/dynamic-tree.md にあります。
