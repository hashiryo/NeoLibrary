---
title: オフラインの動的連結性
documentation_of: ../../neo/base/ds/tree/OfflineDynamicForest.hpp
---

全クエリを先に積み、辺の削除時刻を重みにした最大全域森を 1 本持って再生する、オフラインの動的連結性です。木は親ポインタと部分木の大きさ・部分木和だけで、splay を使いません。2/3 の大きさ平衡 (anti-monopoly) を保ち、期限の短い辺ほど上に来ます。削除は `owner[時刻]` で切る節を引くだけで、代替辺を探しません。

```cpp
template <class G> class OfflineDynamicForest {
public:
 using T = typename G::T;
 OfflineDynamicForest(int n, const std::vector<T>& a);
 void link(int u, int v);                  // 辺を張る (単純グラフ)
 void cut(int u, int v);                   // 生きている辺を消す
 void add(int v, T x);                     // 頂点 v の値を x だけ動かす
 void query(int v);                        // v の成分の和を記録する
 std::vector<T> solve();                   // 積んだ操作を再生し、query の答えを順に返す
};
```

`G` は可換群です (`using T`、`static T ti()`、`static T op(T, T)`、`static T inv(T)`、`using commute`)。和、法つきの和、xor などが載ります。差の取れないモノイド (min、max、非可換の積) は載りません。

オンラインでクエリを処理する必要がある問題には使えません。`link`、`cut`、`add`、`query` を積んでから `solve` を 1 回呼びます。`solve` の中で、辺の寿命の前計算 (どの link がどの cut で消えるか) と、森の再生をします。

## 中身

辺の重みは削除時刻です。森は最大全域森なので、閉路ができる辺は、経路のいちばん早く消える辺と比べて、長く生きる方だけを残します。削除のときは、消える辺が森の経路の最大の重み (いちばん早い期限) になっているので、`owner[時刻]` から切る節を引いて切るだけで済みます。期限の順序を保つために、節を 2/3 の大きさ平衡で上へ回します。成分の和は、根の部分木和です。頂点への加算は、その頂点から根までの和に足します。

Library Checker の dynamic_graph_vertex_add_component_sum の首位 (nandhagk の 402761、57 ms) の核を写したものです。procon-judge の提出の CI (run 46103525567、x64-gcc、Xeon 6975P-C) は、最大 59 ms、合計 373 ms、メモリ 25 MB でした。記録は algo-notes の notes/dynamic-tree.md にあります。
