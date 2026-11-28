---
title: オンラインの動的連結性
documentation_of: ../../neo/base/ds/tree/OnlineDynamicConnectivity.hpp
---

辺の追加と削除、連結性、連結成分の積をオンラインで答える動的連結性です。段ごとに Euler Tour Tree を持ち、森の辺と非木辺を段に振り分けます (Holm–de Lichtenberg–Thorup)。今の Library の `OnlineDynamicConnectivity` と `EulerTourTree` の核を写しました。

```cpp
template <class M = void, std::size_t NODE_SIZE = 4'000'000> class OnlineDynamicConnectivity {
public:
 explicit OnlineDynamicConnectivity(int n);
 void link(int x, int y);                  // 辺を張る (単純グラフ)
 void cut(int x, int y);                   // 生きている辺を消す
 bool connected(int x, int y);
 int size(int x);                          // x の連結成分の頂点数
 const T& get(int x); void set(int x, T val);
 T prod(int x);                            // x の連結成分の積
 void apply(int x, E e);                   // 連結成分への作用
};
```

`M` はモノイド (`using T`、`static T ti()`、`static T op(T, T)`) か作用 (`using T`、`using E`、`mp`、`cp`) を持つ型です。`M= void` のときは連結性と辺の切り貼りだけを使えます。段の ETT は型ごとに 1 つの置き場 (`NODE_SIZE` 個) を共有するので、同時に 2 つ以上の `OnlineDynamicConnectivity` を作らないでください。頂点数と操作数から `NODE_SIZE` を見積もります。

クエリを先に読めるなら、`OfflineDynamicForest` の方がずっと速いです。このクラスは、クエリを順に処理する必要があるとき (前の答えで次の入力を決める、インタラクティブなど) に使います。

## 中身

各段の ETT は、辺の両端、部分木の頂点数、高い段の辺、非木辺の印を 64 bit の flag に詰めます。辺を消したときは、段を上から下へ下がりながら、小さい側の成分の辺を列挙します。同じ段でつながっている非木辺は 1 段上へ移し、切れている非木辺を見つけたら、それを段 0 からその段までの木の辺にして、代替辺にします。

`M= void` のときは節の値を持たないので、`NODE_SIZE` 個の節が 24 byte になります。

今の Library の提出 (`mylib/data_structure/OnlineDynamicConnectivity.hpp`) と同じ作りです。Library Checker の dynamic_graph_vertex_add_component_sum では、オフラインの `OfflineDynamicForest` が 59 ms、このクラスの提出は CI で測っています。記録は algo-notes の notes/dynamic-tree.md にあります。
