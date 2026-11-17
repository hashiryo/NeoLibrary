---
title: 無向グラフの頂点の部分集合ごとの数え上げ
documentation_of: ../../neo/set_power_series/graph.hpp
---

頂点 $n$ 個の無向グラフ (多重辺と自己ループも可) について、頂点の部分集合 $S$ ごとに、$S$ が誘導する部分グラフの辺の部分集合などを数えます。辺は `sps::Edges` (`std::vector<std::pair<int, int>>`) で渡します。今の Library の UndirectedGraphSetPowerSeries のメンバ関数を、グラフを受ける関数に分けたものです。

```cpp
sps::Edges es = {{0, 1}, {1, 2}};
auto c = sps::connected_graph<Z>(n, es);      // 頂点の集合ごとの、連結な辺の部分集合の数
auto b = sps::biconnected_graph<Z>(n, es);    // 点双連結
auto p = sps::chromatic_polynomial<Z>(n, es); // 色多項式の係数
auto x = sps::only_biconnected(c);            // 連結なものの重みから、点双連結な塊の重みへ
```

| 関数 | 返すもの |
| --- | --- |
| `edge_num(n, es)` | 頂点の集合の中の辺の数 (自己ループは数えない) |
| `connected_component_num(n, es)` | 連結成分の数 |
| `cycle_space_rank(n, es)`、`odd_deg_num(n, es)`、`selfloop_num(n, es)` | 閉路空間の次元、次数が奇数の頂点の数、自己ループの数 |
| `graph<T>(n, es)` | $2^{\text{辺の数}}$ |
| `connected_graph<T>(n, es)` | 連結な辺の部分集合の数 |
| `cycle_space_size<T>(n, es)`、`eulerian_graph<T>(n, es)` | $2^{\text{閉路空間の次元}}$、連結なオイラーグラフの数 |
| `connected_bipartite_graph<T>(n, es)` | 連結な二部グラフの数 (T は 2 と $2^{\text{辺の数}}$ で割れること) |
| `tree<T>(n, es)`、`forest<T>(n, es)` | 全域木の数、全域森の数 |
| `cycle_graph<T>(n, es)` | 頂点の集合をちょうど通る閉路の数 (T は 2 で割れること) |
| `biconnected_graph<T>(n, es)`、`two_edge_connected_graph<T>(n, es)` | 点双連結、二辺連結な辺の部分集合の数 |
| `cactus_graph<T>(n, es)` | サボテンの数 |
| `acyclic_orientations<T>(n, es)` | 非巡回な向き付けの数 |
| `colorings_using_exactly_k_colors_num<T>(n, es)` | 全体をちょうど k 色で塗る数 (色の名前を区別しない、k = 0..n) |
| `chromatic_polynomial<T>(n, es)` | 色多項式の係数 |
| `tutte_polynomial<T>(n, es, x, y)` | Tutte 多項式の値 |
| `all_matching<T>(n, es)`、`k_matching<T>(n, es)` | マッチングの数、辺 k 本のマッチングの数 (k = 0..n/2) |
| `only_biconnected(f)`、`articulation_union(f)` | 連結なものの重みと、点双連結な塊の重みを互いに変換する |
| `only_2edge_connected(f, n, es)`、`bridge_union(f, n, es)` | 連結なものの重みと、二辺連結な塊の重みを互いに変換する |

完全マッチングの数 (辺の重みの積の和) は [hafnian](hafnian.md) を使います。集合冪級数の演算 (log、exp、部分集合畳み込み、逆元、$f^k/k!$ の表) は `sps` の関数を呼ぶので、法が奇数の `ZMod` なら速い核を通ります。
