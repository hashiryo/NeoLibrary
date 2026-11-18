---
title: グラフの頂点の部分集合ごとの数え上げ
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

## 有向グラフ

辺を `(a, b)` で a → b と読む有向グラフ (多重辺と自己ループも可) の数え上げもあります。辺ごとの重み `w` (長さ m) を渡すと、条件を満たす辺の部分集合 $E'$ についての $\prod_{e \in E'} w_e$ の和を返します。`w` を省くと数 (重みが全部 1) です。

```cpp
sps::Edges es = {{0, 1}, {1, 2}, {2, 0}};
auto f = sps::strongly_connected_digraph<Z>(n, es);  // 頂点の集合ごとの、強連結な辺の部分集合の数
cout << f.back() << '\n';                            // 主旋律 (UOJ 37)
std::vector<Z> w = {Z(2), Z(3), Z(5)};
auto d = sps::acyclic_digraph<Z>(n, es, w);          // DAG になる辺の部分集合の重み
```

| 関数 | 返すもの |
| --- | --- |
| `strongly_connected_digraph<T>(n, es[, w])` | 強連結になる辺の部分集合 ($[\emptyset]$ は 0) |
| `acyclic_digraph<T>(n, es[, w])` | DAG になる辺の部分集合 ($[\emptyset]$ は 1、自己ループは選べない) |

どちらも $O(3^n)$ です。源になる強連結成分 (入ってくる辺の無い成分) の頂点の集合で包除すると、$g = 1 - \exp(-f)$ は

$$g[S] = X(S) - \sum_{\emptyset \ne T \subsetneq S} g[T]\, X(S \setminus T)\, F(T \to S \setminus T)$$

を満たします。$X(U)$ は $U$ の中の辺の $(1 + w_e)$ の積、$F(T \to U)$ は $T$ から $U$ へ向かう辺の $(1 + w_e)$ の積です。この和で $g$ を求め、強連結は $f = -\log(1 - g)$ で出します。DAG は源の頂点の集合で包除し、$b[S] = (-1)^{|S|} a[S]$ と置くと、辺を逆向きにしたグラフで同じ形の和になります。

法が奇数の `ZMod` で $n \ge 4$ なら、添字の下の 3 bit を AVX2 のレーンに置きます。$T$ と $U$ の上の部分の組ごとに、両側でレーンの向きに作った 8 本の列の内積を取り、組ごとの値 $X(U_h) F(T_h \to U_h)$ を掛けます。組ごとの値は、$S$ の上の部分ごとに $2^{|S_h|}$ 個まとめて作ります。$S$ の上の bit を 1 つずつ決め、両側の表の値を掛ける形で、割り算を使いません。同じ形の書き比べの提出 scc_simd は、luogu-P11714 (主旋律、$n = 15$) で 8.98 ms でした (procon-judge の run 44286626441、x64-gcc)。素直な 3^n の DP (68.58 ms) の 7.6 倍の速さです。
