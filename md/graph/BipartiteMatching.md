---
title: 二部マッチング
documentation_of: ../../neo/graph/BipartiteMatching.hpp
---

二部グラフの最大マッチングを求めます。最小点被覆、最大独立集合、最小辺被覆、Dulmage-Mendelsohn 分解も取り出せます。

```cpp
class BipartiteMatching;
BipartiteMatching(int l, int r, const std::vector<std::array<int, 2>>& edges);
BipartiteMatching(const std::vector<int>& color, const std::vector<std::array<int, 2>>& edges);
int size() const;
int mate(int v) const;
std::vector<int> matching() const;
std::vector<int> min_vertex_cover(std::vector<int> order = {});
std::vector<int> max_independent_set(std::vector<int> order = {});
std::vector<int> min_edge_cover() const;
std::pair<int, std::vector<int>> dulmage_mendelsohn();
int edge_kind(int i);
int vertex_kind(int v);
```

入口が 2 つあります。`BipartiteMatching(l, r, edges)` は、左の頂点 l 個と右の頂点 r 個と、辺 (左 a、右 b) の列を受けます。頂点は、左 a を a、右 b を l + b の番号で扱います。`BipartiteMatching(color, edges)` は、一般のグラフの辺の列と、2 色の塗り分け (`bipartite_coloring` で求める) を受けます。頂点は元の番号のまま扱い、色 0 の頂点を左、色 1 の頂点を右とします。辺はどちらの入口でも、渡した列の添字で扱います。多重辺があってもかまいません。

マッチングは構築の中で求めます。`size()` は最大マッチングの辺の数、`mate(v)` は頂点 v の相手 (いなければ -1) です。`matching()` は、最大マッチングに使う辺の番号を昇順に返します。

`min_vertex_cover(order)` と `max_independent_set(order)` の `order` は頂点の順列で、前の頂点ほど被覆 (独立集合) に入れます。order の順で辞書順最小の解になり、返す頂点も order の順に並びます。省くと番号の順にします。`min_edge_cover()` は、最小辺被覆に使う辺の番号を昇順に返します。辺の無い頂点があれば空を返します。

`dulmage_mendelsohn()` は、ブロックの数 K と、頂点ごとのブロックの番号を返します。番号はトポロジカル順で、左 a と右 b を結ぶ辺があれば、a の番号は b の番号以下です。ブロック 0 は空いた右の頂点から交互路で届く頂点を集めたもので、その中の右の頂点はどれも空けられます。ブロック K - 1 は空いた左の頂点から届く頂点で、その中の左の頂点はどれも空けられます。この 2 つは空のこともあります。1 から K - 2 は残りの頂点の強連結成分で、どれも左と右の頂点の数が等しくなります。

`edge_kind(i)` は、辺 i がどの最大マッチングにも入らないなら 0、入るものと入らないものがあるなら 1、どの最大マッチングにも入るなら 2 を返します。`vertex_kind(v)` は、頂点 v がどの最大マッチングでも空いているなら 0 (辺が無い)、空くものと空かないものがあるなら 1、どの最大マッチングでも相手がいるなら 2 を返します。Dulmage-Mendelsohn 分解は、被覆、独立集合、分類のどれかを初めて聞かれたときに求めるので、マッチングだけを使うときには費用がかかりません。

## 中身

初期のマッチングは Karp-Sipser で作ります。空いている隣が 1 つだけになった頂点を先に組ませ、そういう頂点が無くなったら残りから 1 組を選ぶことを繰り返します。そのあと push-relabel (double push) で仕上げます。右の頂点に、空いている右の頂点までの距離の下界を持たせ、空いている左の頂点を FIFO で取り出して、下界が最小の隣と組ませます。(L + R) / 2 回組ませるごとに、空いている右の頂点からの逆向きの BFS で、下界を正確な距離に直します。

procon-judge の yosupo-bipartitematching ($L, R \le 10^5$、$M \le 2 \times 10^5$) で書き比べました。EPYC 7763 の x64-gcc で、最も遅いケースが 56 ms でした。同じ回に、今の Library の bipartite_matching は 356 ms でした。記録は algo-notes の notes/network-flow.md にあります。
