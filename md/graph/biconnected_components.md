---
title: 二重連結成分分解
documentation_of: ../../neo/graph/biconnected_components.hpp
---

無向グラフを二重連結成分 (関節点で分けた成分) に分けます。

```cpp
struct BiconnectedComponents { int count; std::vector<int> start, vertex; };
BiconnectedComponents biconnected_components(int n, const std::vector<std::array<int, 2>>& edges);
```

頂点は 0 から n - 1 の番号で、`edges` は無向辺 {a, b} の列です。自己ループは無いものとし、多重辺はあってもかまいません。返す `count` は成分の数で、成分 k の頂点は `vertex[start[k]]` から `vertex[start[k + 1] - 1]` です (`start` の長さは `count + 1`)。関節点は複数の成分に入り、辺を持たない頂点はそれだけで 1 つの成分になります。成分と頂点の並びは決めていません。頂点と辺の数の和に比例する時間で求めます。

## 中身

DFS の前に、次数 1 の頂点を待ち行列で剥がします。剥がした頂点 v と、そのとき残っていた唯一の隣接 u は、それだけで 1 つの成分 {v, u} です (辺 vu は橋)。頂点ごとに次数と、残っている隣接の番号の XOR を持てば、u は XOR の値そのものなので、隣接を組まずに辺の列だけで剥がせます。待ち行列に積むかどうかと成分を置くかどうかは分岐にせず、待ち行列の先の頂点を先読みします。

剥がしたあとの 2-core の頂点に 0 から番号を振り直し、2-core の辺だけで CSR を組んで、Hopcroft–Tarjan の反復の DFS で頂点のスタックから成分を切り出します。DFS は頂点ごとの記録 {CSR の始まり, 終わり, 最初の隣接, 2 番目の隣接} で回します。下りた直後は親への辺を番号の比較で飛ばし、次数 2 の頂点では親でない方の隣接 (2 つの隣接と親の XOR) だけを見ます。親への辺を飛ばしても、子 c から帰るときの判定 low[c] >= ord[親] は変わりません。CSR を組む走査では辺の端を先読みし、端が剥がれた辺は、辺の番号で回す 8 つの番兵の頂点に書いて捨てます。作業領域は、2 MB 以上なら 2 MB 境界に置いて huge page を頼みます (Linux のみ)。

procon-judge の yosupo-biconnected-components ($N, M \le 5 \times 10^5$) で書き比べました。x64-gcc の最大ケースは EPYC 7763 で 19.37 ms で、今の Library の block_cut_tree の 3.35 倍の速さにあたります。最大は 47 万頂点の閉路のケースで決まっていて、DFS が 1 頂点ごとに読み込みを待つ時間がそのまま出ています。記録は algo-notes の notes/graph-basics.md にあります。
