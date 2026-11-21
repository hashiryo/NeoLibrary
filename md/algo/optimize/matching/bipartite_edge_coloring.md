---
title: 二部グラフの辺彩色
documentation_of: ../../neo/algo/optimize/matching/bipartite_edge_coloring.hpp
---

二部グラフの辺を、最大次数と同じ数の色で、同じ頂点に接する辺の色が違うように塗ります。

```cpp
std::vector<int> bipartite_edge_coloring(int l, int r, const std::vector<std::array<int, 2>>& edges);
```

左の頂点 l 個、右の頂点 r 個と、辺 {左 a, 右 b} の列を受けます。多重辺があってもかまいません。最大次数を D として、辺ごとの色 (0 から D - 1) を返します。正則二部グラフなら、色ごとの辺の集まりがそれぞれ完全マッチングになります。

## 中身

多重辺を「辺と重複度」にまとめて持ち、D 正則の多重グラフを再帰で分けます。

1. 頂点を、次数の和が D 以下になるように左右それぞれ先頭から詰めて組にまとめ、組の間の辺を重複度付きの辺にまとめる。左右の組の数を揃え、足りない次数は組の間に辺を足して D 正則にする。
2. 次数 d が奇数なら、完全マッチングを 1 つ取って色を 1 つ塗り、d を 1 減らす。
3. d が偶数なら、重複度の半分ずつを両側へ送る。重複度が奇数の辺の残り 1 本ずつは、閉路に沿って交互に両側へ振り分ける (Euler 分割)。両側はどちらも d / 2 正則になる。
4. 辺の種類が 64 以下になったら、完全マッチングを取り、その上の重複度の最小だけ色をまとめて塗る (Birkhoff-von Neumann の分解)。

3 の閉路は、各頂点で辺を隣り合う 2 本ずつ組にし、右の端の組の相手、左の端の組の相手、と交互にたどって作ります。二部グラフなので長さは偶数になります。完全マッチングは `BipartiteMatching` で取ります。最後に、まとめた辺ごとに塗った色を元の辺へ配ります。

procon-judge の yosupo-bipartite-edge-coloring ($L, R, M \le 10^5$) で書き比べました。EPYC 9V45 の x64-gcc で、最も遅いケースが 13.38 ms でした。同じ回に、今の Library の bipartite_edge_coloring は 72.58 ms でした。記録は algo-notes の notes/network-flow.md にあります。
