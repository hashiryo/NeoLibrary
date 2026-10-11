---
title: 強連結成分分解
documentation_of: ../../neo/graph/strongly_connected_components.hpp
---

有向グラフを強連結成分に分けます。

```cpp
struct StronglyConnectedComponents { int count; std::vector<int> comp; };
StronglyConnectedComponents strongly_connected_components(int n, const std::vector<std::array<int, 2>>& edges);
```

頂点は 0 から n - 1 の番号で、`edges` は有向辺 (a, b) の列です。自己ループや多重辺があってもかまいません。返す `count` は成分の数、`comp[v]` は頂点 v の成分の番号 (0 以上 `count` 未満) です。番号は成分を縮めた DAG のトポロジカル順に振り、どの辺 (a, b) でも `comp[a] <= comp[b]` になります。頂点と辺の数の和に比例する時間で求めます。

## 中身

DFS の前に、入次数 0 の頂点を待ち行列で剥がします。剥がした頂点はどれも 1 頂点の成分で、剥がした順がトポロジカル順の先頭に来ます。剥がされない頂点から剥がした頂点への辺は無いので、残りの頂点は、剥がした頂点を見ないまま path-based (Gabow) の反復の DFS で分けます。DFS は読み込みが 1 本の鎖につながりますが、剥がしは待ち行列の頂点ごとに独立なので、読み込みを並べて出せます。辺の両端を一様に選んだ平均の出次数 1 のグラフでは、97% 以上の頂点が剥がれます。

鍵、最初の隣接、残りの隣接の連結リストの先頭は、頂点ごとの別々の 4 byte の配列に置きます。DFS が頂点に下りる 1 歩は、最初の隣接の読み込み 1 段で済みます。隣接は辺の列を 1 回なめるだけで組み、同じ走査で入次数も数えます。作業領域は、2 MB 以上なら 2 MB 境界に置いて huge page を頼みます (Linux のみ)。

procon-judge の yosupo-scc ($N, M \le 5 \times 10^5$) で書き比べました。速い持ち方は CPU で入れ替わりましたが、この形はどの CPU でも首位から 0.5 ms 以内でした。x64-gcc の最大ケースは、EPYC 9V74 で 7.9 ms、Xeon Platinum 8370C で 15.6 ms です。今の Library の StronglyConnectedComponents (Kosaraju) の 3.8 倍から 5.1 倍の速さにあたります。記録は algo-notes の notes/graph-basics.md にあります。
