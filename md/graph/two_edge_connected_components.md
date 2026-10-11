---
title: 二辺連結成分分解
documentation_of: ../../neo/graph/two_edge_connected_components.hpp
---

無向グラフを二辺連結成分に分けます。

```cpp
struct TwoEdgeConnectedComponents { int count; std::vector<int> comp; };
TwoEdgeConnectedComponents two_edge_connected_components(int n, const std::vector<std::array<int, 2>>& edges);
```

頂点は 0 から n - 1 の番号で、`edges` は無向辺 {a, b} の列です。自己ループや多重辺があってもかまいません。返す `count` は成分の数、`comp[v]` は頂点 v の成分の番号 (0 以上 `count` 未満) です。成分にどの番号を振るかは決めていません。頂点と辺の数の和に比例する時間で求めます。

## 中身

DFS の前に、次数 1 以下の頂点を待ち行列で剥がします。剥がした頂点の辺は橋なので、剥がした頂点はどれも 1 頂点の成分です。頂点ごとに次数と、残っている隣接の番号の XOR を持てば、次数が 1 になった頂点の残りの隣接は XOR の値そのものなので、隣接を組まずに辺の列だけで剥がせます。待ち行列に積むかどうかは分岐にせず、待ち行列の先の頂点を先読みします。

剥がしたあとの 2-core の頂点に 0 から番号を振り直し、2-core の辺だけで CSR を組んで、path-based (Gabow) の反復の DFS で分けます。DFS は頂点ごとの記録 {CSR の始まり, 終わり, 最初の隣接, 2 番目の隣接} で回します。次数 2 の頂点に下りたら、親でない方の隣接 (2 つの隣接と親の XOR) だけを見るので、閉路の頂点では、DFS で下りる 1 歩は記録の読み込み 1 段で済みます。親への辺は 1 本だけ飛ばし、2 本目からの平行な辺は後ろ向きの辺として扱います。CSR を組む走査では辺の端を先読みし、端が剥がれた辺は、辺の番号で回す 8 つの番兵の頂点に書いて捨てます。作業領域は、2 MB 以上なら 2 MB 境界に置いて huge page を頼みます (Linux のみ)。

procon-judge の yosupo-two-edge-connected-components ($N, M \le 2 \times 10^5$) で書き比べました。x64-gcc の最大ケースは、EPYC 7763 で 7.45 ms、EPYC 9V45 で 4.33 ms です。今の Library の IncrementalBridgeConnectivity の 1.26 倍から 1.55 倍の速さにあたります。記録は algo-notes の notes/graph-basics.md にあります。
