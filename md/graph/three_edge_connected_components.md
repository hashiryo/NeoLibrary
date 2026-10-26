---
title: 三辺連結成分分解
documentation_of: ../../neo/graph/three_edge_connected_components.hpp
---

無向グラフを三辺連結成分に分けます。2 頂点が同じ成分に入るのは、その 2 頂点を結ぶ辺を共有しない道が 3 本とれるときです。

```cpp
struct ThreeEdgeConnectedComponents { int count; std::vector<int> comp; };
ThreeEdgeConnectedComponents three_edge_connected_components(int n, const std::vector<std::array<int, 2>>& edges);
```

頂点は 0 から n - 1 の番号で、`edges` は無向辺 {a, b} の列です。自己ループや多重辺があってもかまいません。返す `count` は成分の数、`comp[v]` は頂点 v の成分の番号 (0 以上 `count` 未満) です。成分にどの番号を振るかは決めていません。頂点と辺の数の和に比例する時間で求めます。

## 中身

DFS の前に、次数 1 以下の頂点を待ち行列で剥がします。剥がした頂点の辺は橋なので、剥がした頂点はどれも 1 頂点の成分です。頂点ごとに次数と、残っている隣接の番号の XOR を持てば、次数が 1 になった頂点の残りの隣接は XOR の値そのものなので、隣接を組まずに辺の列だけで剥がせます。

剥がしたあとの 2-core を、Tsin の吸収の形の 1 回の DFS で分けます。2-core の次数 2 の頂点はどれも 1 頂点の成分で (2 本の辺を切れば離れる)、両隣をつなぐ 1 本の辺に置き換えてもほかの頂点の成分は変わりません。そこで DFS には次数 3 以上の頂点だけを載せ、次数 2 の頂点が続く鎖に当たったら、鎖の向こう端まで歩いて、向こう端を隣接とみなします。鎖の両端の頂点には両端の頂点の XOR を印として残し、反対の端から当たったときはそれを引きます。鎖の 1 歩は記録 1 つの読み込みで済み、次数 3 以上の頂点を含まない閉路には触りません。DFS の頂点ごとの状態は 32 byte の構造体 1 つにまとめます。作業領域は、2 MB 以上なら 2 MB 境界に置いて huge page を頼みます (Linux のみ)。

procon-judge の yosupo-three-edge-connected-components ($N, M \le 2 \times 10^5$) で書き比べました。x64-gcc の最大ケースは、EPYC 9V74 で 6.66 ms、EPYC 9V45 で 5.32 ms です。Library Checker の参照解と同じ方法 (Tsin の吸収の形を std::vector の隣接で回すもの) の 2.4 倍の速さにあたります。今の Library には実装がありません。記録は algo-notes の notes/graph-basics.md にあります。
