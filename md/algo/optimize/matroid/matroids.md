---
title: マトロイドの型
documentation_of: ../../neo/algo/optimize/matroid/matroids.hpp
---

[マトロイドの交差](./matroid_intersection.md) に渡すマトロイドの型です。どの型も、要素を $0, 1, \dots$ の番号で表します。

```cpp
GraphicMatroid(int n, const std::vector<std::array<int, 2>>& edges);
PartitionMatroid(const std::vector<int>& part, const std::vector<int>& cap);
XorMatroid(const std::vector<unsigned long long>& v);
TransversalMatroid(int l, int r, const std::vector<std::array<int, 2>>& edges);
CographicMatroid(int n, const std::vector<std::array<int, 2>>& edges);
```

| 型 | 要素 | 独立な集合 | 持つ操作 |
| --- | --- | --- | --- |
| GraphicMatroid | n 頂点のグラフの辺 | 閉路を含まない (森) | can_exchange、circuit |
| PartitionMatroid | 組 part[x] に属する要素 (-1 なら制限なし) | 組 c の要素が cap[c] 個以下 | can_exchange、circuit |
| XorMatroid | 64 bit のベクトル | xor で線形独立 | can_exchange、circuit |
| TransversalMatroid | 二部グラフの左の頂点 | マッチングで覆える | circuit |
| CographicMatroid | n 頂点のグラフの辺 | 取り除いても連結成分の数が変わらない | can_exchange |

辺は多重辺と自己ループがあってもかまいません。TransversalMatroid の辺は {左の頂点, 右の頂点} です。

## 中身

- GraphicMatroid: `build` で森を根付き木にし、深さ、親への辺、オイラー順の入りと出を持ちます。`can_exchange(y, x)` は、y の下側の頂点の部分木に x の両端のちょうど一方が入るかを $O(1)$ で見ます。`circuit(x)` は、x の両端から深いほうを親へたどります。
- PartitionMatroid: 閉路は、x と同じ組で I に入っている要素の全部です。
- XorMatroid: `build` で I の基底を掃き出し、各行が I のどの要素の xor かを 64 bit の印で持ちます。x を基底で消して 0 になれば従属で、残った印が閉路です。消した結果は `build` ごとに覚えます。
- TransversalMatroid: `build` は前のマッチングのうち I に残る組を残し、I の空いた頂点から増加路を探して覆い直します。そのあと、空いた右の頂点から交互路を逆にたどり、空いた右の頂点へ届く右の頂点に印を付けます。`can_add(x)` は x の隣に印があるかで、`circuit(x)` は x から交互路をたどって届いた右の頂点の相手を列挙します。
- CographicMatroid: `build` は G から I を取り除いたグラフの DFS の森を作り、lowlink で橋を求めます。`can_add(x)` は x が橋でないことで、`can_exchange(y, x)` は、x が橋なら x の下側の頂点の部分木に y の両端のちょうど一方が入るかを見ます。
