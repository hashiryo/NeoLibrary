---
title: 二部グラフの塗り分け
documentation_of: ../../neo/graph/bipartite_coloring.hpp
---

無向グラフを 2 色に塗り分けます。

```cpp
std::vector<int> bipartite_coloring(int n, const std::vector<std::array<int, 2>>& edges);
```

頂点は 0 から n - 1 の番号で、`edges` は無向辺の列です。返す列の v 番目は頂点 v の色 (0 か 1) で、どの辺も両端の色が違います。連結成分ごとに、番号の最も小さい頂点を 0 にします。二部グラフでなければ、空の列を返します。自己ループがあるときも空です。頂点と辺の数の和に比例する時間で求めます。

`BipartiteMatching` に一般のグラフを渡すときの塗り分けにも使います。
