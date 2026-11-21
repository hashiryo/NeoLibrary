---
title: 一般グラフの重み付きマッチング
documentation_of: ../../neo/algo/optimize/matching/weighted_matching.hpp
---

二部とは限らないグラフで、重みの和が最大のマッチングと、重みの和が最小の完全マッチングを求めます。

```cpp
template <class W> struct WeightedMatching { W weight; std::vector<int> mate; };
template <class W> WeightedMatching<W> max_weight_matching(int n, const std::vector<std::tuple<int, int, W>>& edges);
template <class W> std::optional<WeightedMatching<W>> min_weight_perfect_matching(int n, const std::vector<std::tuple<int, int, W>>& edges);
```

頂点 n 個と、重み付きの無向辺 {u, v, w} の列を受けます。`W` は整数の型です。多重辺と自己ループがあってもかまいません (自己ループは使いません)。返す `weight` は重みの和、`mate[v]` は頂点 v の相手 (いなければ -1) です。重みの和が `W` に収まる必要があります。

`max_weight_matching` は、重みの和が最大のマッチングを求めます。完全でなくてかまいません。重みが 0 以下の辺は使いません。重みの最大を W として、$8W < 2^{62}$ である必要があります。

`min_weight_perfect_matching` は、重みの和が最小の完全マッチングを求めます。完全マッチングが無ければ `nullopt` を返します。重みは負でもかまいません。重みの最大と最小の差を D として、$8 ((n / 2) D + 1) < 2^{62}$ である必要があります。

重み最大の完全マッチングと重み最小のマッチングは、重みの符号を反転して渡してください。頂点の数の 2 乗の大きさの行列を持つので、頂点は数千までが目安です。計算量は $O(n^3)$ です。

## 中身

よく知られた $O(n^3)$ の主双対の花の方法 (頂点 2n の行列で花を持つ形) に、交互森を増加路ごとに作り直さない手を足したものです。増加路を流したら、その 2 本の木の印だけを消し、ほかの木の印と slack はそのまま使い続けます。双対は動かした量の和 T を使って遅延で持ち、slack も縮約費用に T の倍数を足した値で持ちます。そのため、双対を動かす手間は T を足すだけで済みます。木を壊すと、壊した木の頂点を指す slack が古くなります。古い値は本当の値以下なので、双対の動かし方を決めるときに最小を与えたものだけを確かめて付け直します。

元の頂点についての走査 (偶の頂点の行、slack の付け直し、双対の動かし方の最小、花の行の計算など) は AVX2 で回します。重みの最大 W が $8W < 2^{30}$ なら値を int32 で持って 8 列ずつ、そうでなければ i64 で 4 列ずつ回します。

重み最小の完全マッチングは、重みを C - w に替えた重み最大のマッチングで解きます。C は、重みの最大と最小を hi と lo として $C = (n / 2)(\mathrm{hi} - \mathrm{lo}) + \mathrm{lo} + 1$ とします。こうすると、どの完全マッチングも、どの完全でないマッチングより重くなります。

procon-judge の yosupo-general-weighted-matching ($N \le 500$) で書き比べました。EPYC 9V74 の x64-gcc で、最も遅いケースが 35.4 ms でした。同じ回に、今の Library の WeightedMatching は 267 ms でした。重み最小の完全マッチングは、atcoder-arc080-d で確かめています。記録は algo-notes の notes/network-flow.md にあります。
