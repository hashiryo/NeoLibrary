---
title: 支配木
documentation_of: ../../neo/graph/dominator_tree.hpp
---

有向グラフで、頂点 s を根とする支配木を求めます。頂点 d が頂点 v を支配するとは、s から v へのどの道も d を通ることです。

```cpp
std::vector<int> dominator_tree(int n, int s, const std::vector<std::array<int, 2>>& edges);
```

頂点は 0 から n - 1 の番号で、`edges` は有向辺 (a, b) の列です。自己ループや多重辺があってもかまいません。返す列の v 番目は頂点 v の直接支配点 (v を支配する頂点のうち v に最も近いもの) で、v が s なら s、s から届かなければ -1 です。頂点と辺の数の和を m として、ほぼ $O(m \log n)$ の時間で求めます。

## 中身

SEMI-NCA で求めます。s から DFS で届く頂点に行きがけ順の番号を振り、番号の大きい順に、逆向きの辺の始点から森を登って準支配点を求めます。直接支配点は、番号の小さい順に、DFS の木の親から準支配点より上にある直接支配点をたどって決めます。Lengauer–Tarjan のようなバケットは要りません。森を登る eval は、道を 1 回登りながら各点の親を祖父に付け替える経路半減で書いています。

DFS の隣接は、頂点ごとに最初の出辺の行き先と残りの出辺の連結リストを別々の配列に置き、辺の列を 1 回なめて組みます。DFS が下りる 1 歩は読み込み 1 段で済みます。DFS のあとは、届く頂点どうしの辺だけで逆向きの隣接を番号の上で組み直します。作業領域は、2 MB 以上なら 2 MB 境界に置いて huge page を頼みます (Linux のみ)。

procon-judge の yosupo-dominator-tree ($N, M \le 2 \times 10^5$) で書き比べました。x64-gcc の最大ケースは EPYC 9V45 で 6.53 ms で、Library Checker の参照解と同じ方法 (経路圧縮だけの Lengauer–Tarjan) の 2.45 倍の速さにあたります。今の Library には実装がありません。記録は algo-notes の notes/graph-basics.md にあります。
