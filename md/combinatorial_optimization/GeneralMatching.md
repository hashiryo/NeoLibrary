---
title: 一般グラフの最大マッチング
documentation_of: ../../neo/combinatorial_optimization/GeneralMatching.hpp
---

二部とは限らないグラフの最大マッチングを求めます。

```cpp
class GeneralMatching;
GeneralMatching(int n, const std::vector<std::array<int, 2>>& edges);
int size() const;
int mate(int v) const;
std::vector<int> matching() const;
```

頂点 n 個と、無向辺 {u, v} の列を受けます。多重辺と自己ループがあってもかまいません (自己ループは使いません)。辺は、渡した列の添字で扱います。

マッチングは構築の中で求めます。`size()` は最大マッチングの辺の数、`mate(v)` は頂点 v の相手 (いなければ -1) です。`matching()` は、最大マッチングに使う辺の番号を昇順に返します。多重辺はどれか 1 本を返します。

二部グラフと分かっているなら、`bipartite_coloring` で塗り分けて `BipartiteMatching` に渡すほうが 2 倍ほど速くなります。

## 中身

次の順に解きます。

1. 初期のマッチングを作る。疎なグラフは Karp-Sipser で、空いている隣が 1 つだけになった頂点を先に組ませる。隣の数の平均が 16 以上の密なグラフは、頂点の番号の順に最初の空いている隣と組ませる貪欲で作る。密なグラフは貪欲でほぼ最大になり、Karp-Sipser が隣を全部なめる分を省ける。
2. 空いている頂点すべてを根にした交互森を、1 本の待ち行列で同時に育てる。花は Edmonds の方法で縮める。違う木がつながれば増加路なので流し、その 2 本の木はその回の残りで使わない。1 本も流せない回が来れば最大なので終わる。
3. 2 の森を 2 回回しても終わらなければ、二部の二重被覆の最大マッチングを、今のマッチングの各辺を両向きに使った形から push-relabel で求める。
4. 二重被覆のマッチングは元のグラフの道と閉路に分かれるので、隣り合う 2 頂点ずつ組ませて丸め、2 の交互森で最大まで仕上げる。

二重被覆は、頂点 v を左の v と右の v に分け、辺 {u, v} を左 u と右 v、左 v と右 u の 2 本にした二部グラフです。増加路の短いグラフは、2 の森だけで終わります。長い増加路が残るグラフでは、奇閉路が少なければ丸めた時点でほぼ最大になり、長い増加路をたどる探索がほとんど要りません。

procon-judge の自作の self-flow-general-matching (20 万頂点、30 万辺の乱択、3 正則、三角格子、増加路の長い形など) で書き比べました。EPYC 7763 の x64-gcc で、最も遅いケースが 180 ms でした。同じ回に、今の Library の general_matching は 3515 ms でした。頂点 500 以下の yosupo-general-matching では、初期のマッチングを密さで選ぶ前は、Karp-Sipser の分で Library の探索を写した gabow の 2.2〜3.4 倍かかっていました (どれも 2 ms 以下)。EPYC 7763 の x64-gcc で、最も遅いケースが 1.21 ms (前の版、procon-judge da5309a2 の回) から 0.76 ms (run 44319692428) になりました。同じ回の gabow は 0.49 ms です。記録は algo-notes の notes/network-flow.md にあります。
