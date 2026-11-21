---
title: 有向グラフの数え上げの核 (内部)
documentation_of: ../../neo/base/core/set_fps/digraph_simd.hpp
---

[graph.hpp](graph.md) の `strongly_connected_digraph` と `acyclic_digraph` が呼ぶ内部の核です。直接は使いません。どちらも、頂点の集合 $S$ ごとの値を

$$h[S] = c[S] - \sum_{T \subsetneq S} h[T]\, X(S \setminus T)\, F(T \to S \setminus T) \quad (h[\emptyset] = h_0)$$

で $S$ の小さい順に求める $3^n$ の和に帰着します。$F(T \to U)$ は $T$ から $U$ へ向かう辺の $(1 + w_e)$ の積、$X(U)$ は $U$ の中の辺 (自己ループを含む) の $(1 + w_e)$ の積です。強連結は $X$ を掛けて $c[S] = X(S)$、$h_0 = 0$ とし、DAG は辺を逆向きにしたグラフで $X = 1$、$c = 0$、$h_0 = 1$ とします。

## 汎用の道

どの可換環でも動く形です。$S$ ごとに、$S$ の部分集合 $T$ の $F(T \to S \setminus T)$ を $2^{|S|}$ 個まとめて作ります。$S$ の bit を 1 つずつ決め、$T$ に入れるなら決め終えた $U$ 側への辺の積を、$U$ に入れるなら決め終えた $T$ 側からの辺の積を掛けます。表は頂点ごとの $\prod_{b \in X} q(v, b)$ と $\prod_{a \in X} q(a, v)$ で、割り算を使いません。

## AVX2 の道

法が奇数の `ZMod` ($p < 2^{30}$) で $n \ge 4$ のときの形です。値は Montgomery の表現で持ちます。

- 添字の下の 3 bit をレーンに置く。上の bit は、$T$ の上の部分 $T_h$ と $U = S \setminus T$ の上の部分 $U_h$ の組で回す。
- $X(U) F(T \to U)$ は、組ごとの値 $X(U_h) F(T_h \to U_h)$、$T_h$ と $U_l$ の表、$U_h$ とレーンの表、下の 3 bit だけの $8 \times 8$ の表の積に分かれる。$T_h$ の側と $U_h$ の側で、レーンの向きの 8 本の列を先に作る。
- 1 組の寄与は、長さ 8 の内積 (積を u64 にためて 1 回で還元) と、組ごとの値の掛け算になる。組ごとの値は、汎用の道と同じく $S$ の上の bit を 1 つずつ決める形で、$S_h$ ごとに $2^{|S_h|}$ 個まとめて作る。
- 同じ塊の中の項 ($U_h = \emptyset$) は、塊ごとにレーンの小さい順に引く。

procon-judge の luogu-P11714 (主旋律) で書き比べた scc_simd (数だけの形) を、辺の重みと DAG に広げたものです。記録は algo-notes の notes/set-power-series.md にあります。
