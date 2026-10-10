---
title: 最大流
documentation_of: ../../neo/flow/MaxFlow.hpp
---

有向グラフの s から t への最大流を求めます。値だけでなく、最小カット、辺ごとの流量、道への分解も取り出せます。

```cpp
template <class Cap> class MaxFlow;
explicit MaxFlow(int n = 0);
int add_vertex();
int add_edge(int from, int to, Cap cap, Cap rev_cap = 0);
Cap flow(int s, int t);
Cap flow(int s, int t, Cap limit);
std::vector<bool> min_cut(int s);
struct Edge { int from, to; Cap cap, rev_cap, flow; };
Edge get_edge(int i);
std::vector<Edge> edges();
Cap change_cap(int i, Cap cap, int s, int t);
Cap change_cap(int i, Cap cap, Cap rev_cap, int s, int t);
std::vector<std::pair<Cap, std::vector<int>>> decompose(int s, int t);
int size() const;
```

頂点は 0 から n - 1 の番号で、`add_vertex()` で 1 つ足せます。`add_edge` は、from から to へ容量 cap、to から from へ容量 rev_cap の辺を足し、足した順の番号 (0 から) を返します。無向辺は `rev_cap` を `cap` と同じにします。自己ループも足せますが、流れません。名前と意味は、ACL の mf_graph にあるものはそれに合わせています。

`flow(s, t)` は、今の流れから続けて s から t へ流せるだけ流し、増えた量を返します。`flow(s, t, limit)` は、増やす量を limit までにします。何度呼んでもよく、呼ぶたびに終点が違ってもかまいません。flow のあとに辺を足すこともできます。

`min_cut(s)` は、残余グラフで s から届く頂点を true にした列を返します。直前に上限なしの `flow(s, t)` を呼んだあとなら、true の側が s を含む最小カットになります。`get_edge(i)` と `edges()` の `flow` は、from から to へ流れていれば正、to から from へ流れていれば負です。

`change_cap(i, cap, s, t)` は、辺 i の容量を cap に変えます (`rev_cap` も変える版もあります)。今の流れは s から t への流れとし、流量が新しい容量に収まらない分は、まず辺の両端の間で迂回させ、迂回できない分は流れを戻して直します。返すのは、s から t への流量が減った分です。容量を下げたあとの流れは、変更後のグラフの最大流になっています。容量を上げても流量は増えないので、続けて `flow(s, t)` を呼びます。

`decompose(s, t)` は、今の流れ (s から t への流れとする) を s から t への道に分け、道ごとの流量と、通る辺の番号の列を返します。逆向きの容量のある辺は、to から from へ通ることもあります。流れの中の循環は道に含めません。

`Cap` は整数の型で、辺ごとの `cap + rev_cap` と、s から t への流量の答えが `Cap` に収まる必要があります。s から出る辺の容量の和は収まらなくてかまわないので、無限大の代わりに大きな容量の辺を何本張っても溢れません。

## 中身

最初の `flow` は push-relabel で流します。高さが最大の頂点から余りを流し、同じ高さの中では先に入れたものから取り出します。gap と global relabel も使います。前流の段 (t へ流し、届かなかった余りは途中に残す) で値を返し、余りを s へ戻す段は、辺ごとの流量か最小カットを初めて聞かれたときまで遅らせます。値だけの使い方では、この段の費用がかかりません。2 回目からの `flow` は Dinic で流します。容量を少し変えて少しだけ流し直す使い方では、呼ぶたびにグラフ全体を BFS でなめる push-relabel より速いためです。

方式は procon-judge の問題で比べて選びました。self-flow-maxflow には、性質の違う 7 つのグラフの族が並んでいます。同じ方式で値だけを求める版 (容量を int にしたもの) は、最大のケースが EPYC 9V45 の x64-gcc で 264 ms でした。今の Library の push-relabel だと 1057 ms、Dinic だと 2.6 秒です。密なグラフの loj-127 では 80 ms ほどでした。容量を変えながら流し直す aoj-2803、aoj-2313、aoj-2835 では、Dinic で流す今の Library と同じくらいの速さです。記録は algo-notes の notes/network-flow.md にあります。
