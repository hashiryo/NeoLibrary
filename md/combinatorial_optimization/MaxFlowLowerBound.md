---
title: 下限付きの最大流
documentation_of: ../../neo/combinatorial_optimization/MaxFlowLowerBound.hpp
---

流量に下限と上限のある網で、条件を満たす流れがあるかを判定し、s から t への流量の最大と最小を求めます。辺ごとの流量も取り出せます。

```cpp
template <class Cap> class MaxFlowLowerBound;
explicit MaxFlowLowerBound(int n = 0);
int add_vertex();
int add_edge(int from, int to, Cap lower, Cap upper);
void add_supply(int v, Cap x);
bool solve();
std::optional<Cap> max_flow(int s, int t);
std::optional<Cap> min_flow(int s, int t);
struct Edge { int from, to; Cap lower, upper, flow; };
Edge get_edge(int i);
std::vector<Edge> edges();
int size() const;
```

頂点は 0 から n - 1 の番号で、`add_vertex()` で 1 つ足せます。`add_edge` は from から to へ、流量が lower 以上 upper 以下の辺を足し、足した順の番号 (0 から) を返します。lower は負でもよく、負の流量は to から from へ流れる量を表します。無向辺は `add_edge(u, v, -c, c)` です。`add_supply(v, x)` は、頂点 v の湧き出しに x を足します。負の値なら吸い込みです。

`solve()` は、次の流れがあるかを返します。どの辺の流量も下限と上限の間にあり、どの頂点でも、出る量から入る量を引いたものが湧き出しに等しい流れです。湧き出しを足していなければ、循環があるかの判定になります。循環は、どの頂点でも入る量と出る量の等しい流れです。

`max_flow(s, t)` と `min_flow(s, t)` は、s と t だけはこの釣り合いを外し、s から t へ F だけ流す流れのうち、F の最大と最小を返します。s では出る量が湧き出しより F 多く、t では入る量が F 多い流れです。F は 0 以上に限ります。そういう流れが無ければ `nullopt` を返します。s から t への流れがあるかだけを知りたいときは、`max_flow(s, t).has_value()` で判定できます。t から s へ下限 0、上限が十分大きい辺を足して `solve()` を呼んでも同じで、こちらは続けて流す分の時間がかかりません。

解いたあとは、`get_edge(i)` と `edges()` で、直前の `solve`、`max_flow`、`min_flow` が求めた流れの辺ごとの流量を読めます。流れが無かったときや、そのあとに辺、頂点、湧き出しを足したときは、`flow` は下限になります。何度問い合わせてもよく、問い合わせの間に辺を足してもかまいません。

`Cap` は符号付きの整数の型です。辺ごとの upper - lower、下限の絶対値の和と湧き出しの絶対値の和、答えの流量は、どれも `Cap` で表せる大きさにしてください。

## 中身

辺ごとに、下限と上限の間で 0 にいちばん近い量を先に流しておき、残りの幅を容量にした辺を [最大流](MaxFlow.md) に張ります。先に流した分で頂点に出る過不足は、超頂点 S から過剰な頂点へ、不足の頂点から超頂点 T への辺で埋め、S から T へ全部流れるかを見ます。s と t があるときは、t から s へ容量が無限大の辺 (補助の辺) を足して同じことをします。

最大は、続けて s から t へ流せるだけ流します。補助の辺の逆向きの弧は s から t へ直接つながる弧で、最大流は必ずこれを使い切るので、流した量がそのまま F の最大になります。最小は、補助の辺の容量を今の流量に下げて止め、t から s へその量まで戻します。下限も湧き出しも無ければ、超頂点を使わずに s から t へ流すだけです。S から T へは push-relabel で流し、続けて流す 2 回目は Dinic になります (最大流の 2 回目からの flow の決まり)。

辺は最大流に直接溜めます。2 回目からの問い合わせと、問い合わせのあとに辺や頂点を足したときは、溜めた辺の元の容量から最大流を組み直します。

小さい網では辺ごとの流量を全部試した答えと、中くらいの網では別に書いた素朴な Dinic の答えと、100 万回ほど突き合わせました。procon-judge では、LOJ 115 (循環と辺ごとの流量)、LOJ 116 (最大)、LOJ 117 (最小)、aoj-1615 で確かめています。記録は algo-notes の notes/network-flow.md にあります。
