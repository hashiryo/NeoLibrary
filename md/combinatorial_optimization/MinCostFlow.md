---
title: 最小費用流
documentation_of: ../../neo/combinatorial_optimization/MinCostFlow.hpp
---

頂点ごとの湧き出しと、辺ごとの下限、上限、費用を受け、費用の和が最小の流れ (最小費用の b-flow) を求めます。s から t へ流せるだけ流す形 (最小費用最大流) も解けます。

```cpp
template <class Flow, class Cost, bool UseNetworkSimplex = false> class MinCostFlow;
explicit MinCostFlow(int n = 0);
int add_vertex();
int add_edge(int from, int to, Flow upper, Cost cost);
int add_edge(int from, int to, Flow lower, Flow upper, Cost cost);
void add_supply(int v, Flow x);
bool solve();
Flow flow(int s, int t);
Flow flow(int s, int t, Flow limit);
bool feasible() const;
__int128 cost() const;  // Cost が浮動小数点数の型なら Cost
std::vector<Cost> potential() const;
struct Edge { int from, to; Flow lower, upper, flow; Cost cost; };
Edge get_edge(int i) const;
std::vector<Edge> edges() const;
int size() const;
```

頂点は 0 から n - 1 の番号で、`add_vertex()` で 1 つ足せます。`add_edge` は from から to への辺を足し、足した順の番号 (0 から) を返します。流量は lower 以上 upper 以下 (lower を省くと 0) で、1 単位あたりの費用は cost です。費用は負でもかまいません。自己ループも足せます。`add_supply(v, x)` は、頂点 v の湧き出しに x を足します。負の値なら吸い込みです。

`solve()` は、流量が下限と上限の間にあり、どの頂点でも出る量から入る量を引いたものが湧き出しに等しい流れ (b-flow) のうち、費用の和が最小のものを求めます。そういう流れがあれば true を返します。費用が負の閉路があれば、上限まで回して費用を下げます。

`flow(s, t, limit)` は、下限と湧き出しを満たしたうえで s から t へ limit まで流せるだけ流し、その中で費用を最小にして、流した量を返します。`flow(s, t)` は上限なしで、ACL の mcf_graph の flow と同じく最小費用最大流になります。limit までのどの量を流しても下限と湧き出しを満たせなければ、-1 を返します。流す量ごとの費用 (slope) が欲しいときや、流したあとに辺を足して続けて流すときは、[逐次最短路の最小費用流](PrimalDual.md) を使ってください。

解いたあとは、`cost()` で費用の和、`get_edge(i)` と `edges()` で辺ごとの流量、`potential()` でポテンシャルを取り出せます。費用の和は 2^63 を超えることがあるので、`__int128` で返します。ポテンシャル p は相補性を満たし、縮約費用 cost + p[from] - p[to] が正の辺は下限まで、負の辺は上限まで流れています。`feasible()` は、直前の `solve` か `flow` で流れがあったかを返します。費用の和を最大にしたいときは、費用の符号を反転して渡し、`cost()` の符号を戻します。

`Flow` は 64 bit までの整数の型、`Cost` は 64 bit までの整数の型か浮動小数点数の型です。辺ごとの (上限 - 下限) の和と、湧き出しの絶対値の和が long long に収まる必要があります。整数の費用では、費用の絶対値の和を $2^{60}$ 以下にしてください。

`Cost` が浮動小数点数の型なら、第 3 引数によらずネットワーク単体法で解き、`cost()` は `Cost` の型で返します。cost scaling は費用が整数であることを使うためです。距離や確率の対数を費用にする問題で使えます。丸めの誤差で終わらなくなるのを避けるため、縮約費用がごくわずかに負の辺は入れないので、解は誤差の分だけ最適からずれることがあります。

## 中身

既定は cost scaling で解き、第 3 引数を true にするとネットワーク単体法で解きます。

cost scaling は、Goldberg と Tarjan の push-relabel 版に 3 つの手を足したものです。Dial のバケツで距離を測る global price update、LEMON の partial augment (長さ 4 まで)、段ごとに元の費用で最適かを確かめて打ち切る手です。費用を (頂点の数 + 1) 倍して整数で扱い、ε を 16 ずつ割っていきます。始める前に、最大流 ([最大流](MaxFlow.md)) で b-flow があるかを確かめます。費用の絶対値の最大に (頂点の数 + 1) を掛けたものが $2^{61}$ を超えるときと、途中でポテンシャルが $2^{61}$ を超えて動いたときは、ネットワーク単体法で解き直します。長い道に大きな費用が並ぶ網では、最適なポテンシャルの幅そのものが整数に収まらなくなるためです。

ネットワーク単体法は、木を LEMON と同じ形 (先行順の thread と部分木の大きさ) で持ち、動かした部分木のポテンシャルを thread をたどって付け直します。入る辺は block search で選びます。浮動小数点数の費用は long double で計算します。人工の辺の費用を LEMON と同じく (費用の絶対値の最大 + 1)(n + 1) にしてポテンシャルを小さく保ち、pivot 1024 回ごとと最後に、ポテンシャルを木の辺の費用から付け直して誤差をなくします。入れる辺は、縮約費用がポテンシャルの大きさの $64 \cdot 1024 \varepsilon$ 倍 (ε は long double の計算機イプシロン) より小さいものだけにします。

cost scaling は一般の網では頂点の数にほぼ比例して伸びますが、鎖や星の形の網では 2 乗で伸びます。保序回帰、株の売買、締切付きの仕事、直線の上の輸送のような網です。費用の細かい段で、前の段の流れを長い道に沿って流し直すことになるためです。global price update はいちばん近い吸い込みまでの距離しか測らないので、その回数が鎖の長さに比例します。そういう網では、第 3 引数を true にしてネットワーク単体法を選んでください。ネットワーク単体法も、頂点の数を 2 倍にすると 2.6 倍から 5.4 倍かかるので、20 万頂点ほどの網はどの形でも数秒かかります。

方式は procon-judge の問題で比べて選びました。次の時間は、どれも同じ回の同じ CPU (EPYC 7763) の x64-gcc で測った最大のケースです。self-flow-mincostflow には、性質の違う 5 つの族が並んでいます。NETGEN 風、疎な割り当て、格子の輸送、s から t への流れ、負の費用と下限を混ぜた網で、頂点は 1 万から 10 万です。cost scaling は 1.30 秒、ネットワーク単体法は 4.24 秒で、今の Library の NetworkSimplex は 5.46 秒でした。self-flow-mincostflow-chain には、上の鎖や星の形の網が並んでいます (頂点 1 万から 30 万)。ネットワーク単体法は 4.66 秒で、cost scaling は 10 秒を超えます。今の Library の NetworkSimplex は 6.79 秒でした。頂点の少ない yosupo-min-cost-b-flow ($N \le 100$) では、cost scaling は始める前の最大流の分だけ遅くなります。cost scaling が 4.0 ms、ネットワーク単体法が 0.6 ms、今の Library が 0.9 ms でした。yuki-1615 は、疎な二部グラフの重み最大マッチングを解く問題です (頂点 5 万、辺 7.5 万)。同じ CPU でプロセス全体の時間を測ると、cost scaling が 303 ms、ネットワーク単体法が 231 ms でした。今の Library の NetworkSimplex とは同じ回で比べていません。入力を読んだあとの計算の部分を手元の M2 で測ると、cost scaling が 225 ms、ネットワーク単体法が 185 ms、Library が 180 ms でした。浮動小数点数の費用は、確率の積の最大を -log を費用にした最小費用流で解く LOJ 2300 で確かめました。記録は algo-notes の notes/network-flow.md にあります。
