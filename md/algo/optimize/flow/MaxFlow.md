---
title: 最大流
documentation_of: ../../neo/algo/optimize/flow/MaxFlow.hpp
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
std::vector<int> cut_kinds(int s, int t);
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

`cut_kinds(s, t)` は、辺ごとに、s-t の最小カットに入るかを返します。0 はどの最小カットにも入らない、1 はどれかの最小カットに入る (すべてには入らない)、2 はすべての最小カットに入る、です。直前に上限なしの `flow(s, t)` を呼んだあとに使います。辺が最小カット (S, T) に入るとは、容量が正の向きで S から T へ渡ることです。from が S で to が T にあるときと、逆向きの容量が正なら to が S で from が T にあるときです。容量 0 の向きは数えません。最大流のあとの残余グラフの強連結成分で決めるので、手間は O(n + m) です。最小カットの S は、残余の弧で閉じていて s を含み t を含まない頂点集合です。そのため、飽和した向きの弧 u → v は、u と v の成分が違えばどれかの最小カットに入り、u が s の成分、v が t の成分にあればすべての最小カットに入ります。洛谷 P4126 ([AHOI2009] 最小割) はこれで答えられます。

2 頂点 c, d が最小カットで同じ側になりうるか、違う側になりうるかを問う問題 (AtCoder CODE FESTIVAL 2015 エキシビション A) には専用の口を置いていません。残余グラフで s から届く頂点の集合を A、t へ届く頂点の集合を B とすると、次のように決まります。

- 同じ側になりうるのは、c も d も B に無いときか、c も d も A に無いときです。
- c が s の側、d が t の側になりうるのは、c が B に無く、d が A に無く、残余グラフで c から d へ届かないときです。

残余の弧は `edges()` から作れます。`flow < cap` なら from → to、`flow > -rev_cap` なら to → from です。この弧で `strongly_connected_components` を求め、成分を縮めた DAG で成分ごとに届く成分の集合を bitset で持てば、問いに定数時間で答えられます。

`change_cap(i, cap, s, t)` は、辺 i の容量を cap に変えます (`rev_cap` も変える版もあります)。今の流れは s から t への流れとし、流量が新しい容量に収まらない分は、まず辺の両端の間で迂回させ、迂回できない分は流れを戻して直します。返すのは、s から t への流量が減った分です。容量を下げたあとの流れは、変更後のグラフの最大流になっています。容量を上げても流量は増えないので、続けて `flow(s, t)` を呼びます。

`decompose(s, t)` は、今の流れ (s から t への流れとする) を s から t への道に分け、道ごとの流量と、通る辺の番号の列を返します。逆向きの容量のある辺は、to から from へ通ることもあります。流れの中の循環は道に含めません。

`Cap` は整数の型で、辺ごとの `cap + rev_cap` と、s から t への流量の答えが `Cap` に収まる必要があります。s から出る辺の容量の和は収まらなくてかまわないので、無限大の代わりに大きな容量の辺を何本張っても溢れません。

## 中身

最初の `flow` は push-relabel で流します。高さが最大の頂点から余りを流し、同じ高さの中では先に入れたものから取り出します。gap と global relabel も使います。前流の段 (t へ流し、届かなかった余りは途中に残す) で値を返し、余りを s へ戻す段は、辺ごとの流量か最小カットを初めて聞かれたときまで遅らせます。値だけの使い方では、この段の費用がかかりません。2 回目からの `flow` は Dinic で流します。容量を少し変えて少しだけ流し直す使い方では、呼ぶたびにグラフ全体を BFS でなめる push-relabel より速いためです。ただし上限の無い `flow` では、Dinic の仕事量 (BFS と DFS で見た頂点と弧の数) が網を 4 回なめる分を超えると、残りを push-relabel で流します。大きく流し直す使い方では、Dinic の段が 100 回を超えて続くことがあるためです。下限付きの最大流の 2 段目がその例です。

t に届かない余りは、残余グラフではまだ t へたどれる間、高い頂点を行き来し続けます。push-relabel はこの行き来で遅くなることがあります。人が椅子より多い割り当てを累積の鎖で表した網 (ARC 076 F の形) では、人 20 万、椅子 10 万で 10 秒を超えていました。そこで global relabel のたびに、t に届かない余りの下限と、t に届いた量の進みを見ます。進まなくなっていれば、余りのある頂点から t への増加路を少しだけ試し、残りの半分以上を流せたら、そのまま増加路で最後まで流します。流せなければ push-relabel を続けます。push-relabel が得意なグラフでは、試すことはまずありません。

メモリは、組んだあとは辺 1 本につき 2 × (8 + sizeof(Cap)) + 5 バイトと、頂点 1 個につき int 6 つと Cap 1 つほどです。元の容量は、逆向きの容量が 0 の辺と、容量と逆向きの容量が等しい辺では持たずに、2 本の弧の残余容量の和から求めます。add_edge で溜めた辺は、CSR に組みながら捨てます。流している間の待ち行列や道の弧の列は、そのとき使っていない頂点ごとの配列を借ります。Linux で、1000 × 1000 の格子 (頂点 100 万、無向辺 300 万) のピークは、容量が int で 112 MB、long long で 139 MB でした (前の版は 195 MB と 249 MB)。

2 MB 以上の配列 (弧、辺ごとの弧の位置と容量の形、頂点ごとの配列) は、Linux では 2 MB の区画 (huge page) で渡すよう OS に頼みます。OS は普段 4 KB の区画でメモリを渡し、プログラムが区画を初めて使うときに、そのつど割り込んで用意します。区画の場所を覚えておく CPU の表 (TLB) にも、大きい網の配列は収まりきりません。さらに、配列の先頭を 2 KB + 64 byte の倍数ずつずらします。同じ添字で引く頂点ごとの配列がどれも 2 MB の境目から始まると、同じ添字の要素がキャッシュの同じ場所を取り合い、Xeon 6975P-C では huge page で縮んだ分を上回って遅くなりました。高さごとの表は、書いた高さの分だけ常駐させたいので普通に確保します。

まだ遅い形もあります。人がそれぞれ R_i 以上の椅子にだけ座れる割り当てを 1 本の累積の鎖で表すと、増加路が鎖を長く往復するので、push-relabel も Dinic も遅くなります。人と椅子が 20 万ずつで 5 秒を超えました (procon-judge の atcoder-arc076-f の seed 8)。

方式は procon-judge の問題で比べて選びました。次の時間は、どれも同じ回の同じ CPU の x64-gcc で測った最大のケースです。self-flow-maxflow には、性質の違う 8 つのグラフの族が並んでいます。この `MaxFlow` は EPYC 7763 で 589 ms で、今の Library の push-relabel は 1.7 秒、Dinic は 2.9 秒でした。self-flow-maxflow-shapes には、想定解は別にあるが網を作れば最大流で解ける問題の形 (格子、累積の鎖、セグメント木、閉包) が並んでいます。EPYC 9V45 で 570 ms で、今の Library の Dinic は 2.0 秒でした (push-relabel は累積の鎖で 10 秒を超える)。密なグラフの loj-127 では、EPYC 9V45 で 85 ms でした。ここだけは今の Library の push-relabel (68 ms) より 2 割ほど遅くなります。容量を変えながら流し直す問題で、Dinic で流す今の Library と比べた時間は次のとおりです。aoj-2803 は 6.3 ms と 4.1 ms (EPYC 9V45)、aoj-2313 は 71 ms と 50 ms、aoj-2835 は 37 ms と 51 ms (どちらも EPYC 9V74) でした。記録は algo-notes の notes/network-flow.md にあります。
