#pragma once
#include <algorithm>
#include <cassert>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>
// 逐次最短路の最小費用流 (primal-dual)。ACL の mcf_graph と同じく、add_edge で辺を足し、flow(s, t, limit) か
// slope(s, t, limit) で s から t へ流す。slope は、流す量 x ごとの費用の最小の折れ線 (下に凸) の折れ点を返す。ACL と違い、
// 呼ぶたびに今の流れから続けて流し (s と t は呼ぶたびに変えてよい)、流したあとにも頂点と辺を足せて、費用が負の辺も張れる。
// 流す時点の残余グラフに負閉路が無いこと (あとから足した辺を含めて)。
// 中身は、Dijkstra (radix heap、ポテンシャルで縮約した費用) で s からの距離を測り、縮約費用が 0 になった弧だけの網で
// 流せるだけ流すことを繰り返す。Dijkstra は t を取り出したあとも、t と同じ距離の頂点を取り出し終えるまで続け (t から
// 出る弧は緩和しない)、取り出した頂点だけポテンシャルを動かす。そのあと、最短路木の s から t への道に流してから、縮約費用
// 0 の網で t から逆向きの BFS で段数を付け、段数が 1 ずつ減る弧だけを使う DFS (Dinic) で流し、s に段数が付かなくなるまで
// 繰り返す。その距離の道がすべて流れるので、Dijkstra 1 回が折れ線のちょうど 1 区間になる。費用の種類が少ない網や容量の
// 大きい密な網では、道を 1 本ずつ流す形より桁で速い。逆向きの BFS は t へ届く頂点しか歩かず、縮約費用 0 の弧だけ相手の
// 弧の容量を読むので、距離がみな違う網でも、同じ距離の道がもう無いことを確かめる 1 回の BFS の分しか増えない。
// 初めの負の費用と、あとから足した辺で縮約費用が負になったときは、流す前にポテンシャルを直す (残余グラフが DAG なら位相
// 順に、そうでなければキューを使う Bellman-Ford で)。あとから足した頂点のポテンシャルは、古い頂点から入る弧の端の値で
// 決める。辺は初めの flow か slope で CSR に組み、そのあと足した辺の弧は端の頂点の区画の後ろに置く (区画に空きが無ければ、
// その頂点の弧を倍の大きさの区画にして配列の末尾へ移す)。辺を足しながら少しずつ流す使い方でも、呼ぶたびに組み直さない。
// 方式は procon-judge の self-flow-mincostflow-slope で比べて選んだ。記録は algo-notes の notes/network-flow.md。
// Flow と Cost は整数の型。(頂点の数) × (費用の絶対値の最大) が 2^62 に収まり、流した量と費用の和が Flow と Cost に収まること。
template <class Flow, class Cost> class PrimalDual {
 static_assert(std::is_integral_v<Flow> && std::is_integral_v<Cost>, "流量と費用は整数の型");
 using i64= long long;
 using u64= unsigned long long;
public:
 struct Edge {
  int from, to;
  Flow cap, flow;
  Cost cost;
 };
 explicit PrimalDual(int n= 0): n_(n) {}
 int add_vertex() { return n_++; }
 // from から to へ、容量 cap、1 単位あたり cost の辺。番号を返す。費用は負でもよい。
 int add_edge(int from, int to, Flow cap, Cost cost) {
  assert(0 <= from && from < n_ && 0 <= to && to < n_ && cap >= 0);
  es_.push_back({from, to, cap, cost});
  return (int)es_.size() - 1;
 }
 // 今の流れから続けて、s から t へ流せるだけ (limit まで) 流す。増えた流量と費用を返す。
 std::pair<Flow, Cost> flow(int s, int t) { return flow(s, t, std::numeric_limits<Flow>::max()); }
 std::pair<Flow, Cost> flow(int s, int t, Flow limit) { return run(s, t, limit, nullptr); }
 // flow と同じく流し、増やした量 x ごとの費用の増え方の最小 g(x) の折れ点を、(0, 0) から x の昇順に返す。傾きの同じ点は
 // 含めない。最後の点が増えた流量と費用。
 std::vector<std::pair<Flow, Cost>> slope(int s, int t) { return slope(s, t, std::numeric_limits<Flow>::max()); }
 std::vector<std::pair<Flow, Cost>> slope(int s, int t, Flow limit) {
  std::vector<std::pair<Flow, Cost>> res{{Flow(0), Cost(0)}};
  run(s, t, limit, &res);
  return res;
 }
 Edge get_edge(int i) const {
  const E& e= es_[i];
  return {e.from, e.to, e.cap, i < bm_ ? Flow(e.cap - g_[pos_[i]].cap) : Flow(0), e.cost};
 }
 std::vector<Edge> edges() const {
  std::vector<Edge> ret(es_.size());
  for(int i= 0; i < (int)ret.size(); ++i) ret[i]= get_edge(i);
  return ret;
 }
 // ポテンシャル p。残余のある弧 (容量の残る辺と、流れのある辺の逆向き) で、縮約費用 cost + p[from] - p[to] が 0 以上。
 std::vector<i64> potential() {
  build();
  std::vector<i64> ret(n_);
  for(int v= 0; v < n_; ++v) ret[v]= vs_[v].p;
  return ret;
 }
 int size() const { return n_; }
private:
 static constexpr i64 INF= std::numeric_limits<i64>::max();
 struct E {
  int from, to;
  Flow cap;
  Cost cost;
 };
 struct Arc {
  int to, rev;
  Flow cap;
  Cost cost;
 };
 // 頂点ごとのポテンシャルと、今の Dijkstra の距離 (16 byte に収めて、Dijkstra で辺 1 本ごとに読む場所を 1 か所にする)。
 struct V {
  i64 p, d;
 };
 int n_, bm_= 0, bn_= 0;  // 頂点の数と、CSR に入れた辺と頂点の数
 std::vector<E> es_;
 std::vector<Arc> g_;
 std::vector<int> bg_, en_, room_;  // 頂点ごとの弧の区画 [bg_, en_) と区画の大きさ
 std::vector<int> pos_, aid_;       // 辺ごとの順向きの弧の位置と、弧ごとの「辺の番号の 2 倍と向き」(辺を足すときに作る)
 std::vector<V> vs_;
 std::vector<char> vis_;                           // 今の Dijkstra で取り出したか
 std::vector<int> prv_, lv_, lr_, cur_;           // 最短路木の弧、段数、段数を付けた回、DFS の今の弧
 std::vector<int> touched_, visited_, path_, bq_;  // 距離を付けた頂点、取り出した頂点、DFS の道、BFS の待ち行列
 int round_id_= 0;
 std::vector<std::pair<u64, int>> bucket_[65];
 u64 last_= 0;
 int heap_size_= 0;

 // 頂点 v の弧を大きさ r の区画にして配列の末尾へ移し、移した弧の相手の逆向きの添字を付け直す。
 void relocate(int v, int r) {
  const int ns= g_.size(), ob= bg_[v], oe= en_[v];
  g_.resize(ns + r), aid_.resize(ns + r);
  for(int k= 0; k < oe - ob; ++k) {
   const int a= ns + k;
   g_[a]= g_[ob + k], aid_[a]= aid_[ob + k];
   if(!(aid_[a] & 1)) pos_[aid_[a] >> 1]= a;
  }
  for(int k= 0; k < oe - ob; ++k) {
   const int a= ns + k, b= g_[a].rev;
   if(b < 0) continue;
   if(ob <= b && b < oe) g_[a].rev= b - ob + ns;
   else g_[b].rev= a;
  }
  bg_[v]= ns, en_[v]= ns + (oe - ob), room_[v]= r;
 }
 int place(int v) {
  if(en_[v] - bg_[v] == room_[v]) relocate(v, std::max(4, 2 * room_[v]));
  return en_[v]++;
 }
 // 足した頂点と辺を入れる。足した辺の縮約費用が負なら、ポテンシャルを直す。
 void build() {
  const int m= es_.size();
  if(m == bm_ && n_ == bn_) return;
  pos_.resize(m);
  if(bm_ == 0) {
   bg_.assign(n_ + 1, 0);
   for(const E& e: es_) ++bg_[e.from + 1], ++bg_[e.to + 1];
   for(int v= 0; v < n_; ++v) bg_[v + 1]+= bg_[v];
   room_.resize(n_);
   for(int v= 0; v < n_; ++v) room_[v]= bg_[v + 1] - bg_[v];
   bg_.pop_back(), en_= bg_;
   g_.resize(2 * m);
   for(int i= 0; i < m; ++i) {
    const E& e= es_[i];
    const int a= en_[e.from]++, b= en_[e.to]++;
    g_[a]= {e.to, b, e.cap, e.cost}, g_[b]= {e.from, a, Flow(0), Cost(-e.cost)}, pos_[i]= a;
   }
  } else {
   if(aid_.size() != g_.size()) {
    aid_.assign(g_.size(), -1);
    for(int i= 0; i < bm_; ++i) aid_[pos_[i]]= 2 * i, aid_[g_[pos_[i]].rev]= 2 * i + 1;
   }
   bg_.resize(n_, (int)g_.size()), en_.resize(n_, (int)g_.size()), room_.resize(n_, 0);
   for(int i= bm_; i < m; ++i) {
    const E& e= es_[i];
    const int a= place(e.from);
    g_[a]= {e.to, -1, e.cap, e.cost}, aid_[a]= 2 * i, pos_[i]= a;
    const int b= place(e.to);
    g_[b]= {e.from, pos_[i], Flow(0), Cost(-e.cost)}, aid_[b]= 2 * i + 1, g_[pos_[i]].rev= b;
   }
  }
  // 足した頂点のポテンシャル。古い頂点から入る弧があれば、その端の値に費用を足したものの最小にし (入る弧の縮約費用が
  // 0 以上になる)、無ければ、古い頂点へ出る弧の端の値から費用を引いたものの最大にする。どちらも無ければ今の最大。
  if(n_ > bn_) {
   i64 mx= 0;
   for(int v= 0; v < bn_; ++v) mx= v ? std::max(mx, vs_[v].p) : vs_[v].p;
   const int old= bn_;
   std::vector<i64> lo(n_ - old, INF), hi(n_ - old, -INF);
   for(int i= bm_; i < m; ++i) {
    const E& e= es_[i];
    if(e.cap == 0) continue;
    if(e.to >= old && e.from < old) lo[e.to - old]= std::min(lo[e.to - old], vs_[e.from].p + (i64)e.cost);
    if(e.from >= old && e.to < old) hi[e.from - old]= std::max(hi[e.from - old], vs_[e.to].p - (i64)e.cost);
   }
   vs_.resize(n_, V{0, INF});
   for(int v= old; v < n_; ++v) vs_[v].p= lo[v - old] != INF ? lo[v - old] : hi[v - old] != -INF ? hi[v - old] : mx;
   vis_.resize(n_, 0), prv_.resize(n_), lv_.resize(n_), lr_.resize(n_, 0), cur_.resize(n_);
  }
  bool bad= false;
  for(int i= bm_; i < m && !bad; ++i) {
   const E& e= es_[i];
   if(e.cap > 0 && vs_[e.from].p + (i64)e.cost < vs_[e.to].p) bad= true;
  }
  bm_= m, bn_= n_;
  if(bad) repair();
 }
 // 残余のある弧すべてで縮約費用を 0 以上にする。今のポテンシャルから始め、残余グラフが DAG なら位相順に、そうでなければ
 // キューを使う Bellman-Ford で下げる。
 void repair() {
  std::vector<int> indeg(n_, 0), order;
  order.reserve(n_);
  for(int u= 0; u < n_; ++u)
   for(int j= bg_[u]; j < en_[u]; ++j)
    if(g_[j].cap > 0) ++indeg[g_[j].to];
  for(int v= 0; v < n_; ++v)
   if(!indeg[v]) order.push_back(v);
  for(size_t k= 0; k < order.size(); ++k) {
   const int u= order[k];
   for(int j= bg_[u]; j < en_[u]; ++j)
    if(g_[j].cap > 0 && !--indeg[g_[j].to]) order.push_back(g_[j].to);
  }
  if((int)order.size() == n_) {
   for(int u: order)
    for(int j= bg_[u]; j < en_[u]; ++j)
     if(const Arc& a= g_[j]; a.cap > 0) vs_[a.to].p= std::min(vs_[a.to].p, vs_[u].p + (i64)a.cost);
   return;
  }
  std::vector<int> q(n_), cnt(n_, 0);
  std::vector<char> inq(n_, 1);
  for(int v= 0; v < n_; ++v) q[v]= v;
  for(size_t h= 0; h < q.size(); ++h) {
   const int u= q[h];
   inq[u]= 0;
   for(int j= bg_[u]; j < en_[u]; ++j) {
    const Arc& a= g_[j];
    if(a.cap > 0 && vs_[u].p + (i64)a.cost < vs_[a.to].p) {
     vs_[a.to].p= vs_[u].p + (i64)a.cost;
     if(!inq[a.to]) {
      if(++cnt[a.to] > n_) {
       assert(!"残余グラフに負閉路がある");
       return;
      }
      inq[a.to]= 1, q.push_back(a.to);
     }
    }
   }
  }
 }
 // radix heap。桶 k (k >= 1) は、最後に取り出した鍵 last_ との排他的論理和の最上位の bit が k - 1 の鍵。今の鍵と同じ鍵は
 // 桶 0 に積まれ、次に取り出される。
 static int bucket_of(u64 x) { return x ? 64 - __builtin_clzll(x) : 0; }
 [[gnu::always_inline]] void hpush(u64 key, int v) { bucket_[bucket_of(key ^ last_)].push_back({key, v}), ++heap_size_; }
 std::pair<u64, int> hpop() {
  if(bucket_[0].empty()) {
   int i= 1;
   while(bucket_[i].empty()) ++i;
   u64 mn= bucket_[i][0].first;
   for(const auto& e: bucket_[i]) mn= std::min(mn, e.first);
   last_= mn;
   for(const auto& e: bucket_[i]) bucket_[bucket_of(e.first ^ last_)].push_back(e);
   bucket_[i].clear();
  }
  --heap_size_;
  const auto r= bucket_[0].back();
  bucket_[0].pop_back();
  return r;
 }
 // s からの縮約費用の距離を、t と同じ距離の頂点を取り出し終えるまで測り、取り出した頂点のポテンシャルを動かす。t に
 // 届かなければ false。prv_ に最短路木の弧を残す。
 bool dijkstra(int s, int t) {
  for(int v: touched_) vs_[v].d= INF;
  for(int v: visited_) vis_[v]= 0;
  touched_.clear(), visited_.clear();
  for(auto& b: bucket_) b.clear();
  last_= 0, heap_size_= 0;
  i64 dt= INF;
  vs_[s].d= 0, touched_.push_back(s), hpush(0, s);
  while(heap_size_) {
   const auto [key, v]= hpop();
   if(vis_[v] || (i64)key != vs_[v].d) continue;
   const i64 dv= vs_[v].d, pv= vs_[v].p;
   if(dv > dt) break;
   vis_[v]= 1, visited_.push_back(v);
   if(v == t) {
    dt= dv;
    continue;
   }
   for(int j= bg_[v], e= en_[v]; j < e; ++j) {
    const Arc& a= g_[j];
    if(!a.cap) continue;
    V& w= vs_[a.to];
    const i64 nd= dv + (i64)a.cost + pv - w.p;
    if(nd < w.d) {
     if(w.d == INF) touched_.push_back(a.to);
     w.d= nd, prv_[a.to]= j, hpush(nd, a.to);
    }
   }
  }
  if(!vis_[t]) return false;
  for(int v: visited_) vs_[v].p+= vs_[v].d - dt;
  return true;
 }
 // 縮約費用 0 の弧だけの網で、t から逆向きに段数を付ける。s に付けば true。
 bool levels(int s, int t) {
  ++round_id_;
  bq_.clear(), bq_.push_back(t), lv_[t]= 0, lr_[t]= round_id_;
  for(size_t h= 0; h < bq_.size(); ++h) {
   const int v= bq_[h];
   if(v == s) return true;
   const i64 pv= vs_[v].p;
   const int L= lv_[v] + 1;
   for(int j= bg_[v], e= en_[v]; j < e; ++j) {
    const Arc& b= g_[j];
    const int u= b.to;
    // 縮約費用が 0 かを先に見て、0 の弧だけ印と相手の弧の容量を読む。
    if(vs_[u].p - (i64)b.cost != pv || lr_[u] == round_id_ || !vis_[u]) continue;
    if(g_[b.rev].cap > 0) lv_[u]= L, lr_[u]= round_id_, bq_.push_back(u);
   }
  }
  return false;
 }
 // 段数が 1 ずつ減る縮約費用 0 の弧だけで、s から t へ limit まで流す (Dinic の DFS を再帰せずに書いたもの)。
 Flow block(int s, int t, Flow limit) {
  for(int v: bq_) cur_[v]= bg_[v];
  Flow pushed= 0;
  path_.clear();
  int u= s;
  for(;;) {
   if(u == t) {
    Flow f= limit - pushed;
    for(int j: path_) f= std::min(f, g_[j].cap);
    for(int j: path_) g_[j].cap-= f, g_[g_[j].rev].cap+= f;
    pushed+= f;
    if(pushed == limit) break;
    size_t k= 0;
    while(g_[path_[k]].cap > 0) ++k;
    path_.resize(k);
    u= k ? g_[path_[k - 1]].to : s;
    continue;
   }
   const i64 pu= vs_[u].p;
   const int want= lv_[u] - 1, e= en_[u];
   int j= cur_[u];
   for(; j < e; ++j) {
    const Arc& a= g_[j];
    if(a.cap > 0 && (i64)a.cost + pu == vs_[a.to].p && lr_[a.to] == round_id_ && lv_[a.to] == want) break;
   }
   cur_[u]= j;
   if(j < e) {
    path_.push_back(j), u= g_[j].to;
    continue;
   }
   lv_[u]= -1;
   if(u == s) break;
   const int b= path_.back();
   path_.pop_back();
   u= g_[g_[b].rev].to;
   ++cur_[u];
  }
  return pushed;
 }
 std::pair<Flow, Cost> run(int s, int t, Flow limit, std::vector<std::pair<Flow, Cost>>* out) {
  assert(0 <= s && s < n_ && 0 <= t && t < n_ && limit >= 0);
  build();
  Flow fl= 0;
  Cost total= 0;
  i64 prev_d= 0;
  if(s == t) return {fl, total};
  while(fl < limit && dijkstra(s, t)) {
   const i64 d= vs_[t].p - vs_[s].p;
   // 最短路木の道に先に流す。距離がみな違う網では、そのあとの BFS は、同じ距離の道がもう無いことを確かめる 1 回で済む。
   Flow f= limit - fl;
   for(int v= t; v != s; v= g_[g_[prv_[v]].rev].to) f= std::min(f, g_[prv_[v]].cap);
   for(int v= t; v != s; v= g_[g_[prv_[v]].rev].to) g_[prv_[v]].cap-= f, g_[g_[prv_[v]].rev].cap+= f;
   while(f < limit - fl && levels(s, t)) f+= block(s, t, limit - fl - f);
   fl+= f, total+= Cost(f) * Cost(d);
   if(out) {
    if(out->size() >= 2 && prev_d == d) out->pop_back();
    out->push_back({fl, total});
   }
   prev_d= d;
  }
  return {fl, total};
 }
};
