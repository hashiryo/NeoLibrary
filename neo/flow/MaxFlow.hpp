#pragma once
#include <algorithm>
#include <cassert>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>
// 最大流。ACL の mf_graph と同じく、add_edge で辺を足して番号を受け取り、flow(s, t) か flow(s, t, limit) で流し、
// min_cut(s)、get_edge(i)、edges() で解を取り出す。ACL に無いものとして、逆向きの容量 (無向辺)、流れを保ったまま
// 容量を変える change_cap、流れを道に分ける decompose がある。
// 最初の flow は push-relabel (highest label、同じ高さの中は FIFO、gap、global relabel) で流す。前流の段 (t へ流し、
// 届かなかった余りは途中に残す) だけで値を返し、余りを s へ戻す段は、辺ごとの流量か最小カットを初めて聞かれたとき
// (と、辺を足すときや次の flow の前) まで遅らせる。値だけの使い方では、この段の費用がかからない。始めるときは、s から
// 出る容量の和が Cap に収まれば s の辺をすべて飽和させる。和が溢れるときだけ、s を普通の頂点として、s から出る和と
// t へ入る和の小さい方の余りを s に持たせる。どちらでも余りの合計が Cap に収まるので、大きな容量の辺を張っても溢れない。
// 2 回目からの flow と、上限が効く最初の flow は Dinic で流す (DFS は再帰しない)。容量を少し変えて少しだけ流し直す
// 使い方では、呼ぶたびにグラフ全体の BFS がかかる push-relabel より速い。辺は flow を呼んだときに CSR に組み、そのあと
// 辺を足したら、流量を保ったまま組み直す。今の流れから続けて流すので、終点の違う flow を何度呼んでもよい。
// 方式は procon-judge の loj-127、self-flow-maxflow、aoj-2803、aoj-2313、aoj-2835 で比べて選んだ。記録は algo-notes の
// notes/network-flow.md と notes/network-flow-api.md。
// Cap は整数の型。辺ごとの cap + rev_cap と、s から t への流量の答えが Cap に収まること。
template <class Cap> class MaxFlow {
 static_assert(std::is_integral_v<Cap>, "容量は整数の型 (浮動小数点数は丸めの誤差での振る舞いをまだ確かめていない)");
public:
 struct Edge {
  int from, to;
  Cap cap, rev_cap, flow;  // flow は from から to へ正、逆向きに流れていれば負
 };
 explicit MaxFlow(int n= 0): n_(n) {}
 int add_vertex() { return sync(), built_= false, n_++; }
 // from から to へ容量 cap、逆向きに容量 rev_cap の辺。番号を返す。
 int add_edge(int from, int to, Cap cap, Cap rev_cap= 0) {
  assert(0 <= from && from < n_ && 0 <= to && to < n_ && cap >= 0 && rev_cap >= 0);
  sync();
  eu_.push_back(from), ev_.push_back(to), ecap_.push_back(cap), erev_.push_back(rev_cap), eflow_.push_back(0);
  built_= false;
  return int(eu_.size()) - 1;
 }
 // s から t へ流せるだけ流し、増えた量を返す。今の流れから続ける。
 Cap flow(int s, int t) { return run(s, t, false, 0); }
 // 増やす量を limit までにする。
 Cap flow(int s, int t, Cap limit) { return run(s, t, true, limit); }
 // 直前の flow のあとの残余グラフで s から届く頂点。
 std::vector<bool> min_cut(int s) {
  build(), finish();
  std::vector<bool> vis(n_);
  std::vector<int> q{s};
  vis[s]= true;
  for(size_t i= 0; i < q.size(); ++i)
   for(int k= st_[q[i]]; k < st_[q[i] + 1]; ++k)
    if(const A& a= as_[k]; a.cap > 0 && !vis[a.to]) vis[a.to]= true, q.push_back(a.to);
  return vis;
 }
 Edge get_edge(int i) {
  sync();
  return {eu_[i], ev_[i], ecap_[i], erev_[i], eflow_[i]};
 }
 std::vector<Edge> edges() {
  sync();
  std::vector<Edge> ret(eu_.size());
  for(size_t i= 0; i < eu_.size(); ++i) ret[i]= {eu_[i], ev_[i], ecap_[i], erev_[i], eflow_[i]};
  return ret;
 }
 // 辺 i の容量を cap、逆向きの容量を rev_cap に変える。今の流れは s から t への流れとする。流量が新しい容量に収まらない分は、
 // まず辺の両端の間で迂回させ、迂回できない分は t の側と s の側から流れを戻して直す。s から t への流量が減った分を返す。
 Cap change_cap(int i, Cap cap, Cap rev_cap, int s, int t) {
  assert(cap >= 0 && rev_cap >= 0);
  build(), finish();
  if(pos_[i] < 0) return ecap_[i]= cap, erev_[i]= rev_cap, 0;
  A& a= as_[pos_[i]];
  int u= eu_[i], v= ev_[i];
  Cap f= cap_of_forward(i), ext= 0;
  ecap_[i]= cap, erev_[i]= rev_cap;
  if(f > cap) ext= f - cap, f= cap;
  else if(f < -rev_cap) ext= -rev_cap - f, f= -rev_cap, std::swap(u, v);
  a.cap= cap - f, as_[a.rev].cap= rev_cap + f, eflow_[i]= f;
  if(ext == 0) return 0;
  // u に ext だけ余り、v に ext だけ足りない。
  ext-= flow(u, v, ext);
  if(ext > 0) {
   if(v != t) flow(t, v, ext);
   if(u != s) flow(u, s, ext);
  }
  finish();
  return ext;
 }
 Cap change_cap(int i, Cap cap, int s, int t) { return change_cap(i, cap, erev_[i], s, t); }
 // 今の流れ (s から t への流れとする) を、s から t への道に分ける。道ごとの流量と、通る辺の番号の列を返す。逆向きの容量
 // のある辺は、to から from へ通ることもある。流れの中の循環は道に含めない。
 std::vector<std::pair<Cap, std::vector<int>>> decompose(int s, int t) {
  sync();
  const int m= eu_.size();
  std::vector<Cap> rem(m);
  std::vector<int> head(n_ + 1, 0), arc;
  auto tail= [&](int i) { return eflow_[i] > 0 ? eu_[i] : ev_[i]; };
  auto tip= [&](int i) { return eflow_[i] > 0 ? ev_[i] : eu_[i]; };
  for(int i= 0; i < m; ++i)
   if(eflow_[i] != 0 && eu_[i] != ev_[i]) rem[i]= eflow_[i] > 0 ? eflow_[i] : -eflow_[i], ++head[tail(i) + 1];
  for(int x= 0; x < n_; ++x) head[x + 1]+= head[x];
  arc.resize(head[n_]);
  {
   std::vector<int> p(head.begin(), head.end() - 1);
   for(int i= 0; i < m; ++i)
    if(rem[i] > 0) arc[p[tail(i)]++]= i;
  }
  std::vector<int> it(head.begin(), head.end() - 1), on(n_, -1), pv{s}, pe;
  std::vector<std::pair<Cap, std::vector<int>>> ret;
  on[s]= 0;
  // pe[from..] の辺から bottleneck を引き、頂点 pv[from] まで戻る。
  auto cancel= [&](int from) {
   Cap d= rem[pe[from]];
   for(size_t k= from; k < pe.size(); ++k) d= std::min(d, rem[pe[k]]);
   for(size_t k= from; k < pe.size(); ++k) rem[pe[k]]-= d;
   return d;
  };
  auto back_to= [&](int k) {
   while((int)pv.size() > k + 1) on[pv.back()]= -1, pv.pop_back(), pe.pop_back();
  };
  for(;;) {
   const int x= pv.back();
   if(x == t) {
    const Cap d= cancel(0);
    ret.emplace_back(d, pe);
    back_to(0);
    continue;
   }
   while(it[x] < head[x + 1] && rem[arc[it[x]]] == 0) ++it[x];
   if(it[x] == head[x + 1]) {
    if(x == s) break;
    rem[pe.back()]= 0;  // 正しい流れなら起きない行き止まり。入ってきた辺を捨てて戻る
    back_to((int)pv.size() - 2);
    continue;
   }
   const int e= arc[it[x]], y= tip(e);
   if(on[y] >= 0) {
    pe.push_back(e);
    cancel(on[y]);
    pe.pop_back();
    back_to(on[y]);
    continue;
   }
   on[y]= pv.size(), pv.push_back(y), pe.push_back(e);
  }
  return ret;
 }
 int size() const { return n_; }
private:
 struct A {
  int to, rev;
  Cap cap;
 };
 int n_;
 bool built_= false, pending_= false, used_= false;  // used_: 一度 flow を呼んだ
 int ps_= -1, pt_= -1;                               // 余りを戻す段が残っている flow の s と t
 std::vector<int> eu_, ev_;
 std::vector<Cap> ecap_, erev_, eflow_;
 std::vector<int> st_, pos_;  // pos_[i] = arc index of edge i (forward)
 std::vector<A> as_;
 std::vector<int> h_, cur_, anext_, atail_, ahead_, bnext_, bprev_, bhead_, q_;
 std::vector<Cap> ex_;
 Cap cap_of_forward(int i) const { return pos_[i] < 0 ? Cap(0) : ecap_[i] - as_[pos_[i]].cap; }
 void sync() {
  if(!built_) return;
  finish();
  for(size_t i= 0; i < pos_.size(); ++i)
   if(pos_[i] >= 0) eflow_[i]= ecap_[i] - as_[pos_[i]].cap;
 }
 // 保留している、余りを s へ戻す段を回す。
 void finish() {
  if(!pending_) return;
  discharge_all(ps_, pt_);
  ex_[ps_]= 0, pending_= false;
 }
 void build() {
  if(built_) return;
  const int n= n_, m= eu_.size();
  st_.assign(n + 1, 0), pos_.assign(m, -1);
  for(int i= 0; i < m; ++i)
   if(eu_[i] != ev_[i]) ++st_[eu_[i] + 1], ++st_[ev_[i] + 1];
  for(int v= 0; v < n; ++v) st_[v + 1]+= st_[v];
  as_.resize(st_[n]);
  std::vector<int> p(st_.begin(), st_.end() - 1);
  for(int i= 0; i < m; ++i)
   if(eu_[i] != ev_[i]) {
    const int a= p[eu_[i]]++, b= p[ev_[i]]++;
    as_[a]= {ev_[i], b, ecap_[i] - eflow_[i]}, as_[b]= {eu_[i], a, erev_[i] + eflow_[i]};
    pos_[i]= a;
   }
  h_.assign(n, 0), cur_.assign(n, 0), anext_.assign(n, -1), atail_.assign(n + 1, -1), ahead_.assign(n + 1, -1);
  bnext_.assign(n, -1), bprev_.assign(n, -1), bhead_.assign(n + 1, -1), q_.assign(n, 0), ex_.assign(n, 0);
  built_= true;
 }
 // 終点 sink に向けて、ex_ の余りを highest label で流す。excl は使わない頂点 (高さ n に置く)。
 void discharge_all(int sink, int excl) {
  const int n= n_;
  int amax= -1, bmax= -1;
  auto ins_b= [&](int v, int d) {
   bnext_[v]= bhead_[d], bprev_[v]= -1;
   if(bhead_[d] >= 0) bprev_[bhead_[d]]= v;
   bhead_[d]= v;
  };
  auto del_b= [&](int v, int d) {
   if(bprev_[v] >= 0) bnext_[bprev_[v]]= bnext_[v];
   else bhead_[d]= bnext_[v];
   if(bnext_[v] >= 0) bprev_[bnext_[v]]= bprev_[v];
  };
  auto ins_a= [&](int v, int d) {
   anext_[v]= -1;
   if(ahead_[d] < 0) ahead_[d]= atail_[d]= v;
   else anext_[atail_[d]]= v, atail_[d]= v;
   amax= std::max(amax, d);
  };
  auto global_relabel= [&] {
   std::fill(h_.begin(), h_.end(), n);
   std::fill(ahead_.begin(), ahead_.end(), -1), std::fill(bhead_.begin(), bhead_.end(), -1);
   amax= bmax= -1;
   int qh= 0, qt= 0;
   h_[sink]= 0, q_[qt++]= sink;
   while(qh < qt) {
    const int u= q_[qh++], d= h_[u] + 1;
    for(int k= st_[u]; k < st_[u + 1]; ++k)
     if(const A& a= as_[k]; h_[a.to] == n && a.to != excl && as_[a.rev].cap > 0) {
      h_[a.to]= d, q_[qt++]= a.to, cur_[a.to]= st_[a.to];
      ins_b(a.to, d), bmax= d;
      if(ex_[a.to] > 0) ins_a(a.to, d);
     }
   }
  };
  global_relabel();
  const long long freq= 2 * (6LL * n + st_[n] / 2);
  long long work= 0;
  while(amax >= 0) {
   const int v= ahead_[amax];
   if(v < 0) {
    --amax;
    continue;
   }
   ahead_[amax]= anext_[v];
   if(h_[v] != amax || ex_[v] == 0) continue;
   for(int hv= amax;;) {
    int k= cur_[v];
    for(const int ed= st_[v + 1]; k < ed; ++k)
     if(A& a= as_[k]; a.cap > 0 && h_[a.to] == hv - 1) {
      const int w= a.to;
      const Cap d= std::min(ex_[v], a.cap);
      if(w != sink && ex_[w] == 0) ins_a(w, hv - 1);
      a.cap-= d, as_[a.rev].cap+= d, ex_[v]-= d, ex_[w]+= d;
      if(ex_[v] == 0) break;
     }
    if(ex_[v] == 0) {
     cur_[v]= k;
     break;
    }
    int nh= n, nk= st_[v];
    for(int j= st_[v]; j < st_[v + 1]; ++j)
     if(const A& a= as_[j]; a.cap > 0 && h_[a.to] + 1 < nh) nh= h_[a.to] + 1, nk= j;
    work+= 12 + st_[v + 1] - st_[v];
    del_b(v, hv);
    if(bhead_[hv] < 0) {
     for(int d= hv; d <= bmax; ++d) {
      for(int u= bhead_[d]; u >= 0; u= bnext_[u]) h_[u]= n;
      bhead_[d]= -1;
     }
     bmax= hv - 1, h_[v]= n;
     break;
    }
    if(nh >= n) {
     h_[v]= n;
     break;
    }
    h_[v]= hv= nh, cur_[v]= nk;
    ins_b(v, nh), bmax= std::max(bmax, nh);
   }
   if(work > freq) work= 0, global_relabel();
  }
 }
 // Dinic で s から t へ limit まで流す。今の流れから続け、流し終えても正しい流れのまま。DFS は再帰しない。
 Cap augment(int s, int t, Cap limit) {
  Cap total= 0;
  std::vector<int>& lv= h_;
  std::vector<int> stk;
  while(total < limit) {
   std::fill(lv.begin(), lv.end(), -1);
   int qh= 0, qt= 0;
   lv[s]= 0, q_[qt++]= s;
   while(qh < qt && lv[t] < 0) {
    const int u= q_[qh++];
    for(int k= st_[u]; k < st_[u + 1]; ++k)
     if(const A& a= as_[k]; a.cap > 0 && lv[a.to] < 0) lv[a.to]= lv[u] + 1, q_[qt++]= a.to;
   }
   if(lv[t] < 0) break;
   std::copy(st_.begin(), st_.end() - 1, cur_.begin());
   // stk は辿っている弧の列。t に着いたら bottleneck を流し、bottleneck の弧の手前まで戻る。
   stk.clear();
   int v= s;
   while(total < limit) {
    if(v == t) {
     Cap d= limit - total;
     int cut= 0;
     for(int i= 0; i < (int)stk.size(); ++i)
      if(as_[stk[i]].cap < d) d= as_[stk[i]].cap, cut= i;
     for(int k: stk) as_[k].cap-= d, as_[as_[k].rev].cap+= d;
     total+= d;
     stk.resize(cut);
     v= stk.empty() ? s : as_[stk.back()].to;
     continue;
    }
    int& k= cur_[v];
    for(; k < st_[v + 1]; ++k)
     if(const A& a= as_[k]; a.cap > 0 && lv[a.to] == lv[v] + 1) break;
    if(k == st_[v + 1]) {
     if(v == s) break;
     lv[v]= -1;  // 行き止まり
     stk.pop_back();
     v= stk.empty() ? s : as_[stk.back()].to;
     ++cur_[v];
     continue;
    }
    stk.push_back(k), v= as_[k].to;
   }
  }
  return total;
 }
 static Cap add_sat(Cap x, Cap y) { return x > std::numeric_limits<Cap>::max() - y ? std::numeric_limits<Cap>::max() : x + y; }
 Cap run(int s, int t, bool limited, Cap limit) {
  assert(0 <= s && s < n_ && 0 <= t && t < n_ && s != t);
  build(), finish();
  // 2 回目からは、今の流れから少しだけ増やす使い方 (容量を変えて流し直すなど) が多いので Dinic で流す。
  if(used_) return augment(s, t, limited ? limit : std::numeric_limits<Cap>::max());
  used_= true;
  Cap out= 0, in= 0;
  for(int k= st_[s]; k < st_[s + 1]; ++k) out= add_sat(out, as_[k].cap);
  for(int k= st_[t]; k < st_[t + 1]; ++k) in= add_sat(in, as_[as_[k].rev].cap);
  const Cap u= std::min(out, in);  // 流量はこれを超えない
  if(!(u > 0)) return 0;
  // 上限が効くときは Dinic で流す。効かないなら上限なしと同じ。
  if(limited && limit < u) return augment(s, t, limit);
  std::fill(ex_.begin(), ex_.end(), Cap(0));
  if(out < std::numeric_limits<Cap>::max()) {
   // s から出る容量の和が Cap に収まる。s の辺をすべて飽和させ、s は高さ n に置いて使わない。余りの合計はこの和を超えない。
   for(int k= st_[s]; k < st_[s + 1]; ++k)
    if(A& a= as_[k]; a.cap > 0) ex_[a.to]+= a.cap, as_[a.rev].cap+= a.cap, a.cap= 0;
   discharge_all(t, s);
  } else {
   // 和が溢れる。s を普通の頂点として u だけの余りを持たせる (仮想の始点から容量 u の辺で s に流し込むのと同じ)。
   // 余りの合計は u を超えない。s の隣から余りが押し戻されて配り直すことがあり、上の形より遅くなりうる。
   ex_[s]= u;
   discharge_all(t, -1);
  }
  const Cap got= ex_[t];
  ex_[t]= 0, ex_[s]= 0;
  pending_= true, ps_= s, pt_= t;
  return got;
 }
};
