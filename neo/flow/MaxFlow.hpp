#pragma once
#include <algorithm>
#include <cassert>
#include <limits>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>
// 最大流。ACL の mf_graph と同じく、add_edge で辺を足して番号を受け取り、flow(s, t) か flow(s, t, limit) で流し、
// min_cut(s)、get_edge(i)、edges() で解を取り出す。ACL に無いものとして、逆向きの容量 (無向辺)、流れを保ったまま
// 容量を変える change_cap、流れを道に分ける decompose がある。
// 最初の flow は push-relabel (highest label、同じ高さの中は FIFO、gap、global relabel) で流す。t に届かない余りが
// 高い頂点を行き来して進まなくなったら (人が椅子より多い累積の鎖など)、余りのある頂点から t への増加路で流す形に移る
// (discharge_all の説明)。前流の段 (t へ流し、届かなかった余りは途中に残す) だけで値を返し、余りを s へ戻す段は、
// 辺ごとの流量か最小カットを初めて聞かれたとき (と、辺を足すときや次の flow の前) まで遅らせる。値だけの使い方では、
// この段の費用がかからない。始めるときは、s から出る容量の和が Cap に収まれば s の辺をすべて飽和させる。和が溢れる
// ときだけ、s を普通の頂点として、s から出る和と t へ入る和の小さい方の余りを s に持たせる。どちらでも余りの合計が
// Cap に収まるので、大きな容量の辺を張っても溢れない。
// 2 回目からの flow と、上限が効く最初の flow は Dinic で流す (DFS は再帰しない)。容量を少し変えて少しだけ流し直す
// 使い方では、呼ぶたびにグラフ全体の BFS がかかる push-relabel より速い。辺は flow を呼んだときに CSR に組み、そのあと
// 辺を足したら、流量を保ったまま組み直す。今の流れから続けて流すので、終点の違う flow を何度呼んでもよい。
// メモリは、組んだあとに、弧 1 本 (辺ごとに 2 本) につき 8 + sizeof(Cap) バイト、辺ごとに 5 バイト (順向きの弧の位置と容量の
// 形)、頂点ごとに int 6 つと Cap 1 つ (と、使った高さの分の表) を使う。逆向きの容量が 0 でも容量と同じでもない辺があると、
// 辺ごとに Cap 1 つが増える。add_edge で溜めた辺は、組みながら捨てる。流している間の待ち行列や道の弧の列は、そのとき
// 使っていない頂点ごとの配列を借りる。
// 方式は procon-judge の loj-127、self-flow-maxflow、self-flow-maxflow-shapes、aoj-2803、aoj-2313、aoj-2835 で比べて
// 選んだ。記録は algo-notes の notes/network-flow.md と notes/network-flow-api.md。
// Cap は整数の型。辺ごとの cap + rev_cap と、s から t への流量の答えが Cap に収まること。
template <class Cap> class MaxFlow {
 static_assert(std::is_integral_v<Cap>, "容量は整数の型 (浮動小数点数は丸めの誤差での振る舞いをまだ確かめていない)");
public:
 struct Edge {
  int from, to;
  Cap cap, rev_cap, flow;  // flow は from から to へ正、逆向きに流れていれば負
 };
 explicit MaxFlow(int n= 0): n_(n) {}
 int add_vertex() { return finish(), built_= false, n_++; }
 // from から to へ容量 cap、逆向きに容量 rev_cap の辺。番号を返す。
 int add_edge(int from, int to, Cap cap, Cap rev_cap= 0) {
  assert(0 <= from && from < n_ && 0 <= to && to < n_ && cap >= 0 && rev_cap >= 0);
  finish();
  if((npe_ & (CH - 1)) == 0 && npe_ > 0) pe_.emplace_back().reserve(CH);
  else if(npe_ == 0 && pe_.empty()) pe_.emplace_back();
  pe_.back().push_back({from, to, cap, rev_cap});
  built_= false;
  return m_ + npe_++;
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
 Edge get_edge(int i) { return finish(), edge_at(i); }
 std::vector<Edge> edges() {
  finish();
  std::vector<Edge> ret(m_ + npe_);
  for(int i= 0; i < (int)ret.size(); ++i) ret[i]= edge_at(i);
  return ret;
 }
 // 辺 i の容量を cap、逆向きの容量を rev_cap に変える。今の流れは s から t への流れとする。流量が新しい容量に収まらない分は、
 // まず辺の両端の間で迂回させ、迂回できない分は t の側と s の側から流れを戻して直す。s から t への流量が減った分を返す。
 Cap change_cap(int i, Cap cap, Cap rev_cap, int s, int t) {
  assert(cap >= 0 && rev_cap >= 0);
  build(), finish();
  if(pos_[i] < 0) {
   L& l= loops_[-1 - pos_[i]];
   return l.c= cap, l.r= rev_cap, 0;
  }
  A &a= as_[pos_[i]], &b= as_[a.rev];
  int u= b.to, v= a.to;
  Cap f= cap_of(i) - a.cap, ext= 0;
  set_kind(i, cap, rev_cap);
  if(f > cap) ext= f - cap, f= cap;
  else if(f < -rev_cap) ext= -rev_cap - f, f= -rev_cap, std::swap(u, v);
  a.cap= cap - f, b.cap= rev_cap + f;
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
 Cap change_cap(int i, Cap cap, int s, int t) { return build(), finish(), change_cap(i, cap, edge_at(i).rev_cap, s, t); }
 // 今の流れ (s から t への流れとする) を、s から t への道に分ける。道ごとの流量と、通る辺の番号の列を返す。逆向きの容量
 // のある辺は、to から from へ通ることもある。流れの中の循環は道に含めない。
 std::vector<std::pair<Cap, std::vector<int>>> decompose(int s, int t) {
  finish();
  const int m= m_ + npe_;
  std::vector<Cap> rem(m);
  std::vector<int> head(n_ + 1, 0), arc;
  // 流れのある辺の、流れの向きでの始点と終点。
  auto tail= [&](int i) { const Edge e= edge_at(i); return e.flow > 0 ? e.from : e.to; };
  auto tip= [&](int i) { const Edge e= edge_at(i); return e.flow > 0 ? e.to : e.from; };
  for(int i= 0; i < m; ++i)
   if(const Edge e= edge_at(i); e.flow != 0 && e.from != e.to) rem[i]= e.flow > 0 ? e.flow : -e.flow, ++head[tail(i) + 1];
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
  // resize で 0 を書かない (組むときに全部の位置へ書く)。書くまでページを触らないので、溜めた辺を捨てながら組める。
  A() {}
  A(int to, int rev, Cap cap): to(to), rev(rev), cap(cap) {}
 };
 struct P {  // まだ CSR に組んでいない辺
  int u, v;
  Cap c, r;
 };
 struct L {  // 自己ループ (流れないので CSR に入れない)
  int u;
  Cap c, r;
 };
 // 組んでいない辺 (番号が m_ 以上) は、CH 本ずつの塊に分けて pe_ に溜める。組むときに読み終えた塊から捨てるので、塊と
 // CSR の弧の両方を全部持つ時間を短くできる (弧へ書く場所が前から順に進むグラフでは、ピークがほぼ弧の分だけになる)。
 // 最初の塊は普通に伸ばし、小さいグラフで大きな塊を確保しない。
 static constexpr int CH= 1 << 14;
 int n_, m_= 0, npe_= 0;
 bool built_= false, pending_= false, used_= false;  // used_: 一度 flow を呼んだ
 int ps_= -1, pt_= -1;                               // 余りを戻す段が残っている flow の s と t
 std::vector<std::vector<P>> pe_;
 std::vector<L> loops_;
 // 組んだ辺 i について、pos_[i] は順向きの弧の位置 (自己ループなら -1 - loops_ の番号)。元の容量は、2 本の弧の残余容量の
 // 和 (cap + rev_cap のまま変わらない) から、kind_[i] が 0 (rev_cap = 0) なら和、1 (rev_cap = cap) なら和の半分として
 // 求める。2 (どちらでもない) の辺があるときだけ、全部の辺の元の容量を ecap_ に持つ。端点と流量も弧から求める。
 std::vector<int> st_, pos_;
 std::vector<unsigned char> kind_;
 std::vector<Cap> ecap_;
 std::vector<A> as_;
 // 高さごとの表 (ahead_、atail_、bhead_) は頂点の数だけ確保するが、0 で埋めず、使った高さ (hinit_ 未満) までだけ -1 で
 // 初期化する。書いていないページは常駐しないので、高さが低いままのグラフでは表の分のメモリがほとんど要らない。
 std::vector<int> h_, cur_, anext_, bnext_, bprev_;
 std::unique_ptr<int[]> ahead_, atail_, bhead_;
 int hinit_= 0;
 std::vector<Cap> ex_;
 int amax_= -1, bmax_= -1;  // global relabel が決める、活性な頂点と頂点の高さの最大
 // 組んだ辺 i (自己ループでない) の元の容量。
 Cap cap_of(int i) const {
  if(kind_[i] == 2) return ecap_[i];
  const A& a= as_[pos_[i]];
  const Cap sum= a.cap + as_[a.rev].cap;
  return kind_[i] ? sum / 2 : sum;
 }
 static unsigned char kind_of(Cap cap, Cap rev_cap) { return rev_cap == 0 ? 0 : rev_cap == cap ? 1 : 2; }
 // 組んだ辺 i の容量の形を、新しい容量に合わせて直す。ecap_ を持っていなければ、今の容量から作ってから直す。
 void set_kind(int i, Cap cap, Cap rev_cap) {
  const unsigned char k= kind_of(cap, rev_cap);
  if(k == 2 && ecap_.empty()) {
   ecap_.resize(m_);
   for(int j= 0; j < m_; ++j)
    if(pos_[j] >= 0) ecap_[j]= cap_of(j);
  }
  kind_[i]= k;
  if(!ecap_.empty()) ecap_[i]= cap;
 }
 Edge edge_at(int i) const {
  if(i >= m_) {
   const P& e= pe_[(i - m_) / CH][(i - m_) % CH];
   return {e.u, e.v, e.c, e.r, 0};
  }
  if(pos_[i] < 0) {
   const L& l= loops_[-1 - pos_[i]];
   return {l.u, l.u, l.c, l.r, 0};
  }
  const A &a= as_[pos_[i]], &b= as_[a.rev];
  const Cap c= cap_of(i);
  return {b.to, a.to, c, Cap(a.cap + b.cap - c), Cap(c - a.cap)};
 }
 // 保留している、余りを s へ戻す段を回す。
 void finish() {
  if(!pending_) return;
  discharge_all(ps_, pt_);
  ex_[ps_]= 0, pending_= false;
 }
 void init_heights(int d) {  // 高さ d までの表を使えるようにする
  for(; hinit_ <= d; ++hinit_) ahead_[hinit_]= bhead_[hinit_]= -1;
 }
 // 組んでいない辺を CSR に組む。前に組んだ辺があれば、今の残余容量のまま組んでいない辺の前に戻してから組み直す。CSR の
 // 弧は頂点ごとに辺の番号の順に並ぶ。組んだあとに、元の容量の形 kind_ と頂点ごとの配列を確保する。
 void build() {
  if(built_) return;
  const int n= n_, m0= m_;
  if(m0 > 0) {
   std::vector<std::vector<P>> all;
   int cnt= 0;
   auto push= [&](const P& e) {
    if(cnt % CH == 0) all.emplace_back().reserve(CH);
    all.back().push_back(e), ++cnt;
   };
   for(int i= 0; i < m0; ++i)
    if(pos_[i] >= 0) {
     const A &a= as_[pos_[i]], &b= as_[a.rev];
     push({b.to, a.to, a.cap, b.cap});
    } else {
     const L& l= loops_[-1 - pos_[i]];
     push({l.u, l.u, l.c, l.r});
    }
   for(const auto& ch: pe_)
    for(const P& e: ch) push(e);
   pe_.swap(all), npe_= cnt;
   std::vector<A>().swap(as_);
  }
  const int m= npe_;
  st_.assign(n + 1, 0), pos_.resize(m), kind_.resize(m), loops_.clear();
  for(const auto& ch: pe_)
   for(const P& e: ch)
    if(e.u != e.v) ++st_[e.u + 1], ++st_[e.v + 1];
  for(int v= 0; v < n; ++v) st_[v + 1]+= st_[v];
  as_.resize(st_[n]);
  bool general= !ecap_.empty();
  {
   std::vector<int> p(st_.begin(), st_.end() - 1);
   int i= 0;
   for(auto& ch: pe_) {
    for(const P& e: ch) {
     if(e.u == e.v) pos_[i]= -1 - (int)loops_.size(), loops_.push_back({e.u, e.c, e.r});
     else {
      const int a= p[e.u]++, b= p[e.v]++;
      as_[a]= {e.v, b, e.c}, as_[b]= {e.u, a, e.r};
      pos_[i]= a;
     }
     // 前に組んだ辺は元の容量の形をそのまま使う (残余容量は元の容量ではない)。
     if(i >= m0) general|= (kind_[i]= e.u == e.v ? 0 : kind_of(e.c, e.r)) == 2;
     ++i;
    }
    std::vector<P>().swap(ch);  // 読み終えた塊は捨てる
   }
  }
  std::vector<std::vector<P>>().swap(pe_), npe_= 0;
  if(general) {
   if(ecap_.empty()) {
    ecap_.resize(m0);
    for(int j= 0; j < m0; ++j)
     if(pos_[j] >= 0) ecap_[j]= cap_of(j);
   }
   ecap_.resize(m);
   for(int j= m0; j < m; ++j)
    if(pos_[j] >= 0) ecap_[j]= as_[pos_[j]].cap;
  }
  m_= m;
  h_.assign(n, 0), cur_.assign(n, 0), anext_.assign(n, -1), bnext_.assign(n, -1), bprev_.assign(n, -1), ex_.assign(n, 0);
  ahead_.reset(new int[n + 1]), atail_.reset(new int[n + 1]), bhead_.reset(new int[n + 1]), hinit_= 0;
  built_= true;
 }
 // 高さを sink からの距離に直す (global relabel)。excl は使わない頂点 (高さ n に置く)。活性な頂点の表と高さごとの頂点の
 // 表を作り直し、amax_ と bmax_ を決める。alive には生きている (sink へたどれる) 余りの和を入れる。CUT が 0 でなければ
 // sink に届かない余りの下限を返す。高さ d 以上の頂点の集合から出る残余容量は、高さ d から d - 1 への弧の容量の和 C_d
 // なので、sink に届く余りは min_d (E_{<d} + C_d) を超えない (E_{<d} は高さ d 未満の余り)。CUT が 1 のときは、高さ d の
 // 頂点を BFS で取り出したときに、もう読んでいる弧の容量と相手の高さから C_d を数える。CUT が 2 のときは、余りのある
 // 頂点から残余グラフでたどれる頂点だけで数える (余りの届かない頂点の容量は使えないので、下限が締まる)。その分、前向き
 // の BFS が 1 回増える。どちらも判定にしか使わないので double で持つ (和が Cap を溢れても困らない)。
 template <int CUT> double global_relabel(int sink, int excl, double& alive) {
  const int n= n_;
  std::fill(h_.begin(), h_.end(), n);
  std::fill(ahead_.get(), ahead_.get() + hinit_, -1), std::fill(bhead_.get(), bhead_.get() + hinit_, -1);
  // 待ち行列は cur_ を借りる (最後に今の弧を先頭に書き直す)。
  std::vector<int>& q= cur_;
  int amax= -1, bmax= -1, qh= 0, qt= 0, lay= 0;
  double e= 0, reach= 1e300, c= 0, e0= 0;  // e0: 取り出している層より下の余り
  h_[sink]= 0, q[qt++]= sink;
  while(qh < qt) {
   const int u= q[qh++], du= h_[u], d= du + 1;
   if(du != lay) {
    if constexpr(CUT == 1)
     if(lay > 0) reach= std::min(reach, e0 + c);
    lay= du, c= 0, e0= e;
   }
   if(d >= hinit_) init_heights(d);
   for(int k= st_[u]; k < st_[u + 1]; ++k) {
    const A& a= as_[k];
    if(const int w= a.to, hw= h_[w]; hw == n) {
     if(w != excl && as_[a.rev].cap > 0) {
      h_[w]= d, q[qt++]= w;
      bnext_[w]= bhead_[d], bprev_[w]= -1;
      if(bhead_[d] >= 0) bprev_[bhead_[d]]= w;
      bhead_[d]= w, bmax= d;
      if(ex_[w] > 0) {
       anext_[w]= -1;
       if(ahead_[d] < 0) ahead_[d]= atail_[d]= w;
       else anext_[atail_[d]]= w, atail_[d]= w;
       amax= d;
      }
     }
    } else if constexpr(CUT == 1) {
     if(hw == du - 1) c+= double(a.cap);
    }
   }
   if(du > 0) e+= double(ex_[u]);
  }
  if constexpr(CUT == 1)
   if(lay > 0) reach= std::min(reach, e0 + c);
  amax_= amax, bmax_= bmax, alive= e;
  if constexpr(CUT == 2) {
   // 前向きの BFS で余りからたどれる頂点を回りながら、層ごとの余り el と、層 d から d - 1 への弧の容量 cl を足す。
   // 待ち行列は、続けて cur_ を借りる。
   std::vector<char> mk(n);
   std::vector<double> cl(bmax + 2), el(bmax + 2);
   int fn= 0;
   for(int u= 0; u < n; ++u)
    if(h_[u] < n && u != sink && ex_[u] > 0) mk[u]= 1, q[fn++]= u;
   for(int i= 0; i < fn; ++i) {
    const int u= q[i], du= h_[u];
    el[du]+= double(ex_[u]);
    for(int k= st_[u]; k < st_[u + 1]; ++k)
     if(const A& a= as_[k]; a.cap > 0) {
      const int w= a.to, hw= h_[w];
      if(hw == du - 1) cl[du]+= double(a.cap);
      if(!mk[w] && hw < n && w != sink) mk[w]= 1, q[fn++]= w;
     }
   }
   double below= 0;
   for(int d= 1; d <= bmax; ++d) reach= std::min(reach, below + cl[d]), below+= el[d];
   reach= std::min(reach, below);
  }
  std::copy(st_.begin(), st_.end() - 1, cur_.begin());
  if constexpr(CUT == 0) return 0;
  else return std::max(0.0, e - reach);
 }
 // global relabel のあとの状態から、余りを highest label で流す (同じ高さの中は FIFO、gap)。高さを上げた仕事量 (1 回に
 // つき 12 と弧の数) が freq を超えたら false を返し (高さはもう正確でないので、呼ぶ側が global relabel する)、活性な頂点
 // がなくなったら true を返す。global relabel のたびに作り直す状態だけを持つ関数に分けて、ループの中の変数を減らす。
 bool pump(int sink, long long freq) {
  const int n= n_;
  int amax= amax_, bmax= bmax_, hinit= hinit_;
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
    if(nh >= hinit) {
     for(; hinit <= nh; ++hinit) ahead_[hinit]= bhead_[hinit]= -1;
    }
    ins_b(v, nh), bmax= std::max(bmax, nh);
   }
   if(work > freq) return hinit_= hinit, false;
  }
  return hinit_= hinit, true;
 }
 // 余りのある頂点から sink への増加路で流す (始点の多い Dinic)。層は sink からの距離で、余りを途中の頂点に残さない。
 // 余りのある頂点から sink へ残余グラフでたどれなくなれば true。仕事量 (BFS 1 回につき n と弧の数、DFS は見た弧の数、
 // 流すたびに道の長さ) が budget を超えたら、そこで止めて false を返す (流れは前流のまま)。
 bool augment_multi(int sink, int excl, long long budget) {
  const int n= n_;
  long long spent= 0;
  // 待ち行列 q、余りのある頂点の並び src、辿っている弧の列 stk は、ここでは使わない層の表と活性な頂点の表を借りる
  // (このあと global relabel で作り直す)。
  std::vector<int>&q= bnext_, &src= anext_, &stk= bprev_;
  for(;;) {
   if((spent+= n + st_[n]) > budget) return false;
   std::fill(h_.begin(), h_.end(), n);
   int qh= 0, qt= 0, ns= 0;
   h_[sink]= 0, q[qt++]= sink;
   while(qh < qt) {
    const int u= q[qh++], d= h_[u] + 1;
    for(int k= st_[u]; k < st_[u + 1]; ++k)
     if(const A& a= as_[k]; h_[a.to] == n && a.to != excl && as_[a.rev].cap > 0) {
      h_[a.to]= d, q[qt++]= a.to, cur_[a.to]= st_[a.to];
      if(ex_[a.to] > 0) src[ns++]= a.to;
     }
   }
   if(ns == 0) return true;
   for(int i0= 0; i0 < ns; ++i0) {
    // v0 から sink へ、高さが 1 ずつ下がる弧をたどる。sink に着いたら流し、詰まった弧の手前まで戻る。
    const int v0= src[i0];
    int v= v0, sp= 0;
    while(ex_[v0] > 0) {
     if(v == sink) {
      Cap d= ex_[v0];
      int cut= 0;
      for(int i= 0; i < sp; ++i)
       if(as_[stk[i]].cap < d) d= as_[stk[i]].cap, cut= i;
      for(int i= 0; i < sp; ++i) as_[stk[i]].cap-= d, as_[as_[stk[i]].rev].cap+= d;
      ex_[v0]-= d, ex_[sink]+= d;
      spent+= sp;
      sp= cut;
      v= sp == 0 ? v0 : as_[stk[sp - 1]].to;
      continue;
     }
     int& k= cur_[v];
     const int hv= h_[v], k0= k;
     for(; k < st_[v + 1]; ++k)
      if(const A& a= as_[k]; a.cap > 0 && h_[a.to] == hv - 1) break;
     if((spent+= 1 + k - k0) > budget) return false;
     if(k == st_[v + 1]) {
      h_[v]= n;  // 行き止まり
      if(v == v0) break;
      --sp;
      v= sp == 0 ? v0 : as_[stk[sp - 1]].to;
      ++cur_[v];
      continue;
     }
     stk[sp++]= k, v= as_[k].to;
    }
   }
  }
 }
 // 終点 sink に向けて、ex_ の余りを highest label で流す。excl は使わない頂点 (高さ n に置く)。
 // HLPP は、sink に届かない余りが、残余グラフではまだ sink へたどれる間 (人が椅子より多い累積の鎖、真ん中が詰まる 2 段の
 // 鎖など) に、高い頂点を行き来して 2 乗の時間がかかる。highest label なので、その間は届く余りが後回しになり、sink に
 // 何も届かない区間が続くこともある。そこで global relabel のたびに進み具合を見て、詰まっていそうなら、余りのある頂点
 // から sink への増加路 (augment_multi) を BFS 4 回分の仕事量だけ試す。試すのは、直前の区間で sink に届いた量 got が
 // 生きている余りより少なく、届かない余りの下限 stuck が初めに生きていた余りの 1/1024 以上あり、さらに残り (生きている
 // 余り - stuck) に今の進みで 4 区間より多くかかりそうか、stuck が got 以上のとき。試して、残りの半分以上を流せたら
 // 増加路で最後まで流し、流せなかったら増加路が向かないグラフとみて、もう試さずに HLPP を続ける。下限を求めて試さな
 // かったときは、次に求めるのを global relabel の回数が 1.5 倍になってからにする。下限は、got が生きている余りの 1/16
 // より少ないときだけ、余りからたどれる頂点に絞って数える (global_relabel<2>。間隔は 4 倍ずつ空ける)。余りの届かない
 // 頂点が sink の隣に多いと (両端の椅子にしか座れない人が多い網など)、絞らない下限はそれを数えて 0 になる。
 void discharge_all(int sink, int excl) {
  const int n= n_;
  double alive= 0;
  global_relabel<0>(sink, excl, alive);
  const double alive0= alive;
  const long long freq= 2 * (6LL * n + st_[n] / 2);
  // next: 次に下限を求める global relabel の番号 (試したあとは INT_MAX)。nextx: 次に絞った下限を求めてよい番号
  int grs= 0, next= 0, nextx= 0;
  Cap last= ex_[sink];
  while(!pump(sink, freq)) {
   ++grs;
   const double got= double(ex_[sink] - last);
   last= ex_[sink];
   if(grs >= next && got < alive) {
    const bool narrow= got * 16 < alive && grs >= nextx;
    if(narrow) nextx= 4 * grs;
    const double stuck= narrow ? global_relabel<2>(sink, excl, alive) : global_relabel<1>(sink, excl, alive), rest= alive - stuck;
    if(stuck * 1024 >= alive0 && (got * 4 < rest || stuck >= got)) {
     bool done= augment_multi(sink, excl, 4LL * (n + st_[n]));
     if(!done && double(ex_[sink] - last) * 2 >= rest) done= augment_multi(sink, excl, std::numeric_limits<long long>::max());
     if(done) return;
     next= std::numeric_limits<int>::max();
     global_relabel<0>(sink, excl, alive), last= ex_[sink];
    } else next= grs + (grs + 1) / 2;
   } else global_relabel<0>(sink, excl, alive);
  }
 }
 // Dinic で s から t へ limit まで流す。今の流れから続け、流し終えても正しい流れのまま。DFS は再帰しない。
 Cap augment(int s, int t, Cap limit) {
  Cap total= 0;
  // 段は h_、待ち行列は bnext_、辿っている弧の列は bprev_ を借りる (余りを戻す段は global relabel で作り直す)。
  std::vector<int>&lv= h_, &q= bnext_, &stk= bprev_;
  while(total < limit) {
   std::fill(lv.begin(), lv.end(), -1);
   int qh= 0, qt= 0;
   lv[s]= 0, q[qt++]= s;
   while(qh < qt && lv[t] < 0) {
    const int u= q[qh++];
    for(int k= st_[u]; k < st_[u + 1]; ++k)
     if(const A& a= as_[k]; a.cap > 0 && lv[a.to] < 0) lv[a.to]= lv[u] + 1, q[qt++]= a.to;
   }
   if(lv[t] < 0) break;
   std::copy(st_.begin(), st_.end() - 1, cur_.begin());
   // stk の先頭 sp 本が辿っている弧の列。t に着いたら bottleneck を流し、bottleneck の弧の手前まで戻る。
   int v= s, sp= 0;
   while(total < limit) {
    if(v == t) {
     Cap d= limit - total;
     int cut= 0;
     for(int i= 0; i < sp; ++i)
      if(as_[stk[i]].cap < d) d= as_[stk[i]].cap, cut= i;
     for(int i= 0; i < sp; ++i) as_[stk[i]].cap-= d, as_[as_[stk[i]].rev].cap+= d;
     total+= d;
     sp= cut;
     v= sp == 0 ? s : as_[stk[sp - 1]].to;
     continue;
    }
    int& k= cur_[v];
    for(; k < st_[v + 1]; ++k)
     if(const A& a= as_[k]; a.cap > 0 && lv[a.to] == lv[v] + 1) break;
    if(k == st_[v + 1]) {
     if(v == s) break;
     lv[v]= -1;  // 行き止まり
     --sp;
     v= sp == 0 ? s : as_[stk[sp - 1]].to;
     ++cur_[v];
     continue;
    }
    stk[sp++]= k, v= as_[k].to;
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
