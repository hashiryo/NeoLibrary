#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>
#include "neo/flow/MaxFlow.hpp"
// 最小費用の b-flow。頂点ごとの湧き出し b (負なら吸い込み) と、辺ごとの下限、上限、費用を受け、流量が下限と上限の間にあり、
// 各頂点の出る量から入る量を引いたものが b に等しい流れのうち、費用の和が最小のものを求める。費用は負でもよく、自己ループも
// 受ける。solve() で解き、費用の和 cost()、辺ごとの流量 get_edge(i).flow、相補性を満たすポテンシャル potential() を取り出す。
// flow(s, t, limit) は、下限と b を満たしたうえで s から t へ limit まで流せるだけ流し、その中で費用を最小にする。
//
// 既定の中身は cost scaling。Goldberg と Tarjan の push-relabel 版に、Dial のバケツで距離を測る global price update、LEMON の
// partial augment (長さ 4 まで)、段ごとに元の費用で最適かを確かめて打ち切る手を足したもの。第 3 引数を true にすると、
// ネットワーク単体法 (LEMON 型の木と block search) で解く。cost scaling は一般の網では頂点の数にほぼ比例して伸びるが、鎖や
// 星の形の網 (保序回帰、株の売買、締切付きの仕事など) では段ごとに流し直す量が鎖の長さに比例し、2 乗で伸びる。そういう網では
// ネットワーク単体法を選ぶ。Cost が浮動小数点数の型なら、第 3 引数によらずネットワーク単体法で解き、費用の和を Cost の型で
// 返す (cost scaling は費用が整数であることを使う)。浮動小数点数の費用の解は、丸めの誤差の分だけ最適からずれうる。
//
// cost scaling は費用を (頂点の数 + 1) 倍して整数で扱う。費用の絶対値の最大に (頂点の数 + 1) を掛けたものが 2^61 を超える
// ときと、途中でポテンシャルが 2^61 を超えて動いたときは、ネットワーク単体法で解き直す。ネットワーク単体法は、費用の絶対値の
// 和が 2^60 以下であること (整数の費用のとき)。どちらも、辺ごとの (上限 - 下限) の和と b の絶対値の和が long long に収まること。
// 方式は procon-judge の yosupo-min-cost-b-flow、self-flow-mincostflow、self-flow-mincostflow-chain で比べて選んだ。記録は
// algo-notes の notes/network-flow.md。
namespace min_cost_flow_internal {
using i64= long long;
// 辺。C は中で使う費用の型 (整数の費用なら long long、浮動小数点数の費用なら long double)。
template <class C> struct Arc {
 int s, t;
 i64 l, u;
 C c;
};
template <class C> struct Result {
 bool ok= false, fallback= false;
 std::vector<C> pot;
 std::vector<i64> flow;
};
// cost scaling。解けなければ ok が false、整数があふれそうならそこで止めて fallback を立てる。
inline Result<i64> cost_scaling(int n, const std::vector<i64>& b, const std::vector<Arc<i64>>& es) {
 constexpr int ALPHA= 16, MAXLEN= 4, ET= 100;
 constexpr i64 LIM= i64(1) << 61;
 const int m= es.size();
 Result<i64> res;
 std::vector<i64> ex(b);
 for(const Arc<i64>& e: es) ex[e.s]-= e.l, ex[e.t]+= e.l;
 {
  i64 sum= 0, need= 0;
  for(int v= 0; v < n; ++v) sum+= ex[v], need+= ex[v] > 0 ? ex[v] : 0;
  if(sum != 0) return res;
  if(need > 0) {
   MaxFlow<i64> mf(n + 2);
   for(const Arc<i64>& e: es)
    if(e.s != e.t && e.u > e.l) mf.add_edge(e.s, e.t, e.u - e.l);
   for(int v= 0; v < n; ++v)
    if(ex[v] > 0) mf.add_edge(n, v, ex[v]);
    else if(ex[v] < 0) mf.add_edge(v, n + 1, -ex[v]);
   if(mf.flow(n, n + 1) != need) return res;
  }
 }
 const i64 S= n + 1;
 {
  i64 mx= 0;
  for(const Arc<i64>& e: es)
   if(e.s != e.t) mx= std::max(mx, std::abs(e.c));
  if(mx > LIM / S) return res.fallback= true, res;
 }
 // 弧 k の容量と、逆向きの弧 rev[k] の容量を cp[k] に並べて持つ。逆向きの弧の費用は -cost[k]。
 std::vector<int> st(n + 1), pos(m, -1);
 for(const Arc<i64>& e: es)
  if(e.s != e.t) ++st[e.s + 1], ++st[e.t + 1];
 for(int v= 0; v < n; ++v) st[v + 1]+= st[v];
 const int A= st[n];
 std::vector<int> to(A), rev(A), from(A);
 std::vector<std::array<i64, 2>> cp(A);
 std::vector<i64> cost(A);
 i64 eps= 1;
 {
  std::vector<int> p(st.begin(), st.end() - 1);
  for(int i= 0; i < m; ++i) {
   const Arc<i64>& e= es[i];
   if(e.s == e.t) continue;
   const int a= p[e.s]++, r= p[e.t]++;
   to[a]= e.t, from[a]= e.s, rev[a]= r, cp[a]= {e.u - e.l, 0}, cost[a]= e.c * S;
   to[r]= e.s, from[r]= e.t, rev[r]= a, cp[r]= {0, e.u - e.l}, cost[r]= -e.c * S;
   pos[i]= a;
   eps= std::max(eps, std::abs(e.c) * S);
  }
 }
 // ポテンシャル p は 0 から下がるだけで、-LIM 以上に保つ (縮約費用の計算があふれないように)。
 std::vector<i64> p(n, 0);
 std::vector<int> cur(n), q(n), rank(n), bhead(n + 2), bnext, bnode;
 std::vector<char> inq(n), onpath(A);
 bool over= false;
 auto push_arc= [&](int k, i64 d) {
  cp[k][0]-= d, cp[k][1]+= d;
  auto& r= cp[rev[k]];
  r[0]+= d, r[1]-= d;
 };
 long long relabels= 0;
 const long long pu_every= std::max<long long>(1, 2LL * n);
 // extra が 0 以上なら、その弧も (容量が無くても) 残余の辺として扱う。partial augment で道の先を relabel するときに、来た
 // 辺の逆向きを入れて、来た辺が admissible でなくなるようにする (行き止まりの頂点で同じ辺を行き来しないため)。
 auto relabel= [&](int v, int extra) {
  i64 best= -(i64(1) << 62);
  bool any= false;
  for(int j= st[v]; j < st[v + 1]; ++j)
   if(cp[j][0] > 0) best= std::max(best, p[to[j]] - cost[j]), any= true;
  if(extra >= 0) best= std::max(best, p[to[extra]] - cost[extra]), any= true;
  if(!any) return;
  p[v]= best - eps, cur[v]= st[v], ++relabels;
  if(p[v] < -LIM) over= true;
 };
 // global price update。余りの負の頂点から残余の辺を逆にたどる最短路 (辺の長さは floor(縮約費用 / ε) + 1) を、距離ごとの
 // バケツで求める。距離が n を超える頂点はたどらず、余りのある頂点をすべて取り出したら打ち切り、残りは最後の距離にする。
 auto price_update= [&] {
  std::fill(bhead.begin(), bhead.end(), -1);
  bnext.clear(), bnode.clear();
  auto push= [&](int r, int v) { bnode.push_back(v), bnext.push_back(bhead[r]), bhead[r]= (int)bnode.size() - 1; };
  int left= 0;
  for(int v= 0; v < n; ++v) {
   rank[v]= n + 1, left+= ex[v] > 0;
   if(ex[v] < 0) rank[v]= 0, push(0, v);
  }
  int last= 0;
  for(int r= 0; r <= n && left > 0; ++r) {
   while(bhead[r] >= 0 && left > 0) {
    const int id= bhead[r];
    bhead[r]= bnext[id];
    const int w= bnode[id];
    if(rank[w] != r) continue;  // 取り出し済みか、もっと近いバケツに積み直した
    rank[w]= ~r, last= r, left-= ex[w] > 0;
    const i64 pw= p[w];
    for(int k= st[w]; k < st[w + 1]; ++k) {
     if(cp[k][1] == 0) continue;
     const int v= to[k];
     if(rank[v] < 0) continue;
     const i64 rc= -cost[k] + p[v] - pw;
     const i64 nd= r + (rc >= 0 ? rc / eps : -((-rc + eps - 1) / eps)) + 1;
     if(nd < rank[v] && nd <= n) rank[v]= (int)nd, push((int)nd, v);
    }
   }
  }
  for(int v= 0; v < n; ++v) {
   const i64 r= rank[v] < 0 ? ~rank[v] : last;
   if(r > (p[v] + LIM) / eps) return void(over= true);
   p[v]-= eps * r, cur[v]= st[v];
  }
  relabels= 0;
 };
 // ポテンシャルを S で割った値から、残余グラフで元の費用の label-correcting をかけて相補性を満たすように直す。関係を直す
 // 回数が budget (負なら無制限) を超えたら false を返す (負の閉路があり、まだ最適でない)。
 std::vector<i64> pi(n);
 bool settled= false;
 auto settle= [&](long long budget) {
  for(int v= 0; v < n; ++v) pi[v]= p[v] >= 0 ? p[v] / S : -((-p[v] + S - 1) / S);
  int qh= 0, qn= 0;
  for(int v= 0; v < n; ++v) q[qn++]= v, inq[v]= 1;
  long long relax= 0;
  while(qn) {
   const int v= q[qh];
   qh= qh + 1 == n ? 0 : qh + 1, --qn, inq[v]= 0;
   for(int k= st[v]; k < st[v + 1]; ++k)
    if(cp[k][0] > 0 && pi[v] + cost[k] / S < pi[to[k]]) {
     pi[to[k]]= pi[v] + cost[k] / S;
     if(!inq[to[k]]) inq[to[k]]= 1, q[(qh + qn++) % n]= to[k];
     if(++relax == budget) {
      std::fill(inq.begin(), inq.end(), 0);
      return false;
     }
    }
  }
  return settled= true;
 };
 std::vector<int> path;
 do {
  eps= std::max<i64>(1, eps / ALPHA);
  for(int v= 0; v < n; ++v)
   for(int k= st[v]; k < st[v + 1]; ++k)
    if(cp[k][0] > 0 && cost[k] + p[v] - p[to[k]] < 0) {
     const i64 d= cp[k][0];
     push_arc(k, d), ex[v]-= d, ex[to[k]]+= d;
    }
  price_update();
  if(over) return res.fallback= true, res;
  int qh= 0, qn= 0;
  for(int v= 0; v < n; ++v) {
   inq[v]= ex[v] > 0;
   if(inq[v]) q[(qh + qn++) % n]= v;
  }
  while(qn) {
   const int v= q[qh];
   qh= qh + 1 == n ? 0 : qh + 1, --qn, inq[v]= 0;
   while(ex[v] > 0) {
    // v から admissible な辺 (縮約費用が負の残余の辺) をたどって道を伸ばし、余りの負の頂点に着くか、長さが MAXLEN に
    // なったら、道に沿って流せるだけ流す。行き詰まったら、その頂点を relabel して 1 歩戻る。
    int tip= v;
    path.clear();
    while((int)path.size() < MAXLEN && ex[tip] >= 0) {
     int& k= cur[tip];
     while(k < st[tip + 1] && !(cp[k][0] > 0 && cost[k] + p[tip] - p[to[k]] < 0)) ++k;
     if(k < st[tip + 1]) {
      const int a= k;
      path.push_back(a);
      if(onpath[a]) break;  // 輪になった
      onpath[a]= 1, tip= to[a];
      continue;
     }
     relabel(tip, tip != v ? rev[path.back()] : -1);
     if(over) return res.fallback= true, res;
     if(tip != v) {
      const int a= path.back();
      path.pop_back(), onpath[a]= 0, tip= from[a];
     }
    }
    for(const int a: path) {
     onpath[a]= 0;
     const int u= from[a], w= to[a];
     const i64 d= std::min(cp[a][0], ex[u]);
     push_arc(a, d), ex[u]-= d, ex[w]+= d;
     if(w != v && ex[w] > 0 && !inq[w]) inq[w]= 1, q[(qh + qn++) % n]= w;
    }
    if(relabels >= pu_every) {
     price_update();
     if(over) return res.fallback= true, res;
    }
    if(path.empty()) break;  // 起きない (余りのある頂点には残余の辺がある)
   }
  }
  // 費用を S 倍しているので、ε が S を下回った段で元の費用では 1 未満の最適性まで来ている。多くの網では ε が S / 100 ほどの
  // 段で最適になるので、そこからは段を終えるたびに確かめ、最適なら終える。
  if(eps > 1 && eps <= S / ET && settle(A / 2 + n)) break;
 } while(eps > 1);
 if(!settled) settle(-1);
 res.ok= true;
 res.pot= std::move(pi);
 res.flow.resize(m);
 for(int i= 0; i < m; ++i) {
  const Arc<i64>& e= es[i];
  res.flow[i]= pos[i] < 0 ? (e.c < 0 ? e.u : e.l) : e.l + cp[pos[i]][1];
 }
 return res;
}
// ネットワーク単体法。木を LEMON の NetworkSimplex と同じ形で持つ。木の頂点を先行順に並べた thread (と逆向きの rev_thread)、
// 部分木の頂点の数 succ_num、部分木の thread の上で最後の頂点 last_succ を持ち、動かした部分木のポテンシャルは thread を
// たどるだけで付け直す。入る辺は block search (大きさは ceil(sqrt(m + n))) で選ぶ。根 n と、各頂点から根へ (または
// 根から) の人工の辺から始める。人工の辺の費用は、整数なら費用の絶対値の和 + 1、浮動小数点数なら (費用の絶対値の最大 + 1)
// (n + 1) (LEMON と同じく、ポテンシャルを小さく保って丸めの誤差を抑える)。
// 浮動小数点数の費用では、丸めの誤差で同じ pivot を繰り返さないように、縮約費用が -tol より小さい辺だけを入れる。tol は
// ポテンシャルの大きさ (人工の辺の費用) に、1024 回の更新でたまる丸めの誤差の 64 倍を掛けたもの。ポテンシャルは pivot
// 1024 回ごとと最後に、木の辺の費用から付け直して誤差をなくす。
template <class C> Result<C> network_simplex(int n, const std::vector<i64>& b, const std::vector<Arc<C>>& es) {
 constexpr bool FP= std::is_floating_point_v<C>;
 constexpr int RECOMPUTE= 1024;
 Result<C> res;
 const int m= es.size(), M= m + n, root= n;
 {
  i64 sum= 0;
  for(int v= 0; v < n; ++v) sum+= b[v];
  if(sum != 0) return res;
 }
 std::vector<int> src(M), dst(M);
 std::vector<i64> cap(M), flow(M, 0);
 std::vector<C> cost(M);
 std::vector<int8_t> state(M, 1);  // 1 が下限、-1 が上限、0 が木
 std::vector<i64> sup(b);
 C art= 1;
 [[maybe_unused]] C mx= 0;
 for(int e= 0; e < m; ++e) {
  const Arc<C>& a= es[e];
  src[e]= a.s, dst[e]= a.t, cap[e]= a.u - a.l, cost[e]= a.c;
  sup[a.s]-= a.l, sup[a.t]+= a.l;
  const C ac= a.c < 0 ? -a.c : a.c;
  if constexpr(FP) mx= std::max(mx, ac);
  else art+= ac, assert(art <= (i64(1) << 60));
 }
 if constexpr(FP) art= (mx + 1) * (n + 1);
 const C tol= FP ? art * (C)(64 * RECOMPUTE) * std::numeric_limits<C>::epsilon() : C(0);
 // 木。pred_dir は木の辺が子から親へ向いていれば 1、親から子なら -1。
 std::vector<int> parent(n + 1), pred(n + 1), thread(n + 1), rev_thread(n + 1), succ_num(n + 1), last_succ(n + 1);
 std::vector<int8_t> pred_dir(n + 1);
 std::vector<C> pi(n + 1);
 parent[root]= -1, pred[root]= -1, thread[root]= 0, rev_thread[0]= root, succ_num[root]= n + 1, last_succ[root]= root - 1, pi[root]= 0;
 for(int u= 0, e= m; u < n; ++u, ++e) {
  parent[u]= root, pred[u]= e, thread[u]= u + 1, rev_thread[u + 1]= u, succ_num[u]= 1, last_succ[u]= u;
  state[e]= 0, cost[e]= art;
  if(sup[u] >= 0) pred_dir[u]= 1, pi[u]= -art, src[e]= u, dst[e]= root, cap[e]= flow[e]= sup[u];
  else pred_dir[u]= -1, pi[u]= art, src[e]= root, dst[e]= u, cap[e]= flow[e]= -sup[u];
 }
 if(n > 0) thread[n - 1]= root, rev_thread[root]= n - 1;
 std::vector<int> dirty;
 const int block_size= std::max(1, std::max(int(std::ceil(std::sqrt((double)M))), std::min(10, n + 1)));
 // 木の辺の縮約費用が 0 になるように、ポテンシャルを根から thread の順 (先行順) に付け直す。
 [[maybe_unused]] auto recompute= [&] {
  for(int u= thread[root]; u != root; u= thread[u]) pi[u]= pi[parent[u]] - pred_dir[u] * cost[pred[u]];
 };
 int next_arc= 0;
 [[maybe_unused]] int pivots= 0;
 for(;;) {
  if constexpr(FP)
   if(++pivots % RECOMPUTE == 0) recompute();
  // 入る辺 (block search)。
  int in_arc= -1;
  {
   C mn= -tol;
   int cnt= block_size, e= next_arc;
   for(int seen= M; seen--; e= e + 1 == M ? 0 : e + 1) {
    if(const C c= state[e] * (cost[e] + pi[src[e]] - pi[dst[e]]); c < mn) mn= c, in_arc= e;
    if(--cnt == 0) {
     if(in_arc >= 0) break;
     cnt= block_size;
    }
   }
   if(in_arc < 0) break;
   next_arc= e;
  }
  // 最小の共通祖先。
  int join;
  {
   int u= src[in_arc], v= dst[in_arc];
   while(u != v) {
    if(succ_num[u] < succ_num[v]) u= parent[u];
    else v= parent[v];
   }
   join= u;
  }
  // 出る辺。first から join と、second から join をたどる (後で見つけたものを優先して、木を strongly feasible に保つ)。
  int first, second, u_out= -1, result= 0;
  if(state[in_arc] == 1) first= src[in_arc], second= dst[in_arc];
  else first= dst[in_arc], second= src[in_arc];
  i64 delta= cap[in_arc];
  for(int u= first; u != join; u= parent[u]) {
   const int e= pred[u];
   const i64 d= pred_dir[u] == -1 ? cap[e] - flow[e] : flow[e];
   if(d < delta) delta= d, u_out= u, result= 1;
  }
  for(int u= second; u != join; u= parent[u]) {
   const int e= pred[u];
   const i64 d= pred_dir[u] == 1 ? cap[e] - flow[e] : flow[e];
   if(d <= delta) delta= d, u_out= u, result= 2;
  }
  int u_in, v_in;
  if(result == 1) u_in= first, v_in= second;
  else u_in= second, v_in= first;
  // 流す。
  if(delta > 0) {
   const i64 val= state[in_arc] * delta;
   flow[in_arc]+= val;
   for(int u= src[in_arc]; u != join; u= parent[u]) flow[pred[u]]-= pred_dir[u] * val;
   for(int u= dst[in_arc]; u != join; u= parent[u]) flow[pred[u]]+= pred_dir[u] * val;
  }
  if(result == 0) {
   state[in_arc]= -state[in_arc];
   continue;
  }
  state[in_arc]= 0;
  state[pred[u_out]]= flow[pred[u_out]] == 0 ? 1 : -1;
  // 木を付け替える (LEMON の updateTreeStructure)。u_out の部分木を切り離し、u_in から u_out までの幹の頂点の親を逆向きに
  // して、u_in を v_in の子として付け直す。
  {
   const int old_rev_thread= rev_thread[u_out], old_succ_num= succ_num[u_out], old_last_succ= last_succ[u_out];
   const int v_out= parent[u_out];
   if(u_in == u_out) {
    parent[u_in]= v_in, pred[u_in]= in_arc, pred_dir[u_in]= u_in == src[in_arc] ? 1 : -1;
    if(thread[v_in] != u_out) {
     int after= thread[old_last_succ];
     thread[old_rev_thread]= after, rev_thread[after]= old_rev_thread;
     after= thread[v_in];
     thread[v_in]= u_out, rev_thread[u_out]= v_in;
     thread[old_last_succ]= after, rev_thread[after]= old_last_succ;
    }
   } else {
    const int thread_continue= old_rev_thread == v_in ? thread[old_last_succ] : thread[v_in];
    int stem= u_in, par_stem= v_in, last= last_succ[u_in], after= thread[last];
    thread[v_in]= u_in;
    dirty.clear(), dirty.push_back(v_in);
    while(stem != u_out) {
     const int next_stem= parent[stem];
     thread[last]= next_stem, dirty.push_back(last);
     const int before= rev_thread[stem];
     thread[before]= after, rev_thread[after]= before;
     parent[stem]= par_stem, par_stem= stem, stem= next_stem;
     last= last_succ[stem] == last_succ[par_stem] ? rev_thread[par_stem] : last_succ[stem];
     after= thread[last];
    }
    parent[u_out]= par_stem;
    thread[last]= thread_continue, rev_thread[thread_continue]= last;
    last_succ[u_out]= last;
    if(old_rev_thread != v_in) thread[old_rev_thread]= after, rev_thread[after]= old_rev_thread;
    for(int u: dirty) rev_thread[thread[u]]= u;
    int tmp_sc= 0;
    const int tmp_ls= last_succ[u_out];
    for(int u= u_out, p= parent[u]; u != u_in; u= p, p= parent[u]) {
     pred[u]= pred[p], pred_dir[u]= -pred_dir[p];
     tmp_sc+= succ_num[u] - succ_num[p];
     succ_num[u]= tmp_sc, last_succ[p]= tmp_ls;
    }
    pred[u_in]= in_arc, pred_dir[u_in]= u_in == src[in_arc] ? 1 : -1;
    succ_num[u_in]= old_succ_num;
   }
   const int up_limit_out= last_succ[join] == v_in ? join : -1;
   const int last_succ_out= last_succ[u_out];
   for(int u= v_in; u != -1 && last_succ[u] == v_in; u= parent[u]) last_succ[u]= last_succ_out;
   if(join != old_rev_thread && v_in != old_rev_thread) {
    for(int u= v_out; u != up_limit_out && last_succ[u] == old_last_succ; u= parent[u]) last_succ[u]= old_rev_thread;
   } else if(last_succ_out != old_last_succ) {
    for(int u= v_out; u != up_limit_out && last_succ[u] == old_last_succ; u= parent[u]) last_succ[u]= last_succ_out;
   }
   for(int u= v_in; u != join; u= parent[u]) succ_num[u]+= old_succ_num;
   for(int u= v_out; u != join; u= parent[u]) succ_num[u]-= old_succ_num;
  }
  // 動かした部分木のポテンシャルを付け直す。
  {
   const C sigma= pi[v_in] - pi[u_in] - pred_dir[u_in] * cost[in_arc];
   const int end= thread[last_succ[u_in]];
   for(int u= u_in; u != end; u= thread[u]) pi[u]+= sigma;
  }
 }
 for(int e= m; e < M; ++e)
  if(flow[e] != 0) return res;
 if constexpr(FP) recompute();
 res.ok= true;
 res.pot.assign(pi.begin(), pi.begin() + n);
 res.flow.resize(m);
 for(int e= 0; e < m; ++e) res.flow[e]= flow[e] + es[e].l;
 return res;
}
}
template <class Flow, class Cost, bool UseNetworkSimplex= false> class MinCostFlow {
 static_assert(std::is_integral_v<Flow> && sizeof(Flow) <= 8, "流量は 64 bit までの整数の型");
 static_assert((std::is_integral_v<Cost> && sizeof(Cost) <= 8) || std::is_floating_point_v<Cost>, "費用は 64 bit までの整数の型か浮動小数点数の型");
 static constexpr bool FP= std::is_floating_point_v<Cost>;
 using i64= long long;
 using C= std::conditional_t<FP, long double, i64>;  // 中で使う費用の型
 using Arc= min_cost_flow_internal::Arc<C>;
 using Sum= std::conditional_t<FP, Cost, __int128>;  // 費用の和の型
public:
 struct Edge {
  int from, to;
  Flow lower, upper, flow;
  Cost cost;
 };
 explicit MinCostFlow(int n= 0): n_(n), b_(n) {}
 int add_vertex() { return b_.push_back(0), n_++; }
 // from から to へ、流量が 0 以上 upper 以下、1 単位あたり cost の辺。番号を返す。
 int add_edge(int from, int to, Flow upper, Cost cost) { return add_edge(from, to, 0, upper, cost); }
 // 流量が lower 以上 upper 以下の辺。
 int add_edge(int from, int to, Flow lower, Flow upper, Cost cost) {
  assert(0 <= from && from < n_ && 0 <= to && to < n_ && lower <= upper);
  es_.push_back({from, to, (i64)lower, (i64)upper, (C)cost});
  return (int)es_.size() - 1;
 }
 // 頂点 v の湧き出しに x を足す (負なら吸い込み)。
 void add_supply(int v, Flow x) { b_[v]+= x; }
 // b-flow を解く。下限と上限と b を満たす流れがあれば true。
 bool solve() { return run(b_); }
 // 下限と b を満たしたうえで、s から t へ limit まで流せるだけ流し、その中で費用を最小にする。流した量を返す。limit までの
 // どの量を流しても下限と b を満たせなければ -1 を返す。
 Flow flow(int s, int t) { return flow(s, t, std::numeric_limits<Flow>::max()); }
 Flow flow(int s, int t, Flow limit) {
  assert(0 <= s && s < n_ && 0 <= t && t < n_ && limit >= 0);
  if(s == t) return run(b_) ? 0 : -1;
  const i64 f= max_extra(s, t, limit);
  if(f < 0) return ok_= false, cost_= Sum(0), pot_.clear(), flow_.clear(), -1;
  std::vector<i64> b(b_);
  b[s]+= f, b[t]-= f;
  run(b);
  return (Flow)f;
 }
 // 直前の solve か flow の結果。
 bool feasible() const { return ok_; }
 Sum cost() const { return cost_; }
 // 相補性を満たすポテンシャル p。縮約費用 cost + p[from] - p[to] が正の辺は下限まで、負の辺は上限まで流れている。
 std::vector<Cost> potential() const { return std::vector<Cost>(pot_.begin(), pot_.end()); }
 Edge get_edge(int i) const {
  const Arc& e= es_[i];
  return {e.s, e.t, (Flow)e.l, (Flow)e.u, (Flow)(ok_ ? flow_[i] : e.l), (Cost)e.c};
 }
 std::vector<Edge> edges() const {
  std::vector<Edge> ret(es_.size());
  for(int i= 0; i < (int)ret.size(); ++i) ret[i]= get_edge(i);
  return ret;
 }
 int size() const { return n_; }
private:
 int n_;
 std::vector<i64> b_;
 std::vector<Arc> es_;
 bool ok_= false;
 Sum cost_= Sum(0);
 std::vector<C> pot_;
 std::vector<i64> flow_;
 bool run(const std::vector<i64>& b) {
  min_cost_flow_internal::Result<C> r;
  if constexpr(!UseNetworkSimplex && !FP) r= min_cost_flow_internal::cost_scaling(n_, b, es_);
  if(UseNetworkSimplex || FP || r.fallback) r= min_cost_flow_internal::network_simplex<C>(n_, b, es_);
  ok_= r.ok, cost_= Sum(0);
  if(!ok_) return pot_.clear(), flow_.clear(), false;
  pot_= std::move(r.pot), flow_= std::move(r.flow);
  if constexpr(FP) {
   long double sum= 0;
   for(int i= 0; i < (int)es_.size(); ++i) sum+= (long double)flow_[i] * es_[i].c;
   cost_= (Sum)sum;
  } else {
   for(int i= 0; i < (int)es_.size(); ++i) cost_+= (__int128)flow_[i] * es_[i].c;
  }
  return true;
 }
 // 下限を流したあとの余り ex から、超頂点 S = n から余りの正の頂点へ、余りの負の頂点から超頂点 T = n + 1 へ辺を張った網。
 MaxFlow<i64> net(const std::vector<i64>& ex, i64& need) const {
  MaxFlow<i64> mf(n_ + 2);
  for(const Arc& e: es_)
   if(e.s != e.t && e.u > e.l) mf.add_edge(e.s, e.t, e.u - e.l);
  need= 0;
  for(int v= 0; v < n_; ++v)
   if(ex[v] > 0) mf.add_edge(n_, v, ex[v]), need+= ex[v];
   else if(ex[v] < 0) mf.add_edge(v, n_ + 1, -ex[v]);
  return mf;
 }
 // 下限と b を満たしたうえで s から t へ流せる量の最大 (limit まで)。どの量でも満たせなければ -1。まず S から T への最大流で、
 // s から t へ流さずに満たせるかを見る。満たせれば、続けて S から s へ、t から T へ容量 limit の辺を足して流し直すと、増えた量
 // が答えになる (S から出る辺と T へ入る辺は、増加路で流しても減らないので、下限と b は満たしたまま)。満たせなければ、t から
 // s へ容量 limit の辺を足した網で満たせる量 f0 を 1 つ見つけ、b を f0 だけずらしてから同じことをする。満たせる量は区間を
 // なすので、これで最大が求まる。
 i64 max_extra(int s, int t, i64 limit) {
  std::vector<i64> ex(b_);
  for(const Arc& e: es_) ex[e.s]-= e.l, ex[e.t]+= e.l;
  i64 sum= 0, need;
  for(int v= 0; v < n_; ++v) sum+= ex[v];
  if(sum != 0) return -1;
  const int S= n_, T= n_ + 1;
  i64 f0= 0;
  {
   MaxFlow<i64> mf= net(ex, need);
   if(need == 0 || mf.flow(S, T) == need) return mf.add_edge(S, s, limit), mf.add_edge(t, T, limit), mf.flow(S, T);
  }
  {
   MaxFlow<i64> mf= net(ex, need);
   const int k= mf.add_edge(t, s, limit);
   if(mf.flow(S, T) != need) return -1;
   f0= mf.get_edge(k).flow;
  }
  ex[s]+= f0, ex[t]-= f0;
  MaxFlow<i64> mf= net(ex, need);
  if(need > 0) mf.flow(S, T);
  return mf.add_edge(S, s, limit - f0), mf.add_edge(t, T, limit - f0), f0 + mf.flow(S, T);
 }
};
