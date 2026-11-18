#pragma once
#include <algorithm>
#include <cassert>
#include <limits>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
// 無向グラフの Gomory-Hu 木 (全点対の最小カットを表す木)。gomory_hu(n, edges) は、頂点 0 を根とする木を、各頂点の親
// parent と、辺 (v, parent[v]) の重み weight で返す。どの 2 頂点 u, v についても、u-v の最小カットの値は木の u-v の道の
// 上の重みの最小に等しく (min_cut(u, v))、木の辺 (v, parent[v]) を切って分かれる v の側の頂点集合のカットの容量は
// weight[v] に等しい。辺は {u, v, 容量} の列で、多重辺、自己ループ (使わない)、容量 0 の辺を受け、連結でなくてもよい
// (つながっていない組の最小カットは 0)。
//
// 作り方は Gusfield の方法 (グラフを縮めずに最大流を n - 1 回流す)。処理済みの頂点も付け替え、t の親が s の側なら s を
// 間に入れる形 (後ろの頂点だけを付け替える簡単な形は、全組の値は合うが、木の辺で切った集合が最小カットにならないことが
// ある)。頂点を処理する順は、同じ親を持つ未処理の頂点が続くように選ぶ (どの順でも正しく作れる)。
//
// 最大流は、網を 1 回だけ CSR に組み、flow ごとに変えた弧の容量だけを戻す push-relabel (highest label、同じ高さの中は
// 後入れ先出し、gap、global relabel) の前流の段で流す。最初の高さには、元のグラフで t から測った距離を使う。s の弧を
// 飽和させて s を高さ n に置くと、残りの残余の弧はどれも元のグラフにある弧なので、条件 h[u] ≤ h[w] + 1 を満たす。t からの
// BFS は途中で止めて状態を残し、余りのある頂点 (s の隣) に距離が付くところまでだけ進め、t が同じ間は続きから使う (印の
// 無い頂点の高さは、待ち行列の先頭の距離に 1 を足したもの)。Gusfield の t は、カットが {s} のままなら同じ頂点が続く。
// global relabel の BFS も、余りのある頂点に全部高さが付いた段で止める。
// カットは、流量が s の重み付きの次数に等しければ {s}、t の次数に等しければ t 以外の全部、どちらでもなければ s の側の
// 最小の集合 (余りを s へ戻してから s から BFS) にする。
// 方式は procon-judge の self-flow-gomory-hu、loj-2042、luogu-P4897、luogu-P3329 で比べて選んだ。記録は algo-notes の
// notes/network-flow.md。
// Cap は 64 bit までの整数の型。各頂点に接する辺の容量の和が Cap と long long に収まること (最小カットの値はこれを超えない)。
namespace gomory_hu_internal {
using i64= long long;
class Flow {
 int n, m2= 0, s_= 0, t_= 0, cut_= 0, amax= -1, bmax= -1, dist_t= -1, dqh= 0, dqt= 0;
 std::vector<int> st, to, rv, h, cur, anext, bnext, bprev, ahead, bhead, q, touched, dist, dq;
 std::vector<i64> cap0, cap, ex, deg;
 std::vector<char> mark, in_s;
 void touch(int k) {
  if(!mark[k]) mark[k]= 1, touched.push_back(k);
 }
 void b_insert(int v, int d) {
  bprev[v]= -1, bnext[v]= bhead[d];
  if(bhead[d] >= 0) bprev[bhead[d]]= v;
  bhead[d]= v;
 }
 void b_erase(int v, int d) {
  if(bprev[v] >= 0) bnext[bprev[v]]= bnext[v];
  else bhead[d]= bnext[v];
  if(bnext[v] >= 0) bprev[bnext[v]]= bprev[v];
 }
 void a_push(int v, int d) { anext[v]= ahead[d], ahead[d]= v, amax= std::max(amax, d); }
 // 高さを sink からの距離に付け直す。excl は高さ n に置いて使わない。余りのある頂点に全部高さが付いた段で BFS を止め、
 // まだ付いていない頂点を次の段の高さにする (止めた段より手前の頂点の隣は BFS が付け終えているので、条件が崩れない)。
 void relabel_all(int sink, int excl) {
  std::fill(h.begin(), h.end(), -1);
  for(int d= 0; d <= n; ++d) ahead[d]= bhead[d]= -1;
  amax= bmax= -1;
  int need= 0;
  for(int v= 0; v < n; ++v) need+= ex[v] > 0 && v != sink && v != excl;
  int qh= 0, qt= 0, stop= need == 0 ? 0 : n;  // この高さの頂点を取り出したら止める
  h[sink]= 0, q[qt++]= sink, h[excl]= n;
  while(qh < qt) {
   const int v= q[qh];
   if(h[v] >= stop) break;
   ++qh;
   const int d= h[v] + 1;
   for(int k= st[v]; k < st[v + 1]; ++k)
    if(const int w= to[k]; h[w] < 0 && cap[rv[k]] > 0) {
     h[w]= d, q[qt++]= w;
     if(ex[w] > 0 && --need == 0) stop= d;
    }
  }
  const int rest= stop < n ? stop + 1 : n;
  for(int v= 0; v < n; ++v) {
   if(h[v] < 0) h[v]= rest;
   if(v != sink && h[v] < n) {
    b_insert(v, h[v]), bmax= std::max(bmax, h[v]), cur[v]= st[v];
    if(ex[v] > 0) a_push(v, h[v]);
   }
  }
 }
 // 今の高さのまま、余りを sink へ押し出す (excl には押さない)。sink に届かない余りは残る。relabel で見た弧の数が
 // 6n + 弧の数を超えるたびに高さを付け直す。
 void discharge(int sink, int excl) {
  const long long limit= 6LL * n + m2;
  long long work= 0;
  while(amax >= 0) {
   const int u= ahead[amax];
   if(u < 0) {
    --amax;
    continue;
   }
   ahead[amax]= anext[u];
   while(ex[u] > 0) {
    const int hu= h[u], end= st[u + 1];
    int k= cur[u];
    for(; k < end; ++k)
     if(const int w= to[k]; cap[k] > 0 && h[w] == hu - 1) {
      const i64 d= std::min(ex[u], cap[k]);
      if(ex[w] == 0 && w != sink) a_push(w, h[w]);
      touch(k), cap[k]-= d, cap[rv[k]]+= d, ex[u]-= d, ex[w]+= d;
      if(ex[u] == 0) break;
     }
    cur[u]= k;
    if(ex[u] == 0) break;
    b_erase(u, hu);
    if(bhead[hu] < 0) {  // gap: hu より高い頂点は sink に届かない
     for(int d= hu; d <= bmax; ++d) {
      for(int v= bhead[d]; v >= 0; v= bnext[v]) h[v]= n;
      bhead[d]= ahead[d]= -1;
     }
     h[u]= n, bmax= hu - 1, amax= std::min(amax, hu - 1);
     break;
    }
    int nh= n, nk= end;
    for(int j= st[u]; j < end; ++j)
     if(cap[j] > 0 && h[to[j]] + 1 < nh) nh= h[to[j]] + 1, nk= j;
    work+= end - st[u] + 12;
    h[u]= nh, cur[u]= nk;
    if(nh >= n) break;
    b_insert(u, nh), bmax= std::max(bmax, nh);
    if(work > limit) {
     work= 0;
     relabel_all(sink, excl);
     break;
    }
   }
  }
 }
public:
 template <class Cap> Flow(int n_, const std::vector<std::tuple<int, int, Cap>>& edges): n(n_), st(n_ + 1, 0), h(n_), cur(n_), anext(n_), bnext(n_), bprev(n_), ahead(n_ + 2), bhead(n_ + 2), q(n_), dist(n_, -1), dq(n_), ex(n_, 0), deg(n_, 0), in_s(n_, 0) {
  for(const auto& [u, v, c]: edges) {
   assert(0 <= u && u < n && 0 <= v && v < n);
   if constexpr(std::is_signed_v<Cap>) assert(c >= 0);
   if(u != v && c > 0) ++st[u + 1], ++st[v + 1];
  }
  for(int v= 0; v < n; ++v) st[v + 1]+= st[v];
  m2= st[n];
  to.resize(m2), rv.resize(m2), cap0.resize(m2), mark.assign(m2, 0);
  std::vector<int> pos(st.begin(), st.end() - 1);
  for(const auto& [u, v, c]: edges)
   if(u != v && c > 0) {
    const int a= pos[u]++, b= pos[v]++;
    to[a]= v, to[b]= u, rv[a]= b, rv[b]= a, cap0[a]= cap0[b]= (i64)c;
    deg[u]+= (i64)c, deg[v]+= (i64)c;
   }
  cap= cap0;
 }
 // s から t への最大流の値を返し、side(v) が使う最小カットを決める。
 i64 flow(int s, int t) {
  s_= s, t_= t;
  for(int k: touched) cap[k]= cap0[k], cap[rv[k]]= cap0[rv[k]], mark[k]= 0;
  touched.clear();
  std::fill(ex.begin(), ex.end(), 0);
  for(int k= st[s]; k < st[s + 1]; ++k) {
   const i64 c= cap[k];
   touch(k), ex[to[k]]+= c, cap[rv[k]]+= c, cap[k]= 0;
  }
  // 最初の高さ。元のグラフの t からの BFS を、s の隣に距離が付くまで進める (要る頂点に -2 を付けて数える)。
  if(dist_t != t) std::fill(dist.begin(), dist.end(), -1), dist[t]= 0, dq[0]= t, dqh= 0, dqt= 1, dist_t= t;
  int need= 0;
  for(int k= st[s]; k < st[s + 1]; ++k)
   if(dist[to[k]] == -1) dist[to[k]]= -2, ++need;
  while(need > 0 && dqh < dqt) {
   const int v= dq[dqh++];
   for(int k= st[v]; k < st[v + 1]; ++k)
    if(const int w= to[k]; dist[w] < 0) need-= dist[w] == -2, dist[w]= dist[v] + 1, dq[dqt++]= w;
  }
  const int rest= dqh < dqt ? dist[dq[dqh]] + 1 : n;
  for(int d= 0; d <= n; ++d) ahead[d]= bhead[d]= -1;
  amax= bmax= -1;
  for(int v= 0; v < n; ++v) {
   h[v]= v == s ? n : dist[v] >= 0 ? dist[v] : rest;
   if(v != t && h[v] < n) {
    b_insert(v, h[v]), bmax= std::max(bmax, h[v]), cur[v]= st[v];
    if(ex[v] > 0) a_push(v, h[v]);
   }
  }
  discharge(t, s);
  const i64 f= ex[t];
  if(f == deg[s]) return cut_= 1, f;
  if(f == deg[t]) return cut_= 2, f;
  // 余りを s へ戻し、s から届く頂点を求める。
  ex[t]= 0;
  relabel_all(s, t);
  discharge(s, t);
  std::fill(in_s.begin(), in_s.end(), 0);
  int qh= 0, qt= 0;
  in_s[s]= 1, q[qt++]= s;
  while(qh < qt) {
   const int u= q[qh++];
   for(int k= st[u]; k < st[u + 1]; ++k)
    if(cap[k] > 0 && !in_s[to[k]]) in_s[to[k]]= 1, q[qt++]= to[k];
  }
  return cut_= 0, f;
 }
 // 直前の flow の最小カットで、v が s の側にあるか。
 bool side(int v) const { return cut_ == 0 ? in_s[v] : cut_ == 1 ? v == s_ : v != t_; }
};
}  // namespace gomory_hu_internal
template <class Cap> struct GomoryHuTree {
 std::vector<int> parent;  // 根は 0 で parent[0] = -1
 std::vector<Cap> weight;  // weight[v] は辺 (v, parent[v]) の重み (= v と parent[v] の最小カット)。weight[0] = 0
 std::vector<int> depth;   // 根からの辺の数
 // u と v の最小カットの値 (木の道の上の重みの最小)。木の高さに比例する手間。u == v なら Cap の最大値。
 [[gnu::always_inline]] Cap min_cut(int u, int v) const {
  assert(0 <= u && u < (int)parent.size() && 0 <= v && v < (int)parent.size());
  Cap r= std::numeric_limits<Cap>::max();
  while(u != v) {
   if(depth[u] < depth[v]) std::swap(u, v);
   r= std::min(r, weight[u]), u= parent[u];
  }
  return r;
 }
};
template <class Cap> GomoryHuTree<Cap> gomory_hu(int n, const std::vector<std::tuple<int, int, Cap>>& edges) {
 static_assert(std::is_integral_v<Cap> && sizeof(Cap) <= 8, "容量は 64 bit までの整数の型");
 GomoryHuTree<Cap> g;
 std::vector<int>& p= g.parent;
 p.assign(n, 0), g.weight.assign(n, 0), g.depth.assign(n, 0);
 if(n == 0) return g;
 p[0]= -1;
 if(n > 1) {
  gomory_hu_internal::Flow e(n, edges);
  std::vector<long long> w(n, 0);
  // 未処理の頂点を、親ごとの双方向リストに持つ。今の親 t の未処理の子を先に処理し尽くし、途中で子を持った頂点は後に回す。
  std::vector<int> head(n, -1), nxt(n, -1), prv(n, -1), pend{0};
  std::vector<char> done(n, 0);
  auto link= [&](int v, int t) {
   prv[v]= -1, nxt[v]= head[t];
   if(head[t] >= 0) prv[head[t]]= v;
   head[t]= v;
  };
  auto unlink= [&](int v, int t) {
   if(prv[v] >= 0) nxt[prv[v]]= nxt[v];
   else head[t]= nxt[v];
   if(nxt[v] >= 0) prv[nxt[v]]= prv[v];
  };
  for(int v= n - 1; v >= 1; --v) link(v, 0);
  while(!pend.empty()) {
   const int t= pend.back();
   if(head[t] < 0) {
    pend.pop_back();
    continue;
   }
   const int s= head[t];
   unlink(s, t), done[s]= 1;
   const long long f= e.flow(s, t);
   w[s]= f;
   bool got= false;
   for(int i= 1; i < n; ++i)
    if(i != s && p[i] == t && e.side(i)) {
     p[i]= s;
     if(!done[i]) unlink(i, t), link(i, s), got= true;
    }
   if(t != 0 && e.side(p[t])) p[s]= p[t], p[t]= s, w[s]= w[t], w[t]= f;
   if(got) pend.insert(pend.end() - 1, s);
  }
  for(int v= 0; v < n; ++v) g.weight[v]= (Cap)w[v];
  // 深さ。根へ上がって深さの分かっている頂点まで行き、戻りながら付ける。
  std::vector<int> path;
  std::vector<char> known(n, 0);
  known[0]= 1;
  for(int v= 0; v < n; ++v) {
   int x= v;
   for(; !known[x]; x= p[x]) path.push_back(x);
   for(int d= g.depth[x]; !path.empty(); path.pop_back()) g.depth[path.back()]= ++d, known[path.back()]= 1;
  }
 }
 return g;
}
