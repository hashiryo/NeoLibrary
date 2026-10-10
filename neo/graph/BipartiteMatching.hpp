#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <numeric>
#include <utility>
#include <vector>
#include "neo/graph/bipartite_coloring.hpp"
// 二部グラフの最大マッチング。中の表し方は「左の頂点 [0, L)、右の頂点 [0, R)、その間の辺」の 1 つで、入口が 2 つある。
// - BipartiteMatching(l, r, edges): 左右の大きさと、辺 (左 a、右 b) の列を渡す。頂点は左 a を a、右 b を l + b の番号で扱う。
// - BipartiteMatching(color, edges): 一般のグラフの辺の列と 2 色の塗り分け (bipartite_coloring で求める) を渡す。頂点は
//   元の番号で扱う。色 0 の頂点を左、色 1 の頂点を右として、番号の付け替えは型の中で持つ。
// 辺はどちらの入口でも、渡した列の添字で扱う。多重辺があってもよい。
// マッチングは構築の中で求める。Karp-Sipser で初期のマッチングを作り、push-relabel (double push、FIFO、(L + R) / 2 回ごとの
// global relabel) で仕上げる。Dulmage-Mendelsohn 分解は、点被覆、独立集合、辺と頂点の分類のどれかを初めて聞かれたときに
// 求めるので、大きさやマッチングだけの使い方では費用がかからない。
// 方式は procon-judge の yosupo-bipartitematching で比べて選んだ (pr_ks_half_pack)。記録は algo-notes の
// notes/network-flow.md と notes/network-flow-api.md。
class BipartiteMatching {
public:
 BipartiteMatching(int l, int r, const std::vector<std::array<int, 2>>& edges): L(l), R(r), es(edges) { solve(); }
 BipartiteMatching(const std::vector<int>& color, const std::vector<std::array<int, 2>>& edges): L(0), R(0), es(edges.size()) {
  const int n= color.size();
  vid.resize(n), org.resize(n);
  for(int v= 0; v < n; ++v) vid[v]= color[v] ? R++ : L++;
  for(int v= 0; v < n; ++v) org[vid[v]+= color[v] ? L : 0]= v;
  for(size_t i= 0; i < edges.size(); ++i) {
   auto [u, v]= edges[i];
   assert(color[u] != color[v]);
   if(color[u]) std::swap(u, v);
   es[i]= {vid[u], vid[v] - L};
  }
  solve();
 }
 // 最大マッチングの辺の数。
 int size() const { return sz; }
 // 頂点 v の相手 (入口の番号)。いなければ -1。
 int mate(int v) const {
  const int x= vid.empty() ? v : vid[v];
  const int y= x < L ? (ml[x] < 0 ? -1 : L + ml[x]) : mr[x - L];
  return y < 0 || org.empty() ? y : org[y];
 }
 // 最大マッチングに使う辺の番号 (昇順)。
 std::vector<int> matching() const {
  std::vector<int> ret;
  ret.reserve(sz);
  std::vector<char> done(L);
  for(int i= 0; i < (int)es.size(); ++i)
   if(auto [a, b]= es[i]; ml[a] == b && !done[a]) done[a]= 1, ret.push_back(i);
  return ret;
 }
 // 最小点被覆。order は頂点の順列 (入口の番号) で、前の頂点ほど被覆に入れる。省けば番号の順。返す頂点は order の順。
 std::vector<int> min_vertex_cover(std::vector<int> order= {}) { return cover(std::move(order), true); }
 // 最大独立集合。order の前の頂点ほど集合に入れる。返す頂点は order の順。
 std::vector<int> max_independent_set(std::vector<int> order= {}) { return cover(std::move(order), false); }
 // 最小辺被覆に使う辺の番号 (昇順)。辺の無い頂点があれば空。
 std::vector<int> min_edge_cover() const {
  std::vector<char> use(es.size()), ok(L + R);
  for(int i= 0; i < (int)es.size(); ++i)
   if(auto [a, b]= es[i]; ml[a] == b && !ok[a]) use[i]= ok[a]= ok[L + b]= 1;
  for(int i= 0; i < (int)es.size(); ++i)
   if(auto [a, b]= es[i]; !ok[a] || !ok[L + b]) use[i]= ok[a]= ok[L + b]= 1;
  std::vector<int> ret;
  for(int x= 0; x < L + R; ++x)
   if(!ok[x]) return ret;
  for(int i= 0; i < (int)es.size(); ++i)
   if(use[i]) ret.push_back(i);
  return ret;
 }
 // Dulmage-Mendelsohn 分解。ブロックの数 K と、頂点ごとのブロックの番号 (入口の番号で引く) を返す。番号はトポロジカル順で、
 // 左 a と右 b を結ぶ辺があれば block[a] <= block[b]。ブロック 0 は、空いた右の頂点から交互路で届く頂点 (右の頂点は
 // どれも空けられる)。ブロック K - 1 は、空いた左の頂点から交互路で届く頂点 (左の頂点はどれも空けられる)。この 2 つは
 // 空のこともある。1 から K - 2 は、残りを強連結成分に分けたもので、どれも左と右の頂点の数が等しい。
 std::pair<int, std::vector<int>> dulmage_mendelsohn() {
  dm();
  if(org.empty()) return {K, blk};
  std::vector<int> ret(L + R);
  for(int x= 0; x < L + R; ++x) ret[org[x]]= blk[x];
  return {K, ret};
 }
 // 辺 i が、どの最大マッチングにも入らないなら 0、入るものと入らないものがあるなら 1、どの最大マッチングにも入るなら 2。
 int edge_kind(int i) {
  dm();
  const auto [a, b]= es[i];
  const int x= blk[a];
  if(x != blk[L + b]) return 0;
  return x != 0 && x != K - 1 && bl[x] == 1 && be[x] == 1 ? 2 : 1;
 }
 // 頂点 v が、どの最大マッチングでも空いているなら 0 (辺が無い)、空くものと空かないものがあるなら 1、どの最大マッチング
 // でも相手がいるなら 2。
 int vertex_kind(int v) {
  dm();
  const int x= vid.empty() ? v : vid[v];
  const int d= x < L ? sl[x + 1] - sl[x] : sr[x - L + 1] - sr[x - L];
  if(d == 0) return 0;
  return blk[x] == (x < L ? K - 1 : 0) ? 1 : 2;
 }
private:
 int L, R, sz= 0, K= -1;
 std::vector<std::array<int, 2>> es;  // 中の番号の辺 (左 a、右 b)
 std::vector<int> vid, org;           // 2 つ目の入口だけ。元の番号から中の番号 (左 a は a、右 b は L + b) へと、その逆
 std::vector<int> sl, al, sr, ar, ml, mr;
 std::vector<int> blk, bl, be;          // 中の番号の頂点のブロック、ブロックの左の頂点の数、ブロックの中の辺の数
 std::vector<int> fst, fto, bst, bto;  // ブロックの DAG (前向きと後ろ向き)
 void solve() {
  const int M= es.size(), INF= 2 * R;
  sl.assign(L + 1, 0), al.resize(M), sr.assign(R + 1, 0), ar.resize(M);
  for(auto& [a, b]: es) ++sl[a], ++sr[b];
  for(int i= 0; i < L; ++i) sl[i + 1]+= sl[i];
  for(int i= 0; i < R; ++i) sr[i + 1]+= sr[i];
  for(int i= M; i--;) al[--sl[es[i][0]]]= es[i][1], ar[--sr[es[i][1]]]= es[i][0];
  ml.assign(L, -1), mr.assign(R, -1);
  {
   // Karp-Sipser。dl, dr は空いている隣の数。次数が 1 に落ちた頂点を stk に積む (右の頂点 r は L + r で積む)。
   // 次数は減るだけで 1 になるのは各頂点で一度なので、stk は L + R で足りる。
   std::vector<int> dl(L), dr(R), stk(L + R);
   int sp= 0;
   for(int l= 0; l < L; ++l)
    if((dl[l]= sl[l + 1] - sl[l]) == 1) stk[sp++]= l;
   for(int r= 0; r < R; ++r)
    if((dr[r]= sr[r + 1] - sr[r]) == 1) stk[sp++]= L + r;
   auto take= [&](int l, int r) {
    ml[l]= r, mr[r]= l;
    for(int k= sl[l]; k < sl[l + 1]; ++k)
     if(const int r2= al[k]; mr[r2] < 0 && --dr[r2] == 1) stk[sp++]= L + r2;
    for(int k= sr[r]; k < sr[r + 1]; ++k)
     if(const int l2= ar[k]; ml[l2] < 0 && --dl[l2] == 1) stk[sp++]= l2;
   };
   for(int cur= 0;;) {
    while(sp) {
     const int x= stk[--sp];
     if(x < L) {
      if(ml[x] >= 0) continue;
      for(int k= sl[x]; k < sl[x + 1]; ++k)
       if(const int r= al[k]; mr[r] < 0) {
        take(x, r);
        break;
       }
     } else {
      const int r= x - L;
      if(mr[r] >= 0) continue;
      for(int k= sr[r]; k < sr[r + 1]; ++k)
       if(const int l= ar[k]; ml[l] < 0) {
        take(l, r);
        break;
       }
     }
    }
    while(cur < L && (ml[cur] >= 0 || dl[cur] == 0)) ++cur;
    if(cur == L) break;
    for(int k= sl[cur]; k < sl[cur + 1]; ++k)
     if(const int r= al[k]; mr[r] < 0) {
      take(cur, r);
      break;
     }
   }
  }
  // push-relabel。右の頂点に、空いている右の頂点までの残余グラフでの距離の下界 psi を持ち、相手と並べて rv に置く。
  // 空いている左の頂点 u を FIFO で取り出し、psi が最小の隣 r と組ませ、r の元の相手を空けて FIFO に戻す。r の psi は、
  // u の隣のうち r を除いた psi の最小に 2 を足した値に上げる。psi が上限 (2R) の右の頂点からは空いている右の頂点に
  // 届かないので、隣がすべて上限の左の頂点は捨てる。(L + R) / 2 回組ませるごとに、空いている右の頂点から逆向きの BFS で
  // psi を正確な距離に直す。
  struct Rv {
   int psi, mate;
  };
  std::vector<Rv> rv(R);
  std::vector<int> q(L), bq(R);
  for(int r= 0; r < R; ++r) rv[r].mate= mr[r];
  auto global_relabel= [&] {
   int h= 0, t= 0;
   for(int r= 0; r < R; ++r)
    if(rv[r].mate < 0) rv[r].psi= 0, bq[t++]= r;
    else rv[r].psi= INF;
   while(h < t) {
    const int r= bq[h++], d= rv[r].psi + 2;
    for(int k= sr[r]; k < sr[r + 1]; ++k)
     if(const int r2= ml[ar[k]]; r2 >= 0 && rv[r2].psi == INF) rv[r2].psi= d, bq[t++]= r2;
   }
  };
  global_relabel();
  // 空いている左の頂点は同時に L 個を超えないので、長さ L の輪で FIFO を持つ。
  int qh= 0, qn= 0;
  for(int l= 0; l < L; ++l)
   if(ml[l] < 0) q[qn++]= l;
  const int freq= std::max(1, (L + R) / 2);
  for(int cnt= 0; qn;) {
   const int u= q[qh];
   if(++qh == L) qh= 0;
   --qn;
   int m1= INF, m2= INF, r1= -1;
   for(int k= sl[u]; k < sl[u + 1]; ++k) {
    const int r= al[k], p= rv[r].psi;
    if(p < m1) m2= m1, m1= p, r1= r;
    else if(p < m2) m2= p;
   }
   if(m1 >= INF) continue;
   const int w= rv[r1].mate;
   ml[u]= r1, rv[r1]= {std::min(m2 + 2, INF), u};
   if(w >= 0) {
    ml[w]= -1;
    int t= qh + qn;
    if(t >= L) t-= L;
    q[t]= w, ++qn;
   }
   if(++cnt == freq) cnt= 0, global_relabel();
  }
  for(int r= 0; r < R; ++r) mr[r]= rv[r].mate;
  for(int l= 0; l < L; ++l) sz+= ml[l] >= 0;
 }
 void dm() {
  if(K >= 0) return;
  const int n= L + R;
  blk.assign(n, -1);
  std::vector<int> q(n);
  int h= 0, t= 0;
  // ブロック 0: 空いている右の頂点から、右から左へはどの辺でも、左から右へはマッチングの辺で届く頂点。届いた左の頂点は
  // 空いていない (空いていれば増加路がある)。
  for(int b= 0; b < R; ++b)
   if(mr[b] < 0) blk[L + b]= 0, q[t++]= b;
  while(h < t) {
   const int b= q[h++];
   for(int k= sr[b]; k < sr[b + 1]; ++k)
    if(const int a= ar[k]; blk[a] == -1) {
     blk[a]= 0;
     if(const int b2= ml[a]; blk[L + b2] == -1) blk[L + b2]= 0, q[t++]= b2;
    }
  }
  // ブロック K - 1: 空いている左の頂点から、左から右へはどの辺でも、右から左へはマッチングの辺で届く頂点。仮に -2 を付ける。
  h= t= 0;
  for(int a= 0; a < L; ++a)
   if(ml[a] < 0) blk[a]= -2, q[t++]= a;
  while(h < t) {
   const int a= q[h++];
   for(int k= sl[a]; k < sl[a + 1]; ++k)
    if(const int b= al[k]; blk[L + b] == -1) {
     blk[L + b]= -2;
     if(const int a2= mr[b]; blk[a2] == -1) blk[a2]= -2, q[t++]= a2;
    }
  }
  // 残りの左の頂点を、辺 (a, b) ごとに a から b の相手へ向かう有向グラフの強連結成分に分ける (再帰しない Tarjan)。
  // 成分は閉じた順に番号を振るので、後から閉じた成分ほどトポロジカル順で前にある。
  std::vector<int> ord(L, -1), low(L), comp(L, -1), it(L), stk, cs;
  int cnt= 0, C= 0;
  for(int s= 0; s < L; ++s) {
   if(blk[s] != -1 || ord[s] >= 0) continue;
   ord[s]= low[s]= cnt++, it[s]= sl[s], stk.push_back(s), cs.push_back(s);
   while(!cs.empty()) {
    const int v= cs.back();
    if(it[v] < sl[v + 1]) {
     const int w= mr[al[it[v]++]];
     if(blk[w] != -1) continue;  // ブロック K - 1 の頂点。ブロック 0 の右の頂点には、残りの左の頂点から辺が無い
     if(ord[w] < 0) ord[w]= low[w]= cnt++, it[w]= sl[w], stk.push_back(w), cs.push_back(w);
     else if(comp[w] < 0) low[v]= std::min(low[v], ord[w]);
    } else {
     cs.pop_back();
     if(!cs.empty()) low[cs.back()]= std::min(low[cs.back()], low[v]);
     if(low[v] == ord[v]) {
      int x;
      do x= stk.back(), stk.pop_back(), comp[x]= C;
      while(x != v);
      ++C;
     }
    }
   }
  }
  K= C + 2;
  for(int a= 0; a < L; ++a)
   if(blk[a] == -2) blk[a]= K - 1;
   else if(blk[a] == -1) blk[a]= C - comp[a];
  for(int b= 0; b < R; ++b)
   if(blk[L + b] == -2) blk[L + b]= K - 1;
   else if(blk[L + b] == -1) blk[L + b]= blk[mr[b]];
  bl.assign(K, 0), be.assign(K, 0);
  for(int a= 0; a < L; ++a) ++bl[blk[a]];
  for(auto& [a, b]: es)
   if(blk[a] == blk[L + b]) ++be[blk[a]];
 }
 // ブロックの DAG。ブロック 0 と K - 1 は被覆の側が決まっているので、そこに触れる辺は持たない。
 void dag() {
  if(!fst.empty()) return;
  fst.assign(K + 1, 0), bst.assign(K + 1, 0);
  auto inner= [&](int x, int y) { return x != y && x != 0 && y != K - 1; };
  for(auto& [a, b]: es)
   if(const int x= blk[a], y= blk[L + b]; inner(x, y)) ++fst[x + 1], ++bst[y + 1];
  for(int k= 0; k < K; ++k) fst[k + 1]+= fst[k], bst[k + 1]+= bst[k];
  fto.resize(fst[K]), bto.resize(bst[K]);
  std::vector<int> pf(fst.begin(), fst.end() - 1), pb(bst.begin(), bst.end() - 1);
  for(auto& [a, b]: es)
   if(const int x= blk[a], y= blk[L + b]; inner(x, y)) fto[pf[x]++]= y, bto[pb[y]++]= x;
 }
 // order の順に、まだ側が決まっていないブロックの側を決める。z[k] が 1 ならブロック k の左の頂点を、0 なら右の頂点を
 // 被覆に入れる。辺 (a, b) が x = block[a] < y = block[b] を結ぶとき、z[x] = 0 なら z[y] = 0、z[y] = 1 なら z[x] = 1 で
 // なければならないので、0 は DAG の前向きに、1 は後ろ向きに広げる。in なら被覆の頂点を、そうでなければ残りを返す。
 std::vector<int> cover(std::vector<int> order, bool in) {
  dm(), dag();
  const int n= L + R;
  if(order.empty()) order.resize(n), std::iota(order.begin(), order.end(), 0);
  assert((int)order.size() == n);
  std::vector<int> z(K, -1), q(K), ret;
  z[0]= 1, z[K - 1]= 0;
  for(int v: order) {
   const int x= vid.empty() ? v : vid[v], c= x >= L, k= blk[x];
   if(z[k] < 0) {
    const int s= in ? !c : c;
    const std::vector<int>&st= s ? bst : fst, &to= s ? bto : fto;
    int h= 0, t= 0;
    z[k]= s, q[t++]= k;
    while(h < t)
     for(int u= q[h++], j= st[u]; j < st[u + 1]; ++j)
      if(z[to[j]] < 0) z[to[j]]= s, q[t++]= to[j];
   }
   if((c != z[k]) == in) ret.push_back(v);
  }
  return ret;
 }
};
