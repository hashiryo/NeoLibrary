#pragma once
#include <algorithm>
#include <array>
#include <utility>
#include <vector>
// 一般グラフの最大マッチング。GeneralMatching(n, edges) は、頂点 n 個と無向辺 {u, v} の列を受けて、構築の中で解く。
// 多重辺と自己ループも受ける (自己ループは使わない)。辺は渡した列の添字で扱う。size() は最大マッチングの辺の数、
// mate(v) は頂点 v の相手 (いなければ -1)、matching() は使った辺の番号を昇順に返す。
// 二部グラフと分かっているなら、bipartite_coloring で塗り分けて BipartiteMatching に渡すほうが 2 倍ほど速い。
//
// 中身は次の順に解く。
// 1. Karp-Sipser で初期のマッチングを作る (空いている隣が 1 つだけになった頂点を先に組ませる)。
// 2. 空いている頂点すべてを根にした交互森を、1 本の待ち行列で同時に育てる (花は Edmonds の方法で縮める)。違う木が
//    つながれば増加路なので流し、その 2 本の木はその回の残りで使わない。これを 2 回まで回し、1 本も流せない回が来れば
//    最大なので終わる。増加路の短いグラフはここで終わる。
// 3. 終わらなければ、二部の二重被覆 (頂点 v を左の v と右の v に分け、辺 {u, v} を左 u と右 v、左 v と右 u の 2 本に
//    したもの) の最大マッチングを、今のマッチングの各辺を両向きに使った形から push-relabel で求める。二重被覆の
//    マッチングは元のグラフの道と閉路に分かれるので、隣り合う 2 頂点ずつ組ませて丸め、2 の交互森で最大まで仕上げる。
//    奇閉路の少ないグラフでは、丸めた時点でほぼ最大になり、長い増加路をたどる探索がほとんど要らない。
// 方式は procon-judge の self-flow-general-matching で比べて選んだ (hybrid2)。記録は algo-notes の notes/network-flow.md。
namespace general_matching_internal {
// Karp-Sipser。deg は空いている隣の数 (多重辺は重ねて数える)。1 に落ちた頂点を stk に積む。
inline std::vector<int> karp_sipser(int n, const std::vector<int>& off, const std::vector<int>& adj) {
 std::vector<int> mate(n, -1), deg(n), stk(n);
 int sp= 0;
 for(int v= 0; v < n; ++v)
  if((deg[v]= off[v + 1] - off[v]) == 1) stk[sp++]= v;
 auto take= [&](int a, int b) {
  mate[a]= b, mate[b]= a;
  for(int x: {a, b})
   for(int k= off[x]; k < off[x + 1]; ++k)
    if(const int w= adj[k]; mate[w] < 0 && --deg[w] == 1) stk[sp++]= w;
 };
 auto pick= [&](int v) {
  for(int k= off[v]; k < off[v + 1]; ++k)
   if(const int u= adj[k]; mate[u] < 0) return take(v, u);
 };
 for(int cur= 0;;) {
  while(sp)
   if(const int v= stk[--sp]; mate[v] < 0) pick(v);
  while(cur < n && (mate[cur] >= 0 || deg[cur] == 0)) ++cur;
  if(cur == n) break;
  pick(cur);
 }
 return mate;
}
// mate から始めて、交互森を max_pass 回まで回す。最大になれば true、max_pass 回で終わらなければ false (mate はそこまで
// 大きくしたもの)。z[v] == t なら v はこの回の偶、zo[v] == t なら奇 (花で偶になった頂点は両方が t)。tr[v] は木の根。
// fs[v] は偶の頂点 v から根への交互路の作り方で、奇の頂点の相手として偶になったなら {その奇の頂点を見つけた偶の頂点, -1}、
// 花で偶になったなら花を閉じた辺 {x, y} (x が v の側)。p は花の union-find の親。死んだ木の頂点は、その回のうちに別の木が
// 拾い直してよい (印を付け直す)。
inline bool forest(int n, const std::vector<int>& off, const std::vector<int>& adj, std::vector<int>& mate, int max_pass) {
 std::vector<int> z(n), zo(n), p(n), tr(n), q, stk;
 std::vector<char> alive(n);
 std::vector<std::array<int, 2>> fs(n);
 q.reserve(n);
 int t= 0;
 auto find= [&](int x) {
  int r= x;
  while(z[r] == t && p[r] >= 0) r= p[r];
  while(x != r) {
   const int nx= p[x];
   p[x]= r, x= nx;
  }
  return r;
 };
 // u を v と組ませ、u から根への交互路を反転する。花を通るところでは 2 本に分かれるので、後の 1 本を積んでおく。
 auto rematch= [&](int u, int v) {
  stk.clear();
  for(;;) {
   for(;;) {
    const int w= mate[u];
    mate[u]= v;
    if(w < 0 || mate[w] != u) break;
    if(auto [x, y]= fs[u]; y < 0) mate[w]= x, u= x, v= w;
    else stk.push_back(y), stk.push_back(x), u= x, v= y;
   }
   if(stk.empty()) break;
   v= stk.back(), stk.pop_back(), u= stk.back(), stk.pop_back();
  }
 };
 for(;;) {
  ++t, q.clear();
  for(int v= 0; v < n; ++v)
   if(mate[v] < 0 && off[v] < off[v + 1]) z[v]= t, p[v]= -1, tr[v]= v, alive[v]= 1, fs[v]= {-1, -1}, q.push_back(v);
  int aug= 0;
  for(size_t i= 0; i < q.size(); ++i) {
   const int x= q[i];
   if(z[x] != t || !alive[tr[x]]) continue;
   const int r= tr[x];
   for(int k= off[x]; k < off[x + 1]; ++k) {
    const int y= adj[k];
    const bool labeled= (z[y] == t || zo[y] == t);
    if(labeled && alive[tr[y]]) {
     if(z[y] != t) continue;  // 奇
     if(tr[y] != r) {
      rematch(x, y), rematch(y, x), alive[r]= alive[tr[y]]= 0, ++aug;
      break;
     }
     int u= find(x), v= find(y), w= r;
     if(u == v) continue;
     for(; u != r || v != r; fs[u]= {x, y}, u= find(fs[mate[u]][0])) {
      if(v != r) std::swap(u, v);
      if(fs[u][0] == x && fs[u][1] == y) {
       w= u;
       break;
      }
     }
     for(int a: {find(x), find(y)})
      for(; a != w; a= find(fs[mate[a]][0])) z[a]= t, p[a]= w, q.push_back(a);
    } else {
     const int my= mate[y];
     if(my < 0) continue;  // 起きない (空いている頂点はどれも生きた木の根か、死んだ木で組まれた)
     z[y]= 0, zo[y]= t, tr[y]= r, fs[y]= {-1, -1};
     z[my]= t, p[my]= y, tr[my]= r, fs[my]= {x, -1}, q.push_back(my);
    }
   }
  }
  if(!aug) return true;
  if(t >= max_pass) return false;
 }
}
// 二重被覆の最大マッチング。左 u の相手 (右の頂点、いなければ -1) を返す。もとのグラフのマッチング init の各辺を両向きに
// 使った形から始める。左の u の隣は右の N(u)、右の v の隣は左の N(v) なので、どちらも adj[off[x]..off[x + 1]) になり、
// 二重被覆を組まずに解ける。右の頂点に、空いている右の頂点までの距離の下界 psi を持たせ、空いている左の頂点を FIFO で
// 取り出して、下界が最小の隣と組ませる (double push)。n 回組ませるごとに、空いている右の頂点からの逆向きの BFS で下界を
// 正確な距離に直す (global relabel)。
inline std::vector<int> dcover(int n, const std::vector<int>& off, const std::vector<int>& adj, const std::vector<int>& init) {
 const int INF= 2 * n;
 std::vector<int> ml(n, -1), mr(n, -1), q(n), bq(n);
 for(int u= 0; u < n; ++u)
  if(const int v= init[u]; v >= 0) ml[u]= v, mr[v]= u;
 struct Rv {
  int psi, mate;
 };
 std::vector<Rv> rv(n);
 for(int r= 0; r < n; ++r) rv[r].mate= mr[r];
 auto global_relabel= [&] {
  int h= 0, t= 0;
  for(int r= 0; r < n; ++r)
   if(rv[r].mate < 0) rv[r].psi= 0, bq[t++]= r;
   else rv[r].psi= INF;
  while(h < t) {
   const int r= bq[h++], d= rv[r].psi + 2;
   for(int k= off[r]; k < off[r + 1]; ++k)
    if(const int r2= ml[adj[k]]; r2 >= 0 && rv[r2].psi == INF) rv[r2].psi= d, bq[t++]= r2;
  }
 };
 global_relabel();
 int qh= 0, qn= 0;
 for(int l= 0; l < n; ++l)
  if(ml[l] < 0) q[qn++]= l;
 const int freq= std::max(1, n);
 for(int cnt= 0; qn;) {
  const int u= q[qh];
  if(++qh == n) qh= 0;
  --qn;
  int m1= INF, m2= INF, r1= -1;
  for(int k= off[u]; k < off[u + 1]; ++k) {
   const int r= adj[k], p= rv[r].psi;
   if(p < m1) m2= m1, m1= p, r1= r;
   else if(p < m2) m2= p;
  }
  if(m1 >= INF) continue;
  const int w= rv[r1].mate;
  ml[u]= r1, rv[r1]= {std::min(m2 + 2, INF), u};
  if(w >= 0) {
   ml[w]= -1;
   int t= qh + qn;
   if(t >= n) t-= n;
   q[t]= w, ++qn;
  }
  if(++cnt == freq) cnt= 0, global_relabel();
 }
 return ml;
}
// 二重被覆のマッチング (左 u の相手 out[u]) を道と閉路に分け、隣り合う 2 頂点ずつ組ませたマッチングを返す。
inline std::vector<int> round_dcover(int n, const std::vector<int>& out) {
 std::vector<int> in(n, -1), mate(n, -1);
 for(int u= 0; u < n; ++u)
  if(out[u] >= 0) in[out[u]]= u;
 std::vector<char> seen(n);
 auto walk= [&](int x) {
  while(x >= 0 && !seen[x]) {
   seen[x]= 1;
   const int y= out[x];
   if(y < 0 || seen[y]) break;
   seen[y]= 1, mate[x]= y, mate[y]= x, x= out[y];
  }
 };
 for(int v= 0; v < n; ++v)
  if(in[v] < 0) walk(v);
 for(int v= 0; v < n; ++v) walk(v);
 return mate;
}
}  // namespace general_matching_internal
class GeneralMatching {
public:
 GeneralMatching(int n, const std::vector<std::array<int, 2>>& edges): n(n), es(edges) { solve(); }
 // 最大マッチングの辺の数。
 int size() const { return sz; }
 // 頂点 v の相手。いなければ -1。
 int mate(int v) const { return mt[v]; }
 // 最大マッチングに使う辺の番号 (昇順)。多重辺はどれか 1 本。
 std::vector<int> matching() const {
  std::vector<int> ret;
  ret.reserve(sz);
  std::vector<char> done(n);
  for(int i= 0; i < (int)es.size(); ++i)
   if(auto [u, v]= es[i]; u != v && mt[u] == v && !done[u]) done[u]= done[v]= 1, ret.push_back(i);
  return ret;
 }
private:
 int n, sz= 0;
 std::vector<std::array<int, 2>> es;
 std::vector<int> mt;
 void solve() {
  namespace gi= general_matching_internal;
  // 両向きの CSR。隣は辺の列の順に並べ、自己ループは除く。
  std::vector<int> off(n + 1, 0), adj;
  for(auto [u, v]: es)
   if(u != v) ++off[u], ++off[v];
  for(int i= 0; i < n; ++i) off[i + 1]+= off[i];
  adj.resize(off[n]);
  for(int i= es.size(); i--;)
   if(auto [u, v]= es[i]; u != v) adj[--off[u]]= v, adj[--off[v]]= u;
  mt= gi::karp_sipser(n, off, adj);
  if(!gi::forest(n, off, adj, mt, 2)) {
   mt= gi::round_dcover(n, gi::dcover(n, off, adj, mt));
   gi::forest(n, off, adj, mt, 1 << 30);
  }
  for(int v= 0; v < n; ++v) sz+= mt[v] >= 0;
  sz/= 2;
 }
};
