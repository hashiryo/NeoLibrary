#pragma once
#include <algorithm>
#include <array>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>
// 無向グラフの全域最小カット (頂点を 2 つの空でない集合に分けたとき、間を渡る辺の重みの和の最小) の値。頂点 n 個と、無向辺
// {u, v, 重み} の列を受ける。多重辺、自己ループ (使わない)、重み 0 の辺があってもよく、つながっていなければ 0、n が 1 以下
// なら 0 を返す。Cap は整数の型で、各頂点の重み付きの次数が Cap に収まること。
// Nagamochi-Ono-Ibaraki の方法。上界 λ (はじめは重み付きの次数の最小) を持ち、maximum adjacency の順に頂点をなめる。頂点 v を
// なめたとき、隣の未訪問の w への重みの和 r(w) (なめた頂点から w へ) が λ を越えたら、v と w をまとめる (r(w) は v と w の
// 局所連結度の下界で、λ 以上なら、λ より小さいカットは v と w を分けない)。なめた頂点の集合のカットの値でも λ を下げる。1 回
// なめたら、まとめた頂点を縮めて多重辺を足し、頂点が 1 つになるまで繰り返す。1 回のなめで 1 組もまとまらなければ、Stoer-
// Wagner と同じく最後の 2 頂点をまとめる (その組の局所連結度は最後の頂点の次数で、λ 以上)。なめる前に、Padberg-Rinaldi の
// 2 つの条件で縮める。重みが λ 以上の辺と、2 c(u, v) > min(deg u, deg v) の辺 (等号を除くのは、まとめて縮めても最小カットを
// 壊さないため)。
// 今のグラフが疎なうちは、隣接の列と、古い鍵を残す 2 分ヒープでなめる。辺の数 × 20 が頂点の数の 2 乗以上で、頂点が 4096 以下に
// なったら、隣接行列に移り、最大の r を持つ頂点を毎回なめて選ぶ (1 回のなめが O(k^2))。行列は代表の行と列に足して縮める。
// procon-judge の luogu-P5632、hdu-6081、poj-2914、hdu-3691 で書き比べた noi_hyb (中身は problems/_shared/gmincut/noi.hpp と
// noi_dense.hpp) を写したもの。記録は algo-notes の notes/network-flow.md にある。
namespace global_min_cut_internal {
template <class Cap> struct Edge {
 int u, v;
 Cap c;
};
template <class Cap> Cap dense(int n, const std::vector<Edge<Cap>>& edges, Cap lam) {
 std::vector<Cap> w((std::size_t)n * n, 0), deg(n, 0);
 for(auto& [u, v, c]: edges) w[u * n + v]+= c, w[v * n + u]+= c, deg[u]+= c, deg[v]+= c;
 for(int v= 0; v < n; ++v) lam= std::min(lam, deg[v]);
 std::vector<int> par(n), alive(n), keep, gone;
 for(int v= 0; v < n; ++v) par[v]= v, alive[v]= v;
 auto find= [&](int x) {
  while(par[x] != x) x= par[x]= par[par[x]];
  return x;
 };
 auto unite= [&](int a, int b) {
  a= find(a), b= find(b);
  if(a == b) return false;
  par[a]= b;
  return true;
 };
 // まとめられた頂点 x の行を代表の行に足し、代表の行の中で、まとめられた頂点の列を代表の列へ寄せる。代表 A の次数は、A の
 // 頂点の次数の和から、A の中の重みの 2 倍 (寄せたあとの自分の列) を引いたもの。
 auto contract= [&]() {
  keep.clear(), gone.clear();
  for(int v: alive) (find(v) == v ? keep : gone).push_back(v);
  for(int x: gone) {
   const int a= find(x);
   Cap* ra= &w[(std::size_t)a * n];
   const Cap* rx= &w[(std::size_t)x * n];
   for(int y: alive) ra[y]+= rx[y];
   deg[a]+= deg[x];
  }
  for(int a: keep) {
   Cap* ra= &w[(std::size_t)a * n];
   for(int y: gone) ra[find(y)]+= ra[y], ra[y]= 0;
   deg[a]-= ra[a], ra[a]= 0;
  }
  alive.swap(keep);
  if(alive.size() > 1)
   for(int v: alive) lam= std::min(lam, deg[v]);
 };
 std::vector<Cap> r(n);
 std::vector<char> done(n);
 while(alive.size() > 1 && lam > 0) {
  const int k= alive.size();
  bool any= false;
  for(int i= 0; i < k; ++i)
   for(int j= i + 1; j < k; ++j) {
    const int u= alive[i], v= alive[j];
    if(const Cap c= w[(std::size_t)u * n + v]; c > 0 && (c >= lam || 2 * c > std::min(deg[u], deg[v]))) any|= unite(u, v);
   }
  if(any) {
   contract();
   continue;
  }
  for(int v: alive) r[v]= 0, done[v]= 0;
  Cap cut= 0;
  int prev= -1, last= -1;
  for(int step= 0; step < k; ++step) {
   int v= -1;
   for(int x: alive)
    if(!done[x] && (v < 0 || r[x] > r[v])) v= x;
   if(v < 0) break;  // 起きない (未訪問の頂点は残っている)。gcc の -Wstringop-overflow の誤検知を消す
   done[v]= 1, prev= last, last= v;
   cut+= deg[v] - 2 * r[v];
   if(step + 1 < k && cut < lam) lam= cut;
   const Cap* row= &w[(std::size_t)v * n];
   for(int x: alive) {
    if(done[x] || row[x] == 0) continue;
    if(r[x] < lam && r[x] + row[x] >= lam) any|= unite(v, x);
    r[x]+= row[x];
   }
  }
  if(lam == 0) break;
  if(!any) unite(last, prev);
  contract();
 }
 return lam;
}
}  // namespace global_min_cut_internal
template <class Cap> Cap global_min_cut(int n, const std::vector<std::tuple<int, int, Cap>>& edges) {
 if(n <= 1) return 0;
 int k= n;
 using E= global_min_cut_internal::Edge<Cap>;
 std::vector<E> es, ne, srt;  // 今のグラフの辺。u < v で、同じ組はまとめてある
 es.reserve(edges.size());
 for(auto& [u, v, c]: edges)
  if(u != v && c > 0) es.push_back({std::min(u, v), std::max(u, v), c});
 std::vector<int> par(n), id(n), cnt, seen, st, to, pos;
 std::vector<Cap> deg, r, wt;
 std::vector<char> vis;
 auto find= [&](int x) {
  while(par[x] != x) x= par[x]= par[par[x]];
  return x;
 };
 auto unite= [&](int a, int b) {
  a= find(a), b= find(b);
  if(a == b) return false;
  par[a]= b;
  return true;
 };
 // es の頂点を par でまとめた番号に付け替え、同じ組を足して es を作り直す。左の頂点ごとに並べてから、右の頂点の印でまとめる。
 auto rebuild= [&]() {
  int c= 0;
  for(int v= 0; v < k; ++v)
   if(find(v) == v) id[v]= c++;
  ne.clear();
  for(auto& [u, v, w]: es) {
   const int a= id[find(u)], b= id[find(v)];
   if(a != b) ne.push_back({std::min(a, b), std::max(a, b), w});
  }
  k= c;
  cnt.assign(k + 1, 0);
  for(auto& e: ne) ++cnt[e.u + 1];
  for(int i= 0; i < k; ++i) cnt[i + 1]+= cnt[i];
  srt.resize(ne.size());
  for(auto& e: ne) srt[cnt[e.u]++]= e;
  es.clear(), seen.assign(k, -1);
  for(std::size_t i= 0; i < srt.size();) {
   const int a= srt[i].u, base= es.size();
   for(; i < srt.size() && srt[i].u == a; ++i) {
    const int b= srt[i].v;
    if(seen[b] >= base) es[seen[b]].c+= srt[i].c;
    else seen[b]= es.size(), es.push_back(srt[i]);
   }
  }
  for(int v= 0; v < k; ++v) par[v]= v;
 };
 auto build_adj= [&]() {
  deg.assign(k, 0), st.assign(k + 1, 0);
  for(auto& [u, v, c]: es) deg[u]+= c, deg[v]+= c, ++st[u + 1], ++st[v + 1];
  for(int i= 0; i < k; ++i) st[i + 1]+= st[i];
  to.resize(2 * es.size()), wt.resize(2 * es.size()), pos.assign(st.begin(), st.end() - 1);
  for(auto& [u, v, c]: es) to[pos[u]]= v, wt[pos[u]++]= c, to[pos[v]]= u, wt[pos[v]++]= c;
 };
 for(int v= 0; v < n; ++v) par[v]= v;
 // 密かどうかは、隣接の列を組む前に見る。始めから密なら、多重辺もまとめずに行列で解く。
 auto is_dense= [&](int k) { return k <= 4096 && (long long)es.size() * 20 >= (long long)k * k; };
 if(is_dense(n)) return global_min_cut_internal::dense<Cap>(n, es, std::numeric_limits<Cap>::max());
 rebuild();
 Cap lam= std::numeric_limits<Cap>::max();
 std::vector<std::pair<Cap, int>> heap;
 for(;;) {
  if(k <= 1) return lam;
  if(is_dense(k)) return global_min_cut_internal::dense<Cap>(k, es, lam);
  build_adj();
  lam= std::min(lam, *std::min_element(deg.begin(), deg.end()));
  if(lam == 0) return 0;
  bool any= false;
  for(auto& [u, v, c]: es)
   if(c >= lam || 2 * c > std::min(deg[u], deg[v])) any|= unite(u, v);
  if(any) {
   rebuild();
   continue;
  }
  r.assign(k, 0), vis.assign(k, 0);
  // 古い鍵を残す 2 分ヒープ (領域は回をまたいで使い回す)。
  heap.clear();
  for(int v= 0; v < k; ++v) heap.push_back({0, v});
  Cap cut= 0;
  int scanned= 0, prev= -1, last= -1;
  while(!heap.empty()) {
   std::pop_heap(heap.begin(), heap.end());
   const auto [key, v]= heap.back();
   heap.pop_back();
   if(vis[v] || key != r[v]) continue;
   vis[v]= 1, ++scanned, prev= last, last= v;
   cut+= deg[v] - 2 * r[v];
   if(scanned < k && cut < lam) lam= cut;
   for(int i= st[v]; i < st[v + 1]; ++i) {
    const int x= to[i];
    if(vis[x]) continue;
    if(r[x] < lam && r[x] + wt[i] >= lam) any|= unite(v, x);
    r[x]+= wt[i], heap.push_back({r[x], x}), std::push_heap(heap.begin(), heap.end());
   }
  }
  if(lam == 0) return 0;
  if(!any) unite(last, prev);
  rebuild();
 }
}
