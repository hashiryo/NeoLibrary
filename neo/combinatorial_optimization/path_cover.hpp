#pragma once
#include <array>
#include <cassert>
#include <vector>
#include "neo/combinatorial_optimization/BipartiteMatching.hpp"
#include "neo/combinatorial_optimization/MaxFlow.hpp"
// DAG の最小パス被覆と、到達可能性の順序の鎖と反鎖 (Dilworth の定理)。頂点は 0 から n - 1 の番号で、edges は有向辺
// {u, v} (u から v) の列。閉路と自己ループが無いこと。多重辺はあってもよい。
// - min_path_cover(n, edges): 頂点を共有しない道で全部の頂点を覆う、本数が最小の道の列。道は頂点の列で、隣り合う頂点は
//   辺で結ばれる。頂点を出る側と入る側に分けた二部グラフの最大マッチング (BipartiteMatching) から作る。本数は
//   n - (マッチングの大きさ)。
// - min_chain_cover(n, edges): u から v へ道があるという順序 (到達可能性) で、全部の頂点を本数が最小の鎖に分ける。鎖は頂点の
//   列で、前の頂点から後の頂点へ道がある (辺で結ばれているとは限らない)。どの頂点もちょうど 1 つの鎖に入る。
// - max_antichain(n, edges): 同じ順序の最大反鎖 (どの 2 つの間にも道が無い頂点の集合)。大きさは min_chain_cover の本数と
//   等しい。
// 鎖と反鎖は、推移閉包を作らずに MaxFlow を 1 回流して求める。網は、頂点 v を出る側 v' と入る側 v'' に分け、s から v' へ
// 容量 1、辺 {u, v} ごとに u' から v'' へ、v'' から v' へ (v を通り抜ける) 容量 n、v'' から t へ容量 1 の辺を張る。最大流の
// 値は推移閉包の上の最大マッチングの大きさになる。流れを道に分けたとき、u' から入って w'' から出る 1 単位が、鎖の中で u の
// 次が w であることを表す。反鎖は最小カットで、v' が s の側、v'' が t の側にある頂点 v の集合。
// 関係が推移的 (u から v へ道があれば辺 {u, v} もある) なら、3 つの本数は同じで、min_path_cover の道はそのまま鎖になる。
// 比べられる組を全部辺で渡すなら min_path_cover のほうが 2〜3 倍速く、Hasse 図の辺だけを渡せるなら min_chain_cover と
// max_antichain に渡すほうが速い。鎖が多く深い疎な DAG では MaxFlow が遅い (手元の M2 で、頂点 20 万、辺 50 万に 0.15〜2.5 s)。
// 記録は algo-notes の notes/network-flow.md。
inline std::vector<std::vector<int>> min_path_cover(int n, const std::vector<std::array<int, 2>>& edges) {
 for([[maybe_unused]] auto [u, v]: edges) assert(0 <= u && u < n && 0 <= v && v < n && u != v);
 const BipartiteMatching bm(n, n, edges);
 std::vector<int> nx(n, -1);
 std::vector<char> head(n, 1);
 for(int u= 0; u < n; ++u)
  if(const int w= bm.mate(u); w >= 0) nx[u]= w - n, head[w - n]= 0;
 std::vector<std::vector<int>> ret;
 ret.reserve(n - bm.size());
 for(int v= 0; v < n; ++v)
  if(head[v]) {
   auto& p= ret.emplace_back();
   for(int x= v; x >= 0; x= nx[x]) p.push_back(x);
  }
 return ret;
}
namespace path_cover_internal {
// 頂点 v を出る側を v、入る側を n + v、s を 2n、t を 2n + 1 とする網に流す。辺の番号は、s から u が u、v'' から t が n + v、
// v'' から v' が 2n + v、辺 i が 3n + i。
inline MaxFlow<int> chain_flow(int n, const std::vector<std::array<int, 2>>& edges) {
 for([[maybe_unused]] auto [u, v]: edges) assert(0 <= u && u < n && 0 <= v && v < n && u != v);
 MaxFlow<int> mf(2 * n + 2);
 const int s= 2 * n, t= 2 * n + 1;
 for(int u= 0; u < n; ++u) mf.add_edge(s, u, 1);
 for(int v= 0; v < n; ++v) mf.add_edge(n + v, t, 1);
 for(int v= 0; v < n; ++v) mf.add_edge(n + v, v, n);
 for(auto [u, v]: edges) mf.add_edge(u, n + v, n);
 mf.flow(s, t);
 return mf;
}
}  // namespace path_cover_internal
inline std::vector<std::vector<int>> min_chain_cover(int n, const std::vector<std::array<int, 2>>& edges) {
 if(n == 0) return {};
 MaxFlow<int> mf= path_cover_internal::chain_flow(n, edges);
 const int m= edges.size();
 // 流れのある辺を始点ごとに並べ、流れを 1 単位ずつ s から t へたどって、u の次の w を決める。網に閉路は無いので、たどる道は
 // 必ず t に着く。v'' に着いたとき、t へ流れが残っていればそこで止め、残っていなければ v' へ通り抜ける。
 std::vector<int> st(n + 1, 0), to, rem, pass(n), sink(n);
 for(int v= 0; v < n; ++v) pass[v]= mf.get_edge(2 * n + v).flow, sink[v]= mf.get_edge(n + v).flow;
 std::vector<int> f(m);
 for(int i= 0; i < m; ++i)
  if((f[i]= mf.get_edge(3 * n + i).flow) > 0) ++st[edges[i][0] + 1];
 for(int v= 0; v < n; ++v) st[v + 1]+= st[v];
 to.resize(st[n]), rem.resize(st[n]);
 {
  std::vector<int> p(st.begin(), st.end() - 1);
  for(int i= 0; i < m; ++i)
   if(f[i] > 0) to[p[edges[i][0]]]= edges[i][1], rem[p[edges[i][0]]++]= f[i];
 }
 std::vector<int> nx(n, -1), it(st.begin(), st.end() - 1);
 std::vector<char> head(n, 1);
 for(int u= 0; u < n; ++u) {
  if(mf.get_edge(u).flow == 0) continue;
  for(int x= u;;) {
   while(rem[it[x]] == 0) ++it[x];
   const int y= to[it[x]];
   --rem[it[x]];
   if(sink[y] > 0) {
    --sink[y], nx[u]= y, head[y]= 0;
    break;
   }
   --pass[y], x= y;
  }
 }
 std::vector<std::vector<int>> ret;
 for(int v= 0; v < n; ++v)
  if(head[v]) {
   auto& p= ret.emplace_back();
   for(int x= v; x >= 0; x= nx[x]) p.push_back(x);
  }
 return ret;
}
inline std::vector<int> max_antichain(int n, const std::vector<std::array<int, 2>>& edges) {
 if(n == 0) return {};
 MaxFlow<int> mf= path_cover_internal::chain_flow(n, edges);
 const std::vector<bool> vis= mf.min_cut(2 * n);
 std::vector<int> ret;
 for(int v= 0; v < n; ++v)
  if(vis[v] && !vis[n + v]) ret.push_back(v);
 return ret;
}
