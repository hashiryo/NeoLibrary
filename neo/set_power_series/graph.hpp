#pragma once
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <utility>
#include <vector>
#include "neo/set_power_series/composite.hpp"
#include "neo/set_power_series/convolve.hpp"
#include "neo/set_power_series/exp.hpp"
#include "neo/set_power_series/inv.hpp"
#include "neo/set_power_series/log.hpp"
#include "neo/set_power_series/power_projection.hpp"
#include "neo/set_power_series/zeta.hpp"
// 無向グラフ (頂点 n 個、多重辺と自己ループも可) の、頂点の部分集合 S ごとの数え上げ。返す列の S 番目は、S が誘導する部分グラフ
// (の辺の部分集合) についての値。辺は es に (u, v) の組で渡す。今の Library の UndirectedGraphSetPowerSeries のメンバ関数を、
// グラフを受ける関数に分けたもの。集合冪級数の演算は sps の関数 (法が奇数の ZMod なら速い核) を使う。
//   edge_num (自己ループは数えない)、connected_component_num、cycle_space_rank、odd_deg_num、selfloop_num: 部分集合ごとの整数の列
//   graph<T>: 2^{辺の数} (辺の部分集合の数)、connected_graph<T>: 連結な辺の部分集合の数、cycle_space_size<T>: 2^{閉路空間の次元}、
//   eulerian_graph<T>: 連結なオイラーグラフ、connected_bipartite_graph<T>: 連結な二部グラフ (辺の部分集合と 2 色の塗り分けの組の
//   半分。T は 2 と 2^{辺の数} で割れること)、tree<T>: 全域木、forest<T>: 全域森、cycle_graph<T>: 閉路 (長さ 3 以上、多重辺の
//   2 本の閉路を含む。T は 2 で割れること)、biconnected_graph<T>: 点双連結、two_edge_connected_graph<T>: 二辺連結、
//   cactus_graph<T>: サボテン、acyclic_orientations<T>: 非巡回な向き付け
//   colorings_using_exactly_k_colors_num<T>: 全体をちょうど k 色で塗る数 (k = 0..n、色の名前を区別しない)
//   chromatic_polynomial<T>: 色多項式の係数、tutte_polynomial<T>(x, y): Tutte 多項式の値
//   all_matching<T>: マッチングの数、k_matching<T>: 辺 k 本のマッチングの数 (k = 0..n/2)
//   only_biconnected(f)、articulation_union(f): 連結なものの重み f と、点双連結な塊の重みを互いに変換する
//   only_2edge_connected(f, n, es)、bridge_union(f, n, es): 連結なものの重みと、二辺連結な塊の重みを互いに変換する (橋の本数を使う)
// 完全マッチングの数は hafnian (neo/set_power_series/hafnian.hpp) を使う。
namespace sps {
namespace sps_graph_internal {
using u64= unsigned long long;
struct G {
 int n, N, m, o;
 std::vector<u64> adj;
 std::vector<int> es;
 G(int n_, const std::vector<std::pair<int, int>>& edges): n(n_), N(1 << n_), m(n_ + (n_ & 1)), o(m / 2), adj(std::size_t(m) * m), es(n_) {
  for(auto [u, v]: edges) {
   assert(0 <= u && u < n && 0 <= v && v < n);
   adj[u * m + v]= (adj[v * m + u]+= 1), es[u]|= 1 << v, es[v]|= 1 << u;
  }
 }
 template <class F> void bfs(int s, const F& f) const {
  for(int t= s, u, j; t;)
   for(f(u= 1 << __builtin_ctz(t)); u;) j= __builtin_ctz(u), t^= 1 << j, u^= 1 << j, u|= es[j] & t;
 }
 std::vector<u64> edge_num() const {
  std::vector<u64> r(N, 0);
  for(int i= n; i--;)
   for(int j= i; j--;) r[(1 << i) | (1 << j)]= adj[i * m + j];
  sps::subset_zeta(r);  // 自己ループは数えない
  return r;
 }
 std::vector<int> connected_component_num() const {
  std::vector<int> r(N, 0);
  for(int s= N; s--;) bfs(s, [&](int) { r[s]++; });
  return r;
 }
 template <class T> std::vector<T> cyc_only() const {
  std::vector<T> cyc(1 << o);
  for(int i= 0; i < o; ++i) {
   int a= i + i, b= a + 1, K= a + 2, I= 1 << i;
   std::vector<T> dp0(std::size_t(K) << i);
   dp0[a]= T(1);
   for(int s= 0; s < I; ++s) {
    T* dp0s= dp0.data() + (s * K);
    for(int u= s | I, S= u, j, j0, j1; S; S^= 1 << j) {
     j= __builtin_ctz(S), j0= j + j, j1= j0 + 1;
     const u64 *A0= adj.data() + (j0 * m), *A1= A0 + m;
     T x0= dp0s[j0], x1= dp0s[j1];
     cyc[u]+= x0 * T(A0[b]) + x1 * T(A1[b]);
     for(int U= I - 1 - s, k, k0, k1; U; U^= 1 << k) {
      k= __builtin_ctz(U), k0= k + k, k1= k0 + 1;
      dp0s[(K << k) + k0]+= x0 * T(A0[k1]) + x1 * T(A1[k1]), dp0s[(K << k) + k1]+= x0 * T(A0[k0]) + x1 * T(A1[k0]);
     }
    }
   }
  }
  return cyc;
 }
 template <class T> std::pair<std::vector<T>, std::vector<T>> cyc_pth() const {
  std::vector<T> cyc(1 << o), pth(1 << o);
  for(int i= 0; i < o; ++i) {
   int a= i + i, b= a + 1, K= a + 2, I= 1 << i;
   std::vector<T> dp0(std::size_t(K) << i), dp1(std::size_t(K) << i);
   dp0[a]= T(1);
   for(int s= 0; s < I; ++s) {
    T *dp0s= dp0.data() + (s * K), *dp1s= dp1.data() + (s * K);
    for(int j= 0; j < K; ++j) dp1s[b]+= dp0s[j];
    for(int u= s | I, S= u, j, j0, j1; S; S^= 1 << j) {
     j= __builtin_ctz(S), j0= j + j, j1= j0 + 1;
     const u64 *A0= adj.data() + (j0 * m), *A1= A0 + m;
     T x0= dp0s[j0], x1= dp0s[j1], y0= dp1s[j0], y1= dp1s[j1];
     cyc[u]+= x0 * T(A0[b]) + x1 * T(A1[b]), pth[u]+= y0 + y1;
     for(int U= I - 1 - s, k, k0, k1; U; U^= 1 << k) {
      k= __builtin_ctz(U), k0= k + k, k1= k0 + 1;
      dp0s[(K << k) + k0]+= x0 * T(A0[k1]) + x1 * T(A1[k1]), dp0s[(K << k) + k1]+= x0 * T(A0[k0]) + x1 * T(A1[k0]);
      dp1s[(K << k) + k0]+= y0 * T(A0[k1]) + y1 * T(A1[k1]), dp1s[(K << k) + k1]+= y0 * T(A0[k0]) + y1 * T(A1[k0]);
     }
    }
   }
  }
  return {cyc, pth};
 }
};
template <class T> T pow2(u64 k) {
 T r= T(1), x= T(2);
 for(; k; k>>= 1, x*= x)
  if(k & 1) r*= x;
 return r;
}
template <class T, class Fn> void transform_articulation(std::vector<T>& f, const Fn& g) {
 const int M= f.size() / 2;
 std::vector<T> tmp(M);
 for(int I= M; I; I>>= 1) {
  for(int t= 0; t < M; t+= I)
   for(int u= I, t2= t << 1; u--;) tmp[t | u]= f[t2 | I | u];
  tmp= g(tmp);
  for(int t= 0; t < M; t+= I)
   for(int u= I, t2= t << 1; u--;) f[t2 | I | u]= tmp[t | u];
 }
}
template <class T, bool B> void transform_bridge(std::vector<T>& f, const G& g) {
 const int M= g.N / 2;
 std::vector<T> tmp(M), tmp2;
 for(int i= g.n, I= M; --i; I>>= 1) {
  for(int t= 0; t < M; t+= I)
   for(int u= I, t2= t << 1; u--;) tmp[t | u]= f[t2 | I | u];
  tmp2.assign(M, T());
  for(int t= 0; t < M; t+= I)
   for(int j= i, J= I, t2= t << 1; J>>= 1, j--;)
    for(int s= J, J2= J * 2; s < I; s+= J2)
     for(int u= s + J; u-- > s;) {
      if constexpr(B) tmp2[t | u]+= f[t2 | u] * T(g.adj[i * g.m + j]);
      else tmp2[t | u]-= f[t2 | u] * T(g.adj[i * g.m + j]);
     }
  tmp= sps::convolve(tmp, sps::exp(tmp2));
  for(int t= 0; t < M; t+= I)
   for(int u= I, t2= t << 1; u--;) f[t2 | I | u]= tmp[t | u];
 }
}
}  // namespace sps_graph_internal
using Edges= std::vector<std::pair<int, int>>;
inline std::vector<unsigned long long> edge_num(int n, const Edges& es) { return sps_graph_internal::G(n, es).edge_num(); }
inline std::vector<int> connected_component_num(int n, const Edges& es) { return sps_graph_internal::G(n, es).connected_component_num(); }
inline std::vector<unsigned long long> cycle_space_rank(int n, const Edges& es) {
 const sps_graph_internal::G g(n, es);
 auto e= g.edge_num();
 auto k= g.connected_component_num();
 for(int s= g.N; s--;) e[s]= e[s] + k[s] - __builtin_popcount(s);
 return e;
}
inline std::vector<unsigned long long> odd_deg_num(int n, const Edges& es) {
 const sps_graph_internal::G g(n, es);
 std::vector<unsigned long long> r(g.N, 0);
 for(int i= n, I= g.N; I>>= 1, i--;)
  for(int t= 0, I2= I << 1; t < g.N; t+= I2)
   for(int u= I, v, j; u--;) {
    unsigned long long cnt= 0;
    for(v= t | u; v; v^= 1 << j) cnt+= g.adj[i * g.m + (j= __builtin_ctz(v))];
    r[t | I | u]+= cnt & 1;  // S に入る頂点 i ごとに、S の中の次数の偶奇を足す
   }
 return r;
}
inline std::vector<unsigned long long> selfloop_num(int n, const Edges& es) {
 const sps_graph_internal::G g(n, es);
 std::vector<unsigned long long> r(g.N, 0);
 for(int i= 0, I= 1; i < n; ++i, I<<= 1)
  for(int u= I; u--;) r[I | u]= r[u] + g.adj[i * g.m + i];
 return r;
}
template <class T> std::vector<T> graph(int n, const Edges& es) {
 auto e= edge_num(n, es);
 std::vector<T> r(e.size());
 for(std::size_t s= 0; s < e.size(); ++s) r[s]= sps_graph_internal::pow2<T>(e[s]);
 return r;
}
template <class T> std::vector<T> cycle_space_size(int n, const Edges& es) {
 auto e= cycle_space_rank(n, es);
 std::vector<T> r(e.size());
 for(std::size_t s= 0; s < e.size(); ++s) r[s]= sps_graph_internal::pow2<T>(e[s]);
 return r;
}
template <class T> std::vector<T> connected_graph(int n, const Edges& es) { return sps::log(graph<T>(n, es)); }
template <class T> std::vector<T> eulerian_graph(int n, const Edges& es) { return sps::log(cycle_space_size<T>(n, es)); }
template <class T> std::vector<T> connected_bipartite_graph(int n, const Edges& es) {
 std::vector<T> tmp= graph<T>(n, es), r(tmp.size(), T(1));
 for(std::size_t s= 0; s < r.size(); ++s) r[s]/= tmp[s];
 r= sps::convolve(r, r);
 for(std::size_t s= 0; s < r.size(); ++s) r[s]*= tmp[s];
 r= sps::log(r);
 for(auto& x: r) x/= T(2);
 return r;
}
template <class T> std::vector<T> tree(int n, const Edges& es) {
 const int N= 1 << n;
 auto e= edge_num(n, es);
 std::vector<T> r= {T(), T(1)};
 if(n == 0) return {T()};
 r.reserve(N);
 for(int I= 2; I < N; I<<= 1) {
  std::vector<T> g(r);
  for(int s= I; --s;) g[s]*= T(e[s | I] - e[s] - e[I]);
  g= sps::exp(g), std::copy(g.begin(), g.end(), std::back_inserter(r));
 }
 return r;
}
template <class T> std::vector<T> forest(int n, const Edges& es) { return sps::exp(tree<T>(n, es)); }
template <class T> std::vector<T> cycle_graph(int n, const Edges& es) {
 const sps_graph_internal::G g(n, es);
 const int N= g.N, m= g.m;
 std::vector<T> r(N, T());
 if(n == 0) return r;
 std::vector<T> dp(std::size_t(N) * std::max(n - 1, 1));
 auto D= [&](int s, int j) -> T& { return dp[std::size_t(s) * std::max(n - 1, 1) + j]; };
 for(int i= n, I= N; I>>= 1, --i;) {
  for(int s= I; --s;) std::fill_n(&D(s, 0), i, T());
  for(int j= i; j--;) D(1 << j, j)= T(g.adj[i * m + j]);
  for(int s= 1; s < I; ++s)
   for(int t= s, j, u, rr, k; t; r[s | I]+= D(s, j) * T(g.adj[j * m + i]))
    for(t^= 1 << (j= __builtin_ctz(t)), u= rr= s ^ (1 << j); u; D(s, j)+= D(rr, k) * T(g.adj[k * m + j])) u^= 1 << (k= __builtin_ctz(u));
 }
 for(int i= n; i--;)
  for(int j= i; j--;) r[(1 << i) | (1 << j)]-= T(g.adj[i * m + j]);
 for(int s= N; --s;) r[s]/= T(2);
 return r;
}
template <class T> std::vector<T> only_biconnected(std::vector<T> f) { return sps_graph_internal::transform_articulation(f, [](const std::vector<T>& x) { return sps::log(x); }), f; }
template <class T> std::vector<T> articulation_union(std::vector<T> f) { return sps_graph_internal::transform_articulation(f, [](const std::vector<T>& x) { return sps::exp(x); }), f; }
template <class T> std::vector<T> only_2edge_connected(std::vector<T> f, int n, const Edges& es) { return sps_graph_internal::transform_bridge<T, false>(f, sps_graph_internal::G(n, es)), f; }
template <class T> std::vector<T> bridge_union(std::vector<T> f, int n, const Edges& es) { return sps_graph_internal::transform_bridge<T, true>(f, sps_graph_internal::G(n, es)), f; }
template <class T> std::vector<T> biconnected_graph(int n, const Edges& es) { return only_biconnected(connected_graph<T>(n, es)); }
template <class T> std::vector<T> two_edge_connected_graph(int n, const Edges& es) { return only_2edge_connected(connected_graph<T>(n, es), n, es); }
template <class T> std::vector<T> cactus_graph(int n, const Edges& es) {
 const sps_graph_internal::G g(n, es);
 auto r= cycle_graph<T>(n, es);
 for(int i= n; i--;)
  for(int j= i; j--;) r[(1 << i) | (1 << j)]+= T(g.adj[i * g.m + j]);
 return articulation_union(r);
}
template <class T> std::vector<T> acyclic_orientations(int n, const Edges& es) {
 auto k= connected_component_num(n, es);
 std::vector<T> f(k.size(), T());
 for(std::size_t s= k.size(); --s;)
  if(k[s] == __builtin_popcount(s)) f[s]= k[s] & 1 ? T() - T(1) : T(1);
 f[0]= T(1);
 return sps::inv(f);
}
template <class T> std::vector<T> colorings_using_exactly_k_colors_num(int n, const Edges& es) {
 if(n == 0) return {T()};
 auto k= connected_component_num(n, es);
 std::vector<T> indep(k.size(), T());
 for(std::size_t s= k.size(); --s;) indep[s]= T(k[s] == __builtin_popcount(s));
 return sps::egf_T(indep);
}
template <class T> std::vector<T> chromatic_polynomial(int n, const Edges& es) {
 auto e= colorings_using_exactly_k_colors_num<T>(n, es);
 if(e.back() == T()) return {T()};
 std::vector<T> r(n + 1, T()), tmp(n);
 tmp[0]= T(1);
 for(int i= 1, j; i < n; ++i)
  for(j= i; j--; tmp[j]*= T() - T(i)) r[j + 1]+= tmp[j] * e[i], tmp[j + 1]+= tmp[j];
 for(int j= n; j--;) r[j + 1]+= tmp[j];
 return r;
}
template <class T> T tutte_polynomial(int n, const Edges& es, T x, T y) {
 const sps_graph_internal::G g(n, es);
 const int N= g.N, m= g.m;
 std::vector<int> sum(N);
 std::vector<T> fum= {T(), T(1)}, gg= {T()}, h;
 gg.reserve(N), h.reserve(N);
 for(int i= 0; i < n; h= sps::exp(h), std::copy(h.begin(), h.end(), std::back_inserter(gg)), ++i) {
  sum[0]= 0;
  for(int j= 0; j < i; j++)
   for(int s= 1 << j, t= s; s--;) sum[s | t]= sum[s] + int(g.adj[i * m + j]);
  h.resize(1 << i);
  for(int s= 1 << i; s--;) {
   while(int(fum.size()) <= sum[s]) fum.push_back(fum.back() * y + T(1));
   h[s]= gg[s] * fum[sum[s]];
  }
 }
 x-= T(1);
 int t= ~0, loops= 0;
 for(int i= n; i--;) loops+= int(g.adj[i * m + i]);
 int s= N;
 g.bfs(s - 1, [&](int u) { t^= u; });
 for(--s&= t; s; --s&= t) gg[s]*= x;
 T yp= T(1);
 for(int i= 0; i < loops; ++i) yp*= y;
 return sps::exp(gg)[N - 1] * yp;
}
template <class T> T all_matching(int n, const Edges& es) {
 auto [cyc, pth]= sps_graph_internal::G(n, es).cyc_pth<T>();
 for(std::size_t s= 0; s < cyc.size(); ++s) cyc[s]+= pth[s];
 return sps::exp(cyc).back();
}
template <class T> std::vector<T> k_matching(int n, const Edges& es) {
 auto [cyc, pth]= sps_graph_internal::G(n, es).cyc_pth<T>();
 auto r= sps::egf_T(pth, sps::exp(cyc));
 std::reverse(r.begin(), r.end()), r.resize(n / 2 + 1);
 return r;
}
}  // namespace sps
