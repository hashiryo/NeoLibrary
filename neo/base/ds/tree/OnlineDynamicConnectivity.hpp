#pragma once
// オンラインの動的連結性 (Holm–de Lichtenberg–Thorup の段つき ETT)。辺の追加と削除、連結性、連結成分の積を
// O(log^2 N) 償却で答える。段ごとに Euler Tour Tree を持ち、森の辺と非木辺を段に振り分ける。段の ETT は
// 64 bit の flag に辺の両端、部分木の頂点数、高い段の辺と非木辺の印を詰める。今の Library の
// OnlineDynamicConnectivity と EulerTourTree の核を写した (記録は algo-notes の notes/dynamic-tree.md)。
// M はモノイド (T、ti、op) か作用 (T、E、mp、cp) を持つ。M= void のときは連結性と辺の切り貼りだけ。
// 段の ETT は型ごとに 1 つの置き場を共有するので、同時に 2 つ以上の OnlineDynamicConnectivity を作らないこと。
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include "neo/base/core/algebra/concepts.hpp"

namespace online_dynamic_connectivity_internal {
template <class M, bool = Monoid<M> || Action<M>>
struct ValueOf {
 using type = std::nullptr_t;
};
template <class M>
struct ValueOf<M, true> {
 using type = typename M::T;
};
template <class M, bool = Action<M>>
struct EOf {
 using type = std::nullptr_t;
};
template <class M>
struct EOf<M, true> {
 using type = typename M::E;
};

template <class M, std::size_t NODE_SIZE>
class EulerTourTree {
 using node_id = int;
 using vertex_id = int;
 static constexpr bool mono = Monoid<M>;
 static constexpr bool dual = Action<M>;
 struct NodeBase {
  node_id ch[2] = {0, 0};
  node_id par = 0;
  std::uint64_t flag = 0;
 };
 template <class, bool MO, bool DU>
 struct NodeType: NodeBase {};
 template <class D>
 struct NodeType<D, true, false>: NodeBase {
  typename M::T val= M::ti(), sum= M::ti();
 };
 template <class D>
 struct NodeType<D, false, true>: NodeBase {
  typename M::E laz{};
  bool lazf= false;
 };
 template <class D>
 struct NodeType<D, true, true>: NodeBase {
  typename M::T val= M::ti(), sum= M::ti();
  typename M::E laz{};
  bool lazf= false;
 };
 using Node = NodeType<void, mono, dual>;
 static inline Node* nd = new Node[NODE_SIZE];
 static inline node_id ni = 1;
 node_id n_st;
 [[gnu::always_inline]] static void update(node_id i) {
  nd[i].flag &= 0xffffffffff00000f;
  nd[i].flag |= ((nd[i].flag >> 44) == ((nd[i].flag >> 24) & 0xfffff)) << 4;
  nd[i].flag &= -11ll, nd[i].flag |= (nd[i].flag << 1) & 0b1111;
  if constexpr (mono) nd[i].sum= nd[i].val;
  if (nd[i].ch[0]) {
   nd[i].flag+= (nd[nd[i].ch[0]].flag & 0xfffff0), nd[i].flag|= nd[nd[i].ch[0]].flag & 0b1010;
   if constexpr (mono) nd[i].sum= M::op(nd[nd[i].ch[0]].sum, nd[i].sum);
  }
  if (nd[i].ch[1]) {
   nd[i].flag+= (nd[nd[i].ch[1]].flag & 0xfffff0), nd[i].flag|= nd[nd[i].ch[1]].flag & 0b1010;
   if constexpr (mono) nd[i].sum= M::op(nd[i].sum, nd[nd[i].ch[1]].sum);
  }
 }
 template <class T, class E>
 static void map(T& v, E x, int sz) {
  if constexpr (std::is_invocable_r_v<void, decltype(M::mp), T&, E, int>) M::mp(v, x, sz);
  else M::mp(v, x);
 }
 static void propagate(node_id i, const typename EOf<M>::type& x) {
  if constexpr (dual) {
   if (nd[i].lazf) M::cp(nd[i].laz, x);
   else nd[i].laz= x, nd[i].lazf= true;
   if ((nd[i].flag >> 44) == ((nd[i].flag >> 24) & 0xfffff)) map(nd[i].val, x, 1);
   if constexpr (mono) map(nd[i].sum, x, int((nd[i].flag >> 4) & 0xfffff));
  }
 }
 static void push(node_id i) {
  if constexpr (dual) {
   if (!nd[i].lazf) return;
   if (nd[i].ch[0]) propagate(nd[i].ch[0], nd[i].laz);
   if (nd[i].ch[1]) propagate(nd[i].ch[1], nd[i].laz);
   nd[i].lazf= false;
  }
 }
 [[gnu::always_inline]] static int dir(node_id i) { return nd[nd[i].par].ch[1] == i; }
 [[gnu::always_inline]] static void rot(node_id i) {
  node_id p= nd[i].par;
  const int d= dir(i);
  if ((nd[p].ch[d]= nd[i].ch[!d])) nd[nd[p].ch[d]].par= p;
  nd[i].ch[!d]= p;
  if ((nd[i].par= nd[p].par)) nd[nd[p].par].ch[dir(p)]= i;
  nd[p].par= i, update(p);
 }
 static void splay(node_id i) {
  if constexpr (dual) push(i);
  for (node_id p= nd[i].par; p; rot(i), p= nd[i].par) {
   node_id pp= nd[p].par;
   if constexpr (dual) {
    if (pp) push(pp);
    push(p), push(i);
   }
   if (pp) rot(dir(i) == dir(p) ? p : i);
  }
  update(i);
 }
 [[gnu::always_inline]] static node_id merge_back(node_id l, node_id r) {
  if (!l) return r;
  if (!r) return l;
  while (nd[l].ch[1]) l= nd[l].ch[1];
  return splay(l), nd[nd[r].par= l].ch[1]= r, update(l), l;
 }
 [[gnu::always_inline]] static std::pair<node_id, node_id> split(node_id i) {
  splay(i);
  const node_id l= nd[i].ch[0];
  return nd[i].ch[0]= nd[l].par= 0, update(i), std::make_pair(l, i);
 }
 [[gnu::always_inline]] static void reroot(node_id v) {
  auto p= split(v);
  merge_back(p.second, p.first), splay(v);
 }
 [[gnu::always_inline]] static bool same_root(node_id i, node_id j) {
  if (i) splay(i);
  if (j) splay(j);
  while (nd[i].par) i= nd[i].par;
  while (nd[j].par) j= nd[j].par;
  return i == j;
 }
 [[gnu::always_inline]] node_id new_edge(int s, int d, bool hi) {
  const int i= ni++, ri= ni++;
  nd[i].flag= (std::uint64_t(s) << 44) | (std::uint64_t(d) << 24) | hi;
  nd[ri].flag= (std::uint64_t(d) << 44) | (std::uint64_t(s) << 24);
  return i;
 }
 std::unordered_map<std::uint64_t, node_id> emp;

 public:
 using T= typename ValueOf<M>::type;
 using E= typename EOf<M>::type;
 EulerTourTree()= default;
 explicit EulerTourTree(int N): n_st(ni) {
  ni+= N;
  for (int i= 0; i < N; ++i) nd[i + n_st].flag= 0x100001000000 * (std::uint64_t)i;
 }
 const T& operator[](vertex_id x) { return get(x); }
 const T& get(vertex_id x) {
  static_assert(mono || dual, "\"get\" is not available");
  return nd[x + n_st].val;
 }
 void set(vertex_id x, T val) {
  static_assert(mono || dual, "\"set\" is not available");
  splay(x+= n_st), nd[x].val= val, update(x);
 }
 bool edge_exist(vertex_id x, vertex_id y) {
  if (x > y) std::swap(x, y);
  return emp.count(((long long)x << 32) | y);
 }
 void link(vertex_id x, vertex_id y, bool hi= true) {
  if (x > y) std::swap(x, y);
  const int ei= new_edge(x, y, hi);
  emp.insert(std::make_pair(((long long)x << 32) | y, ei));
  x+= n_st, y+= n_st, reroot(x), reroot(y);
  nd[nd[x].par= ei].ch[0]= x, nd[nd[y].par= ei].ch[1]= y;
  update(ei), merge_back(ei, ei + 1);
 }
 void cut(vertex_id x, vertex_id y) {
  if (x > y) std::swap(x, y);
  int ei= emp[((long long)x << 32) | y]; const int rei= ei + 1;
  emp.erase(((long long)x << 32) | y);
  auto [pl, pr]= split(ei);
  node_id left, center, right;
  if (pl && same_root(pl, rei)) {
   auto [ql, qr]= split(rei);
   left= ql, center= nd[qr].ch[1], right= nd[pr].ch[1], nd[center].par= 0;
  } else {
   splay(ei), nd[ei= nd[ei].ch[1]].par= 0;
   auto [ql, qr]= split(rei);
   splay(pl), left= pl, right= nd[qr].ch[1];
  }
  nd[right].par= 0, merge_back(left, right);
 }
 bool connected(vertex_id x, vertex_id y) { return same_root(x + n_st, y + n_st); }
 void subedge_set(vertex_id x, bool val) {
  const std::uint64_t one= 1;
  splay(x+= n_st);
  if (val) nd[x].flag|= one << 2;
  else nd[x].flag&= -(one << 2) - 1;
  update(x);
 }
 std::size_t tree_size(vertex_id x) { return splay(x+= n_st), ((nd[x].flag >> 4) & 0xfffff); }
 T prod_tree(vertex_id x) {
  static_assert(mono, "\"prod\" is not available");
  return splay(x+= n_st), nd[x].sum;
 }
 T prod_subtree(vertex_id x, vertex_id par= -1) {
  static_assert(mono, "\"prod\" is not available");
  if (par == -1) return prod_tree(x);
  cut(x, par);
  const T ret= prod_tree(x);
  link(x, par);
  return ret;
 }
 void apply_tree(vertex_id x, E v) {
  static_assert(dual, "\"apply\" is not available");
  splay(x+= n_st), propagate(x, v), push(x);
 }
 void apply_subtree(vertex_id x, vertex_id par, E v) {
  static_assert(dual, "\"apply\" is not available");
  cut(x, par), apply_tree(x, v), link(x, par);
 }
 template <class Func>
 void hilevel_edges(vertex_id v, Func f) {
  splay(v+= n_st);
  while (v && (nd[v].flag & 0b0010))
   while (1) {
    if (nd[v].flag & 0b0001) {
     f((nd[v].flag >> 44), ((nd[v].flag >> 24) & 0xfffff)), splay(v), nd[v].flag&= -2ll, update(v);
     break;
    } else v= nd[v].ch[!(nd[v].ch[0] && (nd[nd[v].ch[0]].flag & 0b0010))];
   }
 }
 template <class Func>
 int subedges(vertex_id v, Func f) {
  splay(v+= n_st);
  while (v && (nd[v].flag & 0b1000))
   for (bool loop= true; loop;) {
    if (nd[v].flag & 0b0100) {
     if (f(vertex_id(nd[v].flag >> 44))) return 1;
     splay(v), loop= false;
    } else v= nd[v].ch[!(nd[v].ch[0] && (nd[nd[v].ch[0]].flag & 0b1000))];
   }
  return 0;
 }
};
}  // namespace online_dynamic_connectivity_internal

template <class M= void, std::size_t NODE_SIZE= 4'000'000>
class OnlineDynamicConnectivity {
 using T= typename online_dynamic_connectivity_internal::ValueOf<M>::type;
 using E= typename online_dynamic_connectivity_internal::EOf<M>::type;
 using Ett= online_dynamic_connectivity_internal::EulerTourTree<M, NODE_SIZE>;
 int N;
 std::vector<Ett> ett;
 std::vector<std::vector<std::unordered_set<int>>> adj;
 void replace(int x, int y, int level) {
  for (int k= 0; k < level; ++k) ett[k].cut(x, y);
  for (int k= level, loop= true; k-- && loop;) {
   if (ett[k].tree_size(x) > ett[k].tree_size(y)) std::swap(x, y);
   ett[k].hilevel_edges(x, [&](int s, int d) { ett[k + 1].link(s, d, true); });
   ett[k].subedges(x, [&](int s) {
    for (auto itr= adj[k][s].begin(), end= adj[k][s].end(); itr != end;) {
     auto d= *itr;
     if (adj[k][s].size() == 1) ett[k].subedge_set(s, 0);
     if (adj[k][d].size() == 1) ett[k].subedge_set(d, 0);
     adj[k][d].erase(s), itr= adj[k][s].erase(itr);
     if (ett[k].connected(s, d)) {
      if (adj[k + 1][s].size() == 0) ett[k + 1].subedge_set(s, 1);
      if (adj[k + 1][d].size() == 0) ett[k + 1].subedge_set(d, 1);
      adj[k + 1][s].insert(d), adj[k + 1][d].insert(s);
     } else {
      for (int kk= k + 1; kk--;) ett[kk].link(s, d, kk == k);
      return loop= false, true;
     }
    }
    return false;
   });
  }
 }

 public:
 explicit OnlineDynamicConnectivity(int n): N(n) { ett.emplace_back(N), adj.emplace_back(N); }
 void link(int x, int y) {
  if (ett[0].connected(x, y)) {
   if (adj[0][x].size() == 0) ett[0].subedge_set(x, 1);
   if (adj[0][y].size() == 0) ett[0].subedge_set(y, 1);
   adj[0][x].insert(y), adj[0][y].insert(x);
  } else ett[0].link(x, y, true);
 }
 void cut(int x, int y) {
  const int m= int(ett.size());
  for (int k= 0; k < m; ++k)
   if (adj[k][x].count(y)) {
    adj[k][x].erase(y), adj[k][y].erase(x);
    if (adj[k][x].size() == 0) ett[k].subedge_set(x, 0);
    if (adj[k][y].size() == 0) ett[k].subedge_set(y, 0);
    return;
   }
  for (int k= m; k--;)
   if (ett[k].edge_exist(x, y)) {
    if (k + 1 == m) ett.emplace_back(N), adj.emplace_back(N);
    replace(x, y, k + 1);
   }
 }
 const T& operator[](int x) { return ett[0][x]; }
 const T& get(int x) { return ett[0].get(x); }
 void set(int x, T val) { ett[0].set(x, val); }
 int size(int x) { return int(ett[0].tree_size(x)); }
 T prod(int x) { return ett[0].prod_tree(x); }
 void apply(int x, E v) { ett[0].apply_tree(x, v); }
 bool connected(int x, int y) { return ett[0].connected(x, y); }
};
