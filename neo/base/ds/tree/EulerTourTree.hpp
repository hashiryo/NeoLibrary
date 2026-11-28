#pragma once
// 森の辺を切り貼りしながら、連結成分の積と、辺 (v, p) について p を親とした v の部分木の積を答える (Euler Tour Tree)。
// 木を、頂点の節と辺ごとの 2 つの向きの節を並べた巡回列で持ち、splay 木で分ける。部分木は p->v の節から
// v->p の節までで、その区間に作用を遅延で当てる。M はモノイド (using T、static T ti()、static T op(T, T)) で、
// 作用は区間の頂点の数を受ける mp(T&, E, int) と cp(E&, E) を持つもの (部分木への加算のような使い方)。
// 記録は algo-notes の notes/dynamic-tree.md。
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>
#include "neo/base/core/algebra/concepts.hpp"
#ifdef __linux__
#include <sys/mman.h>
#endif

namespace euler_tour_tree_internal {
struct Empty {};
template <class M, bool = Action<M>>
struct ActionOf {
 struct type {};
};
template <class M>
struct ActionOf<M, true> {
 using type = typename M::E;
};
// 作用の mp が区間の長さを受ける形。
template <class M, bool = Action<M>>
struct ActionLenV: std::false_type {};
template <class M>
struct ActionLenV<M, true>: std::integral_constant<bool, requires(typename M::T& v, const typename M::E& e, int len) { M::mp(v, e, len); }> {};
template <class M>
concept ActionLen = ActionLenV<M>::value;
// 1 MB 以上の確保は 2 MB 境界に寄せ、Linux では huge page を頼む。
template <class T>
struct Alloc {
 using value_type = T;
 Alloc() = default;
 template <class U>
 Alloc(const Alloc<U>&) {}
 T* allocate(std::size_t n) {
  const std::size_t want = n * sizeof(T);
  if (want < (std::size_t(1) << 20)) return static_cast<T*>(::operator new(want, std::align_val_t(64)));
  const std::size_t H = std::size_t(1) << 21, bytes = (want + H - 1) & ~(H - 1);
  T* p = static_cast<T*>(::operator new(bytes, std::align_val_t(H)));
#ifdef __linux__
  madvise(p, bytes, MADV_HUGEPAGE);
#endif
  return p;
 }
 void deallocate(T* p, std::size_t n) {
  if (n * sizeof(T) < (std::size_t(1) << 20)) ::operator delete(p, std::align_val_t(64));
  else ::operator delete(p, std::align_val_t(std::size_t(1) << 21));
 }
 template <class U>
 bool operator==(const Alloc<U>&) const {
  return true;
 }
 template <class U>
 void construct(U* p) noexcept {
  ::new (static_cast<void*>(p)) U;
 }
 template <class U, class... A>
 void construct(U* p, A&&... a) {
  ::new (static_cast<void*>(p)) U(std::forward<A>(a)...);
 }
};
template <class T>
using Vec = std::vector<T, Alloc<T>>;

template <class M>
struct Node {
 using T = typename M::T;
 using E = typename ActionOf<M>::type;
 static constexpr bool AC = ActionLen<M>;
 int ch[2] = {0, 0};
 int par = 0;
 int cnt = 0;  // splay の部分木にある頂点の節の数
 T val = M::ti(), sum = M::ti();
 [[no_unique_address]] std::conditional_t<AC, E, Empty> laz;
 [[no_unique_address]] std::conditional_t<AC, bool, Empty> lazf;
 Node() {
  if constexpr (AC) lazf = false;
 }
};
}  // namespace euler_tour_tree_internal

template <class M>
class EulerTourTree {
 public:
  // a は u->v の向きの節、b は v->u の向きの節。
  struct Edge {
   int u, v, a, b;
  };

 private:
  using T = typename M::T;
  using E = typename euler_tour_tree_internal::ActionOf<M>::type;
  using Node = euler_tour_tree_internal::Node<M>;
  static constexpr bool AC = euler_tour_tree_internal::ActionLen<M>;
  euler_tour_tree_internal::Vec<Node> t;
  std::vector<int> free_nodes;
  std::unordered_map<unsigned long long, Edge> es;
  int n;
  static unsigned long long key(int u, int v) {
   if (u > v) std::swap(u, v);
   return (unsigned long long)(unsigned)u << 32 | (unsigned)v;
  }
  [[gnu::always_inline]] void pull(int x) {
   Node& a = t[x];
   int l = a.ch[0], r = a.ch[1];
   a.cnt = (l ? t[l].cnt : 0) + (x <= n ? 1 : 0) + (r ? t[r].cnt : 0);
   T s = a.val;
   if (l) s = M::op(t[l].sum, s);
   if (r) s = M::op(s, t[r].sum);
   a.sum = s;
  }
  template <bool A = AC>
  [[gnu::always_inline]] void apply_to(int x, const E& e) {
   if constexpr (A) {
    if (!x) return;
    Node& a = t[x];
    if (x <= n) M::mp(a.val, e, 1);
    M::mp(a.sum, e, a.cnt);
    if (a.lazf) M::cp(a.laz, e);
    else a.laz = e, a.lazf = true;
   }
  }
  [[gnu::always_inline]] void push(int x) {
   if constexpr (AC)
    if (t[x].lazf) apply_to(t[x].ch[0], t[x].laz), apply_to(t[x].ch[1], t[x].laz), t[x].lazf = false;
  }
  [[gnu::always_inline]] void rotate(int x) {
   Node& a = t[x];
   int p = a.par, g = t[p].par, d = t[p].ch[1] == x, b = a.ch[d ^ 1];
   t[p].ch[d] = b;
   if (b) t[b].par = p;
   a.ch[d ^ 1] = p, t[p].par = x, a.par = g;
   if (g) t[g].ch[t[g].ch[1] == p] = x;
   a.cnt = t[p].cnt, a.sum = t[p].sum;
   pull(p);
  }
  [[gnu::always_inline]] void splay(int x) {
   for (int p; (p = t[x].par);) {
    int g = t[p].par;
    if (!g) {
     push(p), push(x), rotate(x);
     break;
    }
    push(g), push(p), push(x);
    if ((t[g].ch[1] == p) == (t[p].ch[1] == x)) rotate(p), rotate(x);
    else rotate(x), rotate(x);
   }
   push(x);
  }
  // x より前を切り離す。x の木は x から始まる。
  [[gnu::always_inline]] int split_before(int x) {
   splay(x);
   int l = t[x].ch[0];
   if (l) t[l].par = 0;
   t[x].ch[0] = 0;
   pull(x);
   return l;
  }
  [[gnu::always_inline]] int split_after(int x) {
   splay(x);
   int r = t[x].ch[1];
   if (r) t[r].par = 0;
   t[x].ch[1] = 0;
   pull(x);
   return r;
  }
  // a の列の後ろに b の列を継ぐ。a と b は別の木。
  [[gnu::always_inline]] int merge(int a, int b) {
   if (!a) return b;
   if (!b) return a;
   int x = a;
   push(x);
   while (t[x].ch[1]) x = t[x].ch[1], push(x);
   splay(x);
   t[x].ch[1] = b, t[b].par = x;
   pull(x);
   return x;
  }
  // 巡回列を x から始まる形に回す。
  [[gnu::always_inline]] int reroot(int x) {
   int l = split_before(x);
   return merge(x, l);
  }
  int root_of(int x) {
   splay(x);
   int y = x;
   push(y);
   while (t[y].ch[0]) y = t[y].ch[0], push(y);
   splay(y);
   return y;
  }
  [[gnu::always_inline]] void reset(int x) {
   Node& a = t[x];
   a.ch[0] = a.ch[1] = a.par = 0, a.cnt = 0, a.val = M::ti(), a.sum = M::ti();
   if constexpr (AC) a.lazf = false;
  }
  int alloc_edge_node() {
   if (!free_nodes.empty()) {
    int x = free_nodes.back();
    free_nodes.pop_back();
    reset(x);
    return x;
   }
   int x = int(t.size());
   t.emplace_back();
   reset(x);
   return x;
  }
  void free_edge_node(int x) {
   reset(x);
   free_nodes.push_back(x);
  }
  // v から p へ向かう節と p から v へ向かう節。辺は (u,v) の向きで覚えてある。
  std::pair<int, int> edge_nodes(int v, int p) const {
   auto it = es.find(key(v, p));
   assert(it != es.end());
   const Edge& e = it->second;
   if (e.u == v) return {e.a, e.b};
   return {e.b, e.a};
  }

 public:
  explicit EulerTourTree(int n_): t(std::size_t(3 * n_ + 8)), n(n_) {
   for (int v = 1; v <= n; ++v) t[v].cnt = 1;
  }
  int size() const { return n; }
  bool edge_exists(int u, int v) const { return es.count(key(u, v)); }
  T get(int v) {
   splay(v + 1);
   return t[v + 1].val;
  }
  void set(int v, T x) {
   int y = v + 1;
   splay(y);
   t[y].val = x;
   pull(y);
  }
  void link(int u, int v) {
   assert(u != v && !edge_exists(u, v));
   int a = alloc_edge_node(), b = alloc_edge_node();
   int ru = reroot(u + 1), rv = reroot(v + 1);
   int x = merge(merge(merge(ru, a), rv), b);
   static_cast<void>(x);
   es.emplace(key(u, v), Edge{u, v, a, b});
  }
  void cut(int u, int v) {
   auto it = es.find(key(u, v));
   assert(it != es.end());
   Edge e = it->second;
   es.erase(it);
   int a, b;
   if (e.u == u) a = e.a, b = e.b;
   else a = e.b, b = e.a;
   reroot(a);  // 列を a から始め、a の次に v の側、その次に b が来る形にする
   int rest = split_after(a);
   static_cast<void>(rest);
   splay(b);
   int mid = split_before(b);  // v の側
   int suf = split_after(b);   // u の側 (a より前の部分は回したので空)
   static_cast<void>(mid);
   free_edge_node(a), free_edge_node(b);
   static_cast<void>(suf);
  }
  bool connected(int u, int v) { return root_of(u + 1) == root_of(v + 1); }
  // 頂点 v の木の全体の積。
  T prod_tree(int v) {
   int x = v + 1;
   splay(x);
   return t[x].sum;
  }
  // 木の全体の頂点の数。
  int tree_size(int v) {
   int x = v + 1;
   splay(x);
   return t[x].cnt;
  }
  // 木の全体への作用。
  template <bool A = AC>
  void apply_tree(int v, const E& e) {
   static_assert(A, "\"apply_tree\" is not available");
   if constexpr (A) splay(v + 1), apply_to(v + 1, e);
  }
  // 辺 (v, p) について、p を親とした v の部分木の積。v と p は辺で結ばれていること。
  T prod_subtree(int v, int p) {
   auto [a, b] = edge_nodes(v, p);  // a: v->p, b: p->v
   reroot(b), splay(a);             // p->v から前方が v の側
   int l = t[a].ch[0];
   return l ? t[l].sum : M::ti();
  }
  // 辺 (v, p) について、p を親とした v の部分木への作用。
  template <bool A = AC>
  void apply_subtree(int v, int p, const E& e) {
   static_assert(A, "\"apply_subtree\" is not available");
   if constexpr (A) {
    auto [a, b] = edge_nodes(v, p);  // a: v->p, b: p->v
    reroot(b), splay(a);
    apply_to(t[a].ch[0], e);
    pull(a);
   }
  }
};
