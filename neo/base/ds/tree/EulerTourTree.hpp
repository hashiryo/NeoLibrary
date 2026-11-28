#pragma once
// 森の辺を切り貼りしながら、連結成分の積と、辺 (v, p) について p を親とした v の部分木の積を答える (Euler Tour Tree)。
// 木を、頂点の節と辺ごとの 2 つの向きの節を並べた巡回列で持ち、splay 木で分ける。部分木は p->v の節から
// v->p の節までで、その区間に作用を遅延で当てる。M はモノイド (using T、static T ti()、static T op(T, T)) で、
// 作用は区間の頂点の数を受ける mp(T&, E, int) と cp(E&, E) を持つもの (部分木への加算のような使い方)。
// 記録は algo-notes の notes/dynamic-tree.md。
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>
#include <vector>
#include "neo/base/core/algebra/concepts.hpp"
#ifdef __linux__
#include <sys/mman.h>
#endif

namespace euler_tour_tree_internal {
// 2 の冪の大きさの開番地法の表。V は trivially copyable で、値 0 が空の印になること。
// 埋まりが半分を超えたら広げ、消した跡は墓標にして広げ直すときに落とす。キーに 0 は使えない。
template <class V>
class OpenTable {
 struct Slot {
  unsigned long long k;
  V v;
 };
 std::vector<Slot> t;
 std::size_t mask, shift, live, used;
 static constexpr unsigned long long TOMB = ~0ull;
 std::size_t home(unsigned long long k) const { return std::size_t((k * 0x9E3779B97F4A7C15ull) >> shift); }
 void put(unsigned long long k, const V& v) {
  std::size_t h = home(k);
  while (t[h].k && t[h].k != TOMB) h = (h + 1) & mask;
  if (!t[h].k) ++used;
  t[h].k = k, t[h].v = v, ++live;
 }
 void rehash(std::size_t cap, std::size_t new_shift) {
  std::vector<Slot> old;
  old.swap(t);
  t.assign(cap, Slot{});
  mask = cap - 1, shift = new_shift, live = 0, used = 0;
  for (const Slot& s : old)
   if (s.k && s.k != TOMB) put(s.k, s.v);
 }

 public:
 OpenTable(): t(16), mask(15), shift(60), live(0), used(0) {}
 void insert(unsigned long long k, const V& v) {
  if ((used + 1) * 2 > mask + 1) {
   if ((live + 1) * 2 <= mask + 1) rehash(mask + 1, shift);
   else rehash((mask + 1) * 2, shift - 1);
  }
  put(k, v);
 }
 // 無ければ nullptr。
 const V* find(unsigned long long k) const {
  std::size_t h = home(k);
  while (t[h].k) {
   if (t[h].k == k) return &t[h].v;
   h = (h + 1) & mask;
  }
  return nullptr;
 }
 bool contains(unsigned long long k) const { return find(k) != nullptr; }
 void erase(unsigned long long k) {
  std::size_t h = home(k);
  while (t[h].k != k) h = (h + 1) & mask;
  t[h].k = TOMB, --live;
 }
};
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
  euler_tour_tree_internal::OpenTable<Edge> es;
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
  // x を、親が goal になるまで回す。goal が 0 なら木の根まで。
  [[gnu::always_inline]] void splay_to(int x, int goal) {
   for (int p; (p = t[x].par) != goal;) {
    int g = t[p].par;
    if (g == goal) {
     push(p), push(x), rotate(x);
     break;
    }
    push(g), push(p), push(x);
    if ((t[g].ch[1] == p) == (t[p].ch[1] == x)) rotate(p), rotate(x);
    else rotate(x), rotate(x);
   }
   push(x);
  }
  [[gnu::always_inline]] void splay(int x) { splay_to(x, 0); }
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
  // 列の [l, r) から平衡な splay 木を組む。根を返す。
  int build_balanced(const std::vector<int>& seq, int l, int r) {
   if (l >= r) return 0;
   const int m = (l + r) >> 1, x = seq[m];
   const int L = build_balanced(seq, l, m), R = build_balanced(seq, m + 1, r);
   t[x].ch[0] = L, t[x].ch[1] = R;
   if (L) t[L].par = x;
   if (R) t[R].par = x;
   t[x].par = 0;
   pull(x);
   return x;
  }
  std::pair<int, int> edge_nodes(int v, int p) const {
   const Edge* found = es.find(key(v, p));
   assert(found);
   const Edge& e = *found;
   if (e.u == v) return {e.a, e.b};
   return {e.b, e.a};
  }

 public:
  explicit EulerTourTree(int n_): t(std::size_t(3 * n_ + 8)), n(n_) {
   for (int v = 1; v <= n; ++v) t[v].cnt = 1;
  }
 // 辺の列 (森でもよい) と頂点の値から O(N) で組む。頂点の節は 1..N、辺の節は N+1 から。
 EulerTourTree(int n_, const std::vector<std::array<int, 2>>& edges, const std::vector<T>& a) {
  n = n_;
  const int E = int(edges.size()), base = n + 1;
  t.resize(std::size_t(n + 2 * E + 2));
  for (int v = 1; v <= n; ++v) t[v].cnt = 1, t[v].val = a[v - 1], t[v].sum = a[v - 1];
  std::vector<int> head(n, -1), to(2 * E), nxt(2 * E), eid(2 * E), dir(2 * E);
  int ec = 0;
  auto add_arc = [&](int u, int v, int e, int d) { to[ec] = v, nxt[ec] = head[u], eid[ec] = e, dir[ec] = d, head[u] = ec++; };
  for (int e = 0; e < E; ++e) {
   const int u = edges[e][0], v = edges[e][1];
   add_arc(u, v, e, 0), add_arc(v, u, e, 1);
   es.insert(key(u, v), Edge{u, v, base + 2 * e, base + 2 * e + 1});
  }
  std::vector<char> seen(n, 0);
  std::vector<int> seq;
  seq.reserve(std::size_t(n + 2 * E));
  std::vector<int> su, sp, se, snode;
  for (int root = 0; root < n; ++root) {
   if (seen[root]) continue;
   const int start = int(seq.size());
   seen[root] = 1, seq.push_back(root + 1);
   su.push_back(root), sp.push_back(-1), se.push_back(head[root]), snode.push_back(-1);
   while (!su.empty()) {
    const int u = su.back(), p = sp.back();
    int& e = se.back();
    int v = -1, ee = -1, dd = -1;
    while (e != -1) {
     const int idx = e;
     e = nxt[idx];
     if (to[idx] != p) {
      v = to[idx], ee = eid[idx], dd = dir[idx];
      break;
     }
    }
    if (v == -1) {
     if (const int in = snode.back(); in != -1) seq.push_back(in);
     su.pop_back(), sp.pop_back(), se.pop_back(), snode.pop_back();
    } else {
     seen[v] = 1;
     const int node = base + 2 * ee + dd;
     seq.push_back(node), seq.push_back(v + 1);
     su.push_back(v), sp.push_back(u), se.push_back(head[v]), snode.push_back(base + 2 * ee + (dd ^ 1));
    }
   }
   build_balanced(seq, start, int(seq.size()));
  }
 }
  int size() const { return n; }
  bool edge_exists(int u, int v) const { return es.contains(key(u, v)); }
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
   es.insert(key(u, v), Edge{u, v, a, b});
  }
  void cut(int u, int v) {
   const Edge* hit = es.find(key(u, v));
   assert(hit);
   Edge e = *hit;
   es.erase(key(u, v));
   int a, b;
   if (e.u == u) a = e.a, b = e.b;
   else a = e.b, b = e.a;
   reroot(a);  // 列を a から始め、a の次に v の側、その次に b が来る形にする
   splay(b);
   splay_to(a, b);
   const bool right = (t[b].ch[1] == a);
   if (right) {  // [L] b [v の側] a [Q]
    int l = t[b].ch[0], q = t[a].ch[1], v1 = t[a].ch[0];
    if (v1) t[v1].par = 0;
    if (q) t[q].par = 0;
    if (l) t[l].par = 0;
    t[a].ch[0] = t[a].ch[1] = 0;
    t[b].ch[0] = 0, t[b].ch[1] = a;  // {b, a} だけ残して捨てる
    static_cast<void>(merge(q, l));
   } else {  // [X] a [Y] b [R]: v の側は R + X、u の側は Y
    int x = t[a].ch[0], y = t[a].ch[1], rr = t[b].ch[1];
    if (x) t[x].par = 0;
    if (y) t[y].par = 0;
    if (rr) t[rr].par = 0;
    t[a].ch[0] = t[a].ch[1] = 0;
    t[b].ch[1] = 0, t[b].ch[0] = a;
    static_cast<void>(merge(rr, x));
    static_cast<void>(y);
   }
   free_edge_node(a), free_edge_node(b);
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
   splay(b);
   splay_to(a, b);
   const bool right = (t[b].ch[1] == a);
   if (right) {  // [L] b [v の側] a [Q]
    int l = t[a].ch[0];
    return l ? t[l].sum : M::ti();
   }
   // [X] a [Y] b [R]: v の側は R + X
   T s = M::ti();
   if (int x = t[b].ch[1]; x) s = t[x].sum;
   if (int x = t[a].ch[0]; x) s = M::op(s, t[x].sum);
   return s;
  }
  // 辺 (v, p) について、p を親とした v の部分木への作用。
  template <bool A = AC>
  void apply_subtree(int v, int p, const E& e) {
   static_assert(A, "\"apply_subtree\" is not available");
   if constexpr (A) {
    auto [a, b] = edge_nodes(v, p);  // a: v->p, b: p->v
    splay(b);
    int r = a;
    while (t[r].par != b) r = t[r].par;
    const bool right = (t[b].ch[1] == r);
    splay_to(a, b);
    if (right) {
     apply_to(t[a].ch[0], e);
    } else {
     apply_to(t[b].ch[1], e), apply_to(t[a].ch[0], e);
    }
    pull(a), pull(b);
   }
  }
};
