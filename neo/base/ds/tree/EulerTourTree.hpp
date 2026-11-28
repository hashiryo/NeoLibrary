#pragma once
// 森の辺を切り貼りしながら、連結成分の積と、辺 (v, p) について p を親とした v の部分木の積を答える (Euler Tour Tree)。
// 木を、頂点の節と辺ごとの 2 つの向きの節を並べた巡回列で持ち、大きさで釣り合いを取る木 (WBT) で分ける。優先度は
// 持たず、継ぐときは大きい側の木へ降りる。部分木は p->v の節から v->p の節までで、その区間に作用を遅延で当てる。
// M はモノイド (using T、static T ti()、static T op(T, T)) で、作用は区間の頂点の数を受ける mp(T&, E, int) と
// cp(E&, E) を持つもの (部分木への加算のような使い方)。記録は algo-notes の notes/dynamic-tree.md。
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
struct Empty {};
// 作用があるときだけ遅延の作用素を持つ。旗は節に直接置く (同じ型の空の節を 2 つ並べられるように分けてある)。
template <class M, bool AC>
struct Extra {};
template <class M>
struct Extra<M, true> {
 typename ActionOf<M>::type laz{};
};
template <class M>
struct Node {
 using T = typename M::T;
 using E = typename ActionOf<M>::type;
 static constexpr bool AC = ActionLen<M>;
 T val = M::ti(), sum = M::ti();
 [[no_unique_address]] Extra<M, AC> ex;  // 8 byte ぞろえなので int の前に置く (あとに置くと 4 byte の隙間ができる)
 int l = 0, r = 0, par = 0;
 int size = 0;  // 部分木の要素の数
 int cnt = 0;   // 部分木の頂点の節の数
 [[no_unique_address]] std::conditional_t<AC, bool, Empty> laz_f{};
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
   a.size = t[a.l].size + 1 + t[a.r].size;
   a.cnt = t[a.l].cnt + (x <= n ? 1 : 0) + t[a.r].cnt;
   a.sum = M::op(M::op(t[a.l].sum, a.val), t[a.r].sum);
  }
  template <bool A = AC>
  [[gnu::always_inline]] void apply_to(int x, const E& e) {
   if constexpr (A) {
    if (!x) return;
    Node& a = t[x];
    if (x <= n) M::mp(a.val, e, 1);
    M::mp(a.sum, e, a.cnt);
    if (a.laz_f) M::cp(a.ex.laz, e);
    else a.ex.laz = e, a.laz_f = true;
   }
  }
  [[gnu::always_inline]] void push(int x) {
   if constexpr (AC)
    if (t[x].laz_f) apply_to(t[x].l, t[x].ex.laz), apply_to(t[x].r, t[x].ex.laz), t[x].laz_f = false;
  }
  // 根まで登って、上から遅延を配る。
  void push_path(int x) {
   if (const int p = t[x].par; p) push_path(p);
   push(x);
  }
  // 二つの列を継ぐ。a と b はそれぞれの木の根で、大きい側の木へ降りる。釣り合いは、降りなかった側の子と
  // 継いだ側を比べ、継いだ側が 4 倍を超えたら継いだ側の根を持ち上げる回転 1 回で直す (列の WBT の join と同じ)。
  int merge(int a, int b) {
   if (!a) return b;
   if (!b) return a;
   if (t[a].size > t[b].size) {
    push(a);
    const int r = merge(t[a].r, b);
    t[a].r = r, t[r].par = a;  // b は空でないので r は必ず非 0
    pull(a);
    if (t[a].l && t[t[a].l].size * 4 >= t[r].size) return a;
    push(r);
    const int b2 = t[r].l;
    t[a].r = b2;
    if (b2) t[b2].par = a;
    t[r].par = t[a].par, t[a].par = r;
    pull(a), t[r].l = a, pull(r);
    return r;
   }
   push(b);
   const int l = merge(a, t[b].l);
   t[b].l = l, t[l].par = b;
   pull(b);
   if (t[b].r && t[t[b].r].size * 4 >= t[l].size) return b;
   push(l);
   const int b2 = t[l].r;
   t[b].l = b2;
   if (b2) t[b2].par = b;
   t[l].par = t[b].par, t[b].par = l;
   pull(b), t[l].r = b, pull(l);
   return l;
  }
  // 列 x を、先頭の k 個と残りに分ける。
  void split(int x, int k, int& a, int& b) {
   if (!x) return a = 0, b = 0, void();
   push(x);
   const int ls = t[t[x].l].size;
   if (k <= ls) {
    split(t[x].l, k, a, t[x].l);
    if (t[x].l) t[t[x].l].par = x;
    t[x].par = 0;
    b = x;
    pull(x);
   } else {
    split(t[x].r, k - ls - 1, t[x].r, b);
    if (t[x].r) t[t[x].r].par = x;
    t[x].par = 0;
    a = x;
    pull(x);
   }
  }
  // 列の先頭から x より前の要素の数。
  int rank(int x) const {
   int r = t[t[x].l].size;
   for (int y = x; t[y].par;) {
    const int p = t[y].par;
    if (t[p].r == y) r += t[t[p].l].size + 1;
    y = p;
   }
   return r;
  }
  int root_of(int x) const {
   while (t[x].par) x = t[x].par;
   return x;
  }
  // 部分木 x の中の位置 [l, r) の積。
  T range_sum(int x, int l, int r) {
   if (!x || l >= r) return M::ti();
   if (l == 0 && r == t[x].size) return t[x].sum;
   push(x);
   const int ls = t[t[x].l].size;
   T s = M::ti();
   if (l < ls) s = range_sum(t[x].l, l, std::min(r, ls));
   if (l <= ls && ls < r) s = M::op(s, t[x].val);
   if (r > ls + 1) s = M::op(s, range_sum(t[x].r, std::max(l, ls + 1) - ls - 1, r - ls - 1));
   return s;
  }
  // 部分木 x の中の位置 [l, r) に作用を当てる。
  template <bool A = AC>
  void range_apply(int x, int l, int r, const E& e) {
   if constexpr (A) {
    if (!x || l >= r) return;
    if (l == 0 && r == t[x].size) return apply_to(x, e), void();
    push(x);
    const int ls = t[t[x].l].size;
    if (l < ls) range_apply(t[x].l, l, std::min(r, ls), e);
    if (l <= ls && ls < r)
     if (x <= n) M::mp(t[x].val, e, 1);
    if (r > ls + 1) range_apply(t[x].r, std::max(l, ls + 1) - ls - 1, r - ls - 1, e);
    pull(x);
   }
  }
  // x を先頭に回した列の根を返す。root は x の木の根 (root_of(x))。
  int reroot(int x, int root) {
   const int k = rank(x);
   int a, b;
   split(root, k, a, b);
   const int m = merge(b, a);
   t[m].par = 0;
   return m;
  }
  // v から p へ向かう節と p から v へ向かう節。辺は (u,v) の向きで覚えてある。
  std::pair<int, int> edge_nodes(int v, int p) const {
   const Edge* found = es.find(key(v, p));
   assert(found);
   if (found->u == v) return {found->a, found->b};
   return {found->b, found->a};
  }
  void reset(int x) {
   Node& a = t[x];
   a.l = a.r = a.par = 0, a.size = 1, a.cnt = x <= n ? 1 : 0;
   a.val = M::ti(), a.sum = a.val;
   if constexpr (AC) a.laz_f = false;
  }
  int alloc_edge_node() {
   if (!free_nodes.empty()) {
    const int x = free_nodes.back();
    free_nodes.pop_back();
    reset(x);
    return x;
   }
   const int x = int(t.size());
   t.emplace_back();
   reset(x);
   return x;
  }
  void free_edge_node(int x) {
   reset(x);
   free_nodes.push_back(x);
  }
  // seq の [l, r) から、区間の中央を根にした釣り合った木を組む。根を返す。
  int build_balanced(const std::vector<int>& seq, int l, int r) {
   if (l >= r) return 0;
   const int m = (l + r) >> 1, x = seq[m];
   const int L = build_balanced(seq, l, m), R = build_balanced(seq, m + 1, r);
   t[x].l = L, t[x].r = R;
   if (L) t[L].par = x;
   if (R) t[R].par = x;
   t[x].par = 0;
   pull(x);
   return x;
  }

 public:
  explicit EulerTourTree(int n_): t(std::size_t(n_ + 1)), n(n_) {
   t.reserve(std::size_t(3 * n_ + 8));  // 森なら辺の節は 2(n-1) 個まで
   for (int v = 1; v <= n; ++v) reset(v);
  }
  // 辺の列 (森でもよい) と頂点の値から O(N) で組む。頂点の節は 1..N、辺の節は N+1 から。
  EulerTourTree(int n_, const std::vector<std::array<int, 2>>& edges, const std::vector<T>& a) {
   n = n_;
   const int E = int(edges.size()), base = n + 1;
   t.resize(std::size_t(n + 2 * E + 2));
   for (int v = 1; v <= n; ++v) reset(v), t[v].val = a[v - 1], t[v].sum = a[v - 1];
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
      if (const int in = snode.back(); in != -1) reset(in), seq.push_back(in);
      su.pop_back(), sp.pop_back(), se.pop_back(), snode.pop_back();
     } else {
      assert(!seen[v]);  // 入力は森 (閉路も重複辺も無い)
      seen[v] = 1;
      const int node = base + 2 * ee + dd;
      reset(node);
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
   const int x = v + 1;
   push_path(x);
   return t[x].val;
  }
  void set(int v, T x) {
   const int y = v + 1;
   push_path(y);
   t[y].val = x;
   for (int z = y; z; z = t[z].par) pull(z);
  }
  void link(int u, int v) {
   assert(u != v && !edge_exists(u, v));
   const int r0u = root_of(u + 1), r0v = root_of(v + 1);
   assert(r0u != r0v);  // 同じ成分の 2 頂点をつなぐと巡回列が二重になる
   const int a = alloc_edge_node(), b = alloc_edge_node();
   const int ru = reroot(u + 1, r0u), rv = reroot(v + 1, r0v);
   const int x = merge(merge(merge(ru, a), rv), b);
   t[x].par = 0;
   es.insert(key(u, v), Edge{u, v, a, b});
  }
  void cut(int u, int v) {
   const Edge* found = es.find(key(u, v));
   assert(found);
   const Edge e = *found;
   es.erase(key(u, v));
   int a, b;
   if (e.u == u) a = e.a, b = e.b;
   else a = e.b, b = e.a;
   const int root = root_of(a);
   const int ra = rank(a), rb = rank(b);
   int p1, p2, p3, p4, p5, mid;
   if (ra > rb) {  // [.. b .. u の側 .. a ..] で、v の側は p1 と p5 に分かれる
    split(root, rb, p1, p2);           // p1: v の側の前半, p2: b から
    split(p2, 1, mid, p3);             // mid = b
    split(p3, ra - rb - 1, mid, p4);   // mid = u の側 (切れて 1 つの木になる)
    split(p4, 1, mid, p5);             // mid = a
    static_cast<void>(merge(p5, p1));  // v の側を継ぎ直す
   } else {  // [.. a .. v の側 .. b ..] で、u の側は p1 と p5 に分かれる
    split(root, ra, p1, p2);
    split(p2, 1, mid, p3);             // mid = a
    split(p3, rb - ra - 1, mid, p4);   // mid = v の側 (切れて 1 つの木になる)
    split(p4, 1, mid, p5);             // mid = b
    static_cast<void>(merge(p5, p1));  // u の側を継ぎ直す
   }
   free_edge_node(a), free_edge_node(b);
  }
  bool connected(int u, int v) { return root_of(u + 1) == root_of(v + 1); }
  // 頂点 v の木の全体の積。
  T prod_tree(int v) {
   const int x = root_of(v + 1);
   return t[x].sum;
  }
  // 木の全体の頂点の数。
  int tree_size(int v) {
   const int x = root_of(v + 1);
   return t[x].cnt;
  }
  // 木の全体への作用。
  template <bool A = AC>
  void apply_tree(int v, const E& e) {
   static_assert(A, "\"apply_tree\" is not available");
   if constexpr (A) apply_to(root_of(v + 1), e);
  }
  // 辺 (v, p) について、p を親とした v の部分木の積。v と p は辺で結ばれていること。
  T prod_subtree(int v, int p) {
   auto [a, b] = edge_nodes(v, p);  // a: v->p, b: p->v
   const int root = root_of(b);
   const int ra = rank(a), rb = rank(b);
   if (ra > rb) return range_sum(root, rb + 1, ra);
   T s = range_sum(root, rb + 1, t[root].size);
   return M::op(s, range_sum(root, 0, ra));
  }
  // 辺 (v, p) について、p を親とした v の部分木への作用。
  template <bool A = AC>
  void apply_subtree(int v, int p, const E& e) {
   static_assert(A, "\"apply_subtree\" is not available");
   if constexpr (A) {
    auto [a, b] = edge_nodes(v, p);  // a: v->p, b: p->v
    const int root = root_of(b);
    const int ra = rank(a), rb = rank(b);
    if (ra > rb) range_apply(root, rb + 1, ra, e);
    else range_apply(root, rb + 1, t[root].size, e), range_apply(root, 0, ra, e);
   }
  }
};
