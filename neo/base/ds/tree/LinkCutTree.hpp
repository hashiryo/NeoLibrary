#pragma once
// 動的木 (Link-Cut Tree)。森の辺を切り貼りしながら、道の積と道への作用を O(log N) で答える。
//   LinkCutTree<M>: 道の積 (u から v への順、可換でなくてよい)。M は SegmentTree と同じモノイド
//     (using T、static T ti()、static T op(T, T))。作用 (using E、mp(T&, E)、cp(E&, E)) があれば道への作用も使える。
//   SubtreeLinkCutTree<M>: 可換群 (using commute、static T inv(T)) を受け、辺 (v, p) について p を親とした
//     v の部分木の積を答える。
// 節は配列に置き、親は符号つきの int で持つ (正: splay の親、負: 道の親、0: なし)。頂点の番号は 1 始まりで、
// 0 番を番兵にする。回転では x の部分木が p の前の部分木と同じ区間になるので、x の集約を p から写してから
// p だけ求め直す。Library Checker の dynamic_tree_vertex_add_path_sum の首位 (chaihf の 270626) と
// dynamic_tree_vertex_set_path_composite の首位 (near1111 の 419210) の書き方を、一般のモノイドにした。
// 記録は algo-notes の notes/dynamic-tree.md。
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <new>
#include <type_traits>
#include <utility>
#include <vector>
#include "neo/base/core/algebra/concepts.hpp"
#ifdef __linux__
#include <sys/mman.h>
#endif

namespace link_cut_tree_internal {
struct Empty {};
template <class M, bool = Action<M>>
struct ActionOf {
 struct type {};
};
template <class M>
struct ActionOf<M, true> {
 using type = typename M::E;
};
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
 static constexpr bool NC = !Commutative<M>;
 static constexpr bool AC = Action<M>;
 int ch[2] = {0, 0};
 int par = 0;
 bool rev = false;
 T val = M::ti(), sum = M::ti();
 [[no_unique_address]] std::conditional_t<NC, T, Empty> rsum;
 [[no_unique_address]] std::conditional_t<AC, E, Empty> laz;
 [[no_unique_address]] std::conditional_t<AC, bool, Empty> lazf;
 Node() {
  if constexpr (NC) rsum = M::ti();
  if constexpr (AC) lazf = false;
 }
};
}  // namespace link_cut_tree_internal

template <class M = void>
class LinkCutTree {
 using T = typename M::T;
 using E = typename link_cut_tree_internal::ActionOf<M>::type;
 using Node = link_cut_tree_internal::Node<M>;
 static constexpr bool NC = !Commutative<M>;
 static constexpr bool AC = Action<M>;
 link_cut_tree_internal::Vec<Node> t;
 int n;
 [[gnu::always_inline]] void pull(int x) {
  Node& a = t[x];
  T s = a.val;
  if (int l = a.ch[0]; l) s = M::op(t[l].sum, s);
  if (int r = a.ch[1]; r) s = M::op(s, t[r].sum);
  a.sum = s;
  if constexpr (NC) {
   T rs = a.val;
   if (int r = a.ch[1]; r) rs = M::op(t[r].rsum, rs);
   if (int l = a.ch[0]; l) rs = M::op(rs, t[l].rsum);
   a.rsum = rs;
  }
 }
 [[gnu::always_inline]] void reverse(int x) {
  if (!x) return;
  std::swap(t[x].ch[0], t[x].ch[1]);
  if constexpr (NC) std::swap(t[x].sum, t[x].rsum);
  t[x].rev ^= 1;
 }
 template <bool A = AC>
 [[gnu::always_inline]] void apply_to(int x, const E& e) {
  if constexpr (A) {
   if (!x) return;
   Node& a = t[x];
   M::mp(a.val, e), M::mp(a.sum, e);
   if constexpr (NC) M::mp(a.rsum, e);
   if (a.lazf) M::cp(a.laz, e);
   else a.laz = e, a.lazf = true;
  }
 }
 [[gnu::always_inline]] void push(int x) {
  Node& a = t[x];
  if (a.rev) reverse(a.ch[0]), reverse(a.ch[1]), a.rev = false;
  if constexpr (AC)
   if (a.lazf) apply_to(a.ch[0], a.laz), apply_to(a.ch[1], a.laz), a.lazf = false;
 }
 // x は splay の親を持つ。回転しても x の部分木の区間は p の前の部分木と同じなので、集約は p から写す。
 [[gnu::always_inline]] void rotate(int x) {
  Node& a = t[x];
  int p = a.par, g = t[p].par, d = t[p].ch[1] == x, b = a.ch[d ^ 1];
  t[p].ch[d] = b;
  if (b) t[b].par = p;
  a.ch[d ^ 1] = p, t[p].par = x, a.par = g;
  if (g > 0) t[g].ch[t[g].ch[1] == p] = x;
  a.sum = t[p].sum;
  if constexpr (NC) a.rsum = t[p].rsum;
  pull(p);
 }
 [[gnu::always_inline]] void splay(int x) {
  for (int p; (p = t[x].par) > 0;) {
   int g = t[p].par;
   if (g <= 0) {
    push(p), push(x), rotate(x);
    break;
   }
   push(g), push(p), push(x);
   if ((t[g].ch[1] == p) == (t[p].ch[1] == x)) rotate(p), rotate(x);
   else rotate(x), rotate(x);
  }
  push(x);
 }
 [[gnu::always_inline]] int access(int x) {
  int last = 0;
  for (int y = x; y; last = y, y = -t[y].par) {
   splay(y);
   int r = t[y].ch[1];
   if (r) t[r].par = -y;
   t[y].ch[1] = last;
   if (last) t[last].par = y;
   pull(y);
  }
  return splay(x), last;
 }
 [[gnu::always_inline]] void evert_(int x) { access(x), reverse(x); }
 int root_(int x) {
  access(x);
  for (int y = x;;) {
   push(y);
   if (!t[y].ch[0]) return splay(y), y;
   y = t[y].ch[0];
  }
 }
 int parent_(int v) {
  access(v);
  if (!t[v].ch[0]) return -1;
  int y = t[v].ch[0];
  push(y);
  while (t[y].ch[1]) y = t[y].ch[1], push(y);
  splay(y);
  return y;
 }
 int lca_(int u, int v) {
  if (u == v) return u;
  access(u);
  return access(v);
 }

public:
 explicit LinkCutTree(int n_): t(std::size_t(n_) + 2), n(n_) {}
 int size() const { return n; }
 // 1 点の値。前の作用は splay で配ってから読む。
 T get(int v) {
  int x = v + 1;
  splay(x);
  return t[x].val;
 }
 // 1 点の代入。
 void set(int v, T x) {
  int y = v + 1;
  splay(y);
  t[y].val = x;
  pull(y);
 }
 // 1 点の値を op で動かす (可換でなくてよい)。
 void mul(int v, const T& d) {
  int y = v + 1;
  splay(y);
  t[y].val = M::op(t[y].val, d);
  pull(y);
 }
 // 初期の値を置く (splay しない)。
 void init_value(int v, T x) {
  int y = v + 1;
  t[y].val = x, t[y].sum = x;
  if constexpr (NC) t[y].rsum = x;
 }
 // 初期の木を親の向きで張る。すべて仮想辺なので splay しない (v の親 p)。
 void init_parent(int v, int p) { t[v + 1].par = -(p + 1); }
 void evert(int v) { evert_(v + 1); }
 // u と v を結ぶ (u と v は別の木にあること)。
 void link(int u, int v) {
  u += 1, v += 1;
  evert_(u), t[u].par = -v;
 }
 // 辺 (u, v) を切る。
 void cut(int u, int v) {
  u += 1, v += 1;
  evert_(u), access(v);
  assert(t[v].ch[0] == u && t[u].ch[1] == 0);
  t[v].ch[0] = 0, t[u].par = 0;
  pull(v);
 }
 bool connected(int u, int v) { return root_(u + 1) == root_(v + 1); }
 int root(int v) { return root_(v + 1) - 1; }
 // v の親 (根なら -1)。
 int parent(int v) {
  int r = parent_(v + 1);
  return r < 0 ? -1 : r - 1;
 }
 // 同じ木の 2 頂点の lca。
 int lca(int u, int v) { return lca_(u + 1, v + 1) - 1; }
 // u から v への道の積 (道の順)。u と v は同じ木にあること。
 T prod(int u, int v) {
  u += 1, v += 1;
  evert_(u), access(v);
  return t[v].sum;  // v の右の子は空で、和は u から v への道の順
 }
 // 道 u-v への作用。
 template <bool A = AC>
 void apply(int u, int v, const E& e) {
  static_assert(A, "\"apply\" is not available");
  if constexpr (A) u += 1, v += 1, evert_(u), access(v), apply_to(v, e);
 }
};

// 可換な群 (inv を持つ) のモノイドを受ける。道の積に加えて、辺 (v, p) で p を親とした v の部分木の積を答える。
// 節は splay の部分木の和 (道の上の頂点と、仮想の子の部分木) と、仮想の子の部分木の和 vir を持つ。access で
// 右の子を仮想へ移すとき vir に足し、仮の道の子を仮想から外すとき vir から引く。
template <class M>
class SubtreeLinkCutTree {
 using T = typename M::T;
 static constexpr bool CG = Commutative<M> && requires(T a) { M::inv(a); };
 static_assert(CG, "SubtreeLinkCutTree requires a commutative group (commute and inv)");
 struct Node {
  int ch[2] = {0, 0};
  int par = 0;
  bool rev = false;
  T val = M::ti(), sum = M::ti(), vir = M::ti();
 };
 link_cut_tree_internal::Vec<Node> t;
 int n;
 [[gnu::always_inline]] void pull(int x) {
  Node& a = t[x];
  T s = a.val;
  if (int l = a.ch[0]; l) s = M::op(t[l].sum, s);
  s = M::op(s, a.vir);
  if (int r = a.ch[1]; r) s = M::op(s, t[r].sum);
  a.sum = s;
 }
 [[gnu::always_inline]] void reverse(int x) {
  if (!x) return;
  std::swap(t[x].ch[0], t[x].ch[1]), t[x].rev ^= 1;
 }
 [[gnu::always_inline]] void push(int x) {
  if (t[x].rev) reverse(t[x].ch[0]), reverse(t[x].ch[1]), t[x].rev = false;
 }
 [[gnu::always_inline]] void rotate(int x) {
  Node& a = t[x];
  int p = a.par, g = t[p].par, d = t[p].ch[1] == x, b = a.ch[d ^ 1];
  t[p].ch[d] = b;
  if (b) t[b].par = p;
  a.ch[d ^ 1] = p, t[p].par = x, a.par = g;
  if (g > 0) t[g].ch[t[g].ch[1] == p] = x;
  a.sum = t[p].sum;
  pull(p);
 }
 [[gnu::always_inline]] void splay(int x) {
  for (int p; (p = t[x].par) > 0;) {
   int g = t[p].par;
   if (g <= 0) {
    push(p), push(x), rotate(x);
    break;
   }
   push(g), push(p), push(x);
   if ((t[g].ch[1] == p) == (t[p].ch[1] == x)) rotate(p), rotate(x);
   else rotate(x), rotate(x);
  }
  push(x);
 }
 [[gnu::always_inline]] int access(int x) {
  int last = 0;
  for (int y = x; y; last = y, y = -t[y].par) {
   splay(y);
   int r = t[y].ch[1];
   if (r != last) {
    if (r) t[r].par = -y, t[y].vir = M::op(t[y].vir, t[r].sum);
    if (last) t[y].vir = M::op(t[y].vir, M::inv(t[last].sum));
   }
   t[y].ch[1] = last;
   if (last) t[last].par = y;
   pull(y);
  }
  return splay(x), last;
 }
 [[gnu::always_inline]] void evert_(int x) { access(x), reverse(x); }
 int root_(int x) {
  access(x);
  for (int y = x;;) {
   push(y);
   if (!t[y].ch[0]) return splay(y), y;
   y = t[y].ch[0];
  }
 }
 int parent_(int v) {
  access(v);
  if (!t[v].ch[0]) return -1;
  int y = t[v].ch[0];
  push(y);
  while (t[y].ch[1]) y = t[y].ch[1], push(y);
  splay(y);
  return y;
 }
 int lca_(int u, int v) {
  if (u == v) return u;
  access(u);
  return access(v);
 }

public:
 explicit SubtreeLinkCutTree(int n_): t(std::size_t(n_) + 2), n(n_) {}
 int size() const { return n; }
 T get(int v) {
  splay(v + 1);
  return t[v + 1].val;
 }
 // 1 点の値を d だけ動かす。
 void add(int v, T d) {
  int x = v + 1;
  access(x);
  t[x].val = M::op(t[x].val, d);
  pull(x);
 }
 void set(int v, T x) {
  int y = v + 1;
  access(y);
  t[y].val = x;
  pull(y);
 }
 void evert(int v) { evert_(v + 1); }
 void link(int u, int v) {
  u += 1, v += 1;
  evert_(u), access(v);
  t[v].vir = M::op(t[v].vir, t[u].sum);
  t[u].par = -v;
  pull(v);
 }
 void cut(int u, int v) {
  u += 1, v += 1;
  evert_(u), access(v);
  assert(t[v].ch[0] == u && t[u].ch[1] == 0);
  t[v].ch[0] = 0, t[u].par = 0;
  pull(v);
 }
 bool connected(int u, int v) { return root_(u + 1) == root_(v + 1); }
 int root(int v) { return root_(v + 1) - 1; }
 int parent(int v) {
  int r = parent_(v + 1);
  return r < 0 ? -1 : r - 1;
 }
 int lca(int u, int v) { return lca_(u + 1, v + 1) - 1; }
 // v の木の全体の積。
 T prod_tree(int v) {
  int x = v + 1;
  access(x);
  return t[x].sum;
 }
 // 辺 (v, p) について、p を親とした v の部分木の積。v と p は辺で結ばれていること。
 T prod_subtree(int v, int p) {
  int x = v + 1, y = p + 1;
  evert_(x), access(y);
  return M::op(t[x].val, t[x].vir);
 }
};

// 値を持たない森 (M = void)。連結性、根、親、lca と辺の切り貼りだけを使う。
template <>
class LinkCutTree<void> {
 struct Node {
  int ch[2] = {0, 0};
  int par = 0;
  bool rev = false;
 };
 link_cut_tree_internal::Vec<Node> t;
 int n;
 [[gnu::always_inline]] void reverse(int x) {
  if (!x) return;
  std::swap(t[x].ch[0], t[x].ch[1]), t[x].rev ^= 1;
 }
 [[gnu::always_inline]] void push(int x) {
  if (t[x].rev) reverse(t[x].ch[0]), reverse(t[x].ch[1]), t[x].rev = false;
 }
 [[gnu::always_inline]] void rotate(int x) {
  Node& a = t[x];
  int p = a.par, g = t[p].par, d = t[p].ch[1] == x, b = a.ch[d ^ 1];
  t[p].ch[d] = b;
  if (b) t[b].par = p;
  a.ch[d ^ 1] = p, t[p].par = x, a.par = g;
  if (g > 0) t[g].ch[t[g].ch[1] == p] = x;
 }
 [[gnu::always_inline]] void splay(int x) {
  for (int p; (p = t[x].par) > 0;) {
   int g = t[p].par;
   if (g <= 0) {
    push(p), push(x), rotate(x);
    break;
   }
   push(g), push(p), push(x);
   if ((t[g].ch[1] == p) == (t[p].ch[1] == x)) rotate(p), rotate(x);
   else rotate(x), rotate(x);
  }
  push(x);
 }
 [[gnu::always_inline]] int access(int x) {
  int last = 0;
  for (int y = x; y; last = y, y = -t[y].par) {
   splay(y);
   int r = t[y].ch[1];
   if (r) t[r].par = -y;
   t[y].ch[1] = last;
   if (last) t[last].par = y;
  }
  return splay(x), last;
 }
 [[gnu::always_inline]] void evert_(int x) { access(x), reverse(x); }
 int root_(int x) {
  access(x);
  for (int y = x;;) {
   push(y);
   if (!t[y].ch[0]) return splay(y), y;
   y = t[y].ch[0];
  }
 }
 int parent_(int v) {
  access(v);
  if (!t[v].ch[0]) return -1;
  int y = t[v].ch[0];
  push(y);
  while (t[y].ch[1]) y = t[y].ch[1], push(y);
  splay(y);
  return y;
 }
 int lca_(int u, int v) {
  if (u == v) return u;
  access(u);
  return access(v);
 }

public:
 explicit LinkCutTree(int n_): t(std::size_t(n_) + 2), n(n_) {}
 int size() const { return n; }
 // 初期の木を親の向きで張る (v の親 p)。すべて仮想辺なので splay しない。
 void init_parent(int v, int p) { t[v + 1].par = -(p + 1); }
 void evert(int v) { evert_(v + 1); }
 void link(int u, int v) {
  u += 1, v += 1;
  evert_(u), t[u].par = -v;
 }
 void cut(int u, int v) {
  u += 1, v += 1;
  evert_(u), access(v);
  assert(t[v].ch[0] == u && t[u].ch[1] == 0);
  t[v].ch[0] = 0, t[u].par = 0;
 }
 bool connected(int u, int v) { return root_(u + 1) == root_(v + 1); }
 int root(int v) { return root_(v + 1) - 1; }
 int parent(int v) {
  int r = parent_(v + 1);
  return r < 0 ? -1 : r - 1;
 }
 int lca(int u, int v) { return lca_(u + 1, v + 1) - 1; }
};
