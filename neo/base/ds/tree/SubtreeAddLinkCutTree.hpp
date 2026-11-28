#pragma once
// 部分木への作用と部分木の積を O(log N) で扱う Link-Cut Tree (416793 の形)。
// 森の辺を切り貼りしながら、辺 (v, p) について p を親とした v の部分木の積を答え、部分木に作用を当てる。
// 節は道の和 sum と頂点数 sz に加えて、仮想の子の和 vir と頂点数 vir_sz、道の遅延 lz、仮想の子への遅延 vir_lz を持つ。
// access で右の子を仮想へ移すときは vir に足して -vir_lz を当て、仮想から道へ取るときは +vir_lz を当ててから
// vir から引く (仮想の子の値は、自分に掛かった vir_lz を含まない形で持つ)。
// M は可換群 (using commute、static T inv(T)) で、作用は長さを受ける mp(T&, E, int) と cp(E&, E)、
// それに作用の逆元 inve(E) を持つ。記録は algo-notes の notes/heavy-algorithms.md。
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

namespace subtree_add_lct_internal {
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
}  // namespace subtree_add_lct_internal

template <class M>
class SubtreeAddLinkCutTree {
 using T = typename M::T;
 using E = typename M::E;
 static constexpr bool CG = Commutative<M> && ActionLen<M> && requires(T a, E e) { M::inv(a); M::inve(e); };
 static_assert(CG, "SubtreeAddLinkCutTree requires a commutative group, a length action and its inverse");

 struct Node {
  int ch[2] = {0, 0};
  int par = 0;
  bool rev = false;
  T val = M::ti(), sum = M::ti(), vir = M::ti();
  E lz{}, vir_lz{};
  bool lazf = false, virlzf = false;
  int sz = 1, vir_sz = 0;
 };
 subtree_add_lct_internal::Vec<Node> t;
 int n;

 // 節 x の補助木の全体 (道の頂点と、その仮想の子の部分木) に e を当てる。
 [[gnu::always_inline]] void apply_all(int x, const E& e) {
  if (!x) return;
  Node& a = t[x];
  M::mp(a.val, e, 1);
  M::mp(a.sum, e, a.sz);
  M::mp(a.vir, e, a.vir_sz);
  if (a.lazf) M::cp(a.lz, e);
  else a.lz = e, a.lazf = true;
  if (a.virlzf) M::cp(a.vir_lz, e);
  else a.vir_lz = e, a.virlzf = true;
 }
 [[gnu::always_inline]] void pull(int x) {
  Node& a = t[x];
  a.sum = M::op(M::op(t[a.ch[0]].sum, a.val), M::op(a.vir, t[a.ch[1]].sum));
  a.sz = t[a.ch[0]].sz + 1 + a.vir_sz + t[a.ch[1]].sz;
 }
 [[gnu::always_inline]] void reverse(int x) {
  if (!x) return;
  std::swap(t[x].ch[0], t[x].ch[1]), t[x].rev ^= 1;
 }
 [[gnu::always_inline]] void push(int x) {
  Node& a = t[x];
  if (a.rev) reverse(a.ch[0]), reverse(a.ch[1]), a.rev = false;
  if (a.lazf) apply_all(a.ch[0], a.lz), apply_all(a.ch[1], a.lz), a.lazf = false;
 }
 [[gnu::always_inline]] void rotate(int x) {
  Node& a = t[x];
  int p = a.par, g = t[p].par, d = t[p].ch[1] == x, b = a.ch[d ^ 1];
  t[p].ch[d] = b;
  if (b) t[b].par = p;
  a.ch[d ^ 1] = p, t[p].par = x, a.par = g;
  if (g > 0) t[g].ch[t[g].ch[1] == p] = x;
  a.sum = t[p].sum, a.sz = t[p].sz;
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
 // x の右の子を仮想へ移し (c が新しい右の子)、仮想の遅延を釣り合わせる。
 [[gnu::always_inline]] void set_r(int x, int c) {
  Node& a = t[x];
  if (int r = a.ch[1]; r) {
   t[r].par = -x;
   a.vir = M::op(a.vir, t[r].sum);
   a.vir_sz += t[r].sz;
   if (a.virlzf) apply_all(r, M::inve(a.vir_lz));
  }
  if (c) {
   if (a.virlzf) apply_all(c, a.vir_lz);
   a.vir = M::op(a.vir, M::inv(t[c].sum));
   a.vir_sz -= t[c].sz;
  }
  a.ch[1] = c;
 }
 [[gnu::always_inline]] int access(int x) {
  splay(x);
  set_r(x, 0);
  pull(x);
  int last = x;
  while (int p = -t[x].par) {
   splay(p);
   set_r(p, x);
   t[x].par = p;
   pull(p);
   splay(x);
   last = p;
  }
  pull(x);
  return last;
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
 explicit SubtreeAddLinkCutTree(int n_): t(std::size_t(n_) + 2), n(n_) {
  // 0 番は番兵で、頂点を 1 つも持たない。
  t[0].sz = 0, t[0].vir_sz = 0;
 }
 int size() const { return n; }
 T get(int v) {
  int x = v + 1;
  splay(x);
  return t[x].val;
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
  // u は v の仮想の子になる。v の仮想への遅延は u には掛からないので、u の値に
  // -vir_lz を当てて (遅延の基準をそろえて) から、真の和を vir に足す。
  T s = t[u].sum;
  if (t[v].virlzf) apply_all(u, M::inve(t[v].vir_lz));
  t[v].vir = M::op(t[v].vir, s);
  t[v].vir_sz += t[u].sz;
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
 // 辺 (v, p) について、p を親とした v の部分木の積。
 T prod_subtree(int v, int p) {
  int x = v + 1, y = p + 1;
  evert_(y), access(x);
  return M::op(t[x].val, t[x].vir);
 }
 // 辺 (v, p) について、p を親とした v の部分木への作用。
 void apply_subtree(int v, int p, const E& e) {
  int x = v + 1, y = p + 1;
  evert_(y), access(x);
  M::mp(t[x].val, e, 1);
  M::mp(t[x].vir, e, t[x].vir_sz);
  if (t[x].virlzf) M::cp(t[x].vir_lz, e);
  else t[x].vir_lz = e, t[x].virlzf = true;
  pull(x);
 }
};
