#pragma once
// 自己調整 top tree (nandhagk の 337329 の形)。森の辺を切り貼りしながら、部分木への作用と部分木の積、
// 道への作用と道の積を O(log N) で扱う。道の splay (l, r) と rake の splay (h) を 1 つの節で持ち、
// len/sum が道、cap/sub が部分木 (rake の道も含む)。lzy は道へ、hzy は部分木へ。
// expose で x を根のクラスタにし、x->h が x の子の部分木 (rake) になる。
// M は可換モノイド (using T、static T ti()、static T op(T, T)、using commute) で、作用は長さを受ける
// mp(T&, E, int) と cp(E&, E) を持つ。作用の逆元は要らない (構造を変える前に遅延を降ろすため)。
// 記録は algo-notes の notes/heavy-algorithms.md。
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <deque>
#include <utility>
#include <vector>
#include "neo/base/core/algebra/concepts.hpp"

template <class M>
class TopTree {
 using T = typename M::T;
 using E = typename M::E;
 static constexpr bool OK = Monoid<M> && Commutative<M> && ActionLen<M>;
 static_assert(OK, "TopTree requires a commutative monoid with a length action");

 struct Node {
  Node *l = nullptr, *r = nullptr, *h = nullptr, *p = nullptr;
  T val = M::ti(), sum = M::ti(), sub = M::ti();
  E lz{}, hzy{};
  bool lazf = false, hzyf = false;
  std::uint32_t len = 0, cap = 0;
  bool rev = false;
  Node() = default;
  explicit Node(T x): val(x), sum(x), len(1) {}
  bool is_root() const { return !p || (p->l != this && p->r != this); }
 };

 std::deque<Node> buf;
 std::vector<Node *> fre;
 int n;

 static std::uint32_t len_(const Node *x) { return x ? x->len : 0; }
 static std::uint32_t cap_(const Node *x) { return x ? x->cap : 0; }
 static const T &sum_(const Node *x) { return x ? x->sum : ti_; }
 static const T &sub_(const Node *x) { return x ? x->sub : ti_; }
 static inline T ti_ = M::ti();

 static void setl(Node *x, Node *y) {
  if ((x->l = y)) y->p = x;
 }
 static void setr(Node *x, Node *y) {
  if ((x->r = y)) y->p = x;
 }
 static void seth(Node *x, Node *y) {
  if ((x->h = y)) y->p = x;
 }

 template <bool F>
 static void pull(Node *x) {
  if (!x) return;
  if constexpr (!F) {
   x->len = len_(x->l) + 1 + len_(x->r);
   x->sum = M::op(M::op(sum_(x->l), x->val), sum_(x->r));
   x->cap = cap_(x->l) + cap_(x->r) + cap_(x->h);
   x->sub = M::op(M::op(sub_(x->l), sub_(x->r)), sub_(x->h));
  } else {
   x->cap = cap_(x->l) + cap_(x->r) + cap_(x->h) + len_(x->h);
   x->sub = M::op(M::op(sub_(x->l), sub_(x->r)), M::op(sub_(x->h), sum_(x->h)));
  }
 }
 static void apply_at(Node *x, const E &e) {
  if (x->lazf) M::cp(x->lz, e);
  else x->lz = e, x->lazf = true;
  M::mp(x->val, e, 1);
  M::mp(x->sum, e, x->len);
 }
 static void apply_sub(Node *x, const E &e) {
  if (x->hzyf) M::cp(x->hzy, e);
  else x->hzy = e, x->hzyf = true;
  M::mp(x->sub, e, x->cap);
 }
 static void flip(Node *x) {
  std::swap(x->l, x->r);
  x->rev ^= true;
 }
 template <bool F>
 static void push(Node *x) {
  if (!x) return;
  if constexpr (!F) {
   if (x->rev) {
    if (x->l) flip(x->l);
    if (x->r) flip(x->r);
    x->rev = false;
   }
   if (x->lazf) {
    if (x->l) apply_at(x->l, x->lz);
    if (x->r) apply_at(x->r, x->lz);
    x->lazf = false;
   }
   if (x->hzyf) {
    if (x->l) apply_sub(x->l, x->hzy);
    if (x->r) apply_sub(x->r, x->hzy);
    if (x->h) apply_sub(x->h, x->hzy);
    x->hzyf = false;
   }
  } else {
   if (x->hzyf) {
    if (x->l) apply_sub(x->l, x->hzy);
    if (x->r) apply_sub(x->r, x->hzy);
    if (x->h) {
     apply_sub(x->h, x->hzy);
     apply_at(x->h, x->hzy);
    }
    x->hzyf = false;
   }
  }
 }

 template <bool F>
 static void rotate(Node *x) {
  Node *p = x->p, *g = p->p;
  if (g) {
   if (g->h == p) g->h = x;
   else if (g->r == p) g->r = x;
   else g->l = x;
  }
  if (p->r == x) {
   setr(p, x->l);
   setl(x, p);
  } else {
   setl(p, x->r);
   setr(x, p);
  }
  x->p = g;
  pull<F>(p), pull<F>(x);
 }
 template <bool F>
 static void down(Node *x) {
  if (!x->is_root()) down<F>(x->p);
  push<F>(x);
 }
 template <bool F>
 static void splay(Node *x) {
  down<F>(x);
  while (!x->is_root()) {
   Node *p = x->p, *q = p->p;
   if (!p->is_root()) rotate<F>(((p == q->l) ^ (x == p->l)) ? x : p);
   rotate<F>(x);
  }
 }
 template <bool F>
 static void down(Node *x, Node *g) {
  if (!x->is_root() && x->p != g) down<F>(x->p, g);
  push<F>(x);
 }
 template <bool F>
 static void splay(Node *x, Node *g) {
  down<F>(x, g);
  while (!x->is_root() && x->p != g) {
   Node *p = x->p, *q = p->p;
   if (!p->is_root() && q != g) rotate<F>(((p == q->l) ^ (x == p->l)) ? x : p);
   rotate<F>(x);
  }
 }
 static void splay_rightmost(Node *x) {
  Node *g = x->p;
  for (; x->r; x = x->r) push<true>(x);
  splay<true>(x, g);
 }

 Node *create() {
  Node *v;
  if (fre.empty()) {
   buf.emplace_back();
   v = &buf.back();
  } else {
   v = fre.back();
   fre.pop_back();
  }
  v->len = 1;
  return v;
 }
 void del(Node *x) {
  fre.push_back(x);
  *x = Node();
 }

 void splice(Node *x) {
  splay<true>(x), x = x->p, splay<false>(x);
  Node *y = x->h;
  push<true>(y);
  if (x->r) {
   std::swap(x->r->p, y->h->p);
   std::swap(x->r, y->h);
  } else {
   setr(x, y->h);
   if (y->l) {
    splay_rightmost(y->l);
    setr(y->l, y->r);
    x->h = y->l;
   } else {
    push<true>(y->r);
    x->h = y->r;
   }
   del(y);
   if ((y = x->h)) y->p = x;
  }
  pull<true>(y), pull<false>(x);
  push<false>(x->r);
  rotate<false>(x->r);
 }

 void expose_(Node *x) {
  splay<false>(x);
  if (x->r) {
   Node *y = create();
   setl(y, x->h);
   seth(y, x->r);
   pull<true>(y);
   x->r = nullptr;
   seth(x, y);
   pull<false>(x);
  }
  while (x->p) splice(x->p);
 }
 void evert_(Node *x) {
  expose_(x);
  flip(x);
 }
 void link_(Node *x, Node *y) {
  evert_(y), expose_(x);
  setr(x, y);
  pull<false>(x);
 }
 void cut_(Node *x) {
  expose_(x);
  push<false>(x);
  if (x->l) {
   x->l->p = nullptr;
   x->l = nullptr;
   pull<false>(x);
  }
 }
 void apply_subtree_(Node *x, const E &e) {
  expose_(x);
  M::mp(x->val, e, 1);
  pull<false>(x);
  if (x->h) apply_sub(x->h, e);
 }
 T fold_subtree_(Node *x) {
  expose_(x);
  return M::op(x->val, sub_(x->h));
 }
 void apply_path_(Node *x, const E &e) {
  expose_(x);
  apply_at(x, e);
 }
 T fold_path_(Node *x) {
  expose_(x);
  return x->sum;
 }

 public:
 explicit TopTree(const std::vector<T> &a): n((int)a.size()) {
  buf.resize(n);
  for (int i = 0; i < n; ++i) buf[i] = Node(a[i]);
 }
 // 辺 (u, v) を張る。u と v は別の木にあること。
 void link(int u, int v) { link_(&buf[u], &buf[v]); }
 // 辺 (u, v) を切る。
 void cut(int u, int v) {
  evert_(&buf[v]);
  cut_(&buf[u]);
 }
 void evert(int v) { evert_(&buf[v]); }
 // 辺 (v, p) について、p を親とした v の部分木に e を当てる。
 void apply_subtree(int v, int p, const E &e) {
  evert_(&buf[p]);
  apply_subtree_(&buf[v], e);
 }
 // 辺 (v, p) について、p を親とした v の部分木の積。
 T prod_subtree(int v, int p) {
  evert_(&buf[p]);
  return fold_subtree_(&buf[v]);
 }
 // u から v への道に e を当てる。
 void apply_path(int u, int v, const E &e) {
  evert_(&buf[u]);
  apply_path_(&buf[v], e);
 }
 // u から v への道の積。
 T prod_path(int u, int v) {
  evert_(&buf[u]);
  return fold_path_(&buf[v]);
 }
 void set(int v, T x) {
  Node *t = &buf[v];
  evert_(t), expose_(t);
  t->val = x;
  pull<false>(t);
 }
 T get(int v) {
  Node *t = &buf[v];
  evert_(t), expose_(t);
  return t->val;
 }
};
