#pragma once
#include <algorithm>
#include <bit>
#include <cassert>
#include <cstddef>
#include <type_traits>
#include <vector>
#include "neo/base/core/algebra/concepts.hpp"
#include "neo/base/ds/range/SegmentTree.hpp"
// 永続遅延セグメント木。構造体 M でモノイドと作用を渡す (SegmentTree と同じ決まり)。
//   using T、static T ti()、static T op(T a, T b)
//   using E、static void mp(T& v, E x) か mp(T& v, E x, int len)、static void cp(E& pre, E suf)
// apply(l, r, x) と copy_range(l, r, src) が新しい版を返す。版の写しは根 1 個 (O(1))。prod / get はその版の値を読む。
// copy_range は、src の [l, r) をこの版の [l, r) に写す。Library Checker の persistent_range_affine_range_sum
// (範囲コピーつき) を解くために要る。src の道の節は写してから遅延を配るので、src は変えない。
// 節は型ごとに 1 本の置き場に取り、0 番を空の節 (子が自分、値は ti()、遅延なし) にする。apply / copy は根から
// 葉への道の節を写すので、1 回 O(log n) 個の節が増える。pool_empty() が真を返したら、生きている版を渡して gc を呼ぶ。
namespace persistent_lazy_segment_tree_internal {
template <class M> struct Node {
 int l, r;
 typename M::T v;
 typename M::E lz;
 bool flag;
};
}
template <class M> class PersistentLazySegmentTree {
 static_assert(Monoid<M>, "M にはモノイド (T, ti, op) が要る");
 static_assert(Action<M>, "M には作用 (T, E, mp, cp) が要る");
 static_assert(!std::is_invocable_r_v<bool, decltype(M::mp), typename M::T&, const typename M::E&, int> && !std::is_invocable_r_v<bool, decltype(M::mp), typename M::T&, const typename M::E&>, "永続 lazy は Beats の作用 (bool を返す mp) に未対応");
 using T= typename M::T;
 using E= typename M::E;
 using Node= persistent_lazy_segment_tree_internal::Node<M>;
 using Pool= std::vector<Node, segment_tree_internal::Alloc<Node>>;
 static inline Pool pool_= [] {
  Pool p;
  p.push_back(Node{0, 0, M::ti(), E{}, false});
  return p;
 }();
 static inline std::size_t gc_at_= std::size_t(1) << 22;
 int n_= 0, sz_= 1, root_= 0;
 static int new_node(const Node& x) {
  const int i= pool_.size();
  pool_.push_back(x);
  return i;
 }
 // mp の 2 つの形をそろえる。Beats の bool を返す形は受けない (static_assert)。
 static void map(T& v, const E& x, int len) {
  if constexpr(std::is_invocable_v<decltype(M::mp), T&, const E&, int>) M::mp(v, x, len);
  else M::mp(v, x);
 }
 static void all_apply(int q, const E& x, int len) {
  if(pool_[q].flag) M::cp(pool_[q].lz, x);
  else pool_[q].lz= x, pool_[q].flag= true;
  map(pool_[q].v, x, len);
 }
 static void push(int q, int len) {
  if(!pool_[q].flag) return;
  const E x= pool_[q].lz;
  const int h= len >> 1;
  const int lc= new_node(pool_[pool_[q].l]), rc= new_node(pool_[pool_[q].r]);
  all_apply(lc, x, h), all_apply(rc, x, h);
  pool_[q].l= lc, pool_[q].r= rc, pool_[q].flag= false;
 }
 static int build_rec(int l, int r, const std::vector<T>& a, int n) {
  int p= new_node(pool_[0]);
  if(r - l == 1) {
   pool_[p].v= l < n ? a[l] : M::ti();
   return p;
  }
  const int m= (l + r) >> 1;
  const int lc= build_rec(l, m, a, n), rc= build_rec(m, r, a, n);
  pool_[p].l= lc, pool_[p].r= rc;
  pool_[p].v= M::op(pool_[lc].v, pool_[rc].v);
  return p;
 }
 int apply_rec(int p, int nl, int nr, int ql, int qr, const E& x) const {
  if(qr <= nl || nr <= ql) return p;
  int q= new_node(pool_[p]);
  if(ql <= nl && nr <= qr) {
   all_apply(q, x, nr - nl);
   return q;
  }
  push(q, nr - nl);
  const int m= (nl + nr) >> 1;
  pool_[q].l= apply_rec(pool_[q].l, nl, m, ql, qr, x);
  pool_[q].r= apply_rec(pool_[q].r, m, nr, ql, qr, x);
  pool_[q].v= M::op(pool_[pool_[q].l].v, pool_[pool_[q].r].v);
  return q;
 }
 int copy_rec(int pt, int ps, int nl, int nr, int ql, int qr) const {
  if(qr <= nl || nr <= ql) return pt;
  if(ql <= nl && nr <= qr) return ps;
  int q= new_node(pool_[pt]);
  push(q, nr - nl);
  int s= new_node(pool_[ps]);
  push(s, nr - nl);
  const int m= (nl + nr) >> 1;
  pool_[q].l= copy_rec(pool_[q].l, pool_[s].l, nl, m, ql, qr);
  pool_[q].r= copy_rec(pool_[q].r, pool_[s].r, m, nr, ql, qr);
  pool_[q].v= M::op(pool_[pool_[q].l].v, pool_[pool_[q].r].v);
  return q;
 }
 // acc はここより上 (根に近い側) の遅延を、下 (葉に近い側) から上へ当てる順に合成したもの。has が偽なら空。
 T prod_rec(int p, int nl, int nr, int ql, int qr, const E& acc, bool has) const {
  if(qr <= nl || nr <= ql) return M::ti();
  const int len= nr - nl;
  if(ql <= nl && nr <= qr) {
   T v= pool_[p].v;
   if(has) map(v, acc, len);
   return v;
  }
  E nacc= acc;
  bool nh= has;
  if(pool_[p].flag) {
   if(nh) {
    E t= pool_[p].lz;
    M::cp(t, nacc), nacc= t;
   } else nacc= pool_[p].lz, nh= true;
  }
  const int m= (nl + nr) >> 1;
  return M::op(prod_rec(pool_[p].l, nl, m, ql, qr, nacc, nh), prod_rec(pool_[p].r, m, nr, ql, qr, nacc, nh));
 }
 static void gc_roots(const std::vector<int*>& roots) {
  const std::size_t before= pool_.size();
  std::vector<char> alive(before, 0);
  std::vector<int> st;
  for(int* r: roots) st.push_back(*r);
  while(!st.empty()) {
   const int c= st.back();
   st.pop_back();
   if(c == 0 || alive[c]) continue;
   alive[c]= 1;
   st.push_back(pool_[c].l);
   st.push_back(pool_[c].r);
  }
  std::vector<int> map(before, 0);
  int k= 1;
  for(std::size_t i= 1; i < before; ++i)
   if(alive[i]) map[i]= k++;
  for(std::size_t i= 1; i < before; ++i)
   if(alive[i]) {
    Node x= pool_[i];
    x.l= map[x.l], x.r= map[x.r];
    pool_[map[i]]= x;
   }
  pool_.resize(k);
  for(int* r: roots) *r= map[*r];
  const std::size_t live= k;
  gc_at_= std::max<std::size_t>(std::size_t(1) << 22, (2 * live > before ? 4 : 2) * live);
 }
public:
 PersistentLazySegmentTree()= default;
 explicit PersistentLazySegmentTree(int n): n_(n), sz_(std::bit_ceil(unsigned(std::max(n, 1)))) {}
 explicit PersistentLazySegmentTree(const std::vector<T>& a): n_(a.size()), sz_(std::bit_ceil(unsigned(std::max<std::size_t>(a.size(), 1)))) {
  if(n_) root_= build_rec(0, sz_, a, n_);
 }
 int size() const { return n_; }
 T get(int i) const {
  int p= root_, l= 0, r= sz_;
  E acc;
  bool has= false;
  while(r - l > 1) {
   if(pool_[p].flag) {
    if(has) {
     E t= pool_[p].lz;
     M::cp(t, acc), acc= t;
    } else acc= pool_[p].lz, has= true;
   }
   const int m= (l + r) >> 1;
   if(i < m) p= pool_[p].l, r= m;
   else p= pool_[p].r, l= m;
  }
  T v= pool_[p].v;
  if(has) map(v, acc, 1);
  return v;
 }
 T prod(int l, int r) const { return prod_rec(root_, 0, sz_, l, r, E{}, false); }
 T all_prod() const { return pool_[root_].v; }
 PersistentLazySegmentTree apply(int l, int r, E x) const {
  PersistentLazySegmentTree res= *this;
  res.root_= apply_rec(root_, 0, sz_, l, r, x);
  return res;
 }
 PersistentLazySegmentTree copy_range(int l, int r, const PersistentLazySegmentTree& src) const {
  assert(sz_ == src.sz_);
  PersistentLazySegmentTree res= *this;
  res.root_= copy_rec(root_, src.root_, 0, sz_, l, r);
  return res;
 }
 static void reserve(std::size_t n) {
  if(pool_.capacity() < n + 1) pool_.reserve(std::max(n + 1, 2 * pool_.capacity()));
 }
 static bool pool_empty() { return pool_.size() >= gc_at_; }
 static void gc(std::vector<PersistentLazySegmentTree>& vs) {
  std::vector<int*> roots;
  roots.reserve(vs.size());
  for(auto& v: vs) roots.push_back(&v.root_);
  gc_roots(roots);
 }
 template <class... Ps> static void gc(PersistentLazySegmentTree& v, Ps&... vs) {
  std::vector<int*> roots{&v.root_, &vs.root_...};
  gc_roots(roots);
 }
};
