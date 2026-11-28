#pragma once
#include <algorithm>
#include <bit>
#include <cassert>
#include <cstddef>
#include <type_traits>
#include <vector>
#include "neo/base/core/algebra/concepts.hpp"
#include "neo/base/ds/range/SegmentTree.hpp"
// 永続セグメント木。構造体 M でモノイドを渡す (SegmentTree と同じ決まり)。
//   using T、static T ti()、static T op(T a, T b)
// set / mul は新しい版を返す。版の写しは根を写すだけ (O(1))。get、prod、all_prod はその版の値を読む。
// 節は型ごとに 1 本の置き場 (pool) に取り、0 番を空の節 (子が自分、値は ti()) にする。set / mul は根から葉への道の
// 節を写すので、1 回 O(log n) 個の節が増える。置き場が伸びたら pool_empty() が真を返すので、生きている版を渡して
// gc を呼ぶ。gc は渡した版から届く節だけを前へ詰める。reserve で先に取ると、配列を作り直す分を省ける。
namespace persistent_segment_tree_internal {
template <class M> struct Node {
 int l, r;
 typename M::T v;
};
}
template <class M> class PersistentSegmentTree {
 static_assert(Monoid<M>, "M にはモノイド (T, ti, op) が要る");
 using T= typename M::T;
 using Node= persistent_segment_tree_internal::Node<M>;
 using Pool= std::vector<Node, segment_tree_internal::Alloc<Node>>;
 static inline Pool pool_= [] {
  Pool p;
  p.push_back(Node{0, 0, M::ti()});
  return p;
 }();
 static inline std::size_t gc_at_= std::size_t(1) << 22;
 int n_= 0, sz_= 1, root_= 0;
 static int new_node(Node x) {
  const int i= pool_.size();
  pool_.push_back(x);
  return i;
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
 int set_rec(int p, int l, int r, int i, T x) const {
  const int q= new_node(pool_[p]);
  if(r - l == 1) {
   pool_[q].v= x;
  } else {
   const int m= (l + r) >> 1;
   if(i < m) pool_[q].l= set_rec(pool_[p].l, l, m, i, x);
   else pool_[q].r= set_rec(pool_[p].r, m, r, i, x);
   pool_[q].v= M::op(pool_[pool_[q].l].v, pool_[pool_[q].r].v);
  }
  return q;
 }
 int mul_rec(int p, int l, int r, int i, T x) const {
  const int q= new_node(pool_[p]);
  if(r - l == 1) {
   pool_[q].v= M::op(pool_[p].v, x);
  } else {
   const int m= (l + r) >> 1;
   if(i < m) pool_[q].l= mul_rec(pool_[p].l, l, m, i, x);
   else pool_[q].r= mul_rec(pool_[p].r, m, r, i, x);
   pool_[q].v= M::op(pool_[pool_[q].l].v, pool_[pool_[q].r].v);
  }
  return q;
 }
 T prod_rec(int p, int nl, int nr, int ql, int qr) const {
  if(qr <= nl || nr <= ql || p == 0) return M::ti();
  if(ql <= nl && nr <= qr) return pool_[p].v;
  const int m= (nl + nr) >> 1;
  return M::op(prod_rec(pool_[p].l, nl, m, ql, qr), prod_rec(pool_[p].r, m, nr, ql, qr));
 }
 // 版 b と版 a の差。覆う節ごとに sub(b の値, a の値) を求め、op でまとめる。
 template <class Sub> static T prod_diff_rec(int pa, int pb, int nl, int nr, int ql, int qr, const Sub& sub) {
  if(qr <= nl || nr <= ql) return M::ti();
  if(ql <= nl && nr <= qr) return sub(pool_[pb].v, pool_[pa].v);
  const int m= (nl + nr) >> 1;
  return M::op(prod_diff_rec(pool_[pa].l, pool_[pb].l, nl, m, ql, qr, sub), prod_diff_rec(pool_[pa].r, pool_[pb].r, m, nr, ql, qr, sub));
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
 PersistentSegmentTree()= default;
 explicit PersistentSegmentTree(int n): n_(n), sz_(std::bit_ceil(unsigned(std::max(n, 1)))) {}
 explicit PersistentSegmentTree(const std::vector<T>& a): n_(a.size()), sz_(std::bit_ceil(unsigned(std::max<std::size_t>(a.size(), 1)))) {
  if(n_) root_= build_rec(0, sz_, a, n_);
 }
 int size() const { return n_; }
 T get(int i) const {
  int p= root_, l= 0, r= sz_;
  while(r - l > 1) {
   const int m= (l + r) >> 1;
   if(i < m) p= pool_[p].l, r= m;
   else p= pool_[p].r, l= m;
  }
  return pool_[p].v;
 }
 T prod(int l, int r) const { return prod_rec(root_, 0, sz_, l, r); }
 T all_prod() const { return pool_[root_].v; }
 // 版 b と版 a の [l, r) の差。sub(版 b の節の値, 版 a の節の値) を覆う節ごとに求め、op でまとめる。
 // sub は T の引き算 (和や数え上げ)。a と b は同じ sz_ を持つこと。
 template <class Sub> static T prod_diff(const PersistentSegmentTree& a, const PersistentSegmentTree& b, int l, int r, const Sub& sub) {
  assert(a.sz_ == b.sz_);
  return prod_diff_rec(a.root_, b.root_, 0, b.sz_, l, r, sub);
 }
 // T が整数 (数え上げ) のとき、版 b と版 a の差の値のうち、x との xor が最大になる値の xor を返す。
 // a と b は同じ sz_ を持ち、x は [0, sz_) にあること。差が空なら 0 を返す。
 static T max_xor(const PersistentSegmentTree& a, const PersistentSegmentTree& b, T x) {
  static_assert(std::is_integral_v<T>, "max_xor は T が整数のときだけ");
  assert(a.sz_ == b.sz_);
  int pa= a.root_, pb= b.root_;
  const int lg= std::countr_zero(unsigned(b.sz_));
  T res= 0;
  for(int h= lg - 1; h >= 0; --h) {
   const int xb= int((x >> h) & 1), want= xb ^ 1;
   const int c0= int(pool_[pool_[pb].l].v - pool_[pool_[pa].l].v);
   const int c1= int(pool_[pool_[pb].r].v - pool_[pool_[pa].r].v);
   const int choose= (want == 0 ? (c0 > 0 ? 0 : 1) : (c1 > 0 ? 1 : 0));
   if(choose != xb) res|= T(1) << h;
   if(choose == 0) pa= pool_[pa].l, pb= pool_[pb].l;
   else pa= pool_[pa].r, pb= pool_[pb].r;
  }
  return res;
 }
 // T が整数 (数え上げ) のとき、版 b と版 a の差の k 番目 (0 始まり) の位置。a と b は同じ sz_ を持つこと。
 // 値の順位を葉に持つ木で、版 i を先頭 i 個を入れた木にすると、range-kth-smallest になる。
 static int kth(const PersistentSegmentTree& a, const PersistentSegmentTree& b, int k) {
  static_assert(std::is_integral_v<T>, "kth は T が整数のときだけ");
  assert(a.sz_ == b.sz_);
  int pa= a.root_, pb= b.root_, l= 0, r= b.sz_;
  while(r - l > 1) {
   const int m= (l + r) >> 1;
   const int cl= int(pool_[pool_[pb].l].v - pool_[pool_[pa].l].v);
   if(k < cl) pa= pool_[pa].l, pb= pool_[pb].l, r= m;
   else k-= cl, pa= pool_[pa].r, pb= pool_[pb].r, l= m;
  }
  return l;
 }
 PersistentSegmentTree set(int i, T x) const {
  PersistentSegmentTree res= *this;
  res.root_= set_rec(root_, 0, sz_, i, x);
  return res;
 }
 PersistentSegmentTree mul(int i, T x) const {
  PersistentSegmentTree res= *this;
  res.root_= mul_rec(root_, 0, sz_, i, x);
  return res;
 }
 static void reserve(std::size_t n) {
  if(pool_.capacity() < n + 1) pool_.reserve(std::max(n + 1, 2 * pool_.capacity()));
 }
 static bool pool_empty() { return pool_.size() >= gc_at_; }
 static void gc(std::vector<PersistentSegmentTree>& vs) {
  std::vector<int*> roots;
  roots.reserve(vs.size());
  for(auto& v: vs) roots.push_back(&v.root_);
  gc_roots(roots);
 }
 template <class... Ps> static void gc(PersistentSegmentTree& v, Ps&... vs) {
  std::vector<int*> roots{&v.root_, &vs.root_...};
  gc_roots(roots);
 }
};
