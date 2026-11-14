#pragma once
#include <algorithm>
#include <bit>
#include <cassert>
#include <limits>
#include <type_traits>
#include <vector>
#include "neo/data_structure/SegmentTree.hpp"
// 区間 chmin、区間 chmax、区間加算と、区間の和、最大、最小 monoid::ChminChmaxAddSum<V> (Segment Tree Beats)。V は 64 bit 以下の
// 符号付き整数。作用 {lb, ub, add} は x -> min(max(x, lb), ub) + add で、chmin(x)、chmax(x)、add(x) で作る。値は
// {sum, max, min} で、sum は long long。
// SegmentTree に渡すと、この 3 つの操作に特化した Beats になる (procon-judge の yosupo-range-chmin-chmax-add-range-sum の
// beats)。節は最大と 2 番目の最大、最小と 2 番目の最小、和、まだ配っていない加算、最大と最小の個数を 64 byte に並べる。
// 根から再帰で下り、chmin は最大だけを、chmax は最小だけを動かせば済む節で止める。配るときは加算を先に当て、子の最大が
// 親の最大を超えていれば子にも chmin を当てる (chmax も同じ)。chmin と chmax の作用は節に持たず、親の最大と最小で表す。
// 値は絶対値 2^61 未満で使う (2^61 を番兵にする)。
namespace monoid {
template <class V> struct ChminChmaxAddSum {
 static_assert(std::is_integral_v<V> && std::is_signed_v<V> && sizeof(V) <= 8, "V は 64 bit 以下の符号付き整数");
 struct T {
  long long sum;
  V max, min;
  T()= default;
  constexpr T(long long s, V mx, V mn): sum(s), max(mx), min(mn) {}
  // 1 要素の値。
  constexpr T(V x): sum(x), max(x), min(x) {}
  constexpr bool operator==(const T&) const= default;
 };
 struct E {
  V lb, ub, add;
 };
 static constexpr E chmin(V x) { return {std::numeric_limits<V>::lowest(), x, 0}; }
 static constexpr E chmax(V x) { return {x, std::numeric_limits<V>::max(), 0}; }
 static constexpr E add(V x) { return {std::numeric_limits<V>::lowest(), std::numeric_limits<V>::max(), x}; }
 static constexpr T ti() { return {0, std::numeric_limits<V>::lowest(), std::numeric_limits<V>::max()}; }
 static constexpr T op(const T& a, const T& b) { return {a.sum + b.sum, a.max < b.max ? b.max : a.max, b.min < a.min ? b.min : a.min}; }
};
}
namespace segtree_beats_internal {
using i64= long long;
template <class V> class Beats {
 using M= monoid::ChminChmaxAddSum<V>;
 using T= typename M::T;
 using E= typename M::E;
 static constexpr i64 INF= i64(1) << 61;
 struct alignas(64) Node {
  i64 mx, mx2, mn, mn2, sum, lz;
  int mxc, mnc;
 };
 int n= 0, lg= 0, sz= 1;
 segment_tree_internal::Vec<Node> t;
 static Node leaf(i64 x) { return {x, -INF, x, INF, x, 0, 1, 1}; }
 void pull(int k) {
  const Node &L= t[k << 1], &R= t[k << 1 | 1];
  Node& X= t[k];
  X.sum= L.sum + R.sum;
  if(L.mx > R.mx) X.mx= L.mx, X.mxc= L.mxc, X.mx2= std::max(L.mx2, R.mx);
  else if(L.mx < R.mx) X.mx= R.mx, X.mxc= R.mxc, X.mx2= std::max(L.mx, R.mx2);
  else X.mx= L.mx, X.mxc= L.mxc + R.mxc, X.mx2= std::max(L.mx2, R.mx2);
  if(L.mn < R.mn) X.mn= L.mn, X.mnc= L.mnc, X.mn2= std::min(L.mn2, R.mn);
  else if(L.mn > R.mn) X.mn= R.mn, X.mnc= R.mnc, X.mn2= std::min(L.mn, R.mn2);
  else X.mn= L.mn, X.mnc= L.mnc + R.mnc, X.mn2= std::min(L.mn2, R.mn2);
 }
 void apply_add(int k, int len, i64 x) {
  Node& X= t[k];
  X.sum+= x * len, X.mx+= x, X.mn+= x, X.lz+= x;
  if(X.mx2 != -INF) X.mx2+= x;
  if(X.mn2 != INF) X.mn2+= x;
 }
 // X.mx2 < x < X.mx のときだけ呼ぶ。
 void apply_chmin(int k, i64 x) {
  Node& X= t[k];
  X.sum-= (X.mx - x) * X.mxc;
  if(X.mn == X.mx) X.mn= x;
  else if(X.mn2 == X.mx) X.mn2= x;
  X.mx= x;
 }
 // X.mn < x < X.mn2 のときだけ呼ぶ。
 void apply_chmax(int k, i64 x) {
  Node& X= t[k];
  X.sum+= (x - X.mn) * X.mnc;
  if(X.mx == X.mn) X.mx= x;
  else if(X.mx2 == X.mn) X.mx2= x;
  X.mn= x;
 }
 void push(int k, int len) {
  Node& X= t[k];
  if(X.lz) apply_add(k << 1, len >> 1, X.lz), apply_add(k << 1 | 1, len >> 1, X.lz), X.lz= 0;
  for(int c= k << 1; c <= (k << 1 | 1); ++c) {
   if(t[c].mx > X.mx) apply_chmin(c, X.mx);
   if(t[c].mn < X.mn) apply_chmax(c, X.mn);
  }
 }
 void chmin_rec(int k, int nl, int nr, int l, int r, i64 x) {
  if(r <= nl || nr <= l || t[k].mx <= x) return;
  if(l <= nl && nr <= r && t[k].mx2 < x) return apply_chmin(k, x);
  push(k, nr - nl);
  const int m= (nl + nr) >> 1;
  chmin_rec(k << 1, nl, m, l, r, x), chmin_rec(k << 1 | 1, m, nr, l, r, x);
  pull(k);
 }
 void chmax_rec(int k, int nl, int nr, int l, int r, i64 x) {
  if(r <= nl || nr <= l || t[k].mn >= x) return;
  if(l <= nl && nr <= r && t[k].mn2 > x) return apply_chmax(k, x);
  push(k, nr - nl);
  const int m= (nl + nr) >> 1;
  chmax_rec(k << 1, nl, m, l, r, x), chmax_rec(k << 1 | 1, m, nr, l, r, x);
  pull(k);
 }
 void add_rec(int k, int nl, int nr, int l, int r, i64 x) {
  if(r <= nl || nr <= l) return;
  if(l <= nl && nr <= r) return apply_add(k, nr - nl, x);
  push(k, nr - nl);
  const int m= (nl + nr) >> 1;
  add_rec(k << 1, nl, m, l, r, x), add_rec(k << 1 | 1, m, nr, l, r, x);
  pull(k);
 }
 T prod_rec(int k, int nl, int nr, int l, int r) {
  if(l <= nl && nr <= r) return val(k);
  push(k, nr - nl);
  const int m= (nl + nr) >> 1;
  if(r <= m) return prod_rec(k << 1, nl, m, l, r);
  if(m <= l) return prod_rec(k << 1 | 1, m, nr, l, r);
  return M::op(prod_rec(k << 1, nl, m, l, r), prod_rec(k << 1 | 1, m, nr, l, r));
 }
 // 節の値。要素の無い節 (大きさを 2 の冪に揃えた余り) は ti() にする。
 T val(int k) const { return t[k].mxc ? T{t[k].sum, V(t[k].mx), V(t[k].mn)} : M::ti(); }
 void push_path(int p) {
  for(int j= lg; j >= 1; --j) push(p >> j, 1 << j);
 }
 template <class F> void build(const F& init) {
  sz= std::bit_ceil(unsigned(std::max(n, 1))), lg= std::countr_zero(unsigned(sz));
  t.resize(2 * sz);
  for(int i= 0; i < sz; ++i) t[sz + i]= i < n ? leaf(i64(T(init(i)).sum)) : Node{-INF, -INF, INF, INF, 0, 0, 0, 0};
  for(int k= sz - 1; k >= 1; --k) t[k].lz= 0, pull(k);
 }
public:
 Beats()= default;
 explicit Beats(int n): n(n) {
  build([](int) { return T(V(0)); });
 }
 Beats(int n, T x): n(n) {
  build([&](int) { return x; });
 }
 Beats(const std::vector<T>& a): n(a.size()) {
  build([&](int i) { return a[i]; });
 }
 template <class F, std::enable_if_t<std::is_invocable_r_v<T, const F&, int>, int> = 0> Beats(int n, const F& init): n(n) { build(init); }
 int size() const { return n; }
 T get(int i) {
  const int p= sz + i;
  push_path(p);
  return val(p);
 }
 T operator[](int i) { return get(i); }
 // a[i] = x (x は 1 要素の値)。
 void set(int i, T x) {
  const int p= sz + i;
  push_path(p);
  t[p]= leaf(x.sum);
  for(int k= p >> 1; k >= 1; k>>= 1) pull(k);
 }
 T prod(int l, int r) { return l < r ? prod_rec(1, 0, sz, l, r) : M::ti(); }
 T all_prod() const { return n ? val(1) : M::ti(); }
 // a[l..r) に x -> min(max(x, lb), ub) + add を当てる。
 void apply(int l, int r, const E& f) {
  if(l >= r) return;
  if(f.lb != std::numeric_limits<V>::lowest()) chmax_rec(1, 0, sz, l, r, f.lb);
  if(f.ub != std::numeric_limits<V>::max()) chmin_rec(1, 0, sz, l, r, f.ub);
  if(f.add) add_rec(1, 0, sz, l, r, f.add);
 }
 // check(prod(l, r)) が真になる最大の r (ACL と同じ。check(ti()) は真、r について単調)。
 template <class C> int max_right(int l, const C& check) {
  assert(0 <= l && l <= n && check(M::ti()));
  if(l == n) return n;
  l+= sz;
  push_path(l);
  T s= M::ti();
  int len= 1;
  do {
   while(!(l & 1)) l>>= 1, len<<= 1;
   if(!check(M::op(s, val(l)))) {
    while(l < sz) {
     push(l, len);
     l<<= 1, len>>= 1;
     if(T u= M::op(s, val(l)); check(u)) s= u, ++l;
    }
    return std::min(l - sz, n);
   }
   s= M::op(s, val(l)), ++l;
  } while((l & -l) != l);
  return n;
 }
 // check(prod(l, r)) が真になる最小の l。
 template <class C> int min_left(int r, const C& check) {
  assert(0 <= r && r <= n && check(M::ti()));
  if(r == 0) return 0;
  r+= sz;
  push_path(r - 1);
  T s= M::ti();
  int len= 1;
  do {
   --r;
   while(r > 1 && (r & 1)) r>>= 1, len<<= 1;
   if(!check(M::op(val(r), s))) {
    while(r < sz) {
     push(r, len);
     r= 2 * r + 1, len>>= 1;
     if(T u= M::op(val(r), s); check(u)) s= u, --r;
    }
    return r + 1 - sz;
   }
   s= M::op(val(r), s);
  } while((r & -r) != r);
  return 0;
 }
};
}  // namespace segtree_beats_internal
template <class V> class SegmentTree<monoid::ChminChmaxAddSum<V>>: public segtree_beats_internal::Beats<V> {
public:
 using segtree_beats_internal::Beats<V>::Beats;
};
