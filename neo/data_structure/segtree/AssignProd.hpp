#pragma once
#include <algorithm>
#include <bit>
#include <cassert>
#include <type_traits>
#include <vector>
#include "neo/data_structure/SegmentTree.hpp"
// 区間代入と区間の積 monoid::AssignProd<M>。M は SegmentTree と同じモノイド (T、ti、op) で、可換でなくてよい
// (1 次関数の合成なら monoid::Composite、neo/data_structure/segtree/Composite.hpp)。作用は T の値の代入。
// SegmentTree に渡すと、二分木の遅延セグ木で、長さ 2^j の節に x を代入した値 x^{2^j} を、代入のたびに 2 乗を重ねて表に
// 作り、節の作用には表の位置を持たせる (procon-judge の yosupo-range-set-range-composite の lazy_pow)。表は代入が覆う節の
// いちばん高い段までだけ作り、配るときは子の段の冪を表から引くので、作用を当てるたびに冪を求め直さない。表は代入ごとに
// 段の数だけ伸びるので、32 MB (節の数の 4 倍より小さければそちら) を超えたら、残っている作用が引く分だけを新しい表へ
// 詰め直す (残る分は節の数の 2 倍以下)。
namespace monoid {
template <class M> struct AssignProd {
 using T= typename M::T;
 using E= T;
 static T ti() { return M::ti(); }
 static T op(const T& a, const T& b) { return M::op(a, b); }
 // 汎用の形で使うときの作用 (x を len 個並べた積を二分累乗で)。SegmentTree はこの構造体を受けると下の形に替わる。
 static void mp(T& v, const E& x, int len) {
  T r= M::ti(), b= x;
  for(; len; len>>= 1, b= M::op(b, b))
   if(len & 1) r= M::op(r, b);
  v= r;
 }
 static void cp(E& pre, const E& suf) { pre= suf; }
};
}
namespace segtree_assign_prod_internal {
using u32= unsigned;
template <class M> class Pow {
 using T= typename M::T;
 int n= 0, lg= 0, sz= 1;
 segment_tree_internal::Vec<T> d;  // 節 1..2sz-1 の値
 segment_tree_internal::Vec<u32> lz;  // 節 k (段 i >= 1) に残る代入の表の位置 (0 なら無し)。子孫は表の o + i - 1 以下を引く
 segment_tree_internal::Vec<T> pool;  // 代入 o の x^{2^j} を pool[o + j] に置く。位置 0 は使わない
 std::size_t cap= 0;   // 表がこの長さを超えたら詰め直す
 void update(int k) { d[k]= M::op(d[k << 1], d[k << 1 | 1]); }
 // 段 i の節 k の代入を子へ配る。
 void push(int k, int i) {
  if(const u32 o= lz[k]) {
   const T v= pool[o + i - 1];
   d[k << 1]= v, d[k << 1 | 1]= v;
   if(i > 1) lz[k << 1]= o, lz[k << 1 | 1]= o;
   lz[k]= 0;
  }
 }
 [[gnu::always_inline]] void push_ends(int l, int r) {
  for(int i= lg; i >= 1; --i) {
   if(((l >> i) << i) != l) push(l >> i, i);
   if(((r >> i) << i) != r) push((r - 1) >> i, i);
  }
 }
 void push_path(int p) {
  for(int i= lg; i >= 1; --i) push(p >> i, i);
 }
 // 残っている代入が引く分だけを新しい表へ詰める。節を根から順に見ると、同じ代入を引く節のうちいちばん高い段の節に
 // 先に当たるので、そこで段の数だけ写せば、それより下の節が引く分も含む。
 void compact() {
  segment_tree_internal::Vec<u32> to(pool.size(), 0);
  segment_tree_internal::Vec<T> np;
  np.reserve(2 * sz + 1);
  np.push_back(pool[0]);
  for(int k= 1; k < sz; ++k)
   if(const u32 o= lz[k]) {
    if(!to[o]) {
     const int i= lg - std::bit_width(unsigned(k)) + 1;
     to[o]= u32(np.size());
     np.insert(np.end(), pool.begin() + o, pool.begin() + o + i);
    }
    lz[k]= to[o];
   }
  pool.swap(np);
 }
 template <class F> void build(const F& init) {
  sz= std::bit_ceil(unsigned(std::max(n, 1))), lg= std::countr_zero(unsigned(sz));
  d.assign(2 * sz, M::ti());
  for(int i= 0; i < n; ++i) d[sz + i]= init(i);
  for(int k= sz - 1; k >= 1; --k) update(k);
  lz.assign(sz, 0);
  pool.assign(1, M::ti());
  // 詰め直したあとは節の数の 2 倍以下なので、節の数の 4 倍までは詰めない。T が小さいときは 32 MB まで伸ばす。
  cap= std::max<std::size_t>(4 * std::size_t(sz) + 64, (std::size_t(32) << 20) / sizeof(T));
 }
public:
 Pow()= default;
 explicit Pow(int n): n(n) {
  build([](int) { return M::ti(); });
 }
 Pow(int n, T x): n(n) {
  build([&](int) { return x; });
 }
 Pow(const std::vector<T>& a): n(a.size()) {
  build([&](int i) { return a[i]; });
 }
 template <class F, std::enable_if_t<std::is_invocable_r_v<T, const F&, int>, int> = 0> Pow(int n, const F& init): n(n) { build(init); }
 int size() const { return n; }
 // 根に近い代入ほど新しいので、根から下りて初めに当たった代入の値を返す (書き込みなし)。
 T get(int i) {
  const int p= sz + i;
  for(int j= lg; j >= 1; --j)
   if(const u32 o= lz[p >> j]) return pool[o];
  return d[p];
 }
 T operator[](int i) { return get(i); }
 void set(int i, T x) {
  const int p= sz + i;
  push_path(p);
  d[p]= x;
  for(int k= p >> 1; k >= 1; k>>= 1) update(k);
 }
 void mul(int i, T x) { set(i, M::op(get(i), x)); }
 T prod(int l, int r) {
  l+= sz, r+= sz;
  push_ends(l, r);
  if constexpr(std::is_arithmetic_v<T>) {
   T sl= M::ti(), sr= M::ti();
   for(; l < r; l>>= 1, r>>= 1) {
    if(l & 1) sl= M::op(sl, d[l++]);
    if(r & 1) sr= M::op(d[--r], sr);
   }
   return M::op(sl, sr);
  }
  // 覆う節を分岐なしで集めてから、左と右の 2 本の鎖で畳む (汎用の形と同じ)。
  int L[32], R[32], nl= 0, nr= 0;
  for(; l < r; l>>= 1, r>>= 1) {
   L[nl]= l, nl+= l & 1, l+= l & 1;
   R[nr]= r - 1, nr+= r & 1, r-= r & 1;
  }
  T sl= M::ti(), sr= M::ti();
  for(int i= 0; i < nl; ++i) sl= M::op(sl, d[L[i]]);
  for(int i= 0; i < nr; ++i) sr= M::op(d[R[i]], sr);
  return M::op(sl, sr);
 }
 T all_prod() const { return d[1]; }
 // a[l..r) を x にする。
 void apply(int l, int r, T x) {
  if(l >= r) return;
  l+= sz, r+= sz;
  push_ends(l, r);
  if(pool.size() + lg + 2 > cap) compact();
  const u32 o= u32(pool.size());
  pool.push_back(x);
  for(int a= l, b= r, i= 0; a < b; a>>= 1, b>>= 1, ++i) {
   if(i) {
    const T p= pool.back();
    pool.push_back(M::op(p, p));
   }
   if(a & 1) {
    d[a]= pool[o + i];
    if(i) lz[a]= o;
    ++a;
   }
   if(b & 1) {
    --b;
    d[b]= pool[o + i];
    if(i) lz[b]= o;
   }
  }
  for(int i= 1; i <= lg; ++i) {
   if(((l >> i) << i) != l) update(l >> i);
   if(((r >> i) << i) != r) update((r - 1) >> i);
  }
 }
 // check(prod(l, r)) が真になる最大の r (ACL と同じ。check(ti()) は真、r について単調)。
 template <class C> int max_right(int l, const C& check) {
  assert(0 <= l && l <= n && check(M::ti()));
  if(l == n) return n;
  l+= sz;
  push_path(l);
  T s= M::ti();
  int i= 0;
  do {
   while(!(l & 1)) l>>= 1, ++i;
   if(!check(M::op(s, d[l]))) {
    while(l < sz) {
     push(l, i);
     l<<= 1, --i;
     if(T t= M::op(s, d[l]); check(t)) s= t, ++l;
    }
    return std::min(l - sz, n);
   }
   s= M::op(s, d[l]), ++l;
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
  int i= 0;
  do {
   --r;
   while(r > 1 && (r & 1)) r>>= 1, ++i;
   if(!check(M::op(d[r], s))) {
    while(r < sz) {
     push(r, i);
     r= 2 * r + 1, --i;
     if(T t= M::op(d[r], s); check(t)) s= t, --r;
    }
    return r + 1 - sz;
   }
   s= M::op(d[r], s);
  } while((r & -r) != r);
  return 0;
 }
};
}  // namespace segtree_assign_prod_internal
template <class M> class SegmentTree<monoid::AssignProd<M>>: public segtree_assign_prod_internal::Pow<M> {
public:
 using segtree_assign_prod_internal::Pow<M>::Pow;
};
