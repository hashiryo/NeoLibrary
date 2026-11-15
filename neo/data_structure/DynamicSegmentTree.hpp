#pragma once
#include <algorithm>
#include <bit>
#include <cassert>
#include <cstdint>
#include <type_traits>
#include <vector>
#include "neo/data_structure/SegmentTree.hpp"
// 添字の範囲 [0, n) が大きい (n は long long) セグメント木。触った所だけ節を作る。構造体 M の決まりは SegmentTree と同じで
// (neo/data_structure/SegmentTree.hpp)、はじめはどの要素も ti() (作用だけなら T())。渡す構造体で中身が替わる。
//   モノイドだけ: パトリシア木 (枝分かれの無い道を 1 つの節に畳んだ二分木)。節は [lo, lo + 2^lv) を受け持ち、内部の節は必ず
//     子を 2 つ持つので、節の数は set した添字の数の 2 倍未満。区間の積は、l と r - 1 が同じ子にある間は根から下り、分かれた節
//     から先は左の道で l 以上の側、右の道で r - 1 以下の側を集める。
//   作用があるとき: 深さ ceil(log2 n) の二分木で、子の無い節は「その下の要素は ti() のまま」を表す。子を作るのは、下りる道で
//     子が無いときと、節に残った作用を子へ配るときだけ。作用は値の等しい要素の並びを値の等しい並びへ移すので、区間に作用を
//     当てても、子ができるのは区間の端の道だけになる。区間の積は書き込みをせず、配っていない作用は戻りに当てる。
// 節は配列に置いて子を添字で持ち、1 MB 以上の配列は 2 MB 境界に置いて huge page を頼む。mp が bool を返す形 (Beats) は受けない。
// 区間の長さを受ける mp には長さを long long のまま渡すので、n が 2^31 を越えるなら long long で受けること。
// procon-judge の yosupo-point-set-range-composite-large-array の patricia と、yosupo-range-affine-range-sum-large-array の
// dyn_uniform を、構造体を受ける形にしたもの。
namespace dynamic_segment_tree_internal {
using segment_tree_internal::Flag;
using segment_tree_internal::Vec;
using u32= unsigned;
using u64= unsigned long long;
template <class M> class Patricia {
 using T= typename M::T;
 struct Node {
  u64 lo;
  int lv;
  u32 ch[2];
  T val;
 };
 long long n= 0;
 Vec<Node> t;  // t[0] は無い節
 u32 root= 0;
 static u64 last(const Node& x) { return x.lo + ((u64(1) << x.lv) - 1); }
 void pull(u32 v) { t[v].val= M::op(t[t[v].ch[0]].val, t[t[v].ch[1]].val); }
public:
 Patricia(): Patricia(0) {}
 explicit Patricia(long long n): n(n) { t.reserve(64), t.push_back({0, 0, {0, 0}, M::ti()}); }
 long long size() const { return n; }
 T get(long long i) const {
  const u64 p= i;
  for(u32 v= root; v;) {
   const Node& x= t[v];
   if((p ^ x.lo) >> x.lv) break;
   if(!x.lv) return x.val;
   v= x.ch[p >> (x.lv - 1) & 1];
  }
  return M::ti();
 }
 T operator[](long long i) const { return get(i); }
 void set(long long i, T x) {
  const u64 p= i;
  // 1 回で作る節は 2 個以下。先に場所を取っておけば、下りる途中で配列が動かない。
  if(t.capacity() < t.size() + 2) t.reserve(2 * t.capacity() + 2);
  if(!root) {
   root= u32(t.size()), t.push_back({p, 0, {0, 0}, x});
   return;
  }
  u32 path[64];
  int np= 0;
  u32* slot= &root;
  for(;;) {
   const u32 v= *slot;
   Node& y= t[v];
   if((p ^ y.lo) >> y.lv) {
    // p は v の範囲の外。v と p の葉を子に持つ節を v の場所に置く。
    const int lv= std::bit_width(p ^ y.lo);
    const u32 w= u32(t.size()), b= p >> (lv - 1) & 1;
    Node nw{p & ~((u64(1) << lv) - 1), lv, {0, 0}, M::ti()};
    nw.ch[b]= w + 1, nw.ch[b ^ 1]= v;
    t.push_back(nw), t.push_back({p, 0, {0, 0}, x});
    *slot= w, pull(w);
    break;
   }
   if(!y.lv) {
    y.val= x;
    break;
   }
   path[np++]= v, slot= &y.ch[p >> (y.lv - 1) & 1];
  }
  while(np) pull(path[--np]);
 }
 void mul(long long i, T x) { set(i, M::op(get(i), x)); }
 T prod(long long l_, long long r_) const {
  if(l_ >= r_) return M::ti();
  const u64 l= l_, r= r_ - 1;
  u32 v= root;
  for(;;) {
   if(!v) return M::ti();
   const Node& x= t[v];
   if(r < x.lo || last(x) < l) return M::ti();
   if(l <= x.lo && last(x) <= r) return x.val;
   const u64 mid= x.lo + (u64(1) << (x.lv - 1));
   if(r < mid) v= x.ch[0];
   else if(l >= mid) v= x.ch[1];
   else break;
  }
  // 左の子で l 以上の側。全部入る節を浅い方から積み、深い方から畳む。
  u32 st[64];
  int ns= 0;
  T sl= M::ti();
  for(u32 u= t[v].ch[0];;) {
   const Node& x= t[u];
   if(x.lo >= l) {
    sl= x.val;
    break;
   }
   if(last(x) < l) break;
   if(l < x.lo + (u64(1) << (x.lv - 1))) st[ns++]= x.ch[1], u= x.ch[0];
   else u= x.ch[1];
  }
  while(ns) sl= M::op(sl, t[st[--ns]].val);
  // 右の子で r - 1 以下の側。左から右の順に出会うので、そのまま畳む。
  T sr= M::ti();
  for(u32 u= t[v].ch[1];;) {
   const Node& x= t[u];
   if(last(x) <= r) {
    sr= M::op(sr, x.val);
    break;
   }
   if(x.lo > r) break;
   if(r >= x.lo + (u64(1) << (x.lv - 1))) sr= M::op(sr, t[x.ch[0]].val), u= x.ch[1];
   else u= x.ch[0];
  }
  return M::op(sl, sr);
 }
 T all_prod() const { return root ? t[root].val : M::ti(); }
 // check(prod(l, r)) が真になる最大の r (ACL と同じ。check(ti()) は真、r について単調)。
 template <class C> long long max_right(long long l, const C& check) const {
  assert(0 <= l && l <= n && check(M::ti()));
  T acc= M::ti();
  long long res= n;
  auto rec= [&](auto&& self, u32 v) -> bool {
   if(!v) return true;
   const Node& x= t[v];
   if(last(x) < u64(l)) return true;
   if(x.lo >= u64(l)) {
    if(T s= M::op(acc, x.val); check(s)) return acc= s, true;
    if(!x.lv) return res= x.lo, false;
   }
   return self(self, x.ch[0]) && self(self, x.ch[1]);
  };
  if(l < n) rec(rec, root);
  return std::min(res, n);
 }
 // check(prod(l, r)) が真になる最小の l。
 template <class C> long long min_left(long long r, const C& check) const {
  assert(0 <= r && r <= n && check(M::ti()));
  T acc= M::ti();
  long long res= 0;
  auto rec= [&](auto&& self, u32 v) -> bool {
   if(!v) return true;
   const Node& x= t[v];
   if(x.lo >= u64(r)) return true;
   if(last(x) < u64(r)) {
    if(T s= M::op(x.val, acc); check(s)) return acc= s, true;
    if(!x.lv) return res= x.lo + 1, false;
   }
   return self(self, x.ch[1]) && self(self, x.ch[0]);
  };
  if(r > 0) rec(rec, root);
  return res;
 }
};
template <class M> class Lazy {
 using T= typename M::T;
 using E= typename M::E;
 using Ea= std::conditional_t<std::is_trivially_copyable_v<E> && sizeof(E) <= 32, E, const E&>;
 static constexpr bool MO= segment_tree_internal::monoid<M>::value;
 static_assert(segment_tree_internal::action<M>::value);
 struct Node {
  T val;  // モノイドがあれば部分木の積 (節の作用まで当てた値)、作用だけなら葉の値
  E lz;   // 子へ配っていない作用
  u32 ch[2];
  Flag fl;
 };
 long long n= 0;
 int H= 1;     // 根は [0, 2^H) を受け持つ
 Vec<Node> t;  // t[0] は無い節 (値は ti())、t[1] は根
 static T id() {
  if constexpr(MO) return M::ti();
  else return T();
 }
 static void map(T& v, const E& x, u64 len) {
  static_assert(!std::is_invocable_r_v<bool, decltype(M::mp), T&, const E&, int> && !std::is_invocable_r_v<bool, decltype(M::mp), T&, const E&>, "Beats は受けない");
  if constexpr(std::is_invocable_v<decltype(M::mp), T&, const E&, long long>) M::mp(v, x, (long long)len);
  else M::mp(v, x);
 }
 void reserve(std::size_t k) {
  if(t.capacity() < t.size() + k) t.reserve(2 * t.capacity() + k);
 }
 u32 make() { return t.push_back({id(), E(), {0, 0}, Flag::off}), u32(t.size() - 1); }
 u32 child(u32 v, int c) {
  if(!t[v].ch[c]) {
   const u32 k= make();
   t[v].ch[c]= k;
  }
  return t[v].ch[c];
 }
 // 長さ len の節 v の全体に x を当てる。
 void all_apply(u32 v, Ea x, u64 len) {
  Node& y= t[v];
  if constexpr(MO) map(y.val, x, len);
  else if(len == 1) map(y.val, x, 1);
  if(len > 1) {
   if(y.fl == Flag::on) M::cp(y.lz, x);
   else y.lz= x;
   y.fl= Flag::on;
  }
 }
 void push(u32 v, u64 len) {
  if(t[v].fl == Flag::on) {
   const E x= t[v].lz;
   all_apply(child(v, 0), x, len >> 1), all_apply(child(v, 1), x, len >> 1);
   t[v].fl= Flag::off;
  }
 }
 void pull(u32 v) {
  if constexpr(MO) t[v].val= M::op(t[t[v].ch[0]].val, t[t[v].ch[1]].val);
 }
 void apply_rec(u32 v, u64 lo, int lv, u64 l, u64 r, Ea x) {
  const u64 len= u64(1) << lv;
  if(l <= lo && lo + len <= r) return all_apply(v, x, len);
  push(v, len);
  const u64 mid= lo + (len >> 1);
  if(l < mid) apply_rec(child(v, 0), lo, lv - 1, l, r, x);
  if(r > mid) apply_rec(child(v, 1), mid, lv - 1, l, r, x);
  pull(v);
 }
 T prod_rec(u32 v, u64 lo, int lv, u64 l, u64 r) const {
  if(!v) return M::ti();
  const Node& y= t[v];
  const u64 len= u64(1) << lv;
  if(l <= lo && lo + len <= r) return y.val;
  const u64 mid= lo + (len >> 1);
  T s= M::ti();
  if(l < mid) s= prod_rec(y.ch[0], lo, lv - 1, l, r);
  if(r > mid) s= M::op(s, prod_rec(y.ch[1], mid, lv - 1, l, r));
  if(y.fl == Flag::on) map(s, y.lz, std::min(r, lo + len) - std::max(l, lo));
  return s;
 }
public:
 Lazy(): Lazy(0) {}
 explicit Lazy(long long n): n(n), H(std::max(1, int(std::bit_width(u64(std::max(n, 1LL) - 1))))) {
  t.reserve(64), t.push_back({id(), E(), {0, 0}, Flag::off}), t.push_back({id(), E(), {0, 0}, Flag::off});
 }
 long long size() const { return n; }
 T get(long long i) const {
  u32 path[64];
  int d= 0;
  u32 v= 1;
  for(int lv= H; lv > 0 && v; --lv) path[d++]= v, v= t[v].ch[u64(i) >> (lv - 1) & 1];
  T x= v ? t[v].val : id();
  // 根から葉への道では上の節ほど新しい作用なので、下の節から順に当てる。
  while(d--)
   if(t[path[d]].fl == Flag::on) map(x, t[path[d]].lz, 1);
  return x;
 }
 T operator[](long long i) const { return get(i); }
 void set(long long i, T x) {
  reserve(2 * H + 2);
  u32 path[64];
  int d= 0;
  u32 v= 1;
  for(int lv= H; lv > 0; --lv) path[d++]= v, push(v, u64(1) << lv), v= child(v, u64(i) >> (lv - 1) & 1);
  t[v].val= x;
  while(d) pull(path[--d]);
 }
 void mul(long long i, T x) {
  static_assert(MO, "mul はモノイドがあるときだけ");
  set(i, M::op(get(i), x));
 }
 T prod(long long l, long long r) const {
  static_assert(MO, "prod はモノイドがあるときだけ");
  if(l >= r) return M::ti();
  return prod_rec(1, 0, H, l, r);
 }
 T all_prod() const {
  static_assert(MO, "all_prod はモノイドがあるときだけ");
  return t[1].val;
 }
 void apply(long long l, long long r, E x) {
  if(l >= r) return;
  // 1 回で作る節は、両端の道で段ごとに 2 個ずつまで。先に場所を取っておけば、下りる途中で配列が動かない。
  reserve(4 * H + 4);
  apply_rec(1, 0, H, l, r, x);
 }
 // check(prod(l, r)) が真になる最大の r (ACL と同じ。check(ti()) は真、r について単調)。
 template <class C> long long max_right(long long l, const C& check) {
  static_assert(MO, "max_right はモノイドがあるときだけ");
  assert(0 <= l && l <= n && check(M::ti()));
  reserve(2 * H + 2);
  T acc= M::ti();
  long long res= n;
  auto rec= [&](auto&& self, u32 v, u64 lo, int lv) -> bool {
   if(!v) return true;  // 触っていない部分木の積は ti()
   const u64 len= u64(1) << lv;
   if(lo + len <= u64(l)) return true;
   if(lo >= u64(l)) {
    if(T s= M::op(acc, t[v].val); check(s)) return acc= s, true;
    if(!lv) return res= lo, false;
   }
   push(v, len);
   return self(self, t[v].ch[0], lo, lv - 1) && self(self, t[v].ch[1], lo + (len >> 1), lv - 1);
  };
  if(l < n) rec(rec, 1, 0, H);
  return std::min(res, n);
 }
 // check(prod(l, r)) が真になる最小の l。
 template <class C> long long min_left(long long r, const C& check) {
  static_assert(MO, "min_left はモノイドがあるときだけ");
  assert(0 <= r && r <= n && check(M::ti()));
  reserve(2 * H + 2);
  T acc= M::ti();
  long long res= 0;
  auto rec= [&](auto&& self, u32 v, u64 lo, int lv) -> bool {
   if(!v) return true;
   const u64 len= u64(1) << lv;
   if(lo >= u64(r)) return true;
   if(lo + len <= u64(r)) {
    if(T s= M::op(t[v].val, acc); check(s)) return acc= s, true;
    if(!lv) return res= lo + 1, false;
   }
   push(v, len);
   return self(self, t[v].ch[1], lo + (len >> 1), lv - 1) && self(self, t[v].ch[0], lo, lv - 1);
  };
  if(r > 0) rec(rec, 1, 0, H);
  return res;
 }
};
}  // namespace dynamic_segment_tree_internal
template <class M> class DynamicSegmentTree: public std::conditional_t<segment_tree_internal::action<M>::value, dynamic_segment_tree_internal::Lazy<M>, dynamic_segment_tree_internal::Patricia<M>> {
 using Base= std::conditional_t<segment_tree_internal::action<M>::value, dynamic_segment_tree_internal::Lazy<M>, dynamic_segment_tree_internal::Patricia<M>>;
public:
 using Base::Base;
};
