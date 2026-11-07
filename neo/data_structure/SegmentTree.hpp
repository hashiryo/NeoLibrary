#pragma once
#include <algorithm>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <type_traits>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// セグメント木。構造体 M でモノイドと作用を渡す (今の Library の SegmentTree と同じ決まり)。
//   モノイド: using T、static T ti() (単位元)、static T op(T a, T b)
//   作用:     using E、static void mp(T& v, E x) か mp(T& v, E x, int len) (len は区間の長さ)、
//             static void cp(E& pre, E suf) (pre を先に、suf をあとに当てる合成を pre に書く)
//   mp が bool を返す形は Beats で、偽なら子へ下りて作り直す (葉では必ず真を返すこと)。
// モノイドだけなら点の更新と区間の積、両方あれば遅延セグ木、作用だけなら双対セグ木になる。
// neo/data_structure/segtree/ の標準の構造体 (monoid::Sum など) を渡すと、そのヘッダが中身を速い形に替える。
// ここにあるのは、それ以外の構造体を受ける二分木の汎用の形。大きさは 2 の冪に揃え、ACL と同じ非再帰で書く。
// 区間の積は、値が算術型なら節が決まった段ですぐ畳み、それ以外は覆う節の番号を分岐なしで集めてから左と右の 2 本の
// 鎖で畳む。作用が残っているかは別の配列の印で持つので、作用の型に恒等元は要らない。
namespace segment_tree_internal {
template <class M, class= void> struct monoid: std::false_type {};
template <class M> struct monoid<M, std::void_t<typename M::T, decltype(&M::op), decltype(&M::ti)>>: std::true_type {};
template <class M, class= void> struct action: std::false_type {};
template <class M> struct action<M, std::void_t<typename M::T, typename M::E, decltype(&M::mp), decltype(&M::cp)>>: std::true_type {};
template <class M, class= void> struct action_type {
 using type= char;
};
template <class M> struct action_type<M, std::void_t<typename M::E>> {
 using type= typename M::E;
};
// 1 MB 以上の確保は 2 MB 境界に寄せ、Linux では huge page を頼む。
template <class T> struct Alloc {
 using value_type= T;
 Alloc()= default;
 template <class U> Alloc(const Alloc<U>&) {}
 T* allocate(std::size_t n) {
  const std::size_t want= n * sizeof(T);
  if(want < (std::size_t(1) << 20)) return static_cast<T*>(::operator new(want, std::align_val_t(64)));
  const std::size_t H= std::size_t(1) << 21, bytes= (want + H - 1) & ~(H - 1);
  T* p= static_cast<T*>(::operator new(bytes, std::align_val_t(H)));
#ifdef __linux__
  madvise(p, bytes, MADV_HUGEPAGE);
#endif
  return p;
 }
 void deallocate(T* p, std::size_t n) {
  const std::size_t want= n * sizeof(T);
  if(want < (std::size_t(1) << 20)) ::operator delete(p, std::align_val_t(64));
  else ::operator delete(p, std::align_val_t(std::size_t(1) << 21));
 }
 template <class U> bool operator==(const Alloc<U>&) const { return true; }
};
template <class T> using Vec= std::vector<T, Alloc<T>>;
template <class M> class Generic {
 using T= typename M::T;
 using E= typename action_type<M>::type;
 static constexpr bool MO= monoid<M>::value, AC= action<M>::value;
 static_assert(MO || AC, "M にはモノイド (T, ti, op) か作用 (T, E, mp, cp) が要る");
 int n= 0, lg= 0, sz= 1;
 Vec<T> d;  // モノイドがあれば節 1..2sz-1 の値、作用だけなら葉の値 (sz 個)
 Vec<E> lz;
 std::vector<unsigned char> fl;
 // mp の 4 つの形をそろえる。偽を返したら作用が節に収まらなかった (Beats)。
 static bool map(T& v, const E& x, int len) {
  if constexpr(std::is_invocable_r_v<bool, decltype(M::mp), T&, const E&, int>) return M::mp(v, x, len);
  else if constexpr(std::is_invocable_r_v<bool, decltype(M::mp), T&, const E&>) return M::mp(v, x);
  else if constexpr(std::is_invocable_v<decltype(M::mp), T&, const E&, int>) return M::mp(v, x, len), true;
  else return M::mp(v, x), true;
 }
 // 子の添字は 2k と 2k | 1 と書く。2k + 1 と書くと gcc が 2 つの子を 16 byte の 1 回の読み込みにまとめ、直前の段で
 // 8 byte で書いた子から値を受け取れずに、書き込みが終わるまで段ごとに待つ (procon-judge の self-segtree-point-set-range-min)。
 void update(int k) { d[k]= M::op(d[k << 1], d[k << 1 | 1]); }
 // 節 k (長さ len) に作用 x を当てる。
 void all_apply(int k, const E& x, int len) {
  if(k < sz) {
   if(fl[k]) M::cp(lz[k], x);
   else lz[k]= x, fl[k]= 1;
  }
  if constexpr(MO) {
   if(!map(d[k], x, len)) push(k, len), update(k);
  } else if(k >= sz) map(d[k - sz], x, 1);
 }
 void push(int k, int len) {
  if(fl[k]) all_apply(k << 1, lz[k], len >> 1), all_apply(k << 1 | 1, lz[k], len >> 1), fl[k]= 0;
 }
 // 葉 l と r - 1 (どちらも sz を足した番号) の祖先のうち、区間の途中で切れるものを根から配る。
 void push_ends(int l, int r) {
  for(int i= lg; i >= 1; --i) {
   if(((l >> i) << i) != l) push(l >> i, 1 << i);
   if(((r >> i) << i) != r) push((r - 1) >> i, 1 << i);
  }
 }
 void push_path(int p) {
  for(int i= lg; i >= 1; --i) push(p >> i, 1 << i);
 }
 template <class F> void build(const F& init) {
  sz= std::bit_ceil(unsigned(std::max(n, 1))), lg= std::countr_zero(unsigned(sz));
  if constexpr(MO) {
   d.assign(2 * sz, M::ti());
   for(int i= 0; i < n; ++i) d[sz + i]= init(i);
   for(int k= sz - 1; k >= 1; --k) update(k);
  } else {
   d.assign(sz, T());
   for(int i= 0; i < n; ++i) d[i]= init(i);
  }
  if constexpr(AC) lz.assign(sz, E()), fl.assign(sz, 0);
 }
public:
 Generic()= default;
 explicit Generic(int n): n(n) {
  if constexpr(MO) build([](int) { return M::ti(); });
  else build([](int) { return T(); });
 }
 Generic(int n, T x): n(n) {
  build([&](int) { return x; });
 }
 Generic(const std::vector<T>& a): n(a.size()) {
  build([&](int i) { return a[i]; });
 }
 template <class F, std::enable_if_t<std::is_invocable_r_v<T, const F&, int>, int> = 0> Generic(int n, const F& init): n(n) { build(init); }
 int size() const { return n; }
 T get(int i) {
  if constexpr(AC) {
   // 根から葉への道では上の節ほど新しい作用なので、葉の値に下の節から順に当てる (書き込みなし)。
   T v= MO ? d[sz + i] : d[i];
   for(int k= (sz + i) >> 1; k >= 1; k>>= 1)
    if(fl[k]) map(v, lz[k], 1);
   return v;
  } else return d[sz + i];
 }
 T operator[](int i) { return get(i); }
 void set(int i, T x) {
  const int p= sz + i;
  if constexpr(AC) push_path(p);
  if constexpr(MO) {
   d[p]= x;
   for(int k= p >> 1; k >= 1; k>>= 1) update(k);
  } else d[i]= x;
 }
 void mul(int i, T x) {
  static_assert(MO, "mul はモノイドがあるときだけ");
  set(i, M::op(get(i), x));
 }
 T prod(int l, int r) {
  static_assert(MO, "prod はモノイドがあるときだけ");
  l+= sz, r+= sz;
  if constexpr(AC) push_ends(l, r);
  if constexpr(std::is_arithmetic_v<T>) {
   // 値が算術型なら演算は軽いとみなし、節が決まった段ですぐ値を読む 1 本のループにする。分岐なしで集めてから
   // 読む形は、値の読み込みが集め終わるまで始まらず、x64 の最小で今の Library より遅かった (procon-judge の
   // self-segtree-point-set-range-min)。
   T sl= M::ti(), sr= M::ti();
   for(; l < r; l>>= 1, r>>= 1) {
    if(l & 1) sl= M::op(sl, d[l++]);
    if(r & 1) sr= M::op(d[--r], sr);
   }
   return M::op(sl, sr);
  }
  // 覆う節を分岐なしで集める。左は見つけた順 (左から右)、右は見つけた順 (右から左)。重い演算 (1 次関数の合成や
  // 行列の積) では、分岐の予測の外れが無くなるぶん、こちらが速い (procon-judge の yosupo-point-set-range-composite)。
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
 T all_prod() const {
  static_assert(MO, "all_prod はモノイドがあるときだけ");
  return d[1];
 }
 void apply(int l, int r, E x) {
  static_assert(AC, "apply は作用があるときだけ");
  if(l >= r) return;
  l+= sz, r+= sz;
  push_ends(l, r);
  for(int a= l, b= r, len= 1; a < b; a>>= 1, b>>= 1, len<<= 1) {
   if(a & 1) all_apply(a++, x, len);
   if(b & 1) all_apply(--b, x, len);
  }
  if constexpr(MO)
   for(int i= 1; i <= lg; ++i) {
    if(((l >> i) << i) != l) update(l >> i);
    if(((r >> i) << i) != r) update((r - 1) >> i);
   }
 }
 // check(prod(l, r)) が真になる最大の r (ACL と同じ。check(ti()) は真、r について単調)。
 template <class C> int max_right(int l, const C& check) {
  static_assert(MO, "max_right はモノイドがあるときだけ");
  assert(0 <= l && l <= n && check(M::ti()));
  if(l == n) return n;
  l+= sz;
  if constexpr(AC) push_path(l);
  T s= M::ti();
  int len= 1;
  do {
   while(!(l & 1)) l>>= 1, len<<= 1;
   if(!check(M::op(s, d[l]))) {
    while(l < sz) {
     if constexpr(AC) push(l, len);
     l<<= 1, len>>= 1;
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
  static_assert(MO, "min_left はモノイドがあるときだけ");
  assert(0 <= r && r <= n && check(M::ti()));
  if(r == 0) return 0;
  r+= sz;
  if constexpr(AC) push_path(r - 1);
  T s= M::ti();
  int len= 1;
  do {
   --r;
   while(r > 1 && (r & 1)) r>>= 1, len<<= 1;
   if(!check(M::op(d[r], s))) {
    while(r < sz) {
     if constexpr(AC) push(r, len);
     r= 2 * r + 1, len>>= 1;
     if(T t= M::op(d[r], s); check(t)) s= t, --r;
    }
    return r + 1 - sz;
   }
   s= M::op(d[r], s);
  } while((r & -r) != r);
  return 0;
 }
};
}  // namespace segment_tree_internal
// 標準の構造体のヘッダ (neo/data_structure/segtree/) は、この型を構造体ごとに特殊化して中身を替える。
template <class M> class SegmentTree: public segment_tree_internal::Generic<M> {
public:
 using segment_tree_internal::Generic<M>::Generic;
};
