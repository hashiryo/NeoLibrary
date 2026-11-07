#pragma once
#include "neo/data_structure/SegmentTree.hpp"
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
// 区間加算と区間和 monoid::AddSum<T>。SegmentTree<monoid::AddSum<T>> は、T が 64 bit 以下の整数なら、区間加算を差分の
// 配列 d (d[l] += x、d[r] -= x) に読み替えて、節を 8 個に広げた接頭辞和の木 (monoid::Sum の木の写し) を 2 本持つ形に
// なる。a[0..p) の和は、初めの値の接頭辞和 S0[p] と、b1 = Σ_{j<p} d[j] と b2 = Σ_{j<p} d[j] j を使って
// S0[p] + p b1 - b2 で出す (BIT を 2 本使う定番の形と同じ)。区間の和は 2 つの接頭辞の 4 本の和を同じループで読む。
// 区間の加算は 2 本の木に 2 点ずつ足す。値は u64 の桁あふれのまま計算し、差を T に戻すので、T の範囲に収まる和は
// 正しく出る。max_right と min_left は区間の和を使った二分探索。それ以外の T は二分木の汎用の形になる。
namespace monoid {
template <class V> struct AddSum {
 using T= V;
 using E= V;
 static constexpr T ti() { return T(); }
 static constexpr T op(T a, T b) { return a + b; }
 static constexpr void mp(T& v, E x, int len) { v+= x * T(len); }
 static constexpr void cp(E& pre, E suf) { pre+= suf; }
};
}
namespace segtree_add_sum_internal {
using u64= unsigned long long;
constexpr int LB= 3, B= 1 << LB;
// mask[B - 1 - s + j] は j > s のとき全部 1。s 番目より後ろのレーンを選ぶ窓。
struct Mask {
 alignas(64) u64 v[2 * B];
 constexpr Mask(): v{} {
  for(int j= 0; j < B; ++j) v[B + j]= ~0ull;
 }
};
inline constexpr Mask mask{};
// 位置 [0, n) への 1 点の加算と、[0, r) の和 (0 <= r <= n) を持つ 8 分木。最初はすべて 0。
class Prefix {
 int H= 0;
 std::size_t loff[24]= {};  // 段 h の表は mem の loff[h] から (写しても壊れないよう、ポインタでなくずれで持つ)
 segment_tree_internal::Vec<u64> mem;
 u64* L(int h) { return mem.data() + loff[h]; }
public:
 Prefix()= default;
 explicit Prefix(int n) {
  const std::size_t N= n;
  H= 1;
  while((N >> (H * LB)) != 0) ++H;
  std::size_t total= 0;
  for(int h= 0; h < H; ++h) loff[h]= total, total+= ((N >> ((h + 1) * LB)) + 1) << LB;
  mem.assign(total, 0);
 }
 int height() const { return H; }
 const u64* level(int h) const { return mem.data() + loff[h]; }
 [[gnu::always_inline]] void add(int p, u64 x) {
  const __m256i xv= _mm256_set1_epi64x((long long)x);
  for(int h= 0; h < H; ++h) {
   const std::size_t c= std::size_t(p) >> (h * LB);
   u64* blk= L(h) + ((c >> LB) << LB);
   const u64* m= mask.v + (B - 1 - (c & (B - 1)));
   const __m256i v0= _mm256_load_si256((const __m256i*)blk), v1= _mm256_load_si256((const __m256i*)(blk + 4));
   const __m256i m0= _mm256_loadu_si256((const __m256i*)m), m1= _mm256_loadu_si256((const __m256i*)(m + 4));
   _mm256_store_si256((__m256i*)blk, _mm256_add_epi64(v0, _mm256_and_si256(m0, xv)));
   _mm256_store_si256((__m256i*)(blk + 4), _mm256_add_epi64(v1, _mm256_and_si256(m1, xv)));
  }
 }
 [[gnu::always_inline]] u64 prefix(int r) const {
  u64 s= 0;
  for(int h= 0; h < H; ++h) s+= level(h)[std::size_t(r) >> (h * LB)];
  return s;
 }
};
template <class V> class Wide {
 int n= 0;
 segment_tree_internal::Vec<u64> s0;  // 初めの値の接頭辞和
 Prefix b1, b2;
 template <class F> void build(const F& init) {
  s0.assign(n + 1, 0);
  for(int i= 0; i < n; ++i) s0[i + 1]= s0[i] + u64(V(init(i)));
  b1= Prefix(n), b2= Prefix(n);
 }
 u64 prefix(int p) const { return s0[p] + u64(p) * b1.prefix(p) - b2.prefix(p); }
public:
 Wide()= default;
 explicit Wide(int n): n(n) {
  build([](int) { return V(); });
 }
 Wide(int n, V x): n(n) {
  build([&](int) { return x; });
 }
 Wide(const std::vector<V>& a): n(a.size()) {
  build([&](int i) { return a[i]; });
 }
 template <class F, std::enable_if_t<std::is_invocable_r_v<V, const F&, int>, int> = 0> Wide(int n, const F& init): n(n) { build(init); }
 int size() const { return n; }
 V get(int i) const { return prod(i, i + 1); }
 V operator[](int i) const { return get(i); }
 void set(int i, V x) { apply(i, i + 1, V(u64(x) - u64(get(i)))); }
 void mul(int i, V x) { apply(i, i + 1, x); }
 // [l, r) の和。2 つの接頭辞の b1 と b2 (4 本) を同じループで読む。
 [[gnu::always_inline]] V prod(int l, int r) const {
  const int H= b1.height();
  u64 x1= 0, y1= 0, x2= 0, y2= 0;
  for(int h= 0; h < H; ++h) {
   const std::size_t il= std::size_t(l) >> (h * LB), ir= std::size_t(r) >> (h * LB);
   const u64 *p= b1.level(h), *q= b2.level(h);
   x1+= p[il], y1+= p[ir], x2+= q[il], y2+= q[ir];
  }
  return V(s0[r] - s0[l] + u64(r) * y1 - u64(l) * x1 - (y2 - x2));
 }
 V all_prod() const { return V(prefix(n)); }
 [[gnu::always_inline]] void apply(int l, int r, V x) {
  if(l >= r) return;
  const u64 d= u64(x);
  b1.add(l, d), b2.add(l, d * u64(l));
  if(r < n) b1.add(r, -d), b2.add(r, -d * u64(r));
 }
 template <class C> int max_right(int l, const C& check) const {
  assert(0 <= l && l <= n && check(V()));
  int lo= l, hi= n + 1;
  while(hi - lo > 1) {
   const int mid= (lo + hi) >> 1;
   (check(prod(l, mid)) ? lo : hi)= mid;
  }
  return lo;
 }
 template <class C> int min_left(int r, const C& check) const {
  assert(0 <= r && r <= n && check(V()));
  int lo= -1, hi= r;
  while(hi - lo > 1) {
   const int mid= (lo + hi) >> 1;
   (check(prod(mid, r)) ? hi : lo)= mid;
  }
  return hi;
 }
};
template <class V> constexpr bool use_wide= std::is_integral_v<V> && !std::is_same_v<V, bool> && sizeof(V) <= 8;
template <class V> using Impl= std::conditional_t<use_wide<V>, Wide<V>, segment_tree_internal::Generic<monoid::AddSum<V>>>;
}
template <class V> class SegmentTree<monoid::AddSum<V>>: public segtree_add_sum_internal::Impl<V> {
 using Base= segtree_add_sum_internal::Impl<V>;
public:
 using Base::Base;
};
