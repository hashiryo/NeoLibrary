#pragma once
#include "neo/data_structure/SegmentTree.hpp"
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
// 和のモノイド monoid::Sum<T>。SegmentTree<monoid::Sum<T>> は、T が 64 bit 以下の整数なら、節を 8 個に広げた
// 接頭辞和の木になる (procon-judge の yosupo-point-add-range-sum の wide8)。段 h の表の f 番目には、f と同じ塊
// (8 個ずつ) の中で f より前にある子の和を持ち、[0, r) の和は段ごとに表の (r >> 3h) 番目を 1 つずつ読んで足す。
// 1 点の加算は、段ごとに塊の後ろ側へマスクを取った値を AVX2 で足す。塊の 8 個は u64 でキャッシュの 1 行に収まる。
// 値は u64 の桁あふれのまま足し、差を T に戻すので、T の範囲に収まる和は正しく出る。get のために元の値も持つ。
// max_right と min_left は、段ごとに塊の 8 個の境目の接頭辞和を見て、条件が切り替わる子へ下りる。
// それ以外の T (浮動小数点数など) は二分木の汎用の形になる。
namespace monoid {
template <class V> struct Sum {
 using T= V;
 static constexpr T ti() { return T(); }
 static constexpr T op(T a, T b) { return a + b; }
};
}
namespace segtree_sum_internal {
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
template <class V> class Wide {
 int n= 0, H= 0;
 u64* lv[24]= {};
 segment_tree_internal::Vec<u64> mem;
 segment_tree_internal::Vec<V> a;
 template <class F> void build(const F& init) {
  const std::size_t N= n;
  H= 1;
  while((N >> (H * LB)) != 0) ++H;
  std::size_t off[24], total= 0;
  for(int h= 0; h < H; ++h) off[h]= total, total+= ((N >> ((h + 1) * LB)) + 1) << LB;
  mem.assign(total, 0), a.resize(n);
  for(int h= 0; h < H; ++h) lv[h]= mem.data() + off[h];
  for(int i= 0; i < n; ++i) a[i]= init(i);
  // 段 0 は塊の中の左からの和、段 h + 1 の子の値は段 h の塊の合計。下から順に作る。
  std::vector<u64> cur(a.begin(), a.end()), nxt;
  for(int h= 0; h < H; ++h) {
   const std::size_t m= ((N >> ((h + 1) * LB)) + 1) << LB;
   nxt.assign(m >> LB, 0);
   for(std::size_t f= 0; f < m; ++f) {
    const std::size_t j= f & (B - 1);
    lv[h][f]= j ? lv[h][f - 1] + (f - 1 < cur.size() ? cur[f - 1] : 0) : 0;
    if(j == B - 1) nxt[f >> LB]= lv[h][f] + (f < cur.size() ? cur[f] : 0);
   }
   cur.swap(nxt);
  }
 }
 [[gnu::always_inline]] u64 prefix(int r) const {
  u64 s= 0;
  for(int h= 0; h < H; ++h) s+= lv[h][std::size_t(r) >> (h * LB)];
  return s;
 }
 [[gnu::always_inline]] void add(int p, u64 x) {
  const __m256i xv= _mm256_set1_epi64x((long long)x);
  for(int h= 0; h < H; ++h) {
   const std::size_t c= std::size_t(p) >> (h * LB);
   u64* blk= lv[h] + ((c >> LB) << LB);
   const u64* m= mask.v + (B - 1 - (c & (B - 1)));
   const __m256i v0= _mm256_load_si256((const __m256i*)blk), v1= _mm256_load_si256((const __m256i*)(blk + 4));
   const __m256i m0= _mm256_loadu_si256((const __m256i*)m), m1= _mm256_loadu_si256((const __m256i*)(m + 4));
   _mm256_store_si256((__m256i*)blk, _mm256_add_epi64(v0, _mm256_and_si256(m0, xv)));
   _mm256_store_si256((__m256i*)(blk + 4), _mm256_add_epi64(v1, _mm256_and_si256(m1, xv)));
  }
 }
public:
 Wide()= default;
 explicit Wide(int n): n(n) {
  build([](int) { return V(); });
 }
 Wide(int n, V x): n(n) {
  build([&](int) { return x; });
 }
 Wide(const std::vector<V>& v): n(v.size()) {
  build([&](int i) { return v[i]; });
 }
 template <class F, std::enable_if_t<std::is_invocable_r_v<V, const F&, int>, int> = 0> Wide(int n, const F& init): n(n) { build(init); }
 int size() const { return n; }
 V get(int i) const { return a[i]; }
 V operator[](int i) const { return a[i]; }
 [[gnu::always_inline]] void set(int i, V x) { add(i, u64(x) - u64(a[i])), a[i]= x; }
 [[gnu::always_inline]] void mul(int i, V x) { add(i, u64(x)), a[i]+= x; }
 // [l, r) の和。2 つの接頭辞を同じループで読む。
 [[gnu::always_inline]] V prod(int l, int r) const {
  u64 s= 0;
  for(int h= 0; h < H; ++h) s+= lv[h][std::size_t(r) >> (h * LB)] - lv[h][std::size_t(l) >> (h * LB)];
  return V(s);
 }
 V all_prod() const { return V(prefix(n)); }
 // check(prod(l, r)) が真になる最大の r。境目 x の真偽 g(x) を「x <= l なら真、x > n なら偽、ほかは
 // check(prefix(x) - prefix(l))」とすると、真の並びのあとに偽の並びが来るので、各段で塊の中の最後の真の子へ下りる。
 template <class C> int max_right(int l, const C& check) const {
  assert(0 <= l && l <= n && check(V()));
  const u64 base= prefix(l);
  u64 acc= 0;
  std::size_t b= 0;
  for(int h= H - 1;; --h) {
   const int sh= h * LB;
   const std::size_t start= b << (sh + LB);
   int jb= 0;
   for(int j= 1; j < B; ++j) {
    const std::size_t x= start + (std::size_t(j) << sh);
    if(x > std::size_t(n)) break;
    if(x > std::size_t(l) && !check(V(acc + lv[h][(b << LB) + j] - base))) break;
    jb= j;
   }
   if(h == 0) return int(start + jb);
   acc+= lv[h][(b << LB) + jb], b= (b << LB) + jb;
  }
 }
 // check(prod(l, r)) が真になる最小の l。g(x) を「x >= r なら真、ほかは check(prefix(r) - prefix(x))」とすると、
 // 偽の並びのあとに真の並びが来るので、各段で塊の中の最後の偽の子へ下りる。
 template <class C> int min_left(int r, const C& check) const {
  assert(0 <= r && r <= n && check(V()));
  const u64 top= prefix(r);
  if(check(V(top))) return 0;
  u64 acc= 0;
  std::size_t b= 0;
  for(int h= H - 1;; --h) {
   const int sh= h * LB;
   const std::size_t start= b << (sh + LB);
   int jf= 0;  // 塊の始まり (j = 0) はいつも偽
   for(int j= B - 1; j >= 1; --j) {
    const std::size_t x= start + (std::size_t(j) << sh);
    if(x >= std::size_t(r) || check(V(top - (acc + lv[h][(b << LB) + j])))) continue;
    jf= j;
    break;
   }
   if(h == 0) return int(start + jf + 1);
   acc+= lv[h][(b << LB) + jf], b= (b << LB) + jf;
  }
 }
};
template <class V> constexpr bool use_wide= std::is_integral_v<V> && !std::is_same_v<V, bool> && sizeof(V) <= 8;
}
template <class V> class SegmentTree<monoid::Sum<V>>: public std::conditional_t<segtree_sum_internal::use_wide<V>, segtree_sum_internal::Wide<V>, segment_tree_internal::Generic<monoid::Sum<V>>> {
 using Base= std::conditional_t<segtree_sum_internal::use_wide<V>, segtree_sum_internal::Wide<V>, segment_tree_internal::Generic<monoid::Sum<V>>>;
public:
 using Base::Base;
};
