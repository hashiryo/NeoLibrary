#pragma once
#include <limits>
#include "neo/base/ds/range/SegmentTree.hpp"
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
// 区間代入と区間和 monoid::AssignSum<T>、区間代入と区間最小 monoid::AssignMin<T>、区間代入と区間最大 monoid::AssignMax<T>。
// SegmentTree に渡すと、T が 64 bit 以下の符号付き整数なら、子を 16 個持つ遅延セグ木になる。段 0 は要素を 16 個ずつの
// 塊に並べ、段 h >= 1 の塊は子 16 個ぶんの値 s (和か最小か最大) と、まだ配っていない代入 t を並べて持つ (t が持つ型の
// 最小なら代入なし。代入する値に持つ型の最小は使えない)。最小と最大で T が 32 bit 以下なら i32 で持ち、16 個の塊が
// 64 byte に収まる。和と 64 bit の T は i64 で持つ。区間の代入は、根から両端の道を配り、下から全部掛かるレーンの
// s と t を AVX2 でまとめて書き、両端の道の s を子の 16 個から取り直す。代入は上の段ほど新しく、子の部分木全体を同じ
// 値にするので、区間の取得と 1 点の取得は、道の途中で代入の残った子に当たればその値 (和なら値と個数の積) を返して
// 下りずに済む。区間の取得は書き込みをせず、根から l と r - 1 が同じ子にある間は下り、分かれたら左右の道を下りる。
// 最後の塊だけはレーンごとの長さを使い、要素の無いレーンには代入しない。1 レーンの書き込みも 32 byte をベクトルで読んで混ぜて
// 書き、あとのベクトルの読み込みが書き込みから値を受け取れるようにする。それ以外の T は二分木の汎用の形になる。
namespace monoid {
template <class V> struct AssignSum {
 using T= V;
 using E= V;
 static constexpr T ti() { return T(); }
 static constexpr T op(T a, T b) { return a + b; }
 static constexpr void mp(T& v, E x, int len) { v= x * T(len); }
 static constexpr void cp(E& pre, E suf) { pre= suf; }
};
template <class V> struct AssignMin {
 using T= V;
 using E= V;
 static constexpr T ti() { return std::numeric_limits<V>::max(); }
 static constexpr T op(T a, T b) { return b < a ? b : a; }
 static constexpr void mp(T& v, E x) { v= x; }
 static constexpr void cp(E& pre, E suf) { pre= suf; }
};
template <class V> struct AssignMax {
 using T= V;
 using E= V;
 static constexpr T ti() { return std::numeric_limits<V>::lowest(); }
 static constexpr T op(T a, T b) { return a < b ? b : a; }
 static constexpr void mp(T& v, E x) { v= x; }
 static constexpr void cp(E& pre, E suf) { pre= suf; }
};
}
namespace segtree_assign_internal {
using i64= long long;
using i32= int;
enum Kind { SUM, MIN, MAX };
inline __m256i ld(const void* p) { return _mm256_load_si256((const __m256i*)p); }
inline void st(void* p, __m256i x) { _mm256_store_si256((__m256i*)p, x); }
// レーン [lo, hi) を選ぶマスク (64 bit の 16 レーンを 4 本)。gcc は __m256i の配列をスタックへ退避するので別々の変数で返す。
inline void masks64(int lo, int hi, __m256i& k0, __m256i& k1, __m256i& k2, __m256i& k3) {
 const __m256i L= _mm256_set1_epi64x(lo - 1), R= _mm256_set1_epi64x(hi);
 const __m256i i0= _mm256_setr_epi64x(0, 1, 2, 3), i1= _mm256_setr_epi64x(4, 5, 6, 7);
 const __m256i i2= _mm256_setr_epi64x(8, 9, 10, 11), i3= _mm256_setr_epi64x(12, 13, 14, 15);
 k0= _mm256_and_si256(_mm256_cmpgt_epi64(i0, L), _mm256_cmpgt_epi64(R, i0));
 k1= _mm256_and_si256(_mm256_cmpgt_epi64(i1, L), _mm256_cmpgt_epi64(R, i1));
 k2= _mm256_and_si256(_mm256_cmpgt_epi64(i2, L), _mm256_cmpgt_epi64(R, i2));
 k3= _mm256_and_si256(_mm256_cmpgt_epi64(i3, L), _mm256_cmpgt_epi64(R, i3));
}
// 32 bit の 16 レーンを 2 本。
inline void masks32(int lo, int hi, __m256i& k0, __m256i& k1) {
 const __m256i L= _mm256_set1_epi32(lo - 1), R= _mm256_set1_epi32(hi);
 const __m256i i0= _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7), i1= _mm256_setr_epi32(8, 9, 10, 11, 12, 13, 14, 15);
 k0= _mm256_and_si256(_mm256_cmpgt_epi32(i0, L), _mm256_cmpgt_epi32(R, i0));
 k1= _mm256_and_si256(_mm256_cmpgt_epi32(i1, L), _mm256_cmpgt_epi32(R, i1));
}
// 塊 (16 レーン) の操作を、持つ型 W (i64 か i32) で分ける。p は塊の先頭 (32 byte 境界)。
template <class W, Kind K> struct Lanes {
 static constexpr bool S= sizeof(W) == 4;
 static constexpr W PAD= K == SUM ? 0 : K == MIN ? std::numeric_limits<W>::max() : std::numeric_limits<W>::min();
 static W agg(W a, W b) {
  if constexpr(K == SUM) return a + b;
  else if constexpr(K == MIN) return b < a ? b : a;
  else return a < b ? b : a;
 }
 static __m256i vagg(__m256i a, __m256i b) {
  if constexpr(S) return K == SUM ? _mm256_add_epi32(a, b) : K == MIN ? _mm256_min_epi32(a, b) : _mm256_max_epi32(a, b);
  else if constexpr(K == SUM) return _mm256_add_epi64(a, b);
  else if constexpr(K == MIN) return _mm256_blendv_epi8(a, b, _mm256_cmpgt_epi64(a, b));
  else return _mm256_blendv_epi8(a, b, _mm256_cmpgt_epi64(b, a));
 }
 static W hagg(__m256i m) {
  if constexpr(S) {
   __m128i c= _mm256_castsi256_si128(m);
   const __m128i h= _mm256_extracti128_si256(m, 1);
   auto f= [](__m128i x, __m128i y) { return K == SUM ? _mm_add_epi32(x, y) : K == MIN ? _mm_min_epi32(x, y) : _mm_max_epi32(x, y); };
   c= f(c, h), c= f(c, _mm_shuffle_epi32(c, 0x4e)), c= f(c, _mm_shuffle_epi32(c, 0xb1));
   return _mm_cvtsi128_si32(c);
  } else {
   W v[4];
   _mm256_storeu_si256((__m256i*)v, m);
   return agg(agg(v[0], v[1]), agg(v[2], v[3]));
  }
 }
 static W block(const W* p) {
  if constexpr(S) return hagg(vagg(ld(p), ld(p + 8)));
  else return hagg(vagg(vagg(ld(p), ld(p + 4)), vagg(ld(p + 8), ld(p + 12))));
 }
 // レーン [lo, hi) の値 (塊の外のレーンは PAD とみなす)。
 static W lanes(const W* p, int lo, int hi) {
  if constexpr(S) {
   __m256i k0, k1;
   masks32(lo, hi, k0, k1);
   const __m256i pad= _mm256_set1_epi32(PAD);
   return hagg(vagg(_mm256_blendv_epi8(pad, ld(p), k0), _mm256_blendv_epi8(pad, ld(p + 8), k1)));
  } else {
   __m256i k0, k1, k2, k3;
   masks64(lo, hi, k0, k1, k2, k3);
   const __m256i pad= _mm256_set1_epi64x(PAD);
   const __m256i a= vagg(_mm256_blendv_epi8(pad, ld(p), k0), _mm256_blendv_epi8(pad, ld(p + 4), k1));
   return hagg(vagg(a, vagg(_mm256_blendv_epi8(pad, ld(p + 8), k2), _mm256_blendv_epi8(pad, ld(p + 12), k3))));
  }
 }
 // p のレーン [lo, hi) を v に、q (null でなければ) のレーン [lo, hi) を w にする。マスクは 1 回だけ作る。
 static void fill2(W* p, W v, W* q, W w, int lo, int hi) {
  if constexpr(S) {
   __m256i k0, k1;
   masks32(lo, hi, k0, k1);
   const __m256i x= _mm256_set1_epi32(v), y= _mm256_set1_epi32(w);
   st(p, _mm256_blendv_epi8(ld(p), x, k0)), st(p + 8, _mm256_blendv_epi8(ld(p + 8), x, k1));
   if(q) st(q, _mm256_blendv_epi8(ld(q), y, k0)), st(q + 8, _mm256_blendv_epi8(ld(q + 8), y, k1));
  } else {
   __m256i k0, k1, k2, k3;
   masks64(lo, hi, k0, k1, k2, k3);
   const __m256i x= _mm256_set1_epi64x(v), y= _mm256_set1_epi64x(w);
   st(p, _mm256_blendv_epi8(ld(p), x, k0)), st(p + 4, _mm256_blendv_epi8(ld(p + 4), x, k1));
   st(p + 8, _mm256_blendv_epi8(ld(p + 8), x, k2)), st(p + 12, _mm256_blendv_epi8(ld(p + 12), x, k3));
   if(q) {
    st(q, _mm256_blendv_epi8(ld(q), y, k0)), st(q + 4, _mm256_blendv_epi8(ld(q + 4), y, k1));
    st(q + 8, _mm256_blendv_epi8(ld(q + 8), y, k2)), st(q + 12, _mm256_blendv_epi8(ld(q + 12), y, k3));
   }
  }
 }
 // j 番目のレーンを v にする。スカラーで書くと、あとで同じ塊をベクトルで読むときに書き込みから値を受け取れず、
 // 書き込みが終わるまで待つ。ベクトルで読んで混ぜて書けば受け取れる。
 static void put(W* p, int j, W v) {
  if constexpr(S) {
   __m256i* q= (__m256i*)(p + (j & ~7));
   const __m256i k= _mm256_cmpeq_epi32(_mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7), _mm256_set1_epi32(j & 7));
   _mm256_store_si256(q, _mm256_blendv_epi8(_mm256_load_si256(q), _mm256_set1_epi32(v), k));
  } else {
   __m256i* q= (__m256i*)(p + (j & ~3));
   const __m256i k= _mm256_cmpeq_epi64(_mm256_setr_epi64x(0, 1, 2, 3), _mm256_set1_epi64x(j & 3));
   _mm256_store_si256(q, _mm256_blendv_epi8(_mm256_load_si256(q), _mm256_set1_epi64x(v), k));
  }
 }
};
template <class W> struct alignas(64) Blk {
 W s[16], t[16];
};
template <class V, Kind K> class Wide {
 // 塊に置く型。和は 64 bit が要るので、32 bit にするのは最小と最大で V が 32 bit 以下のときだけ。
 using W= std::conditional_t<K != SUM && sizeof(V) <= 4, i32, i64>;
 using Ln= Lanes<W, K>;
 using B_= Blk<W>;
 static constexpr W NONE= std::numeric_limits<W>::min();
 static constexpr W PAD= Ln::PAD;
 static constexpr V TI= K == SUM ? V() : K == MIN ? std::numeric_limits<V>::max() : std::numeric_limits<V>::lowest();
 int n= 0, H= 0;
 std::size_t loff[8]= {};  // 段 h >= 1 の塊は mem1 の loff[h] から。写しても壊れないよう、ポインタでなくずれで持つ
 int last[8]= {};
 i64 plen[8]= {};                // 段 h の子の長さ 16^h
 alignas(32) i64 tl[8][16]= {};  // 段 h の最後の塊のレーンごとの長さ (要素の無いレーンは 0)
 segment_tree_internal::Vec<W> mem0;
 segment_tree_internal::Vec<B_> mem1;
 W* V0() { return mem0.data(); }
 const W* V0() const { return mem0.data(); }
 B_* LV(int h) { return mem1.data() + loff[h]; }
 const B_* LV(int h) const { return mem1.data() + loff[h]; }
 static W agg(W a, W b) { return Ln::agg(a, b); }
 // 値 x を len 個の区間に代入したときの値。
 static W val(W x, i64 len) {
  if constexpr(K == SUM) return x * len;
  else return x;
 }
 const W* row(int h, int b) const { return h == 0 ? V0() + (std::size_t(b) << 4) : LV(h)[b].s; }
 W block_agg(int h, int b) const { return Ln::block(row(h, b)); }
 W lanes_agg(int h, int b, int lo, int hi) const { return Ln::lanes(row(h, b), lo, hi); }
 W tag(int h, int c) const { return LV(h)[c >> 4].t[c & 15]; }
 // 段 h の塊 b のレーン [lo, hi) に x を代入する。最後の塊ではレーンごとの長さを使い、要素の無いレーンは書かない。
 void assign_lanes(int h, int b, int lo, int hi, W x) {
  if(b == last[h]) hi= std::min(hi, int(std::count_if(tl[h], tl[h] + 16, [](i64 v) { return v > 0; })));
  if(lo >= hi) return;
  if(h == 0) return Ln::fill2(V0() + (std::size_t(b) << 4), x, nullptr, 0, lo, hi);
  B_& B= LV(h)[b];
  if constexpr(K == SUM) {
   if(b == last[h]) {
    // レーンごとの長さ (2^32 未満) との積を、x の下と上の 32 bit に分けた 2 回の 32 x 32 の積で作る。
    __m256i k0, k1, k2, k3;
    masks64(lo, hi, k0, k1, k2, k3);
    const __m256i xv= _mm256_set1_epi64x(x), xh= _mm256_srli_epi64(xv, 32);
    auto mul= [&](__m256i len) { return _mm256_add_epi64(_mm256_mul_epu32(xv, len), _mm256_slli_epi64(_mm256_mul_epu32(xh, len), 32)); };
    const i64* L= tl[h];
    st(B.s, _mm256_blendv_epi8(ld(B.s), mul(ld(L)), k0)), st(B.s + 4, _mm256_blendv_epi8(ld(B.s + 4), mul(ld(L + 4)), k1));
    st(B.s + 8, _mm256_blendv_epi8(ld(B.s + 8), mul(ld(L + 8)), k2)), st(B.s + 12, _mm256_blendv_epi8(ld(B.s + 12), mul(ld(L + 12)), k3));
    return Ln::fill2(B.t, x, nullptr, 0, lo, hi);
   }
  }
  Ln::fill2(B.s, val(x, plen[h]), B.t, x, lo, hi);
 }
 void push(int h, int c) {
  B_& B= LV(h)[c >> 4];
  const W x= B.t[c & 15];
  if(x == NONE) return;
  Ln::put(B.t, c & 15, NONE);
  assign_lanes(h - 1, c, 0, 16, x);
 }
 // 段 h の子 (l を含む) から下の、[l, 子の終わり) の値。
 W left_path(int h, int l) const {
  W acc= PAD;
  for(int g= h;; --g) {
   const int c= l >> (4 * g);
   if(const W t= tag(g, c); t != NONE) return agg(acc, val(t, (i64(c + 1) << (4 * g)) - l));
   const int sh= 4 * (g - 1), b= c, j= (l >> sh) & 15;
   if(g == 1 || !(l & ((1 << sh) - 1))) return agg(acc, lanes_agg(g - 1, b, j, 16));
   if(j < 15) acc= agg(acc, lanes_agg(g - 1, b, j + 1, 16));
  }
 }
 // 段 h の子 (r - 1 を含む) から下の、[子の始まり, r) の値。
 W right_path(int h, int r) const {
  const int x= r - 1;
  W acc= PAD;
  for(int g= h;; --g) {
   const int c= x >> (4 * g);
   if(const W t= tag(g, c); t != NONE) return agg(acc, val(t, r - (i64(c) << (4 * g))));
   const int sh= 4 * (g - 1), b= c, j= (x >> sh) & 15;
   if(g == 1 || !(r & ((1 << sh) - 1))) return agg(acc, lanes_agg(g - 1, b, 0, j + 1));
   if(j > 0) acc= agg(acc, lanes_agg(g - 1, b, 0, j));
  }
 }
 W range_agg(int l, int r) const {
  const int x= r - 1;
  int h= H - 1, b= 0;
  for(;;) {
   const int sh= 4 * h, jl= (l >> sh) & 15, jr= (x >> sh) & 15;
   if(h == 0) return lanes_agg(0, b, jl, jr + 1);
   const int sm= (1 << sh) - 1;
   if(jl != jr) {
    W res= PAD;
    int fl= jl, fr= jr + 1;
    if(l & sm) res= agg(res, left_path(h, l)), fl= jl + 1;
    if(r & sm) res= agg(res, right_path(h, r)), fr= jr;
    if(fl < fr) res= agg(res, lanes_agg(h, b, fl, fr));
    return res;
   }
   const B_& B= LV(h)[b];
   if(B.t[jl] != NONE) return val(B.t[jl], r - l);
   if(!(l & sm) && !(r & sm)) return B.s[jl];
   b= (b << 4) + jl, --h;
  }
 }
 void range_assign(int l, int r, W x) {
  const int y= r - 1;
  for(int h= H - 1; h >= 1; --h) {
   const int sh= 4 * h, sm= (1 << sh) - 1;
   if(l & sm) push(h, l >> sh);
   if(r & sm) push(h, y >> sh);
  }
  for(int h= 0, p= l, q= r; p < q; ++h) {
   const int bp= p >> 4, bq= q >> 4;
   if(bp == bq || h == H - 1) {
    assign_lanes(h, bp, p & 15, q - (bp << 4), x);
    break;
   }
   if(p & 15) assign_lanes(h, bp, p & 15, 16, x), p= bp + 1;
   else p= bp;
   if(q & 15) assign_lanes(h, bq, 0, q & 15, x), q= bq;
   else q= bq;
  }
  for(int h= 1; h < H; ++h) {
   const int sh= 4 * h, sm= (1 << sh) - 1;
   if(l & sm) Ln::put(LV(h)[l >> (sh + 4)].s, (l >> sh) & 15, block_agg(h - 1, l >> sh));
   if(r & sm) Ln::put(LV(h)[y >> (sh + 4)].s, (y >> sh) & 15, block_agg(h - 1, y >> sh));
  }
 }
 template <class F> void build(const F& init) {
  int cnt[8], nb= std::max((n + 15) >> 4, 1);
  H= 0;
  std::size_t total1= 0;
  for(;;) {
   cnt[H++]= nb;
   if(nb == 1) break;
   nb= (nb + 15) >> 4;
   total1+= nb;
  }
  mem0.assign(std::size_t(cnt[0]) << 4, PAD);
  mem1.resize(std::max<std::size_t>(total1, 1));
  std::size_t o= 0;
  for(int h= 1; h < H; ++h) loff[h]= o, o+= cnt[h];
  for(int h= 0; h < H; ++h) {
   last[h]= cnt[h] - 1, plen[h]= i64(1) << (4 * h);
   for(int j= 0; j < 16; ++j) {
    const i64 start= (i64(last[h]) * 16 + j) << (4 * h);
    tl[h][j]= std::clamp<i64>(n - start, 0, plen[h]);
   }
  }
  for(int i= 0; i < n; ++i) V0()[i]= W(init(i));
  for(int h= 1; h < H; ++h)
   for(int b= 0; b < cnt[h]; ++b)
    for(int j= 0; j < 16; ++j) {
     const int c= b * 16 + j;
     LV(h)[b].s[j]= c < cnt[h - 1] ? block_agg(h - 1, c) : PAD, LV(h)[b].t[j]= NONE;
    }
 }
public:
 Wide()= default;
 explicit Wide(int n): n(n) {
  build([](int) { return TI; });
 }
 Wide(int n, V x): n(n) {
  build([&](int) { return x; });
 }
 Wide(const std::vector<V>& a): n(a.size()) {
  build([&](int i) { return a[i]; });
 }
 template <class F, std::enable_if_t<std::is_invocable_r_v<V, const F&, int>, int> = 0> Wide(int n, const F& init): n(n) { build(init); }
 int size() const { return n; }
 // 代入は上の段ほど新しいので、根から見て最初に当たった代入がその要素の値。
 V get(int i) const {
  for(int h= H - 1; h >= 1; --h)
   if(const W t= tag(h, i >> (4 * h)); t != NONE) return V(t);
  return V(V0()[i]);
 }
 V operator[](int i) const { return get(i); }
 void set(int i, V x) { range_assign(i, i + 1, W(x)); }
 void mul(int i, V x) {
  if constexpr(K == SUM) set(i, V(get(i) + x));
  else set(i, V(agg(W(get(i)), W(x))));
 }
 [[gnu::always_inline]] V prod(int l, int r) const { return l < r ? V(range_agg(l, r)) : TI; }
 V all_prod() const { return n ? V(block_agg(H - 1, 0)) : TI; }
 [[gnu::always_inline]] void apply(int l, int r, V x) {
  if(l < r) range_assign(l, r, W(x));
 }
 template <class C> int max_right(int l, const C& check) const {
  assert(0 <= l && l <= n && check(TI));
  int lo= l, hi= n + 1;
  while(hi - lo > 1) {
   const int mid= (lo + hi) >> 1;
   (check(prod(l, mid)) ? lo : hi)= mid;
  }
  return lo;
 }
 template <class C> int min_left(int r, const C& check) const {
  assert(0 <= r && r <= n && check(TI));
  int lo= -1, hi= r;
  while(hi - lo > 1) {
   const int mid= (lo + hi) >> 1;
   (check(prod(mid, r)) ? hi : lo)= mid;
  }
  return hi;
 }
};
template <class V> constexpr bool use_wide= std::is_integral_v<V> && std::is_signed_v<V> && sizeof(V) <= 8;
template <class V, Kind K, class M> using Impl= std::conditional_t<use_wide<V>, Wide<V, K>, segment_tree_internal::Generic<M>>;
}
template <class V> class SegmentTree<monoid::AssignSum<V>>: public segtree_assign_internal::Impl<V, segtree_assign_internal::SUM, monoid::AssignSum<V>> {
 using Base= segtree_assign_internal::Impl<V, segtree_assign_internal::SUM, monoid::AssignSum<V>>;
public:
 using Base::Base;
};
template <class V> class SegmentTree<monoid::AssignMin<V>>: public segtree_assign_internal::Impl<V, segtree_assign_internal::MIN, monoid::AssignMin<V>> {
 using Base= segtree_assign_internal::Impl<V, segtree_assign_internal::MIN, monoid::AssignMin<V>>;
public:
 using Base::Base;
};
template <class V> class SegmentTree<monoid::AssignMax<V>>: public segtree_assign_internal::Impl<V, segtree_assign_internal::MAX, monoid::AssignMax<V>> {
 using Base= segtree_assign_internal::Impl<V, segtree_assign_internal::MAX, monoid::AssignMax<V>>;
public:
 using Base::Base;
};
