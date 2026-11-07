#pragma once
#include <limits>
#include "neo/data_structure/SegmentTree.hpp"
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
// 最小のモノイド monoid::Min<T> と最大のモノイド monoid::Max<T>。SegmentTree<monoid::Min<T>> は、T が 64 bit 以下の
// 符号付き整数なら、子を 16 個持つ木になる。段 0 は要素を 16 個ずつの塊に並べ、段 h >= 1 の塊は子 16 個ぶんの最小を
// 並べて持つ (monoid::AddMin の木から加算の量を除いた形)。1 点の代入は葉を書き換え、段ごとに子の塊の 16 個の最小を
// AVX2 で取り直し、親の値が変わらなくなった段で止める。1 レーンの書き込みも、そのレーンを含む 32 byte をベクトルで読んで混ぜて
// 書く (あとのベクトルの読み込みが書き込みから値を受け取れるように)。mul (a[i] = min(a[i], x)) は値が小さくなる向きにしか動かない
// ので、子を読み直さずに上へ min を取るだけで済ます。区間の最小は、根から l と r - 1 が同じ子にある間は下り、
// 分かれたら左右の道を下りて、全部掛かるレーンの最小を比べる。値は i64 で、16 個を 4 本の __m256i で扱う (gcc が
// __m256i の配列をスタックへ退避するので、マスクは別々の変数で作る)。要素の無いレーンは、最小なら i64 の最大、
// 最大なら i64 の最小。それ以外の T は二分木の汎用の形になる。
namespace monoid {
template <class V> struct Min {
 using T= V;
 static constexpr T ti() { return std::numeric_limits<V>::max(); }
 static constexpr T op(T a, T b) { return b < a ? b : a; }
};
template <class V> struct Max {
 using T= V;
 static constexpr T ti() { return std::numeric_limits<V>::lowest(); }
 static constexpr T op(T a, T b) { return a < b ? b : a; }
};
}
namespace segtree_min_internal {
using i64= long long;
inline __m256i ld(const i64* p) { return _mm256_load_si256((const __m256i*)p); }
// p は 16 レーンの塊の先頭 (32 byte 境界)。j 番目のレーンを v にする。スカラーで書くと、あとで同じ塊をベクトルで読む
// ときに書き込みから値を受け取れず、書き込みが終わるまで待つ。ベクトルで読んで混ぜて書けば受け取れる。
inline void put_lane(i64* p, int j, i64 v) {
 __m256i* q= (__m256i*)(p + (j & ~3));
 const __m256i k= _mm256_cmpeq_epi64(_mm256_setr_epi64x(0, 1, 2, 3), _mm256_set1_epi64x(j & 3));
 _mm256_store_si256(q, _mm256_blendv_epi8(_mm256_load_si256(q), _mm256_set1_epi64x(v), k));
}
// レーン [lo, hi) を選ぶマスク (16 レーンを 4 本)。
inline void masks(int lo, int hi, __m256i& k0, __m256i& k1, __m256i& k2, __m256i& k3) {
 const __m256i L= _mm256_set1_epi64x(lo - 1), R= _mm256_set1_epi64x(hi);
 const __m256i i0= _mm256_setr_epi64x(0, 1, 2, 3), i1= _mm256_setr_epi64x(4, 5, 6, 7);
 const __m256i i2= _mm256_setr_epi64x(8, 9, 10, 11), i3= _mm256_setr_epi64x(12, 13, 14, 15);
 k0= _mm256_and_si256(_mm256_cmpgt_epi64(i0, L), _mm256_cmpgt_epi64(R, i0));
 k1= _mm256_and_si256(_mm256_cmpgt_epi64(i1, L), _mm256_cmpgt_epi64(R, i1));
 k2= _mm256_and_si256(_mm256_cmpgt_epi64(i2, L), _mm256_cmpgt_epi64(R, i2));
 k3= _mm256_and_si256(_mm256_cmpgt_epi64(i3, L), _mm256_cmpgt_epi64(R, i3));
}
template <class V, bool MAX> class Wide {
 static constexpr i64 PAD= MAX ? std::numeric_limits<i64>::min() : std::numeric_limits<i64>::max();
 static constexpr V TI= MAX ? std::numeric_limits<V>::lowest() : std::numeric_limits<V>::max();
 int n= 0, H= 0;
 std::size_t loff[8]= {};  // 段 h (0 は要素、h >= 1 は子の最小) は mem の loff[h] から。写しても壊れないようずれで持つ
 segment_tree_internal::Vec<i64> mem;
 i64* lv(int h) { return mem.data() + loff[h]; }
 const i64* lv(int h) const { return mem.data() + loff[h]; }
 static bool better(i64 a, i64 b) {
  if constexpr(MAX) return a > b;
  else return a < b;
 }
 static i64 best(i64 a, i64 b) { return better(b, a) ? b : a; }
 static __m256i vbest(__m256i a, __m256i b) {
  if constexpr(MAX) return _mm256_blendv_epi8(a, b, _mm256_cmpgt_epi64(b, a));
  else return _mm256_blendv_epi8(a, b, _mm256_cmpgt_epi64(a, b));
 }
 static i64 hbest(__m256i x0, __m256i x1, __m256i x2, __m256i x3) {
  const __m256i m= vbest(vbest(x0, x1), vbest(x2, x3));
  __m128i c= _mm256_castsi256_si128(m);
  const __m128i h= _mm256_extracti128_si256(m, 1);
  c= MAX ? _mm_blendv_epi8(c, h, _mm_cmpgt_epi64(h, c)) : _mm_blendv_epi8(c, h, _mm_cmpgt_epi64(c, h));
  const __m128i s= _mm_unpackhi_epi64(c, c);
  c= MAX ? _mm_blendv_epi8(c, s, _mm_cmpgt_epi64(s, c)) : _mm_blendv_epi8(c, s, _mm_cmpgt_epi64(c, s));
  return _mm_cvtsi128_si64(c);
 }
 i64 block_best(int h, int b) const {
  const i64* p= lv(h) + (std::size_t(b) << 4);
  return hbest(ld(p), ld(p + 4), ld(p + 8), ld(p + 12));
 }
 i64 lanes_best(int h, int b, int lo, int hi) const {
  const i64* p= lv(h) + (std::size_t(b) << 4);
  __m256i k0, k1, k2, k3;
  masks(lo, hi, k0, k1, k2, k3);
  const __m256i pad= _mm256_set1_epi64x(PAD);
  return hbest(_mm256_blendv_epi8(pad, ld(p), k0), _mm256_blendv_epi8(pad, ld(p + 4), k1), _mm256_blendv_epi8(pad, ld(p + 8), k2), _mm256_blendv_epi8(pad, ld(p + 12), k3));
 }
 // 段 h の子 c (l を含む) から下の、[l, 子の終わり) の最小。
 i64 left_path(int h, int l) const {
  i64 acc= PAD;
  for(int g= h - 1;; --g) {
   const int sh= 4 * g, b= l >> (sh + 4), j= (l >> sh) & 15;
   if(g == 0 || !(l & ((1 << sh) - 1))) return best(acc, lanes_best(g, b, j, 16));
   if(j < 15) acc= best(acc, lanes_best(g, b, j + 1, 16));
  }
 }
 // 段 h の子 c (r - 1 を含む) から下の、[子の始まり, r) の最小。
 i64 right_path(int h, int r) const {
  const int x= r - 1;
  i64 acc= PAD;
  for(int g= h - 1;; --g) {
   const int sh= 4 * g, b= x >> (sh + 4), j= (x >> sh) & 15;
   if(g == 0 || !(r & ((1 << sh) - 1))) return best(acc, lanes_best(g, b, 0, j + 1));
   if(j > 0) acc= best(acc, lanes_best(g, b, 0, j));
  }
 }
 i64 range_best(int l, int r) const {
  const int x= r - 1;
  int h= H - 1, b= 0;
  for(;;) {
   const int sh= 4 * h, jl= (l >> sh) & 15, jr= (x >> sh) & 15;
   if(h == 0) return lanes_best(0, b, jl, jr + 1);
   const int sm= (1 << sh) - 1;
   if(jl != jr) {
    i64 res= PAD;
    int fl= jl, fr= jr + 1;
    if(l & sm) res= best(res, left_path(h, l)), fl= jl + 1;
    if(r & sm) res= best(res, right_path(h, r)), fr= jr;
    if(fl < fr) res= best(res, lanes_best(h, b, fl, fr));
    return res;
   }
   if(!(l & sm) && !(r & sm)) return lv(h)[(b << 4) + jl];
   b= (b << 4) + jl, --h;
  }
 }
 // 段 h の子 c の中で、check が最初に偽になる要素を、左から (L) か右から (!L) 探す。m はそれまでの最小。
 template <bool L, class C> int descend(int h, int c, i64 m, const C& check) const {
  for(int g= h - 1;; --g) {
   int j= L ? 0 : 15;
   for(;; j+= L ? 1 : -1) {
    const i64 t= best(m, lv(g)[(c << 4) + j]);
    if(!check(V(t))) break;
    m= t;
   }
   if(g == 0) return L ? (c << 4) + j : (c << 4) + j + 1;
   c= (c << 4) + j;
  }
 }
 template <class F> void build(const F& init) {
  int cnt[8], nb= std::max((n + 15) >> 4, 1);
  H= 0;
  std::size_t total= 0;
  for(;;) {
   cnt[H++]= nb, total+= std::size_t(nb) << 4;
   if(nb == 1) break;
   nb= (nb + 15) >> 4;
  }
  mem.assign(total, PAD);
  std::size_t o= 0;
  for(int h= 0; h < H; ++h) loff[h]= o, o+= std::size_t(cnt[h]) << 4;
  for(int i= 0; i < n; ++i) lv(0)[i]= i64(init(i));
  for(int h= 1; h < H; ++h)
   for(int c= 0; c < cnt[h - 1]; ++c) lv(h)[c]= block_best(h - 1, c);
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
 V get(int i) const { return V(lv(0)[i]); }
 V operator[](int i) const { return get(i); }
 [[gnu::always_inline]] void set(int i, V x) {
  put_lane(lv(0) + ((i >> 4) << 4), i & 15, i64(x));
  for(int h= 1, c= i >> 4; h < H; ++h, c>>= 4) {
   const i64 m= block_best(h - 1, c);
   if(lv(h)[c] == m) break;  // ここから上は変わらない
   put_lane(lv(h) + ((c >> 4) << 4), c & 15, m);
  }
 }
 // a[i] = op(a[i], x)。値はよくなる向きにしか動かないので、上へ best を取るだけでよい。
 [[gnu::always_inline]] void mul(int i, V x) {
  const i64 v= i64(x);
  for(int h= 0, c= i; h < H; ++h, c>>= 4) {
   if(!better(v, lv(h)[c])) break;
   put_lane(lv(h) + ((c >> 4) << 4), c & 15, v);
  }
 }
 [[gnu::always_inline]] V prod(int l, int r) const { return l < r ? V(range_best(l, r)) : TI; }
 V all_prod() const { return n ? V(block_best(H - 1, 0)) : TI; }
 // check(prod(l, r)) が真になる最大の r。l から右へ、塊の残りのレーン、その上の段の残りの子と広げて最小を集め、
 // 偽になる子が見つかったらそこへ下りる。
 template <class C> int max_right(int l, const C& check) const {
  assert(0 <= l && l <= n && check(TI));
  if(l == n) return n;
  i64 m= PAD;
  for(int h= 0; h < H; ++h) {
   const int sh= 4 * h, b= l >> (sh + 4);
   for(int j= ((l >> sh) & 15) + (h > 0); j < 16; ++j) {
    const int c= (b << 4) + j;
    if((i64(c) << sh) >= n) return n;
    const i64 t= best(m, lv(h)[c]);
    if(!check(V(t))) return h == 0 ? c : descend<true>(h, c, m, check);
    m= t;
   }
  }
  return n;
 }
 // check(prod(l, r)) が真になる最小の l。r - 1 から左へ広げる。
 template <class C> int min_left(int r, const C& check) const {
  assert(0 <= r && r <= n && check(TI));
  if(r == 0) return 0;
  const int x= r - 1;
  i64 m= PAD;
  for(int h= 0; h < H; ++h) {
   const int sh= 4 * h, b= x >> (sh + 4);
   for(int j= ((x >> sh) & 15) - (h > 0); j >= 0; --j) {
    const int c= (b << 4) + j;
    const i64 t= best(m, lv(h)[c]);
    if(!check(V(t))) return h == 0 ? c + 1 : descend<false>(h, c, m, check);
    m= t;
   }
  }
  return 0;
 }
};
template <class V> constexpr bool use_wide= std::is_integral_v<V> && std::is_signed_v<V> && sizeof(V) <= 8;
template <class V, bool MAX, class M> using Impl= std::conditional_t<use_wide<V>, Wide<V, MAX>, segment_tree_internal::Generic<M>>;
}
template <class V> class SegmentTree<monoid::Min<V>>: public segtree_min_internal::Impl<V, false, monoid::Min<V>> {
 using Base= segtree_min_internal::Impl<V, false, monoid::Min<V>>;
public:
 using Base::Base;
};
template <class V> class SegmentTree<monoid::Max<V>>: public segtree_min_internal::Impl<V, true, monoid::Max<V>> {
 using Base= segtree_min_internal::Impl<V, true, monoid::Max<V>>;
public:
 using Base::Base;
};
