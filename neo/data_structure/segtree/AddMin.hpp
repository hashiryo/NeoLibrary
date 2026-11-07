#pragma once
#include <limits>
#include "neo/data_structure/SegmentTree.hpp"
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
// 区間加算と区間最小 monoid::AddMin<T>、区間加算と区間最大 monoid::AddMax<T>。SegmentTree<monoid::AddMin<T>> は、
// T が 64 bit 以下の符号付き整数なら、子を 16 個持つ作用を配らない木になる (procon-judge の
// yosupo-range-add-range-min の wide3)。段 0 は要素を 16 個ずつの塊に並べ、段 h >= 1 の塊は、子の部分木全体に
// 足した量 t と、子の最小に t を足した s を子 16 個ぶん並べて持つ。加算は最小と可換なので作用を子へ配らない。
// 区間の加算は下から全部掛かるレーンの s と t に AVX2 でまとめて足し、下から両端の道の s を取り直す。区間の最小は、
// 根から l と r - 1 が同じ子にある間は t を足しながら下り、分かれたら左右の道を下りて、全部掛かるレーンの最小に
// そこまでの t を足して比べる。値は i64 で、16 個を 4 本の __m256i で扱う (gcc が __m256i の配列をスタックへ
// 退避するので、マスクは別々の変数で作る)。要素の無いレーンは、最小なら i64 の最大、最大なら i64 の最小。
// それ以外の T は二分木の汎用の形になる。
namespace monoid {
template <class V> struct AddMin {
 using T= V;
 using E= V;
 static constexpr T ti() { return std::numeric_limits<V>::max(); }
 static constexpr T op(T a, T b) { return b < a ? b : a; }
 static constexpr void mp(T& v, E x) { v+= x; }
 static constexpr void cp(E& pre, E suf) { pre+= suf; }
};
template <class V> struct AddMax {
 using T= V;
 using E= V;
 static constexpr T ti() { return std::numeric_limits<V>::lowest(); }
 static constexpr T op(T a, T b) { return a < b ? b : a; }
 static constexpr void mp(T& v, E x) { v+= x; }
 static constexpr void cp(E& pre, E suf) { pre+= suf; }
};
}
namespace segtree_add_min_internal {
using i64= long long;
inline __m256i ld(const i64* p) { return _mm256_load_si256((const __m256i*)p); }
inline void st(i64* p, __m256i x) { _mm256_store_si256((__m256i*)p, x); }
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
struct alignas(64) Blk {
 i64 s[16], t[16];
};
template <class V, bool MAX> class Wide {
 static constexpr i64 PAD= MAX ? std::numeric_limits<i64>::min() : std::numeric_limits<i64>::max();
 static constexpr V TI= MAX ? std::numeric_limits<V>::lowest() : std::numeric_limits<V>::max();
 int n= 0, H= 0;
 std::size_t loff[8]= {};               // 段 h >= 1 の塊は mem1 の loff[h] から。写しても壊れないよう、ポインタでなくずれで持つ
 segment_tree_internal::Vec<i64> mem0;  // 段 0 の要素
 segment_tree_internal::Vec<Blk> mem1;
 i64* V0() { return mem0.data(); }
 const i64* V0() const { return mem0.data(); }
 Blk* LV(int h) { return mem1.data() + loff[h]; }
 const Blk* LV(int h) const { return mem1.data() + loff[h]; }
 static i64 best(i64 a, i64 b) {
  if constexpr(MAX) return a < b ? b : a;
  else return b < a ? b : a;
 }
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
 const i64* row(int h, int b) const { return h == 0 ? V0() + (std::size_t(b) << 4) : LV(h)[b].s; }
 i64 tag(int h, int c) const { return LV(h)[c >> 4].t[c & 15]; }
 void add_lanes(int h, int b, int lo, int hi, i64 x) {
  __m256i k0, k1, k2, k3;
  masks(lo, hi, k0, k1, k2, k3);
  const __m256i xv= _mm256_set1_epi64x(x);
  const __m256i x0= _mm256_and_si256(k0, xv), x1= _mm256_and_si256(k1, xv), x2= _mm256_and_si256(k2, xv), x3= _mm256_and_si256(k3, xv);
  auto upd= [&](i64* p) {
   st(p, _mm256_add_epi64(ld(p), x0)), st(p + 4, _mm256_add_epi64(ld(p + 4), x1));
   st(p + 8, _mm256_add_epi64(ld(p + 8), x2)), st(p + 12, _mm256_add_epi64(ld(p + 12), x3));
  };
  if(h == 0) return upd(V0() + (std::size_t(b) << 4));
  upd(LV(h)[b].s), upd(LV(h)[b].t);
 }
 i64 block_best(int h, int b) const {
  const i64* p= row(h, b);
  return hbest(ld(p), ld(p + 4), ld(p + 8), ld(p + 12));
 }
 i64 lanes_best(int h, int b, int lo, int hi) const {
  const i64* p= row(h, b);
  __m256i k0, k1, k2, k3;
  masks(lo, hi, k0, k1, k2, k3);
  const __m256i pad= _mm256_set1_epi64x(PAD);
  return hbest(_mm256_blendv_epi8(pad, ld(p), k0), _mm256_blendv_epi8(pad, ld(p + 4), k1), _mm256_blendv_epi8(pad, ld(p + 8), k2), _mm256_blendv_epi8(pad, ld(p + 12), k3));
 }
 // 段 h の子 c (l を含む) から下の、[l, 子の終わり) の最小。子の t まで足した値。
 i64 left_path(int h, int c, int l) const {
  i64 G= tag(h, c), acc= PAD;
  for(int g= h - 1;; --g) {
   const int sh= 4 * g, b= l >> (sh + 4), j= (l >> sh) & 15;
   if(g == 0 || !(l & ((1 << sh) - 1))) return best(acc, G + lanes_best(g, b, j, 16));
   if(j < 15) acc= best(acc, G + lanes_best(g, b, j + 1, 16));
   G+= LV(g)[b].t[j];
  }
 }
 // 段 h の子 c (r - 1 を含む) から下の、[子の始まり, r) の最小。
 i64 right_path(int h, int c, int r) const {
  const int x= r - 1;
  i64 G= tag(h, c), acc= PAD;
  for(int g= h - 1;; --g) {
   const int sh= 4 * g, b= x >> (sh + 4), j= (x >> sh) & 15;
   if(g == 0 || !(r & ((1 << sh) - 1))) return best(acc, G + lanes_best(g, b, 0, j + 1));
   if(j > 0) acc= best(acc, G + lanes_best(g, b, 0, j));
   G+= LV(g)[b].t[j];
  }
 }
 i64 range_best(int l, int r) const {
  const int x= r - 1;
  i64 F= 0;
  int h= H - 1, b= 0;
  for(;;) {
   const int sh= 4 * h, jl= (l >> sh) & 15, jr= (x >> sh) & 15;
   if(h == 0) return F + lanes_best(0, b, jl, jr + 1);
   const int sm= (1 << sh) - 1;
   if(jl != jr) {
    i64 res= PAD;
    int fl= jl, fr= jr + 1;
    if(l & sm) res= best(res, left_path(h, (b << 4) + jl, l)), fl= jl + 1;
    if(r & sm) res= best(res, right_path(h, (b << 4) + jr, r)), fr= jr;
    if(fl < fr) res= best(res, lanes_best(h, b, fl, fr));
    return F + res;
   }
   if(!(l & sm) && !(r & sm)) return F + LV(h)[b].s[jl];
   F+= LV(h)[b].t[jl];
   b= (b << 4) + jl, --h;
  }
 }
 void range_add(int l, int r, i64 x) {
  for(int h= 0, p= l, q= r; p < q; ++h) {
   const int bp= p >> 4, bq= q >> 4;
   if(bp == bq || h == H - 1) {
    add_lanes(h, bp, p & 15, q - (bp << 4), x);
    break;
   }
   if(p & 15) add_lanes(h, bp, p & 15, 16, x), p= bp + 1;
   else p= bp;
   if(q & 15) add_lanes(h, bq, 0, q & 15, x), q= bq;
   else q= bq;
  }
  const int y= r - 1;
  for(int h= 1; h < H; ++h) {
   const int sh= 4 * h, sm= (1 << sh) - 1;
   if(l & sm) {
    Blk& B= LV(h)[l >> (sh + 4)];
    const int j= (l >> sh) & 15;
    B.s[j]= block_best(h - 1, l >> sh) + B.t[j];
   }
   if(r & sm) {
    Blk& B= LV(h)[y >> (sh + 4)];
    const int j= (y >> sh) & 15;
    B.s[j]= block_best(h - 1, y >> sh) + B.t[j];
   }
  }
 }
 // l (または r - 1) を含む道の、段 h の値に足す量 (段 h より上の t の和) を off[h] に書く。
 void path_offsets(int p, i64* off) const {
  off[H - 1]= 0;
  for(int h= H - 1; h >= 1; --h) off[h - 1]= off[h] + tag(h, p >> (4 * h));
 }
 // 段 h の子 c の中で、check が最初に偽になる要素を、左から (L) か右から (!L) 探す。m はそれまでの最小、
 // off は段 h - 1 の値に足す量。子の全体を足すと偽になることが分かっている。
 template <bool L, class C> int descend(int h, int c, i64 off, i64 m, const C& check) const {
  for(int g= h - 1;; --g) {
   int j= L ? 0 : 15;
   for(;; j+= L ? 1 : -1) {
    const int cc= (c << 4) + j;
    const i64 t= best(m, (g == 0 ? V0()[cc] : LV(g)[c].s[j]) + off);
    if(!check(V(t))) break;
    m= t;
   }
   const int cc= (c << 4) + j;
   if(g == 0) return L ? cc : cc + 1;
   off+= LV(g)[c].t[j], c= cc;
  }
 }
 template <class F> void build(const F& init) {
  int cnt[8], nb= (n + 15) >> 4;
  if(nb == 0) nb= 1;
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
  for(int i= 0; i < n; ++i) V0()[i]= i64(init(i));
  for(int h= 1; h < H; ++h)
   for(int b= 0; b < cnt[h]; ++b)
    for(int j= 0; j < 16; ++j) {
     const int c= b * 16 + j;
     LV(h)[b].s[j]= c < cnt[h - 1] ? block_best(h - 1, c) : PAD;
     LV(h)[b].t[j]= 0;
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
 V get(int i) const {
  i64 v= V0()[i];
  for(int h= 1; h < H; ++h) v+= tag(h, i >> (4 * h));
  return V(v);
 }
 V operator[](int i) const { return get(i); }
 void set(int i, V x) { range_add(i, i + 1, i64(x) - i64(get(i))); }
 void mul(int i, V x) { set(i, V(best(i64(get(i)), i64(x)))); }
 [[gnu::always_inline]] V prod(int l, int r) const { return l < r ? V(range_best(l, r)) : TI; }
 V all_prod() const { return n ? V(block_best(H - 1, 0)) : TI; }
 [[gnu::always_inline]] void apply(int l, int r, V x) {
  if(l < r) range_add(l, r, i64(x));
 }
 // check(prod(l, r)) が真になる最大の r。l から右へ、塊の残りのレーン、その上の段の残りの子と広げて最小を集め、
 // 偽になる子が見つかったらそこへ下りる。
 template <class C> int max_right(int l, const C& check) const {
  assert(0 <= l && l <= n && check(TI));
  if(l == n) return n;
  i64 off[8];
  path_offsets(l, off);
  i64 m= PAD;
  for(int j= l & 15, base= l & ~15; j < 16; ++j) {
   if(base + j >= n) return n;
   const i64 t= best(m, V0()[base + j] + off[0]);
   if(!check(V(t))) return base + j;
   m= t;
  }
  for(int h= 1; h < H; ++h) {
   const int sh= 4 * h, b= l >> (sh + 4);
   for(int j= ((l >> sh) & 15) + 1; j < 16; ++j) {
    const int c= (b << 4) + j;
    if((i64(c) << sh) >= n) return n;
    const i64 t= best(m, LV(h)[b].s[j] + off[h]);
    if(!check(V(t))) return descend<true>(h, c, off[h] + LV(h)[b].t[j], m, check);
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
  i64 off[8];
  path_offsets(x, off);
  i64 m= PAD;
  for(int j= x & 15, base= x & ~15; j >= 0; --j) {
   const i64 t= best(m, V0()[base + j] + off[0]);
   if(!check(V(t))) return base + j + 1;
   m= t;
  }
  for(int h= 1; h < H; ++h) {
   const int sh= 4 * h, b= x >> (sh + 4);
   for(int j= ((x >> sh) & 15) - 1; j >= 0; --j) {
    const int c= (b << 4) + j;
    const i64 t= best(m, LV(h)[b].s[j] + off[h]);
    if(!check(V(t))) return descend<false>(h, c, off[h] + LV(h)[b].t[j], m, check);
    m= t;
   }
  }
  return 0;
 }
};
template <class V> constexpr bool use_wide= std::is_integral_v<V> && std::is_signed_v<V> && sizeof(V) <= 8;
template <class V, bool MAX, class M> using Impl= std::conditional_t<use_wide<V>, Wide<V, MAX>, segment_tree_internal::Generic<M>>;
}
template <class V> class SegmentTree<monoid::AddMin<V>>: public segtree_add_min_internal::Impl<V, false, monoid::AddMin<V>> {
 using Base= segtree_add_min_internal::Impl<V, false, monoid::AddMin<V>>;
public:
 using Base::Base;
};
template <class V> class SegmentTree<monoid::AddMax<V>>: public segtree_add_min_internal::Impl<V, true, monoid::AddMax<V>> {
 using Base= segtree_add_min_internal::Impl<V, true, monoid::AddMax<V>>;
public:
 using Base::Base;
};
