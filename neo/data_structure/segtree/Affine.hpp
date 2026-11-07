#pragma once
#include <array>
#include <algorithm>
#include "neo/algebra/ZMod.hpp"
#include "neo/data_structure/SegmentTree.hpp"
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
// 区間の 1 次関数と区間和 monoid::AffineSum<Z> (遅延)、区間の 1 次関数と 1 点取得 monoid::Affine<Z> (双対)。作用 {b, c} は
// x -> b x + c。Z が法を定数にした ZMod<MOD> なら、SegmentTree はどちらも子を 16 個持つ木になる (procon-judge の
// yosupo-range-affine-range-sum の wide_lazy3、yosupo-range-affine-point-get の wide_dual2_pf)。段 0 は要素を 16 個ずつ
// の塊に並べ、段 h >= 1 の塊は子 (段 h - 1 の塊) 16 個ぶんの和 s (遅延のときだけ) と、まだ配っていない作用 (m, a) を
// 並べて持つ。s はその子の作用まで当てた値で、祖先の作用は当てていない。区間の作用は、根から両端の道を配り、下から
// 全部掛かるレーンへ AVX2 でまとめて当て、下から両端の道の和を取り直す。初めにその操作が下りる 2 本の道の塊を先読み
// する。区間の和は書き込みをせず、根から l と r - 1 が同じ子にある間は作用を合成しながら下り、分かれたら左右の道を
// 下りて、全部掛かるレーンの和にそこまでに合成した作用を当てて足す。値は剰余をそのままの値で 2p 未満に持ち、掛け算
// は Shoup の方法 (掛ける数 m に floor(m 2^32 / p) を前計算) で 8 レーンずつ。区間の長さは段ごとに決まり、最後の
// 塊だけレーンごとの長さの表を使う。1 レーンの書き込みも 32 byte をベクトルで読んで混ぜて書き、あとのベクトルの読み込みが
// 書き込みから値を受け取れるようにする。1 点の取得は、葉の値に段 1 から順に作用を当てる (根に近いほど新しい作用)。
// それ以外の Z (実行時の法の ZMod<0> や自前の剰余の型) は二分木の汎用の形になる。
namespace monoid {
template <class Z> struct AffineSum {
 using T= Z;
 using E= std::array<Z, 2>;
 static constexpr T ti() { return T(); }
 static constexpr T op(T a, T b) { return a + b; }
 static constexpr void mp(T& v, const E& f, int len) { v= f[0] * v + f[1] * T(len); }
 static constexpr void cp(E& pre, const E& suf) { pre= {suf[0] * pre[0], suf[0] * pre[1] + suf[1]}; }
};
template <class Z> struct Affine {
 using T= Z;
 using E= std::array<Z, 2>;
 static constexpr void mp(T& v, const E& f) { v= f[0] * v + f[1]; }
 static constexpr void cp(E& pre, const E& suf) { pre= {suf[0] * pre[0], suf[0] * pre[1] + suf[1]}; }
};
}
namespace segtree_affine_internal {
using u32= unsigned;
using u64= unsigned long long;
template <unsigned P> struct Arith {
 static constexpr u32 P2= 2 * P;
 // 掛ける数 x (p 未満に直す) と、Shoup の商の近似 floor(x 2^32 / p) を 8 レーンに配ったもの。
 struct Mul {
  __m256i m, q;
  explicit Mul(u32 x) {
   x= x >= P ? x - P : x;
   m= _mm256_set1_epi32(int(x)), q= _mm256_set1_epi32(int(u32((u64(x) << 32) / P)));
  }
 };
 // x m mod p を 2p 未満で。x は u32 なら何でもよい。
 static __m256i mul(__m256i x, const Mul& f) {
  const __m256i e= _mm256_mul_epu32(x, f.q), o= _mm256_mul_epu32(_mm256_srli_epi64(x, 32), f.q);
  const __m256i q= _mm256_blend_epi32(_mm256_srli_epi64(e, 32), o, 0xaa);
  return _mm256_sub_epi32(_mm256_mullo_epi32(x, f.m), _mm256_mullo_epi32(q, _mm256_set1_epi32(int(P))));
 }
 // a + b を 2p 未満で。a + b < 4p。
 static __m256i add(__m256i a, __m256i b) {
  const __m256i s= _mm256_add_epi32(a, b);
  return _mm256_min_epu32(s, _mm256_sub_epi32(s, _mm256_set1_epi32(int(P2))));
 }
 static u32 add1(u32 a, u32 b) {
  const u32 s= a + b;
  return s >= P2 ? s - P2 : s;
 }
 // 16 レーン (2p 未満) の合計を 2p 未満で。
 static u32 hsum(__m256i x0, __m256i x1) {
  const __m256i s= add(x0, x1);
  const __m128i p2= _mm_set1_epi32(int(P2));
  __m128i t= _mm_add_epi32(_mm256_castsi256_si128(s), _mm256_extracti128_si256(s, 1));
  t= _mm_min_epu32(t, _mm_sub_epi32(t, p2));
  t= _mm_add_epi32(t, _mm_shuffle_epi32(t, 0x4e));
  t= _mm_min_epu32(t, _mm_sub_epi32(t, p2));
  t= _mm_add_epi32(t, _mm_shuffle_epi32(t, 0xb1));
  t= _mm_min_epu32(t, _mm_sub_epi32(t, p2));
  return u32(_mm_cvtsi128_si32(t));
 }
};
// レーン [lo, hi) を選ぶマスク (16 レーンを 2 本)。
inline void masks(int lo, int hi, __m256i& k0, __m256i& k1) {
 const __m256i i0= _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7), i1= _mm256_setr_epi32(8, 9, 10, 11, 12, 13, 14, 15);
 const __m256i L= _mm256_set1_epi32(lo - 1), R= _mm256_set1_epi32(hi);
 k0= _mm256_and_si256(_mm256_cmpgt_epi32(i0, L), _mm256_cmpgt_epi32(R, i0));
 k1= _mm256_and_si256(_mm256_cmpgt_epi32(i1, L), _mm256_cmpgt_epi32(R, i1));
}
inline __m256i ld(const u32* p) { return _mm256_load_si256((const __m256i*)p); }
inline void st(u32* p, __m256i x) { _mm256_store_si256((__m256i*)p, x); }
// p は 16 レーンの塊の先頭 (32 byte 境界)。j 番目のレーンを v にする。スカラーで書くと、あとで同じ塊をベクトルで読む
// ときに書き込みから値を受け取れず、書き込みが終わるまで待つ。ベクトルで読んで混ぜて書けば受け取れる。
inline void put_lane(u32* p, int j, u32 v) {
 __m256i* q= (__m256i*)(p + (j & ~7));
 const __m256i k= _mm256_cmpeq_epi32(_mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7), _mm256_set1_epi32(j & 7));
 _mm256_store_si256(q, _mm256_blendv_epi8(_mm256_load_si256(q), _mm256_set1_epi32(int(v)), k));
}
// 段の数を数え、段 0 の要素の置き場と段 h >= 1 の塊の置き場を作る。cnt[h] は段 h の塊の数。
template <class Blk> int layout(int n, int* cnt, segment_tree_internal::Vec<u32>& mem0, segment_tree_internal::Vec<Blk>& mem1, std::size_t* off) {
 int H= 0, nb= std::max((n + 15) >> 4, 1);
 std::size_t total1= 0;
 for(;;) {
  cnt[H++]= nb;
  if(nb == 1) break;
  nb= (nb + 15) >> 4;
  total1+= nb;
 }
 mem0.assign(std::size_t(cnt[0]) << 4, 0);
 mem1.resize(std::max<std::size_t>(total1, 1));
 std::size_t o= 0;
 for(int h= 1; h < H; ++h) off[h]= o, o+= cnt[h];
 return H;
}
struct alignas(64) LBlk {
 u32 s[16], m[16], a[16];
};
template <class Z, unsigned P> class Lazy {
 using A= Arith<P>;
 using Mul= typename A::Mul;
 using E= std::array<Z, 2>;
 int n= 0, H= 0;
 std::size_t loff[8]= {};  // 段 h >= 1 の塊は mem1 の loff[h] から。写しても壊れないよう、ポインタでなくずれで持つ
 u32* V0() { return mem0.data(); }
 const u32* V0() const { return mem0.data(); }
 LBlk* LV(int h) { return mem1.data() + loff[h]; }
 const LBlk* LV(int h) const { return mem1.data() + loff[h]; }
 int last[8]= {};
 u32 pw[8]= {};
 alignas(32) u32 tl[8][16]= {};
 segment_tree_internal::Vec<u32> mem0;
 segment_tree_internal::Vec<LBlk> mem1;
 // 段 h の塊 b のレーン [lo, hi) に作用 x -> m x + a を当てる。al は長さの揃った塊で足す値 a 16^h mod p。
 template <bool ALL> void apply_lanes(int h, int b, int lo, int hi, const Mul& f, u32 a, u32 al) {
  __m256i add0, add1;
  if(b == last[h]) {
   const Mul g(a);
   add0= A::mul(ld(tl[h]), g), add1= A::mul(ld(tl[h] + 8), g);
  } else add0= add1= _mm256_set1_epi32(int(al));
  __m256i k0= _mm256_setzero_si256(), k1= k0;
  if constexpr(!ALL) masks(lo, hi, k0, k1);
  auto upd= [&](u32* p, __m256i x0, __m256i x1) {
   if constexpr(ALL) st(p, x0), st(p + 8, x1);
   else st(p, _mm256_blendv_epi8(ld(p), x0, k0)), st(p + 8, _mm256_blendv_epi8(ld(p + 8), x1, k1));
  };
  if(h == 0) {
   u32* p= V0() + (std::size_t(b) << 4);
   upd(p, A::add(A::mul(ld(p), f), add0), A::add(A::mul(ld(p + 8), f), add1));
   return;
  }
  LBlk& B= LV(h)[b];
  upd(B.s, A::add(A::mul(ld(B.s), f), add0), A::add(A::mul(ld(B.s + 8), f), add1));
  upd(B.m, A::mul(ld(B.m), f), A::mul(ld(B.m + 8), f));
  const __m256i av= _mm256_set1_epi32(int(a));
  upd(B.a, A::add(A::mul(ld(B.a), f), av), A::add(A::mul(ld(B.a + 8), f), av));
 }
 // 段 h の子 c (段 h - 1 の塊 c) へ作用を配る。
 void push(int h, int c) {
  LBlk& B= LV(h)[c >> 4];
  const int j= c & 15;
  const u32 m= B.m[j], a= B.a[j];
  if(m == 1 && a == 0) return;
  put_lane(B.m, j, 1), put_lane(B.a, j, 0);
  apply_lanes<true>(h - 1, c, 0, 16, Mul(m), a, u32(u64(a) * pw[h - 1] % P));
 }
 const u32* row(int h, int b) const { return h == 0 ? V0() + (std::size_t(b) << 4) : LV(h)[b].s; }
 u32 block_sum(int h, int b) const {
  const u32* p= row(h, b);
  return A::hsum(ld(p), ld(p + 8));
 }
 u32 lanes_sum(int h, int b, int lo, int hi) const {
  const u32* p= row(h, b);
  __m256i k0, k1;
  masks(lo, hi, k0, k1);
  return A::hsum(_mm256_and_si256(ld(p), k0), _mm256_and_si256(ld(p + 8), k1));
 }
 // l と r - 1 を含む段 0 から段 H - 2 の塊を先読みする。根の塊はいつも触るので要らない。
 void prefetch(int l, int r) const {
  const int x= r - 1;
  _mm_prefetch((const char*)(V0() + ((std::size_t(l) >> 4) << 4)), _MM_HINT_T0);
  _mm_prefetch((const char*)(V0() + ((std::size_t(x) >> 4) << 4)), _MM_HINT_T0);
  for(int h= 1; h + 1 < H; ++h) {
   const char *p= (const char*)(LV(h) + (l >> (4 * h + 4))), *q= (const char*)(LV(h) + (x >> (4 * h + 4)));
   _mm_prefetch(p, _MM_HINT_T0), _mm_prefetch(p + 64, _MM_HINT_T0), _mm_prefetch(p + 128, _MM_HINT_T0);
   _mm_prefetch(q, _MM_HINT_T0), _mm_prefetch(q + 64, _MM_HINT_T0), _mm_prefetch(q + 128, _MM_HINT_T0);
  }
 }
 void affine(int l, int r, u32 m, u32 a) {
  prefetch(l, r);
  const int x= r - 1;
  for(int h= H - 1; h >= 1; --h) {
   const int sh= 4 * h, sm= (1 << sh) - 1;
   if(l & sm) push(h, l >> sh);
   if(r & sm) push(h, x >> sh);
  }
  const Mul f(m);
  u32 al= a;
  for(int h= 0, p= l, q= r; p < q; ++h) {
   const int bp= p >> 4, bq= q >> 4;
   if(bp == bq || h == H - 1) {
    apply_lanes<false>(h, bp, p & 15, q - (bp << 4), f, a, al);
    break;
   }
   if(p & 15) apply_lanes<false>(h, bp, p & 15, 16, f, a, al), p= bp + 1;
   else p= bp;
   if(q & 15) apply_lanes<false>(h, bq, 0, q & 15, f, a, al), q= bq;
   else q= bq;
   al= u32(u64(al) * 16 % P);
  }
  for(int h= 1; h < H; ++h) {
   const int sh= 4 * h, sm= (1 << sh) - 1;
   if(l & sm) put_lane(LV(h)[l >> (sh + 4)].s, (l >> sh) & 15, block_sum(h - 1, l >> sh));
   if(r & sm) put_lane(LV(h)[x >> (sh + 4)].s, (x >> sh) & 15, block_sum(h - 1, x >> sh));
  }
 }
 // 区間の和に当てる 1 次関数 (x, 個数) -> m x + a 個数。m と a は p 未満。
 struct Af {
  u64 m, a;
 };
 static Af comp(Af F, u32 tm, u32 ta) { return {F.m * tm % P, (F.m * ta + F.a) % P}; }
 static u64 app(Af F, u64 x, u64 cnt) { return (F.m * x + F.a * cnt) % P; }
 // 段 h の子 c (l を含む) から下の、[l, 子の終わり) の和。子の作用まで当てた値。
 u64 left_path(int h, int c, int l) const {
  Af G{LV(h)[c >> 4].m[c & 15] % P, LV(h)[c >> 4].a[c & 15] % P};
  u64 acc= 0;
  for(int g= h - 1;; --g) {
   const int sh= 4 * g, b= l >> (sh + 4), j= (l >> sh) & 15, end= (b + 1) << (sh + 4);
   if(g == 0 || !(l & ((1 << sh) - 1))) return (acc + app(G, lanes_sum(g, b, j, 16), u64(end - l))) % P;
   if(j < 15) acc+= app(G, lanes_sum(g, b, j + 1, 16), u64(end - (((l >> sh) + 1) << sh)));
   G= comp(G, LV(g)[b].m[j], LV(g)[b].a[j]);
  }
 }
 // 段 h の子 c (r - 1 を含む) から下の、[子の始まり, r) の和。
 u64 right_path(int h, int c, int r) const {
  const int x= r - 1;
  Af G{LV(h)[c >> 4].m[c & 15] % P, LV(h)[c >> 4].a[c & 15] % P};
  u64 acc= 0;
  for(int g= h - 1;; --g) {
   const int sh= 4 * g, b= x >> (sh + 4), j= (x >> sh) & 15, start= b << (sh + 4);
   if(g == 0 || !(r & ((1 << sh) - 1))) return (acc + app(G, lanes_sum(g, b, 0, j + 1), u64(r - start))) % P;
   if(j > 0) acc+= app(G, lanes_sum(g, b, 0, j), u64(((x >> sh) << sh) - start));
   G= comp(G, LV(g)[b].m[j], LV(g)[b].a[j]);
  }
 }
 u32 sum(int l, int r) const {
  prefetch(l, r);
  const int x= r - 1;
  Af F{1, 0};
  int h= H - 1, b= 0;
  for(;;) {
   const int sh= 4 * h, jl= (l >> sh) & 15, jr= (x >> sh) & 15;
   if(h == 0) return u32(app(F, lanes_sum(0, b, jl, jr + 1), u64(r - l)));
   const int sm= (1 << sh) - 1;
   if(jl != jr) {
    u64 res= 0;
    int fl= jl, fr= jr + 1;
    if(l & sm) res+= left_path(h, (b << 4) + jl, l), fl= jl + 1;
    if(r & sm) res+= right_path(h, (b << 4) + jr, r), fr= jr;
    if(fl < fr) res+= lanes_sum(h, b, fl, fr);
    return u32(app(F, res % P, u64(r - l)));
   }
   if(!(l & sm) && !(r & sm)) return u32(app(F, LV(h)[b].s[jl], u64(r - l)));
   F= comp(F, LV(h)[b].m[jl], LV(h)[b].a[jl]);
   b= (b << 4) + jl, --h;
  }
 }
 template <class F> void build(const F& init) {
  int cnt[8];
  H= layout(n, cnt, mem0, mem1, loff);
  for(int h= 0; h < H; ++h) {
   last[h]= cnt[h] - 1, pw[h]= u32((u64(1) << (4 * h)) % P);
   for(int j= 0; j < 16; ++j) {
    const long long start= (static_cast<long long>(last[h]) * 16 + j) << (4 * h);
    tl[h][j]= u32(std::clamp<long long>(n - start, 0, 1LL << (4 * h)));
   }
  }
  for(int i= 0; i < n; ++i) V0()[i]= Z(init(i)).val();
  for(int h= 1; h < H; ++h)
   for(int b= 0; b < cnt[h]; ++b)
    for(int j= 0; j < 16; ++j) {
     const int c= b * 16 + j;
     LV(h)[b].s[j]= c < cnt[h - 1] ? block_sum(h - 1, c) : 0, LV(h)[b].m[j]= 1, LV(h)[b].a[j]= 0;
    }
 }
public:
 Lazy()= default;
 explicit Lazy(int n): n(n) {
  build([](int) { return Z(); });
 }
 Lazy(int n, Z x): n(n) {
  build([&](int) { return x; });
 }
 Lazy(const std::vector<Z>& a): n(a.size()) {
  build([&](int i) { return a[i]; });
 }
 template <class F, std::enable_if_t<std::is_invocable_r_v<Z, const F&, int>, int> = 0> Lazy(int n, const F& init): n(n) { build(init); }
 int size() const { return n; }
 Z get(int i) const {
  u64 v= V0()[i];
  for(int h= 1; h < H; ++h) {
   const int c= i >> (4 * h);
   v= (u64(LV(h)[c >> 4].m[c & 15]) * v + LV(h)[c >> 4].a[c & 15]) % P;
  }
  return Z::raw(u32(v));
 }
 Z operator[](int i) const { return get(i); }
 void set(int i, Z x) {
  for(int h= H - 1; h >= 1; --h) push(h, i >> (4 * h));
  put_lane(V0() + ((i >> 4) << 4), i & 15, x.val());
  for(int h= 1; h < H; ++h) put_lane(LV(h)[i >> (4 * h + 4)].s, (i >> (4 * h)) & 15, block_sum(h - 1, i >> (4 * h)));
 }
 void mul(int i, Z x) { set(i, get(i) + x); }
 [[gnu::always_inline]] Z prod(int l, int r) const { return l < r ? Z::raw(sum(l, r)) : Z(); }
 Z all_prod() const { return Z::raw(block_sum(H - 1, 0)); }
 [[gnu::always_inline]] void apply(int l, int r, const E& f) {
  if(l < r) affine(l, r, f[0].val(), f[1].val());
 }
 // 和に単調な条件はふつう無いので、prod を使った二分探索で済ませる。
 template <class C> int max_right(int l, const C& check) const {
  assert(0 <= l && l <= n && check(Z()));
  int lo= l, hi= n + 1;
  while(hi - lo > 1) {
   const int mid= (lo + hi) >> 1;
   (check(prod(l, mid)) ? lo : hi)= mid;
  }
  return lo;
 }
 template <class C> int min_left(int r, const C& check) const {
  assert(0 <= r && r <= n && check(Z()));
  int lo= -1, hi= r;
  while(hi - lo > 1) {
   const int mid= (lo + hi) >> 1;
   (check(prod(mid, r)) ? hi : lo)= mid;
  }
  return hi;
 }
};
struct alignas(64) DBlk {
 u32 m[16], a[16];
};
template <class Z, unsigned P> class Dual {
 using A= Arith<P>;
 using Mul= typename A::Mul;
 using E= std::array<Z, 2>;
 int n= 0, H= 0;
 std::size_t loff[8]= {};  // 段 h >= 1 の塊は mem1 の loff[h] から。写しても壊れないよう、ポインタでなくずれで持つ
 u32* V0() { return mem0.data(); }
 const u32* V0() const { return mem0.data(); }
 DBlk* LV(int h) { return mem1.data() + loff[h]; }
 const DBlk* LV(int h) const { return mem1.data() + loff[h]; }
 segment_tree_internal::Vec<u32> mem0;
 segment_tree_internal::Vec<DBlk> mem1;
 template <bool ALL> void apply_lanes(int h, int b, int lo, int hi, const Mul& f, u32 a) {
  const __m256i av= _mm256_set1_epi32(int(a));
  __m256i k0= _mm256_setzero_si256(), k1= k0;
  if constexpr(!ALL) masks(lo, hi, k0, k1);
  auto upd= [&](u32* p, __m256i x0, __m256i x1) {
   if constexpr(ALL) st(p, x0), st(p + 8, x1);
   else st(p, _mm256_blendv_epi8(ld(p), x0, k0)), st(p + 8, _mm256_blendv_epi8(ld(p + 8), x1, k1));
  };
  if(h == 0) {
   u32* p= V0() + (std::size_t(b) << 4);
   upd(p, A::add(A::mul(ld(p), f), av), A::add(A::mul(ld(p + 8), f), av));
   return;
  }
  DBlk& B= LV(h)[b];
  upd(B.m, A::mul(ld(B.m), f), A::mul(ld(B.m + 8), f));
  upd(B.a, A::add(A::mul(ld(B.a), f), av), A::add(A::mul(ld(B.a + 8), f), av));
 }
 void push(int h, int c) {
  DBlk& B= LV(h)[c >> 4];
  const int j= c & 15;
  const u32 m= B.m[j], a= B.a[j];
  if(m == 1 && a == 0) return;
  put_lane(B.m, j, 1), put_lane(B.a, j, 0);
  apply_lanes<true>(h - 1, c, 0, 16, Mul(m), a);
 }
 void prefetch(int l, int r) const {
  const int x= r - 1;
  _mm_prefetch((const char*)(V0() + ((std::size_t(l) >> 4) << 4)), _MM_HINT_T0);
  _mm_prefetch((const char*)(V0() + ((std::size_t(x) >> 4) << 4)), _MM_HINT_T0);
  for(int h= 1; h + 1 < H; ++h) {
   const char *p= (const char*)(LV(h) + (l >> (4 * h + 4))), *q= (const char*)(LV(h) + (x >> (4 * h + 4)));
   _mm_prefetch(p, _MM_HINT_T0), _mm_prefetch(p + 64, _MM_HINT_T0);
   _mm_prefetch(q, _MM_HINT_T0), _mm_prefetch(q + 64, _MM_HINT_T0);
  }
 }
 template <class F> void build(const F& init) {
  int cnt[8];
  H= layout(n, cnt, mem0, mem1, loff);
  for(int i= 0; i < n; ++i) V0()[i]= Z(init(i)).val();
  for(int h= 1; h < H; ++h)
   for(int b= 0; b < cnt[h]; ++b)
    for(int j= 0; j < 16; ++j) LV(h)[b].m[j]= 1, LV(h)[b].a[j]= 0;
 }
public:
 Dual()= default;
 explicit Dual(int n): n(n) {
  build([](int) { return Z(); });
 }
 Dual(int n, Z x): n(n) {
  build([&](int) { return x; });
 }
 Dual(const std::vector<Z>& a): n(a.size()) {
  build([&](int i) { return a[i]; });
 }
 template <class F, std::enable_if_t<std::is_invocable_r_v<Z, const F&, int>, int> = 0> Dual(int n, const F& init): n(n) { build(init); }
 int size() const { return n; }
 Z get(int i) const {
  u64 v= V0()[i];
  for(int h= 1; h < H; ++h) {
   const int c= i >> (4 * h);
   v= (u64(LV(h)[c >> 4].m[c & 15]) * v + LV(h)[c >> 4].a[c & 15]) % P;
  }
  return Z::raw(u32(v));
 }
 Z operator[](int i) const { return get(i); }
 void set(int i, Z x) {
  for(int h= H - 1; h >= 1; --h) push(h, i >> (4 * h));
  put_lane(V0() + ((i >> 4) << 4), i & 15, x.val());
 }
 [[gnu::always_inline]] void apply(int l, int r, const E& f) {
  if(l >= r) return;
  prefetch(l, r);
  const u32 m= f[0].val(), a= f[1].val();
  const int x= r - 1;
  for(int h= H - 1; h >= 1; --h) {
   const int sh= 4 * h, sm= (1 << sh) - 1;
   if(l & sm) push(h, l >> sh);
   if(r & sm) push(h, x >> sh);
  }
  const Mul mf(m);
  for(int h= 0, p= l, q= r; p < q; ++h) {
   const int bp= p >> 4, bq= q >> 4;
   if(bp == bq || h == H - 1) {
    apply_lanes<false>(h, bp, p & 15, q - (bp << 4), mf, a);
    break;
   }
   if(p & 15) apply_lanes<false>(h, bp, p & 15, 16, mf, a), p= bp + 1;
   else p= bp;
   if(q & 15) apply_lanes<false>(h, bq, 0, q & 15, mf, a), q= bq;
   else q= bq;
  }
 }
};
template <class Z> struct const_zmod: std::false_type {};
template <unsigned M> struct const_zmod<ZMod<M>>: std::bool_constant<M != 0> {};
template <class Z, class= void> struct Pick {
 using lazy= segment_tree_internal::Generic<monoid::AffineSum<Z>>;
 using dual= segment_tree_internal::Generic<monoid::Affine<Z>>;
};
template <class Z> struct Pick<Z, std::enable_if_t<const_zmod<Z>::value>> {
 using lazy= Lazy<Z, Z::mod()>;
 using dual= Dual<Z, Z::mod()>;
};
}
template <class Z> class SegmentTree<monoid::AffineSum<Z>>: public segtree_affine_internal::Pick<Z>::lazy {
 using Base= typename segtree_affine_internal::Pick<Z>::lazy;
public:
 using Base::Base;
};
template <class Z> class SegmentTree<monoid::Affine<Z>>: public segtree_affine_internal::Pick<Z>::dual {
 using Base= typename segtree_affine_internal::Pick<Z>::dual;
public:
 using Base::Base;
};
