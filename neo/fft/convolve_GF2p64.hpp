#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <memory>
#include <type_traits>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
#include "neo/algebra/GF2p64.hpp"
namespace gf2p64_internal::cantor {
constexpr int BB= 14;  // block の大きさ 2^BB
inline __m128i clm(u64 a, u64 b) { return _mm_clmulepi64_si128(_mm_cvtsi64_si128((long long)a), _mm_cvtsi64_si128((long long)b), 0); }
inline u64 red(const __m128i& v) {
 static constexpr u8 RED[]= {0, 27, 45, 54, 90, 65, 119, 108};
 const u64 h= v[1], d= h ^ (h << 1);
 return u64(v[0]) ^ RED[h >> 60] ^ d ^ (d << 3);
}
// Cantor 基底。β_l = CHAIN[62-l]、β_l^2 + β_l = β_{l-1}、β_0 = 1。
constexpr int CHAIN_LEN= 63;
constexpr std::array<u64, CHAIN_LEN> CHAIN= []() {
 std::array<u64, CHAIN_LEN> c{};
 c[0]= 2;
 for(int k= 1; k < CHAIN_LEN; ++k) c[k]= cmul(c[k - 1], c[k - 1]) ^ c[k - 1];
 return c;
}();
inline int msb(u64 n) { return 63 - __builtin_clzll(n); }
inline int clog2(u64 n) { return n <= 1 ? 0 : msb(n - 1) + 1; }
// DIF の twiddle。M[j] = Σ_{j の bit L が立つ} β_{L+1} は j の bit について線型なので、M[j] = lo[j & 1023] ^ hi[j >> 10] で引く。
struct Twiddle {
 alignas(64) u64 lo[1024], hi[2048];
};
inline constexpr Twiddle TW= []() {
 Twiddle t{};
 for(int j= 1; j < 1024; ++j) t.lo[j]= t.lo[j & (j - 1)] ^ CHAIN[61 - __builtin_ctz(j)];
 for(int h= 1; h < 2048; ++h) t.hi[h]= t.hi[h & (h - 1)] ^ CHAIN[51 - __builtin_ctz(h)];
 return t;
}();
inline __m256i ld(const u64* p) { return _mm256_loadu_si256((const __m256i*)p); }
inline void st(u64* p, const __m256i& v) { _mm256_storeu_si256((__m256i*)p, v); }
inline __m256i bc(u64 s) { return _mm256_set1_epi64x((long long)s); }
inline __m256i x(const __m256i& a, const __m256i& b) { return _mm256_xor_si256(a, b); }
// M[j] を引く口。ld(j) は M[j..j+3] (j は 4 の倍数なので上位の表は 4 つで同じ)。
struct TwView {
 const u64 *lo, *hi;
 u64 operator[](size_t j) const { return lo[j & 1023] ^ hi[j >> 10]; }
 __m256i ld(size_t j) const { return x(_mm256_load_si256((const __m256i*)(lo + (j & 1023))), _mm256_set1_epi64x((long long)hi[j >> 10])); }
};
inline TwView tw_view() { return {TW.lo, TW.hi}; }
// 4 つの積を畳む前の形 (下位 64 bit の 4 つと上位 64 bit の 4 つ) で持つ。
struct U4 {
 __m256i lo, hi;
};
template <bool V> inline U4 clm4(const __m256i& a, const __m256i& b) {
 __m256i p0, p1;
 if constexpr(V) p0= _mm256_clmulepi64_epi128(a, b, 0x00), p1= _mm256_clmulepi64_epi128(a, b, 0x11);
 else {
  const __m128i al= _mm256_castsi256_si128(a), ah= _mm256_extracti128_si256(a, 1), bl= _mm256_castsi256_si128(b), bh= _mm256_extracti128_si256(b, 1);
  p0= _mm256_setr_m128i(_mm_clmulepi64_si128(al, bl, 0x00), _mm_clmulepi64_si128(ah, bh, 0x00));
  p1= _mm256_setr_m128i(_mm_clmulepi64_si128(al, bl, 0x11), _mm_clmulepi64_si128(ah, bh, 0x11));
 }
 return {_mm256_unpacklo_epi64(p0, p1), _mm256_unpackhi_epi64(p0, p1)};
}
inline __m256i red4(const U4& p) {
 const __m256i RED256= _mm256_setr_epi8(0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0, 0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0);
 const __m256i d= x(p.hi, _mm256_slli_epi64(p.hi, 1));
 return x(x(p.lo, _mm256_shuffle_epi8(RED256, _mm256_srli_epi64(p.hi, 60))), x(d, _mm256_slli_epi64(d, 3)));
}
// clmul の結果を unpack する前の形 (p0 は積 0 と 2、p1 は積 1 と 3)。lazy に足し合わせる積は生のまま XOR し、畳む直前に 1 回だけ unpack する。
struct R4 {
 __m256i p0, p1;
};
template <bool V> inline R4 clr4(const __m256i& a, const __m256i& b) {
 if constexpr(V) return {_mm256_clmulepi64_epi128(a, b, 0x00), _mm256_clmulepi64_epi128(a, b, 0x11)};
 else {
  const __m128i al= _mm256_castsi256_si128(a), ah= _mm256_extracti128_si256(a, 1), bl= _mm256_castsi256_si128(b), bh= _mm256_extracti128_si256(b, 1);
  return {_mm256_setr_m128i(_mm_clmulepi64_si128(al, bl, 0x00), _mm_clmulepi64_si128(ah, bh, 0x00)), _mm256_setr_m128i(_mm_clmulepi64_si128(al, bl, 0x11), _mm_clmulepi64_si128(ah, bh, 0x11))};
 }
}
inline R4 operator^(const R4& a, const R4& b) { return {x(a.p0, b.p0), x(a.p1, b.p1)}; }
inline __m256i redr(const R4& r) { return red4({_mm256_unpacklo_epi64(r.p0, r.p1), _mm256_unpackhi_epi64(r.p0, r.p1)}); }
template <bool V> inline __m256i mul4(const __m256i& a, const __m256i& b) { return red4(clm4<V>(a, b)); }
// 小さいときの素朴な積 (畳むのは出力ごとに 1 回)。
template <class T> inline std::vector<T> naive(const std::vector<T>& a, const std::vector<T>& b) {
 const int n= a.size(), m= b.size(), N= n + m - 1;
 std::vector<T> c(N);
 for(int k= 0; k < N; ++k) {
  __m128i t= _mm_setzero_si128();
  for(int i= std::max(0, k - m + 1); i <= std::min(k, n - 1); ++i) t= _mm_xor_si128(t, clm(u64(a[i]), u64(b[k - i])));
  c[k]= T(red(t));
 }
 return c;
}
// level 3, 2 を 8 個ずつ (level 3: i = 3..0 で e[i+1] ^= e[4+i]、level 2: 4 個ずつ e[2] ^= e[3], e[1] ^= e[2])。
inline void lch_low(u64* p, size_t lo, size_t hi) {
 for(size_t i= lo; i < hi; i+= 8) {
  u64* e= p + i;
  u64 e1= e[1], e2= e[2], e3= e[3], e4= e[4], e5= e[5], e6= e[6], e7= e[7];
  e4^= e7, e3^= e6, e2^= e5, e1^= e4;
  e2^= e3, e1^= e2, e6^= e7, e5^= e6;
  e[1]= e1, e[2]= e2, e[3]= e3, e[4]= e4, e[5]= e5, e[6]= e6;
 }
}
inline void mono_low(u64* p, size_t lo, size_t hi) {
 for(size_t i= lo; i < hi; i+= 8) {
  u64* e= p + i;
  u64 e1= e[1], e2= e[2], e3= e[3], e4= e[4], e5= e[5], e6= e[6], e7= e[7];
  e1^= e2, e2^= e3, e5^= e6, e6^= e7;
  e1^= e4, e2^= e5, e3^= e6, e4^= e7;
  e[1]= e1, e[2]= e2, e[3]= e3, e[4]= e4, e[5]= e5, e[6]= e6;
 }
}
// k 段 (2^k 要素、1 要素は e 個の u64) の単項式 → LCH。a は k 未満で最大の 2 の冪。s_{a+b} = s_a ∘ s_b と s_a(x) = x^(2^a) + x から、
// f を y = s_a(x) について Taylor 展開し、y について (要素は 2^a 個ずつの塊) と、各塊の中の x について変換する。出力の並びは自然順のまま。
inline void tay_lch(u64* p, int k, int a, size_t e) {
 const size_t tot= e << k;
 for(int L= k; L > a; --L) {
  const size_t H= e << (L - 1), SH= e << (L - 1 - a);
  for(size_t base= 0; base < tot; base+= 2 * H) {
   u64* b= p + base;
   if(H >= 8)
    for(size_t c= H; c;) c-= 4, st(b + c + SH, x(ld(b + c + SH), ld(b + H + c)));
   else
    for(size_t c= H; c--;) b[c + SH]^= b[H + c];
  }
 }
}
inline void tay_mono(u64* p, int k, int a, size_t e) {
 const size_t tot= e << k;
 for(int L= a + 1; L <= k; ++L) {
  const size_t H= e << (L - 1), SH= e << (L - 1 - a);
  for(size_t base= 0; base < tot; base+= 2 * H) {
   u64* b= p + base;
   if(H >= 8)
    for(size_t c= 0; c < H; c+= 4) st(b + c + SH, x(ld(b + c + SH), ld(b + H + c)));
   else
    for(size_t c= 0; c < H; ++c) b[c + SH]^= b[H + c];
  }
 }
}
// 16 要素の塊 (e = 1、k = 4) をレジスタに載せたまま回す。ずらした位置へ読み書きすると store forwarding が効かないため。
inline void lch16(u64* p) {
 u64 e1= p[1], e2= p[2], e3= p[3], e4= p[4], e5= p[5], e6= p[6], e7= p[7], e8= p[8], e9= p[9], e10= p[10], e11= p[11], e12= p[12], e13= p[13], e14= p[14], e15= p[15];
 e9^= e15, e8^= e14, e7^= e13, e6^= e12, e5^= e11, e4^= e10, e3^= e9, e2^= e8;
 e4^= e7, e3^= e6, e2^= e5, e1^= e4, e12^= e15, e11^= e14, e10^= e13, e9^= e12;
 e11^= e15, e10^= e14, e9^= e13, e8^= e12, e7^= e11, e6^= e10, e5^= e9, e4^= e8;
 e2^= e3, e1^= e2, e6^= e7, e5^= e6, e10^= e11, e9^= e10, e14^= e15, e13^= e14;
 p[1]= e1, p[2]= e2, p[3]= e3, p[4]= e4, p[5]= e5, p[6]= e6, p[7]= e7, p[8]= e8, p[9]= e9, p[10]= e10, p[11]= e11, p[12]= e12, p[13]= e13, p[14]= e14, p[15]= e15;
}
inline void mono16(u64* p) {
 u64 e1= p[1], e2= p[2], e3= p[3], e4= p[4], e5= p[5], e6= p[6], e7= p[7], e8= p[8], e9= p[9], e10= p[10], e11= p[11], e12= p[12], e13= p[13], e14= p[14], e15= p[15];
 e13^= e14, e14^= e15, e9^= e10, e10^= e11, e5^= e6, e6^= e7, e1^= e2, e2^= e3;
 e4^= e8, e5^= e9, e6^= e10, e7^= e11, e8^= e12, e9^= e13, e10^= e14, e11^= e15;
 e9^= e12, e10^= e13, e11^= e14, e12^= e15, e1^= e4, e2^= e5, e3^= e6, e4^= e7;
 e2^= e8, e3^= e9, e4^= e10, e5^= e11, e6^= e12, e7^= e13, e8^= e14, e9^= e15;
 p[1]= e1, p[2]= e2, p[3]= e3, p[4]= e4, p[5]= e5, p[6]= e6, p[7]= e7, p[8]= e8, p[9]= e9, p[10]= e10, p[11]= e11, p[12]= e12, p[13]= e13, p[14]= e14, p[15]= e15;
}
// e = 1 で k <= 4 の塊を [0, n) にまとめて (k = 4 は lch16、k = 3 は 8 個ずつ、k = 2 は 4 個ずつ)。
inline void lch_small(u64* p, size_t n, int k) {
 if(k >= 4)
  for(size_t i= 0; i < n; i+= 16) lch16(p + i);
 else if(k == 3) lch_low(p, 0, n);
 else if(k == 2)
  for(size_t i= 0; i < n; i+= 4) p[i + 2]^= p[i + 3], p[i + 1]^= p[i + 2];
}
inline void mono_small(u64* p, size_t n, int k) {
 if(k >= 4)
  for(size_t i= 0; i < n; i+= 16) mono16(p + i);
 else if(k == 3) mono_low(p, 0, n);
 else if(k == 2)
  for(size_t i= 0; i < n; i+= 4) p[i + 1]^= p[i + 2], p[i + 2]^= p[i + 3];
}
// 行 (e 個ずつ、4 <= e <= 256) が 2^k 個 (k <= 4) の変換を、列を 4 本ずつ取って行をレジスタに載せたまま回す。行は別々の変数で持つ (配列だと gcc がスタックへ退避する)。
inline void lch_rows4(u64* p, size_t e) {
 for(size_t c= 0; c < e; c+= 4) {
  u64* q= p + c;
  __m256i v1= ld(q + e), v2= ld(q + 2 * e), v3= ld(q + 3 * e), v4= ld(q + 4 * e), v5= ld(q + 5 * e), v6= ld(q + 6 * e), v7= ld(q + 7 * e), v8= ld(q + 8 * e), v9= ld(q + 9 * e), v10= ld(q + 10 * e), v11= ld(q + 11 * e), v12= ld(q + 12 * e), v13= ld(q + 13 * e), v14= ld(q + 14 * e), v15= ld(q + 15 * e);
  v9= x(v9, v15), v8= x(v8, v14), v7= x(v7, v13), v6= x(v6, v12), v5= x(v5, v11), v4= x(v4, v10), v3= x(v3, v9), v2= x(v2, v8);
  v4= x(v4, v7), v3= x(v3, v6), v2= x(v2, v5), v1= x(v1, v4), v12= x(v12, v15), v11= x(v11, v14), v10= x(v10, v13), v9= x(v9, v12);
  v11= x(v11, v15), v10= x(v10, v14), v9= x(v9, v13), v8= x(v8, v12), v7= x(v7, v11), v6= x(v6, v10), v5= x(v5, v9), v4= x(v4, v8);
  v2= x(v2, v3), v1= x(v1, v2), v6= x(v6, v7), v5= x(v5, v6), v10= x(v10, v11), v9= x(v9, v10), v14= x(v14, v15), v13= x(v13, v14);
  st(q + e, v1), st(q + 2 * e, v2), st(q + 3 * e, v3), st(q + 4 * e, v4), st(q + 5 * e, v5), st(q + 6 * e, v6), st(q + 7 * e, v7), st(q + 8 * e, v8), st(q + 9 * e, v9), st(q + 10 * e, v10), st(q + 11 * e, v11), st(q + 12 * e, v12), st(q + 13 * e, v13), st(q + 14 * e, v14), st(q + 15 * e, v15);
 }
}
inline void mono_rows4(u64* p, size_t e) {
 for(size_t c= 0; c < e; c+= 4) {
  u64* q= p + c;
  __m256i v1= ld(q + e), v2= ld(q + 2 * e), v3= ld(q + 3 * e), v4= ld(q + 4 * e), v5= ld(q + 5 * e), v6= ld(q + 6 * e), v7= ld(q + 7 * e), v8= ld(q + 8 * e), v9= ld(q + 9 * e), v10= ld(q + 10 * e), v11= ld(q + 11 * e), v12= ld(q + 12 * e), v13= ld(q + 13 * e), v14= ld(q + 14 * e), v15= ld(q + 15 * e);
  v13= x(v13, v14), v14= x(v14, v15), v9= x(v9, v10), v10= x(v10, v11), v5= x(v5, v6), v6= x(v6, v7), v1= x(v1, v2), v2= x(v2, v3);
  v4= x(v4, v8), v5= x(v5, v9), v6= x(v6, v10), v7= x(v7, v11), v8= x(v8, v12), v9= x(v9, v13), v10= x(v10, v14), v11= x(v11, v15);
  v9= x(v9, v12), v10= x(v10, v13), v11= x(v11, v14), v12= x(v12, v15), v1= x(v1, v4), v2= x(v2, v5), v3= x(v3, v6), v4= x(v4, v7);
  v2= x(v2, v8), v3= x(v3, v9), v4= x(v4, v10), v5= x(v5, v11), v6= x(v6, v12), v7= x(v7, v13), v8= x(v8, v14), v9= x(v9, v15);
  st(q + e, v1), st(q + 2 * e, v2), st(q + 3 * e, v3), st(q + 4 * e, v4), st(q + 5 * e, v5), st(q + 6 * e, v6), st(q + 7 * e, v7), st(q + 8 * e, v8), st(q + 9 * e, v9), st(q + 10 * e, v10), st(q + 11 * e, v11), st(q + 12 * e, v12), st(q + 13 * e, v13), st(q + 14 * e, v14), st(q + 15 * e, v15);
 }
}
inline void lch_rows3(u64* p, size_t e) {
 for(size_t c= 0; c < e; c+= 4) {
  u64* q= p + c;
  __m256i v1= ld(q + e), v2= ld(q + 2 * e), v3= ld(q + 3 * e), v4= ld(q + 4 * e), v5= ld(q + 5 * e), v6= ld(q + 6 * e), v7= ld(q + 7 * e);
  v4= x(v4, v7), v3= x(v3, v6), v2= x(v2, v5), v1= x(v1, v4), v2= x(v2, v3), v1= x(v1, v2), v6= x(v6, v7), v5= x(v5, v6);
  st(q + e, v1), st(q + 2 * e, v2), st(q + 3 * e, v3), st(q + 4 * e, v4), st(q + 5 * e, v5), st(q + 6 * e, v6), st(q + 7 * e, v7);
 }
}
inline void mono_rows3(u64* p, size_t e) {
 for(size_t c= 0; c < e; c+= 4) {
  u64* q= p + c;
  __m256i v1= ld(q + e), v2= ld(q + 2 * e), v3= ld(q + 3 * e), v4= ld(q + 4 * e), v5= ld(q + 5 * e), v6= ld(q + 6 * e), v7= ld(q + 7 * e);
  v5= x(v5, v6), v6= x(v6, v7), v1= x(v1, v2), v2= x(v2, v3), v1= x(v1, v4), v2= x(v2, v5), v3= x(v3, v6), v4= x(v4, v7);
  st(q + e, v1), st(q + 2 * e, v2), st(q + 3 * e, v3), st(q + 4 * e, v4), st(q + 5 * e, v5), st(q + 6 * e, v6), st(q + 7 * e, v7);
 }
}
inline void lch_rows2(u64* p, size_t e) {
 for(size_t c= 0; c < e; c+= 4) {
  u64* q= p + c;
  __m256i v1= ld(q + e), v2= ld(q + 2 * e), v3= ld(q + 3 * e);
  v2= x(v2, v3), v1= x(v1, v2);
  st(q + e, v1), st(q + 2 * e, v2), st(q + 3 * e, v3);
 }
}
inline void mono_rows2(u64* p, size_t e) {
 for(size_t c= 0; c < e; c+= 4) {
  u64* q= p + c;
  __m256i v1= ld(q + e), v2= ld(q + 2 * e), v3= ld(q + 3 * e);
  v1= x(v1, v2), v2= x(v2, v3);
  st(q + e, v1), st(q + 2 * e, v2), st(q + 3 * e, v3);
 }
}
inline void lch_rows(u64* p, int k, size_t e) {
 if(k == 4) lch_rows4(p, e);
 else if(k == 3) lch_rows3(p, e);
 else if(k == 2) lch_rows2(p, e);
}
inline void mono_rows(u64* p, int k, size_t e) {
 if(k == 4) mono_rows4(p, e);
 else if(k == 3) mono_rows3(p, e);
 else if(k == 2) mono_rows2(p, e);
}
inline void lch_rec(u64* p, int k, size_t e) {
 if(k <= 1) return;
 if(e >= 4 && e <= 256 && k <= 4) return lch_rows(p, k, e);  // 行の間隔が 2 KB 以下
 if(e == 1 && k <= 4) return lch_small(p, size_t(1) << k, k);
 const int a= 1 << (31 - __builtin_clz(k - 1));
 tay_lch(p, k, a, e);
 lch_rec(p, k - a, e << a);
 if(e == 1 && a <= 4) lch_small(p, size_t(1) << k, a);
 else
  for(size_t b= 0, nb= size_t(1) << (k - a); b < nb; ++b) lch_rec(p + (b << a) * e, a, e);
}
inline void mono_rec(u64* p, int k, size_t e) {
 if(k <= 1) return;
 if(e >= 4 && e <= 256 && k <= 4) return mono_rows(p, k, e);
 if(e == 1 && k <= 4) return mono_small(p, size_t(1) << k, k);
 const int a= 1 << (31 - __builtin_clz(k - 1));
 if(e == 1 && a <= 4) mono_small(p, size_t(1) << k, a);
 else
  for(size_t b= 0, nb= size_t(1) << (k - a); b < nb; ++b) mono_rec(p + (b << a) * e, a, e);
 mono_rec(p, k - a, e << a);
 tay_mono(p, k, a, e);
}
// DIF の radix-4 (段 k と k-1)。block 0 は v = u = 0 なので掛け算が 1 つ。順変換は読む場所 s と書く場所 b を分けられる (上の段を飛ばすと、どの block も base から読む)。
template <bool V> inline void r4f_first(const u64* s, u64* b, int q, u64 uh) {
 const __m256i Uh= bc(uh);
 for(int i= 0; i < q; i+= 4) {
  const __m256i Q0= ld(s + i), Q1= ld(s + q + i), Q2= ld(s + 2 * q + i), Q3= ld(s + 3 * q + i);
  const __m256i A3= x(Q1, Q3), R2= x(x(Q0, Q2), mul4<V>(A3, Uh));
  st(b + i, Q0), st(b + q + i, x(Q0, Q1)), st(b + 2 * q + i, R2), st(b + 3 * q + i, x(R2, A3));
 }
}
template <bool V> inline void r4f_lazy(const u64* s, u64* b, int q, u64 v, u64 ul, u64 uh) {
 const __m256i Vv= bc(v), Ul= bc(ul), Uh= bc(uh);
 for(int i= 0; i < q; i+= 4) {
  const __m256i Q0= ld(s + i), Q1= ld(s + q + i), Q2= ld(s + 2 * q + i), Q3= ld(s + 3 * q + i);
  const R4 P2= clr4<V>(Q2, Vv);
  const __m256i A1= x(Q1, redr(clr4<V>(Q3, Vv))), A3= x(A1, Q3);
  const __m256i R0= x(Q0, redr(P2 ^ clr4<V>(A1, Ul))), R2= x(x(Q0, Q2), redr(P2 ^ clr4<V>(A3, Uh)));
  st(b + i, R0), st(b + q + i, x(R0, A1)), st(b + 2 * q + i, R2), st(b + 3 * q + i, x(R2, A3));
 }
}
template <bool V> inline void r4i_first(u64* b, int q, u64 uh) {
 const __m256i Uh= bc(uh);
 for(int i= 0; i < q; i+= 4) {
  const __m256i R0= ld(b + i), R1= ld(b + q + i), R2= ld(b + 2 * q + i), R3= ld(b + 3 * q + i);
  const __m256i A3= x(R2, R3), A1= x(R0, R1), Q2= x(x(R0, R2), mul4<V>(A3, Uh));
  st(b + q + i, A1), st(b + 2 * q + i, Q2), st(b + 3 * q + i, x(A1, A3));
 }
}
template <bool V> inline void r4i_lazy(u64* b, int q, u64 v, u64 ul, u64 uh) {
 const __m256i Vv= bc(v), Ul= bc(ul), Uh= bc(uh);
 for(int i= 0; i < q; i+= 4) {
  const __m256i R0= ld(b + i), R1= ld(b + q + i), R2= ld(b + 2 * q + i), R3= ld(b + 3 * q + i), A3= x(R2, R3), A1= x(R0, R1);
  const R4 P1= clr4<V>(A1, Ul);
  const __m256i Q2= x(x(R0, R2), redr(P1 ^ clr4<V>(A3, Uh))), Q3= x(A1, A3);
  st(b + i, x(R0, redr(P1 ^ clr4<V>(Q2, Vv)))), st(b + q + i, x(A1, redr(clr4<V>(Q3, Vv)))), st(b + 2 * q + i, Q2), st(b + 3 * q + i, Q3);
 }
}
// radix-2 (段 k、h = 2^(k-1))。
template <bool V, bool FIRST> inline void r2f(const u64* s, u64* b, int h, u64 w) {
 const __m256i W= bc(w);
 for(int i= 0; i < h; i+= 4) {
  __m256i u= ld(s + i);
  const __m256i v= ld(s + h + i);
  if constexpr(!FIRST) u= x(u, mul4<V>(v, W));
  st(b + i, u), st(b + h + i, x(u, v));
 }
}
template <bool V, bool FIRST> inline void r2i(u64* b, int h, u64 w) {
 const __m256i W= bc(w);
 for(int i= 0; i < h; i+= 4) {
  const __m256i u= ld(b + i), v= x(u, ld(b + h + i));
  if constexpr(!FIRST) st(b + i, x(u, mul4<V>(v, W)));
  st(b + h + i, v);
 }
}
// 4x4 の転置 (行 r0..r3 の t 番目を集めて列 t に)。自分自身が逆。
inline void tr4(__m256i& r0, __m256i& r1, __m256i& r2, __m256i& r3) {
 const __m256i t0= _mm256_unpacklo_epi64(r0, r1), t1= _mm256_unpackhi_epi64(r0, r1), t2= _mm256_unpacklo_epi64(r2, r3), t3= _mm256_unpackhi_epi64(r2, r3);
 r0= _mm256_permute2x128_si256(t0, t2, 0x20), r1= _mm256_permute2x128_si256(t1, t3, 0x20), r2= _mm256_permute2x128_si256(t0, t2, 0x31), r3= _mm256_permute2x128_si256(t1, t3, 0x31);
}
// 最下段の radix-4 (段 2, 1、q = 1) を block 4 つずつ。lane l が block j+l で、twiddle も lane ごと (M[0] = 0 なので block 0 も同じ式でよい)。
template <bool V> inline void r4l1f(u64* f, size_t j0, size_t j1) {
 const TwView M= tw_view();
 for(size_t j= j0; j < j1; j+= 4) {
  u64* b= f + 4 * j;
  __m256i Q0= ld(b), Q1= ld(b + 4), Q2= ld(b + 8), Q3= ld(b + 12);
  tr4(Q0, Q1, Q2, Q3);
  const __m256i Vv= M.ld(j), m0= M.ld(2 * j), m1= M.ld(2 * j + 4);
  const __m256i Ul= _mm256_permute4x64_epi64(_mm256_unpacklo_epi64(m0, m1), 0xD8), Uh= _mm256_permute4x64_epi64(_mm256_unpackhi_epi64(m0, m1), 0xD8);
  const R4 P2= clr4<V>(Q2, Vv);
  const __m256i A1= x(Q1, redr(clr4<V>(Q3, Vv))), A3= x(A1, Q3);
  __m256i R0= x(Q0, redr(P2 ^ clr4<V>(A1, Ul))), R2= x(x(Q0, Q2), redr(P2 ^ clr4<V>(A3, Uh))), R1= x(R0, A1), R3= x(R2, A3);
  tr4(R0, R1, R2, R3);
  st(b, R0), st(b + 4, R1), st(b + 8, R2), st(b + 12, R3);
 }
}
template <bool V> inline void r4l1i(u64* f, size_t j0, size_t j1) {
 const TwView M= tw_view();
 for(size_t j= j0; j < j1; j+= 4) {
  u64* b= f + 4 * j;
  __m256i R0= ld(b), R1= ld(b + 4), R2= ld(b + 8), R3= ld(b + 12);
  tr4(R0, R1, R2, R3);
  const __m256i Vv= M.ld(j), m0= M.ld(2 * j), m1= M.ld(2 * j + 4);
  const __m256i Ul= _mm256_permute4x64_epi64(_mm256_unpacklo_epi64(m0, m1), 0xD8), Uh= _mm256_permute4x64_epi64(_mm256_unpackhi_epi64(m0, m1), 0xD8);
  const __m256i A3= x(R2, R3), A1= x(R0, R1);
  const R4 P1= clr4<V>(A1, Ul);
  __m256i Q2= x(x(R0, R2), redr(P1 ^ clr4<V>(A3, Uh))), Q3= x(A1, A3), Q0= x(R0, redr(P1 ^ clr4<V>(Q2, Vv))), Q1= x(A1, redr(clr4<V>(Q3, Vv)));
  tr4(Q0, Q1, Q2, Q3);
  st(b, Q0), st(b + 4, Q1), st(b + 8, Q2), st(b + 12, Q3);
 }
}
struct Pass {
 int k;
 bool r4;  // true: 段 k と k-1、false: 段 k だけ
};
// 段 hi から lo までを上から。段の数が奇数なら一番上を radix-2 にする (最下段は必ず radix-4 の q = 1 になる)。
inline int plan(int lo, int hi, Pass* ps) {
 int n= 0;
 if((hi - lo + 1) & 1) ps[n++]= {hi, false}, --hi;
 for(; hi > lo; hi-= 2) ps[n++]= {hi, true};
 return n;
}
// [lo, hi) の block に 1 pass。fused なら、どの block も f の先頭 (上の段を飛ばした base) から読む。base を最後に上書きするよう後ろから。
template <bool V> inline void fwd_pass(u64* f, size_t lo, size_t hi, const Pass& p, bool fused) {
 const TwView M= tw_view();
 const int k= p.k;
 const size_t j0= lo >> k, j1= hi >> k;
 if(p.r4) {
  const int q= 1 << (k - 2);
  if(q == 1) return r4l1f<V>(f, j0, j1);
  for(size_t j= j1; j-- > j0;) {
   u64* b= f + (j << k);
   const u64* s= fused ? f : b;
   if(!j) r4f_first<V>(s, b, q, M[1]);
   else r4f_lazy<V>(s, b, q, M[j], M[2 * j], M[2 * j + 1]);
  }
 } else {
  const int h= 1 << (k - 1);
  for(size_t j= j1; j-- > j0;) {
   u64* b= f + (j << k);
   const u64* s= fused ? f : b;
   if(!j) r2f<V, true>(s, b, h, 0);
   else r2f<V, false>(s, b, h, M[j]);
  }
 }
}
template <bool V> inline void inv_pass(u64* f, size_t lo, size_t hi, const Pass& p) {
 const TwView M= tw_view();
 const int k= p.k;
 const size_t j0= lo >> k, j1= hi >> k;
 if(p.r4) {
  const int q= 1 << (k - 2);
  if(q == 1) return r4l1i<V>(f, j0, j1);
  for(size_t j= j0; j < j1; ++j) {
   u64* b= f + (j << k);
   if(!j) r4i_first<V>(b, q, M[1]);
   else r4i_lazy<V>(b, q, M[j], M[2 * j], M[2 * j + 1]);
  }
 } else {
  const int h= 1 << (k - 1);
  for(size_t j= j0; j < j1; ++j) {
   u64* b= f + (j << k);
   if(!j) r2i<V, true>(b, h, 0);
   else r2i<V, false>(b, h, M[j]);
  }
 }
}
// 順変換の上の段 (2^bb を超える block の段) を全体に。入力が [0, 2^k0) に収まるので段 k0 から始める。
template <bool V> inline void fwd_top(u64* f, int d, int k0, int bb) {
 if(k0 <= bb) return;
 Pass ps[64];
 const int np= plan(bb + 1, k0, ps);
 for(int i= 0; i < np; ++i) fwd_pass<V>(f, 0, size_t(1) << d, ps[i], i == 0 && k0 < d);
}
// 順変換の下の段を block [lo, hi) に。上の段が無く (k0 <= bb) 飛ばした段があるなら、最初の pass は base から読む。
template <bool V> inline void fwd_bottom(u64* f, size_t lo, size_t hi, int d, int k0, int bb) {
 Pass ps[64];
 const int np= plan(1, std::min(bb, k0), ps);
 for(int i= 0; i < np; ++i) fwd_pass<V>(f, lo, hi, ps[i], i == 0 && k0 <= bb && k0 < d);
}
template <bool V> inline void inv_bottom(u64* f, size_t lo, size_t hi, int bb) {
 Pass ps[64];
 const int np= plan(1, bb, ps);
 for(int i= np; i--;) inv_pass<V>(f, lo, hi, ps[i]);
}
template <bool V> inline void inv_top(u64* f, int d, int bb) {
 if(d <= bb) return;
 Pass ps[64];
 const int np= plan(bb + 1, d, ps);
 for(int i= np; i--;) inv_pass<V>(f, 0, size_t(1) << d, ps[i]);
}
// f と g の最下段の pass と各点の積と逆変換の最下段の pass は同じ block 4 つと twiddle を使うので、転置したまま続けて回し、積は畳まずに逆変換の式に入れる。
template <bool V> inline void r4c_fwd(const __m256i& Q0, const __m256i& Q1, const __m256i& Q2, const __m256i& Q3, const __m256i& Vv, const __m256i& Ul, const __m256i& Uh, __m256i* R) {
 const R4 P2= clr4<V>(Q2, Vv);
 const __m256i A1= x(Q1, redr(clr4<V>(Q3, Vv))), A3= x(A1, Q3);
 R[0]= x(Q0, redr(P2 ^ clr4<V>(A1, Ul))), R[2]= x(x(Q0, Q2), redr(P2 ^ clr4<V>(A3, Uh))), R[1]= x(R[0], A1), R[3]= x(R[2], A3);
}
// radix-4 (段 2, 1)。畳むのは f, g で 3 回ずつと、積と逆変換で 5 回 (分けると 3 + 3 + 4 + 3 回)。
template <bool V> inline void fused_r4l1(u64* F, const u64* G, size_t j0, size_t j1) {
 const TwView M= tw_view();
 for(size_t j= j0; j < j1; j+= 4) {
  u64* b= F + 4 * j;
  const u64* c= G + 4 * j;
  const __m256i Vv= M.ld(j), m0= M.ld(2 * j), m1= M.ld(2 * j + 4);
  const __m256i Ul= _mm256_permute4x64_epi64(_mm256_unpacklo_epi64(m0, m1), 0xD8), Uh= _mm256_permute4x64_epi64(_mm256_unpackhi_epi64(m0, m1), 0xD8);
  __m256i f0= ld(b), f1= ld(b + 4), f2= ld(b + 8), f3= ld(b + 12), g0= ld(c), g1= ld(c + 4), g2= ld(c + 8), g3= ld(c + 12), Fv[4], Gv[4];
  tr4(f0, f1, f2, f3), tr4(g0, g1, g2, g3);
  r4c_fwd<V>(f0, f1, f2, f3, Vv, Ul, Uh, Fv), r4c_fwd<V>(g0, g1, g2, g3, Vv, Ul, Uh, Gv);
  const R4 P0= clr4<V>(Fv[0], Gv[0]), P1= clr4<V>(Fv[1], Gv[1]), P2= clr4<V>(Fv[2], Gv[2]), P3= clr4<V>(Fv[3], Gv[3]);
  // r4l1i の式で R_t = 畳み(P_t) としたもの。
  const __m256i A1= redr(P0 ^ P1), A3= redr(P2 ^ P3);
  const R4 X1= clr4<V>(A1, Ul);
  __m256i Q2= redr(P0 ^ P2 ^ X1 ^ clr4<V>(A3, Uh)), Q3= x(A1, A3), Q0= redr(P0 ^ X1 ^ clr4<V>(Q2, Vv)), Q1= x(A1, redr(clr4<V>(Q3, Vv)));
  tr4(Q0, Q1, Q2, Q3);
  st(b, Q0), st(b + 4, Q1), st(b + 8, Q2), st(b + 12, Q3);
 }
}
inline bool same_pass(const Pass& a, const Pass& b) { return a.k == b.k && a.r4 == b.r4; }
inline bool fusable(const Pass& p) { return p.r4 && p.k == 2; }
template <bool V> inline void fused_bottom(u64* F, const u64* G, size_t lo, size_t hi, const Pass&) { fused_r4l1<V>(F, G, lo >> 2, hi >> 2); }
// 順変換の下の段のうち、最後の pass (fused_bottom で回す) を除いたもの。
template <bool V> inline void fwd_bottom_head(u64* f, size_t lo, size_t hi, int d, int k0, int bb) {
 Pass ps[64];
 const int np= plan(1, std::min(bb, k0), ps);
 for(int i= 0; i + 1 < np; ++i) fwd_pass<V>(f, lo, hi, ps[i], i == 0 && k0 <= bb && k0 < d);
}
// 逆変換の下の段のうち、最初の pass (fused_bottom で回す) を除いたもの。
template <bool V> inline void inv_bottom_tail(u64* f, size_t lo, size_t hi, int bb) {
 Pass ps[64];
 const int np= plan(1, bb, ps);
 for(int i= np - 1; i--;) inv_pass<V>(f, lo, hi, ps[i]);
}
template <bool V> inline void pointwise(u64* f, const u64* g, size_t lo, size_t hi) {
 for(size_t i= lo; i < hi; i+= 4) st(f + i, mul4<V>(ld(f + i), ld(g + i)));
}
// 確保したばかりの領域を、触る前に MADV_POPULATE_WRITE (Linux 5.14+) でまとめて用意させる。古い kernel と Linux 以外では何もしない。
inline void advise(void* p, size_t bytes) {
#ifdef __linux__
 const std::uintptr_t a= (reinterpret_cast<std::uintptr_t>(p) + 4095) & ~std::uintptr_t(4095), e= (reinterpret_cast<std::uintptr_t>(p) + bytes) & ~std::uintptr_t(4095);
 if(e > a) madvise(reinterpret_cast<void*>(a), e - a, 23);  // MADV_POPULATE_WRITE
#else
 (void)p, (void)bytes;
#endif
}
// 大きい確保は glibc では mmap のページの先頭 + 16 byte が返り、32 byte の読み書きの半分が cache line をまたぐので、64 byte 境界から使う。
inline u64* align64(u64* p) { return (u64*)((reinterpret_cast<std::uintptr_t>(p) + 63) & ~std::uintptr_t(63)); }
// T は u64 か GF2p64 (8 byte で trivially copyable)。中では u64 の配列として扱う。
template <bool V, class T> inline std::vector<T> conv(const std::vector<T>& a, const std::vector<T>& b) {
 static_assert(sizeof(T) == 8 && std::is_trivially_copyable_v<T>);
 const size_t n= a.size(), m= b.size();
 if(!n || !m) return {};
 const size_t N= n + m - 1;
 if(N <= 64) return naive(a, b);
 const int d= clog2(N), bb= std::min(BB, d), kf= std::max(clog2(n), 4), kg= std::max(clog2(m), 4);
 const size_t S= size_t(1) << d, BS= size_t(1) << bb;
 std::vector<T> f;  // 出力を兼ねる。[n, 2^kf) の 0 はここで入り、最後に先頭へ詰める
 f.reserve(S + 8), advise(f.data(), (S + 8) * sizeof(u64)), f.resize(S + 8);
 std::unique_ptr<u64[]> gbuf(new u64[S + 8]);  // 0 で埋めない (読む前に書く)
 advise(gbuf.get(), (S + 8) * sizeof(u64));
 u64 *f0= reinterpret_cast<u64*>(f.data()), *F= align64(f0), *G= align64(gbuf.get());
 std::memcpy(F, a.data(), n * sizeof(u64));
 std::memcpy(G, b.data(), m * sizeof(u64)), std::fill(G + m, G + (size_t(1) << kg), 0);
 lch_rec(G, kg, 1), fwd_top<V>(G, d, kg, bb), lch_rec(F, kf, 1), fwd_top<V>(F, d, kf, bb);
 // f, g, 逆変換の最下段の pass がそろっていれば、block ごとに g と f の下の段を進めて最下段をつなぐ (g の block も block 0 を最後に)。
 Pass pf[64], pg[64], pi[64];
 const int nf= plan(1, std::min(bb, kf), pf), ng= plan(1, std::min(bb, kg), pg), ni= plan(1, bb, pi);
 const Pass& last= pi[ni - 1];
 if(nf >= 2 && ng >= 2 && fusable(last) && same_pass(pf[nf - 1], last) && same_pass(pg[ng - 1], last)) {
  for(size_t s= S; s;) s-= BS, fwd_bottom_head<V>(G, s, s + BS, d, kg, bb), fwd_bottom_head<V>(F, s, s + BS, d, kf, bb), fused_bottom<V>(F, G, s, s + BS, last), inv_bottom_tail<V>(F, s, s + BS, bb);
 } else {
  for(size_t s= S; s;) s-= BS, fwd_bottom<V>(G, s, s + BS, d, kg, bb);
  for(size_t s= S; s;) s-= BS, fwd_bottom<V>(F, s, s + BS, d, kf, bb), pointwise<V>(F, G, s, s + BS), inv_bottom<V>(F, s, s + BS, bb);
 }
 inv_top<V>(F, d, bb), mono_rec(F, d, 1);
 if(F != f0) std::memmove(f0, F, N * sizeof(u64));
 f.resize(N);
 return f;
}
template <class T> inline std::vector<T> convolve(const std::vector<T>& a, const std::vector<T>& b) {
#ifdef __x86_64__
 if(__builtin_cpu_supports("vpclmulqdq")) return conv<1>(a, b);
#endif
 return conv<0>(a, b);
}
}
inline std::vector<GF2p64> convolve(const std::vector<GF2p64>& a, const std::vector<GF2p64>& b) { return gf2p64_internal::cantor::convolve(a, b); }
