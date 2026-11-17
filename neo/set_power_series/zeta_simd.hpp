#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <cstddef>
// 集合冪級数のゼータ変換・メビウス変換・アダマール変換の AVX2 の核。sps の zeta.hpp から呼ぶ内部の道具で、直接は使わない。
// 値は 4 byte か 8 byte で、32 byte の塊 (4 byte なら 8 個、8 byte なら 4 個) を 1 本の __m256i に読む。塊の中の bit
// (4 byte なら下の 3 bit、8 byte なら 2 bit) は、レーンをずらした値と混ぜて済ませる。塊をまたぐ bit は、塊の列を半分か 4 分の 1 に
// 分けて再帰してから合わせる (今の Library の汎用の道と同じく、下の段は L1 と L2 の中で済む)。4 つに分けたときは 2 段を 1 回の
// 読み書きで合わせる。bit ごとの 2 × 2 の変換は互いに可換なので、段の順は結果に効かない。
//   向きのある変換 (ゼータとメビウス): f(受け手, 渡し手) を受け手に書く。UP なら bit が 1 の添字が受け手 (部分集合の向き)、
//   そうでなければ bit が 0 の添字が受け手 (上位集合の向き)。
//   蝶 (アダマール): bit が 0 と 1 の組 (x, y) を (add(x, y), sub(x, y)) にする。
namespace sps_zeta {
using V= __m256i;
#define SZ_AI [[gnu::always_inline]] inline
constexpr std::size_t BASE= 64;  // 再帰をやめる塊の数 (2 KB)
SZ_AI V ld(const char* p) { return _mm256_loadu_si256((const V*)p); }
SZ_AI void st(char* p, V x) { _mm256_storeu_si256((V*)p, x); }

// ---- 演算 ----
// 受け手 a に渡し手 b を合わせた値を返す。
struct Add32 {
 SZ_AI V operator()(V a, V b) const { return _mm256_add_epi32(a, b); }
};
struct Add64 {
 SZ_AI V operator()(V a, V b) const { return _mm256_add_epi64(a, b); }
};
struct Sub32 {
 SZ_AI V operator()(V a, V b) const { return _mm256_sub_epi32(a, b); }
};
struct Sub64 {
 SZ_AI V operator()(V a, V b) const { return _mm256_sub_epi64(a, b); }
};
template <bool MAX, bool SIGNED> struct MinMax32 {
 SZ_AI V operator()(V a, V b) const {
  if constexpr(SIGNED) return MAX ? _mm256_max_epi32(a, b) : _mm256_min_epi32(a, b);
  else return MAX ? _mm256_max_epu32(a, b) : _mm256_min_epu32(a, b);
 }
};
// 64 bit の max と min は cmpgt で比べて混ぜる (AVX2 に命令が無い)。符号なしは符号の bit を反転してから比べる。
template <bool MAX, bool SIGNED> struct MinMax64 {
 V sb= _mm256_set1_epi64x((long long)(1ull << 63));
 SZ_AI V operator()(V a, V b) const {
  V x= a, y= b;
  if constexpr(!SIGNED) x= _mm256_xor_si256(a, sb), y= _mm256_xor_si256(b, sb);
  return _mm256_blendv_epi8(a, b, MAX ? _mm256_cmpgt_epi64(y, x) : _mm256_cmpgt_epi64(x, y));
 }
};
// ZMod の値 (2M 未満のそのままの値、M < 2^30)。足し算と引き算の結果も 2M 未満で、ZMod の += と -= と同じ値になる。
struct ZAdd {
 V m2;
 SZ_AI V operator()(V a, V b) const {
  const V s= _mm256_add_epi32(a, b);
  return _mm256_min_epu32(s, _mm256_sub_epi32(s, m2));
 }
};
struct ZSub {
 V m2;
 SZ_AI V operator()(V a, V b) const {
  const V s= _mm256_add_epi32(_mm256_sub_epi32(a, b), m2);
  return _mm256_min_epu32(s, _mm256_sub_epi32(s, m2));
 }
};
// 蝶の add と sub の組。
template <class A, class S> struct Bfly {
 A add;
 S sub;
};

// ---- 塊の中の bit ----
template <int E, bool UP, class F> SZ_AI V dir_lanes(V v, const F& f) {
 if constexpr(E == 4) {
  if constexpr(UP) {
   v= _mm256_blend_epi32(v, f(v, _mm256_slli_epi64(v, 32)), 0xAA);
   v= _mm256_blend_epi32(v, f(v, _mm256_slli_si256(v, 8)), 0xCC);
   v= _mm256_blend_epi32(v, f(v, _mm256_permute2x128_si256(v, v, 0x08)), 0xF0);
  } else {
   v= _mm256_blend_epi32(v, f(v, _mm256_srli_epi64(v, 32)), 0x55);
   v= _mm256_blend_epi32(v, f(v, _mm256_srli_si256(v, 8)), 0x33);
   v= _mm256_blend_epi32(v, f(v, _mm256_permute2x128_si256(v, v, 0x81)), 0x0F);
  }
 } else {
  if constexpr(UP) {
   v= _mm256_blend_epi32(v, f(v, _mm256_slli_si256(v, 8)), 0xCC);
   v= _mm256_blend_epi32(v, f(v, _mm256_permute2x128_si256(v, v, 0x08)), 0xF0);
  } else {
   v= _mm256_blend_epi32(v, f(v, _mm256_srli_si256(v, 8)), 0x33);
   v= _mm256_blend_epi32(v, f(v, _mm256_permute2x128_si256(v, v, 0x81)), 0x0F);
  }
 }
 return v;
}
// 蝶。相手のレーン p を作り、bit が 0 のレーンは add(v, p)、1 のレーンは sub(p, v) にする。
template <int E, class B> SZ_AI V bfly_lanes(V v, const B& b) {
 if constexpr(E == 4) {
  V p= _mm256_shuffle_epi32(v, 0xB1);
  v= _mm256_blend_epi32(b.add(v, p), b.sub(p, v), 0xAA);
  p= _mm256_shuffle_epi32(v, 0x4E);
  v= _mm256_blend_epi32(b.add(v, p), b.sub(p, v), 0xCC);
 } else {
  V p= _mm256_shuffle_epi32(v, 0x4E);
  v= _mm256_blend_epi32(b.add(v, p), b.sub(p, v), 0xCC);
 }
 V p= _mm256_permute2x128_si256(v, v, 0x01);
 return _mm256_blend_epi32(b.add(v, p), b.sub(p, v), 0xF0);
}

// ---- 塊をまたぐ bit ----
// A は L 個の塊 (L は 2 冪)。塊の中の bit も済ませる。演算 (定数を持つものがある) は値で受ける。参照のままだと、__m256i の
// 書き込みのたびに gcc が定数をメモリから読み直す (ranked_xor.hpp の MK と同じ)。
template <int E, bool UP, class F> void dir_rec(char* A, std::size_t L, F f) {
 if(L <= BASE) {
  for(std::size_t i= 0; i < L; ++i) st(A + 32 * i, dir_lanes<E, UP>(ld(A + 32 * i), f));
  for(std::size_t k= 1; k < L; k<<= 1)
   for(std::size_t i= 0; i < L; i+= 2 * k)
    for(std::size_t j= i; j < i + k; ++j) {
     char *x= A + 32 * j, *y= A + 32 * (j + k);
     if constexpr(UP) st(y, f(ld(y), ld(x)));
     else st(x, f(ld(x), ld(y)));
    }
  return;
 }
 if(L >= 4 * BASE) {
  const std::size_t q= L >> 2;
  for(int t= 0; t < 4; ++t) dir_rec<E, UP>(A + 32 * q * t, q, f);
  for(std::size_t j= 0; j < q; ++j) {
   char *p0= A + 32 * j, *p1= p0 + 32 * q, *p2= p1 + 32 * q, *p3= p2 + 32 * q;
   V a0= ld(p0), a1= ld(p1), a2= ld(p2), a3= ld(p3);
   if constexpr(UP) {
    a1= f(a1, a0), a3= f(a3, a2);
    a2= f(a2, a0), a3= f(a3, a1);
    st(p1, a1), st(p2, a2), st(p3, a3);
   } else {
    a0= f(a0, a1), a2= f(a2, a3);
    a0= f(a0, a2), a1= f(a1, a3);
    st(p0, a0), st(p1, a1), st(p2, a2);
   }
  }
 } else {
  const std::size_t h= L >> 1;
  dir_rec<E, UP>(A, h, f), dir_rec<E, UP>(A + 32 * h, h, f);
  for(std::size_t j= 0; j < h; ++j) {
   char *x= A + 32 * j, *y= x + 32 * h;
   if constexpr(UP) st(y, f(ld(y), ld(x)));
   else st(x, f(ld(x), ld(y)));
  }
 }
}
template <int E, class B> void bfly_rec(char* A, std::size_t L, B b) {
 auto two= [&](char* x, char* y) __attribute__((always_inline)) {
  const V u= ld(x), w= ld(y);
  st(x, b.add(u, w)), st(y, b.sub(u, w));
 };
 if(L <= BASE) {
  for(std::size_t i= 0; i < L; ++i) st(A + 32 * i, bfly_lanes<E>(ld(A + 32 * i), b));
  for(std::size_t k= 1; k < L; k<<= 1)
   for(std::size_t i= 0; i < L; i+= 2 * k)
    for(std::size_t j= i; j < i + k; ++j) two(A + 32 * j, A + 32 * (j + k));
  return;
 }
 if(L >= 4 * BASE) {
  const std::size_t q= L >> 2;
  for(int t= 0; t < 4; ++t) bfly_rec<E>(A + 32 * q * t, q, b);
  for(std::size_t j= 0; j < q; ++j) {
   char *p0= A + 32 * j, *p1= p0 + 32 * q, *p2= p1 + 32 * q, *p3= p2 + 32 * q;
   const V a0= ld(p0), a1= ld(p1), a2= ld(p2), a3= ld(p3);
   const V b0= b.add(a0, a1), b1= b.sub(a0, a1), b2= b.add(a2, a3), b3= b.sub(a2, a3);
   st(p0, b.add(b0, b2)), st(p2, b.sub(b0, b2)), st(p1, b.add(b1, b3)), st(p3, b.sub(b1, b3));
  }
 } else {
  const std::size_t h= L >> 1;
  bfly_rec<E>(A, h, b), bfly_rec<E>(A + 32 * h, h, b);
  for(std::size_t j= 0; j < h; ++j) two(A + 32 * j, A + 32 * (j + h));
 }
}
// k 進の変換の 1 つの列。d[r] = f(d[r], s[r]) (r < cnt、値 E byte) を 32 byte ずつ行い、済んだ個数を返す (残りは呼び出し元)。
// d と s は cnt 個ずつで、32 byte の中で重ならないこと (d と s の差が 32 byte 以上)。
template <int E, class F> SZ_AI std::size_t run_vec(char* d, const char* s, std::size_t cnt, F f) {
 constexpr std::size_t L= 32 / E;
 std::size_t r= 0;
 for(; r + L <= cnt; r+= L) st(d + r * E, f(ld(d + r * E), ld(s + r * E)));
 return r;
}
// 値 E byte の n 個 (2 冪) を変換する。32 byte に満たなければ false を返す (呼び出し元が汎用の道で行う)。
template <int E, bool UP, class F> inline bool directed(void* p, std::size_t n, F f) {
 if(n * E < 32) return false;
 dir_rec<E, UP>(static_cast<char*>(p), n * E / 32, f);
 return true;
}
template <int E, class B> inline bool butterfly(void* p, std::size_t n, B b) {
 if(n * E < 32) return false;
 bfly_rec<E>(static_cast<char*>(p), n * E / 32, b);
 return true;
}
#undef SZ_AI
}  // namespace sps_zeta
