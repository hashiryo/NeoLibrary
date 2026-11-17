#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <type_traits>
#include <utility>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
#include "neo/algebra/ZMod.hpp"
// 法が奇数の ZMod (法 p < 2^30) の集合冪級数の速い核。sps の各ヘッダから呼ぶ内部の道具で、直接は使わない。
// rank 付きの XOR 変換を使う。f[S] を f[S] t^{|S|} に置いて XOR 畳み込みをすると、S xor T = U の組は t^{|S|+|T|} = t^{|U| + 2|S∩T|} の
// 項に入るので、t^{|U|} の係数が部分集合畳み込みになる。部分集合の大きさより rank の大きい項 (重なる組から出る項) は、足し算と
// 掛け算のどちらを通しても大きさより rank の大きいままなので、途中で逆変換して rank を取り出さずに済み、i bit の集合冪級数は
// rank i までで切れる。これで exp、log、合成、power projection を、上の要素から 1 つずつ足す変換の領域の各点の計算として書ける。
// 変換の添字 y を補集合に移すと t が -t に変わるので、y の bit 0 が 0 の半分だけを求め (bit 0 の変換は f[2S'] と f[2S'+1] を同じ
// 位置に置くだけ)、逆変換で 2/2^n を掛ける。値は u32 で [0, 2p] に持ち、位置 8 個を 1 本の __m256i に並べる (塊 = 位置 8 個 ×
// rank の数)。下の bit から変換すると、bit b まで済んだ位置の値は rank が [残りの上の bit の数, それ + b + 1] にしかないので、その
// 範囲 (窓) だけを変換する。逆変換は上の bit からにして、最後の数段は答えに要る rank だけを持つ。積は [0, p] に直してから 12 項
// までずつ u64 にためて Montgomery で還元する。塊を 2^9 個ずつの区画に分けて下の段を区画の中で、上の段は列ごとに作業の配列で行う。
// procon-judge の Library Checker の 5 問で書き比べた xw1、xe1、xl1、xc1、xp1 の法を実行時の奇数に広げたもの。
namespace sps_xor {
using u32= std::uint32_t;
using u64= std::uint64_t;
using V= __m256i;
#define SPS_AI [[gnu::always_inline]] inline
constexpr int KC= 9;     // 区画の塊の数は 2^KC
constexpr int NMAX= 22;  // 速い核で受ける集合の大きさの上限
inline constexpr u32 inv32(u32 n) {
 u32 x= n;
 for(int i= 0; i < 5; ++i) x*= 2u - n * x;
 return x;
}
inline u32 pw(u64 a, u64 e, u32 p) {
 u64 r= 1;
 for(a%= p; e; e>>= 1, a= a * a % p)
  if(e & 1) r= r * a % p;
 return u32(r);
}
// 法の定数。p は奇数で 3 <= p < 2^30。
struct MK {
 u32 P, R1, R2;
 V p, p2, p64, np64;
 explicit MK(u32 m): P(m), R1(u32((u64(1) << 32) % m)), R2(u32(u64(R1) * R1 % m)) {
  p= _mm256_set1_epi32(int(m)), p2= _mm256_set1_epi32(int(2 * m)), p64= _mm256_set1_epi64x(m), np64= _mm256_set1_epi64x(u32(0u - inv32(m)));
 }
 u32 inv2pow(int e) const { return pw((P + 1) / 2, e, P); }  // 2^{-e} (法は奇数なので 2^{-1} = (p + 1) / 2)
};
SPS_AI V ld(const u32* p) { return _mm256_load_si256((const V*)p); }
SPS_AI void st(u32* p, V x) { _mm256_store_si256((V*)p, x); }
SPS_AI V cst(u32 x) { return _mm256_set1_epi32(int(x)); }
SPS_AI V red2(V x, const MK& k) { return _mm256_min_epu32(x, _mm256_sub_epi32(x, k.p2)); }  // [0, 4p) -> [0, 2p]
SPS_AI V neg2(V x, const MK& k) { return _mm256_sub_epi32(k.p2, x); }
SPS_AI V add2(V a, V b, const MK& k) { return red2(_mm256_add_epi32(a, b), k); }
SPS_AI V sub2(V a, V b, const MK& k) { return red2(_mm256_add_epi32(a, neg2(b, k)), k); }
SPS_AI V norm1(V x, const MK& k) { return _mm256_min_epu32(x, _mm256_sub_epi32(x, k.p)); }  // [0, 2p] -> [0, p]
SPS_AI V mul64(V a, V b) { return _mm256_mul_epu32(a, b); }
// e は偶数のレーン、o は奇数のレーンの u64 の和。それぞれ x R^{-1} に還元して 8 レーンに並べる。
SPS_AI V redc2(V e, V o, const MK& k) {
 V te= _mm256_add_epi64(e, mul64(mul64(e, k.np64), k.p64)), to= _mm256_add_epi64(o, mul64(mul64(o, k.np64), k.p64));
 return _mm256_blend_epi32(_mm256_srli_epi64(te, 32), to, 0xAA);
}
SPS_AI V mont(V x, V y, const MK& k) { return redc2(mul64(x, y), mul64(_mm256_srli_epi64(x, 32), _mm256_srli_epi64(y, 32)), k); }  // x y R^{-1}
SPS_AI V lane_wht(V x, const MK& k) {
 V t= _mm256_shuffle_epi32(x, 0xB1);
 x= red2(_mm256_add_epi32(t, _mm256_blend_epi32(x, neg2(x, k), 0xAA)), k);
 t= _mm256_shuffle_epi32(x, 0x4E);
 x= red2(_mm256_add_epi32(t, _mm256_blend_epi32(x, neg2(x, k), 0xCC)), k);
 t= _mm256_permute2x128_si256(x, x, 0x01);
 x= red2(_mm256_add_epi32(t, _mm256_blend_epi32(x, neg2(x, k), 0xF0)), k);
 return x;
}
SPS_AI V mpc(int j) {
 auto b= [&](int r) { return __builtin_popcount(r) == j ? -1 : 0; };
 return _mm256_setr_epi32(b(0), b(1), b(2), b(3), b(4), b(5), b(6), b(7));
}
inline u32* huge_alloc(std::size_t bytes) {
 const std::size_t H= std::size_t(1) << 21;
 bytes= (bytes + H - 1) & ~(H - 1);
 void* p= std::aligned_alloc(H, bytes);
#ifdef __linux__
 madvise(p, bytes, MADV_HUGEPAGE);
#endif
 return static_cast<u32*>(p);
}
inline u32* small_alloc(std::size_t bytes) { return static_cast<u32*>(std::aligned_alloc(32, (bytes + 31) & ~std::size_t(31))); }
// 塊には rank r の値を場所 r - sb に置き、場所は [0, ns) だけを持つ。
// 入力 16 個 (rank は lo + popcount) を置いてレーンの WHT まで済ませる。lo は塊の上の bit の数 (+ 全体の rank のずれ)。
SPS_AI void init_block(u32* X, const u32* in, int lo, int sb, int ns, const MK& k) {
 V lo_= _mm256_loadu_si256((const V*)in), hi= _mm256_loadu_si256((const V*)(in + 8));
 __m256 lf= _mm256_castsi256_ps(lo_), hf= _mm256_castsi256_ps(hi);
 V E= _mm256_permute4x64_epi64(_mm256_castps_si256(_mm256_shuffle_ps(lf, hf, 0x88)), 0xD8);
 V O= _mm256_permute4x64_epi64(_mm256_castps_si256(_mm256_shuffle_ps(lf, hf, 0xDD)), 0xD8);
 const V M0= mpc(0), M1= mpc(1), M2= mpc(2), M3= mpc(3);
 V v[5];
 v[0]= _mm256_broadcastd_epi32(_mm256_castsi256_si128(E));
 v[1]= lane_wht(_mm256_or_si256(_mm256_and_si256(E, M1), _mm256_and_si256(O, M0)), k);
 v[2]= lane_wht(_mm256_or_si256(_mm256_and_si256(E, M2), _mm256_and_si256(O, M1)), k);
 v[3]= lane_wht(_mm256_or_si256(_mm256_and_si256(E, M3), _mm256_and_si256(O, M2)), k);
 V o7= _mm256_permutevar8x32_epi32(O, cst(7));
 v[4]= _mm256_blend_epi32(o7, neg2(o7, k), 0x96);
 for(int j= 0; j < 5; ++j)
  if(int s= lo + j - sb; s >= 0 && s < ns) st(X + 8 * s, v[j]);
}
// 塊の bit kk の変換。X (bit が 0) は rank [lo, lo+4+kk]、Y は [lo+1, lo+5+kk] を持ち、どちらも [lo, lo+5+kk] になる。
SPS_AI void fwd_pair(u32* X, u32* Y, int lo, int kk, int sb, int ns, const MK& k) {
 if(int s= lo - sb; s >= 0 && s < ns) st(Y + 8 * s, ld(X + 8 * s));
 const int a= std::max(lo + 1 - sb, 0), b= std::min(lo + 4 + kk - sb, ns - 1);
 for(int s= a; s <= b; ++s) {
  V x= ld(X + 8 * s), y= ld(Y + 8 * s);
  st(X + 8 * s, add2(x, y, k)), st(Y + 8 * s, sub2(x, y, k));
 }
 if(int s= lo + 5 + kk - sb; s >= 0 && s < ns) {
  V y= ld(Y + 8 * s);
  st(X + 8 * s, y), st(Y + 8 * s, neg2(y, k));
 }
}
// 逆変換の塊の bit kk。どちらも rank [lo, lo+5+kk] を持ち、X は [lo, lo+4+kk]、Y は [lo+1, lo+5+kk] を残す。
SPS_AI void inv_pair(u32* X, u32* Y, int lo, int kk, int sb, int ns, const MK& k) {
 if(int s= lo - sb; s >= 0 && s < ns) st(X + 8 * s, add2(ld(X + 8 * s), ld(Y + 8 * s), k));
 const int a= std::max(lo + 1 - sb, 0), b= std::min(lo + 4 + kk - sb, ns - 1);
 for(int s= a; s <= b; ++s) {
  V x= ld(X + 8 * s), y= ld(Y + 8 * s);
  st(X + 8 * s, add2(x, y, k)), st(Y + 8 * s, sub2(x, y, k));
 }
 if(int s= lo + 5 + kk - sb; s >= 0 && s < ns) st(Y + 8 * s, sub2(ld(X + 8 * s), ld(Y + 8 * s), k));
}
// 塊 q の、rank r (場所 r - sb) の値をレーンの WHT で戻し、16 個の位置の rank |U| の値を取り出す (偶数と奇数の位置に並べた 2 本)。
SPS_AI void extract(const u32* X, int c, int sb, int ns, V q[2], const MK& k) {
 const V M0= mpc(0), M1= mpc(1), M2= mpc(2), M3= mpc(3), z= _mm256_setzero_si256();
 V v[5];
 for(int j= 0; j < 5; ++j) {
  int s= c + j - sb;
  v[j]= s >= 0 && s < ns ? lane_wht(ld(X + 8 * s), k) : z;
 }
 V ev= _mm256_or_si256(_mm256_or_si256(_mm256_and_si256(v[0], M0), _mm256_and_si256(v[1], M1)), _mm256_or_si256(_mm256_and_si256(v[2], M2), _mm256_and_si256(v[3], M3)));
 V od= _mm256_or_si256(_mm256_or_si256(_mm256_and_si256(v[1], M0), _mm256_and_si256(v[2], M1)), _mm256_or_si256(_mm256_and_si256(v[3], M2), _mm256_and_si256(v[4], M3)));
 V ul= _mm256_unpacklo_epi32(ev, od), uh= _mm256_unpackhi_epi32(ev, od);
 q[0]= _mm256_permute2x128_si256(ul, uh, 0x20), q[1]= _mm256_permute2x128_si256(ul, uh, 0x31);
}
// 和 sum_{i} A_i B_{d-i} の u64 の和を、12 項までずつ還元して [0, 2p] で返す。Ae/Ao と Be/Bo は [0, p] の値の偶数と奇数のレーン。
// 項は a = a0..a1 (含む) で、B の添字は bo - a。
SPS_AI V dot(const u32* Ae, const u32* Ao, const u32* Be, const u32* Bo, int a0, int a1, int bo, V e0, V o0, const MK& k) {
 const V z= _mm256_setzero_si256();
 V r= z;
 bool first= true;
 for(int s= a0; s <= a1; s+= 12) {
  const int t= std::min(a1, s + 11);
  V e= first ? e0 : z, o= first ? o0 : z;
  for(int a= s; a <= t; ++a) {
   e= _mm256_add_epi64(e, mul64(ld(Ae + 8 * a), ld(Be + 8 * (bo - a))));
   o= _mm256_add_epi64(o, mul64(ld(Ao + 8 * a), ld(Bo + 8 * (bo - a))));
  }
  r= first ? red2(redc2(e, o, k), k) : add2(r, red2(redc2(e, o, k), k), k);
  first= false;
 }
 if(first) r= red2(redc2(e0, o0, k), k);
 return r;
}
SPS_AI void split_lanes(const u32* X, int cnt, u32* E, u32* O, const MK& k, bool zext) {
 const V lo32= _mm256_set1_epi64x(0xffffffffu);
 for(int i= 0; i < cnt; ++i) {
  V x= norm1(ld(X + 8 * i), k);
  st(E + 8 * i, zext ? _mm256_and_si256(x, lo32) : x), st(O + 8 * i, _mm256_srli_epi64(x, 32));
 }
}

// 塊の列を変換する。in は 2^m 個 (m >= 4、rank は popcount + roff)、X は 2^{m-4} 個の塊 (幅 bs)。区画の中の段だけ行う。
inline void fwd_chunks(u32* X, const u32* in, int m, int roff, int sb, int ns, std::size_t bs, const MK& k) {
 const int nb= m - 4, K= std::min(nb, KC);
 const std::size_t cnt= std::size_t(1) << nb, csz= std::size_t(1) << K;
 for(std::size_t q0= 0; q0 < cnt; q0+= csz) {
  for(std::size_t q= q0; q < q0 + csz; ++q) init_block(X + q * bs, in + 16 * q, __builtin_popcountll(q) + roff, sb, ns, k);
  for(int kk= 0; kk < K; ++kk)
   for(std::size_t q= q0; q < q0 + csz; ++q)
    if(!(q >> kk & 1)) fwd_pair(X + q * bs, X + (q | std::size_t(1) << kk) * bs, __builtin_popcountll(q >> (kk + 1)) + roff, kk, sb, ns, k);
 }
}
// 列 j の塊 (h << K) + j を c に写して、上の段 (区画をまたぐ段) を変換する。
inline void fwd_column(u32* c, const u32* X, std::size_t j, int m, int roff, int sb, int ns, std::size_t bs, const MK& k) {
 const int nb= m - 4, K= std::min(nb, KC), H= nb - K;
 const std::size_t hs= std::size_t(1) << H;
 for(std::size_t h= 0; h < hs; ++h) std::memcpy(c + h * bs, X + ((h << K) + j) * bs, std::size_t(ns) * 32);
 for(int kk= 0; kk < H; ++kk)
  for(std::size_t h= 0; h < hs; ++h)
   if(!(h >> kk & 1)) fwd_pair(c + h * bs, c + (h | std::size_t(1) << kk) * bs, __builtin_popcountll(h >> (kk + 1)) + roff, kk + K, sb, ns, k);
}
inline void inv_column(u32* c, u32* X, std::size_t j, int m, int roff, int sb, int ns, std::size_t bs, const MK& k) {
 const int nb= m - 4, K= std::min(nb, KC), H= nb - K;
 const std::size_t hs= std::size_t(1) << H;
 for(int kk= H; kk--;)
  for(std::size_t h= 0; h < hs; ++h)
   if(!(h >> kk & 1)) inv_pair(c + h * bs, c + (h | std::size_t(1) << kk) * bs, __builtin_popcountll(h >> (kk + 1)) + roff, kk + K, sb, ns, k);
 for(std::size_t h= 0; h < hs; ++h) std::memcpy(X + ((h << K) + j) * bs, c + h * bs, std::size_t(ns) * 32);
}
// 逆変換の区画の中の段と、塊ごとの取り出し。f(q, 取り出した 2 本) を呼ぶ。
template <class Out> inline void inv_chunks(u32* X, int m, int roff, int sb, int ns, std::size_t bs, const MK& k, const Out& out) {
 const int nb= m - 4, K= std::min(nb, KC);
 const std::size_t cnt= std::size_t(1) << nb, csz= std::size_t(1) << K;
 for(std::size_t q0= 0; q0 < cnt; q0+= csz) {
  for(int kk= K; kk--;)
   for(std::size_t q= q0; q < q0 + csz; ++q)
    if(!(q >> kk & 1)) inv_pair(X + q * bs, X + (q | std::size_t(1) << kk) * bs, __builtin_popcountll(q >> (kk + 1)) + roff, kk, sb, ns, k);
  for(std::size_t q= q0; q < q0 + csz; ++q) {
   V v[2];
   extract(X + q * bs, __builtin_popcountll(q) + roff, sb, ns, v, k);
   out(q, v);
  }
 }
}
inline std::size_t col_blocks(int m) {
 const int nb= m - 4, K= std::min(nb, KC);
 return std::size_t(1) << (nb - K);
}
// 値 x (2p 未満) を R 倍して [0, p) に (n 個、8 の倍数)。
inline void scale_r(const u32* x, u32* y, std::size_t n, const MK& k) {
 const V r2= cst(k.R2);
 for(std::size_t s= 0; s < n; s+= 8) _mm256_store_si256((V*)(y + s), norm1(mont(_mm256_loadu_si256((const V*)(x + s)), r2, k), k));
}

// ---- 部分集合畳み込み ----
// 塊ごとの rank の積 (長さ nm1 で切り詰める)。入力の rank r は場所 r - 1、積の rank r は場所 r - 2。結果は a に書く。
template <int NM1> SPS_AI void conv_block(u32* a, const u32* b, const MK& k) {
 alignas(32) u32 Ae[NM1 * 8], Ao[NM1 * 8], Be[NM1 * 8], Bo[NM1 * 8];
 split_lanes(a, NM1, Ae, Ao, k, false), split_lanes(b, NM1, Be, Bo, k, false);
 const V z= _mm256_setzero_si256();
 for(int s= 0; s < NM1; ++s) st(a + 8 * s, dot(Ae, Ao, Be, Bo, 0, s, s, z, z, k));
}
template <int N> void convolve_n(const u32* f, const u32* g, u32* out, const MK& k) {
 constexpr int NM1= N - 1, NB= N - 4, BS= NM1 * 8;
 const std::size_t nblk= std::size_t(1) << NB;
 u32 *A= huge_alloc(nblk * BS * 4), *B= huge_alloc(nblk * BS * 4);
 fwd_chunks(A, f, N, 0, 1, NM1, BS, k), fwd_chunks(B, g, N, 0, 1, NM1, BS, k);
 // 空集合は入力から外す (rank 0 は場所 -1 なので置かれない)
 if constexpr(NB <= KC) {
  for(std::size_t q= 0; q < nblk; ++q) conv_block<NM1>(A + q * BS, B + q * BS, k);
 } else {
  const std::size_t hs= col_blocks(N), K= KC;
  u32 *cA= small_alloc(hs * BS * 4), *cB= small_alloc(hs * BS * 4);
  for(std::size_t j= 0; j < (std::size_t(1) << K); ++j) {
   fwd_column(cA, A, j, N, 0, 1, NM1, BS, k), fwd_column(cB, B, j, N, 0, 1, NM1, BS, k);
   for(std::size_t h= 0; h < hs; ++h) conv_block<NM1>(cA + h * BS, cB + h * BS, k);
   inv_column(cA, A, j, N, 0, 2, NM1, BS, k);
  }
  std::free(cA), std::free(cB);
 }
 const V Kc= cst(u32(u64(k.inv2pow(N - 1)) * k.R2 % k.P)), F0= cst(u32(u64(f[0] % k.P) * k.R1 % k.P)), G0= cst(u32(u64(g[0] % k.P) * k.R1 % k.P));
 inv_chunks(A, N, 0, 2, NM1, BS, k, [&](std::size_t q, const V* v) {
  for(int t= 0; t < 2; ++t) {
   const u32 *fp= f + 16 * q + 8 * t, *gp= g + 16 * q + 8 * t;
   V fv= norm1(_mm256_loadu_si256((const V*)fp), k), gv= norm1(_mm256_loadu_si256((const V*)gp), k);
   V e= _mm256_add_epi64(_mm256_add_epi64(mul64(v[t], Kc), mul64(gv, F0)), mul64(fv, G0));
   V o= _mm256_add_epi64(_mm256_add_epi64(mul64(_mm256_srli_epi64(v[t], 32), Kc), mul64(_mm256_srli_epi64(gv, 32), F0)), mul64(_mm256_srli_epi64(fv, 32), G0));
   _mm256_storeu_si256((V*)(out + 16 * q + 8 * t), norm1(redc2(e, o, k), k));
  }
 });
 out[0]= u32(u64(f[0] % k.P) * (g[0] % k.P) % k.P);
 std::free(A), std::free(B);
}

// ---- exp ----
// 1 bit 伸ばす。G (rank 1..i を場所 0..i-1、rank 0 は 1) と F (f̂_i を R 倍、rank 0..i を場所 0..i) から D = Ĝ t f̂ (rank i+1 まで) を作り、
// Glo = G + D、Ghi = G - D を rank 1..i+1 (場所 0..i) に書く。Glo は G と同じ場所でよい。
SPS_AI void exp_step(u32* G, u32* Ghi, const u32* F, int i, const MK& k) {
 alignas(32) u32 Ge[(NMAX + 1) * 8], Go[(NMAX + 1) * 8], Fe[(NMAX + 1) * 8], Fo[(NMAX + 1) * 8];
 split_lanes(G, i, Ge + 8, Go + 8, k, false);
 split_lanes(F, i + 1, Fe, Fo, k, true);
 const V z= _mm256_setzero_si256();
 for(int d= 1; d <= i + 1; ++d) {
  // D_d = F_{d-1} + sum_{a=1}^{d-1} G_a F_{d-1-a}
  V D= dot(Ge, Go, Fe, Fo, 1, d - 1, d - 1, ld(Fe + 8 * (d - 1)), ld(Fo + 8 * (d - 1)), k);
  V g= d <= i ? ld(G + 8 * (d - 1)) : z;
  st(Ghi + 8 * (d - 1), sub2(g, D, k)), st(G + 8 * (d - 1), add2(g, D, k));
 }
}
inline void naive_exp(int m, const u32* b, u32* h, u32 P) {
 h[0]= 1;
 for(int i= 0; i < m; ++i)
  for(int s= 0; s < (1 << i); ++s) {
   u64 t= 0;
   for(int a= s;; a= (a - 1) & s) {
    t= (t + u64(h[a]) * b[(1 << i) | (s ^ a)]) % P;
    if(!a) break;
   }
   h[(1 << i) | s]= u32(t);
  }
}
template <int N> void exp_n(const u32* b, u32* out, const MK& k) {
 constexpr int NS= N, BS= NS * 8, NB= N - 4;
 const std::size_t nblk= std::size_t(1) << NB;
 u32 *G= huge_alloc(nblk * BS * 4), *F= huge_alloc((nblk >> 1) * BS * 4), *bs= huge_alloc(std::size_t(4) << N);
 scale_r(b, bs, std::size_t(1) << N, k);
 {
  alignas(32) u32 h[16], bb[16];
  for(int s= 0; s < 16; ++s) bb[s]= b[s] % k.P;
  naive_exp(4, bb, h, k.P);
  init_block(G, h, 0, 1, NS, k);
 }
 u32* cF= small_alloc(col_blocks(N - 1) * BS * 4);
 for(int i= 4; i < N; ++i) {
  const int K= std::min(i - 4, KC);
  const std::size_t cnt= std::size_t(1) << (i - 4);
  fwd_chunks(F, bs + (std::size_t(1) << i), i, 0, 0, i + 1, BS, k);
  if(i - 4 <= KC) {
   for(std::size_t q= 0; q < cnt; ++q) exp_step(G + q * BS, G + (q + cnt) * BS, F + q * BS, i, k);
  } else {
   const std::size_t hs= col_blocks(i);
   for(std::size_t j= 0; j < (std::size_t(1) << K); ++j) {
    fwd_column(cF, F, j, i, 0, 0, i + 1, BS, k);
    for(std::size_t h= 0; h < hs; ++h) {
     const std::size_t q= (h << K) + j;
     exp_step(G + q * BS, G + (q + cnt) * BS, cF + h * BS, i, k);
    }
   }
  }
 }
 if constexpr(NB > KC) {
  const std::size_t hs= col_blocks(N);
  u32* cG= small_alloc(hs * BS * 4);
  for(std::size_t j= 0; j < (std::size_t(1) << KC); ++j) {
   for(std::size_t h= 0; h < hs; ++h) std::memcpy(cG + h * BS, G + ((h << KC) + j) * BS, BS * 4);
   inv_column(cG, G, j, N, 0, 1, NS, BS, k);
  }
  std::free(cG);
 }
 const V Kc= cst(u32(u64(k.inv2pow(N - 1)) * k.R1 % k.P));
 inv_chunks(G, N, 0, 1, NS, BS, k, [&](std::size_t q, const V* v) {
  for(int t= 0; t < 2; ++t) _mm256_storeu_si256((V*)(out + 16 * q + 8 * t), norm1(mont(v[t], Kc, k), k));
 });
 out[0]= 1;
 std::free(G), std::free(F), std::free(bs), std::free(cF);
}

// ---- 割り算 (log と逆元と割り算の共通) ----
// A = Â' (R 倍、rank 0..i を場所 0..i)、B = B̂' (R 倍、rank 1..i を場所 0..i-1、rank 0 は 1)。Q' = A / B を rank i まで求めて A に書く。
// upd なら B を B̂'[y +] = B + t A に、上の半分 Bhi を B - t A にする (rank 1..i+1)。
SPS_AI void div_step(u32* A, u32* B, u32* Bhi, int i, bool upd, const MK& k) {
 alignas(32) u32 Be[(NMAX + 1) * 8], Bo[(NMAX + 1) * 8], Qe[(NMAX + 1) * 8], Qo[(NMAX + 1) * 8], Q[(NMAX + 1) * 8];
 split_lanes(B, i, Be + 8, Bo + 8, k, false);
 const V z= _mm256_setzero_si256();
 st(Q, ld(A));
 split_lanes(A, 1, Qe, Qo, k, false);
 for(int d= 1; d <= i; ++d) {
  V S= dot(Be, Bo, Qe, Qo, 1, d, d, z, z, k);
  V q= sub2(ld(A + 8 * d), S, k);
  st(Q + 8 * d, q);
  split_lanes(Q + 8 * d, 1, Qe + 8 * d, Qo + 8 * d, k, false);
 }
 if(upd)
  for(int r= 1; r <= i + 1; ++r) {
   V b= r <= i ? ld(B + 8 * (r - 1)) : z, a= ld(A + 8 * (r - 1));
   st(Bhi + 8 * (r - 1), sub2(b, a, k)), st(B + 8 * (r - 1), add2(b, a, k));
  }
 for(int d= 0; d <= i; ++d) st(A + 8 * d, ld(Q + 8 * d));
}
inline void naive_log(int m, const u32* h, u32* f, u32 P) {
 f[0]= 0;
 for(int i= 0; i < m; ++i)
  for(int s= 0; s < (1 << i); ++s) {
   u64 t= h[(1 << i) | s];
   for(int a= s; a; a= (a - 1) & s) t= (t + u64(P - h[a]) * f[(1 << i) | (s ^ a)]) % P;
   f[(1 << i) | s]= u32(t % P);
  }
}
// log (h[0] は 1)。
template <int N> void log_n(const u32* h, u32* out, const MK& k) {
 constexpr int NSB= N - 1, BSB= NSB * 8, NSF= N, BSF= NSF * 8;
 u32 *B= huge_alloc((std::size_t(1) << (N - 5)) * BSB * 4), *F= huge_alloc((std::size_t(1) << (N - 5)) * BSF * 4), *hs= huge_alloc(std::size_t(4) << N);
 scale_r(h, hs, std::size_t(1) << N, k);
 {
  u32 hh[16];
  for(int s= 0; s < 16; ++s) hh[s]= h[s] % k.P;
  naive_log(4, hh, out, k.P);
 }
 init_block(B, hs, 0, 1, NSB, k);
 u32* cF= small_alloc(col_blocks(N - 1) * BSF * 4);
 for(int i= 4; i < N; ++i) {
  const int K= std::min(i - 4, KC), ns= i + 1;
  const std::size_t cnt= std::size_t(1) << (i - 4);
  const bool upd= i + 1 < N;
  fwd_chunks(F, hs + (std::size_t(1) << i), i, 0, 0, ns, BSF, k);
  if(i - 4 <= KC) {
   for(std::size_t q= 0; q < cnt; ++q) div_step(F + q * BSF, B + q * BSB, B + (q + cnt) * BSB, i, upd, k);
  } else {
   const std::size_t hs_= col_blocks(i);
   for(std::size_t j= 0; j < (std::size_t(1) << K); ++j) {
    fwd_column(cF, F, j, i, 0, 0, ns, BSF, k);
    for(std::size_t hh= 0; hh < hs_; ++hh) {
     const std::size_t q= (hh << K) + j;
     div_step(cF + hh * BSF, B + q * BSB, B + (q + cnt) * BSB, i, upd, k);
    }
    inv_column(cF, F, j, i, 0, 0, ns, BSF, k);
   }
  }
  const V Kc= cst(k.inv2pow(i - 1));
  u32* o= out + (std::size_t(1) << i);
  inv_chunks(F, i, 0, 0, ns, BSF, k, [&](std::size_t q, const V* v) {
   for(int t= 0; t < 2; ++t) _mm256_storeu_si256((V*)(o + 16 * q + 8 * t), norm1(mont(v[t], Kc, k), k));
  });
 }
 std::free(B), std::free(F), std::free(hs), std::free(cF);
}
// 割り算 f / g (g[0] は 1)。変換の領域の各点で割って 1 回だけ逆変換する。
template <int N> void div_n(const u32* f, const u32* g, u32* out, const MK& k) {
 constexpr int NS= N + 1, BS= NS * 8, NB= N - 4;
 const std::size_t nblk= std::size_t(1) << NB;
 u32 *A= huge_alloc(nblk * BS * 4), *B= huge_alloc(nblk * BS * 4), *gs= huge_alloc(std::size_t(4) << N), *fs= huge_alloc(std::size_t(4) << N);
 scale_r(g, gs, std::size_t(1) << N, k), scale_r(f, fs, std::size_t(1) << N, k);
 // A は f の R 倍 (rank 0..N)、B は g の R 倍 (rank 1..N を場所 0..N-1)
 fwd_chunks(A, fs, N, 0, 0, NS, BS, k), fwd_chunks(B, gs, N, 0, 1, N, BS, k);
 auto body= [&](u32* a, u32* b) {
  // Q' = A' / B' (B の rank 0 は R、R 倍のままの割り算は Q'_d = A'_d - REDC(sum B'_j Q'_{d-j}))
  alignas(32) u32 Be[(NMAX + 2) * 8], Bo[(NMAX + 2) * 8], Qe[(NMAX + 2) * 8], Qo[(NMAX + 2) * 8];
  split_lanes(b, N, Be + 8, Bo + 8, k, false);
  split_lanes(a, 1, Qe, Qo, k, false);
  const V z= _mm256_setzero_si256();
  for(int d= 1; d <= N; ++d) {
   V q= sub2(ld(a + 8 * d), dot(Be, Bo, Qe, Qo, 1, d, d, z, z, k), k);
   st(a + 8 * d, q);
   split_lanes(a + 8 * d, 1, Qe + 8 * d, Qo + 8 * d, k, false);
  }
 };
 if constexpr(NB <= KC) {
  for(std::size_t q= 0; q < nblk; ++q) body(A + q * BS, B + q * BS);
 } else {
  const std::size_t hs= col_blocks(N);
  u32 *cA= small_alloc(hs * BS * 4), *cB= small_alloc(hs * BS * 4);
  for(std::size_t j= 0; j < (std::size_t(1) << KC); ++j) {
   fwd_column(cA, A, j, N, 0, 0, NS, BS, k), fwd_column(cB, B, j, N, 0, 1, N, BS, k);
   for(std::size_t h= 0; h < hs; ++h) body(cA + h * BS, cB + h * BS);
   inv_column(cA, A, j, N, 0, 0, NS, BS, k);
  }
  std::free(cA), std::free(cB);
 }
 const V Kc= cst(k.inv2pow(N - 1));
 inv_chunks(A, N, 0, 0, NS, BS, k, [&](std::size_t q, const V* v) {
  for(int t= 0; t < 2; ++t) _mm256_storeu_si256((V*)(out + 16 * q + 8 * t), norm1(mont(v[t], Kc, k), k));
 });
 std::free(A), std::free(B), std::free(gs), std::free(fs);
}

// ---- 合成 (EGF) ----
// F = ĝ'_i (R 倍、rank 0..i)。G[j] は Ĝ_j の塊 (rank 0..i)、Gh[j] は上の半分の行き先。j = 0..J-1 について
// D = t Ĝ_{j+1} ĝ_i (rank 1..i+1)、Ĝ_j の下を G_j + D、上を G_j - D (rank 0..i+1) にする。
SPS_AI void comp_step(u32* const* G, u32* const* Gh, const u32* F, int i, int J, const MK& k) {
 alignas(32) u32 Fe[(NMAX + 1) * 8], Fo[(NMAX + 1) * 8], Ae[(NMAX + 1) * 8], Ao[(NMAX + 1) * 8];
 split_lanes(F, i + 1, Fe, Fo, k, false);
 const V z= _mm256_setzero_si256();
 for(int j= 0; j < J; ++j) {
  split_lanes(G[j + 1], i + 1, Ae, Ao, k, false);
  u32 *X= G[j], *Y= Gh[j];
  st(Y, ld(X));
  for(int d= 1; d <= i + 1; ++d) {
   V D= dot(Ae, Ao, Fe, Fo, 0, d - 1, d - 1, z, z, k);
   V g= d <= i ? ld(X + 8 * d) : z;
   st(Y + 8 * d, sub2(g, D, k)), st(X + 8 * d, add2(g, D, k));
  }
 }
}
// sum_k e_k g^k / k! (g[0] は 0、e は n + 1 個)。
template <int N> void egf_comp_n(const u32* g, const u32* e, u32* out, const MK& k) {
 constexpr int NB= N - 4;
 u32* gs= huge_alloc(std::size_t(4) << N);
 scale_r(g, gs, std::size_t(1) << N, k);
 u32* G[N - 3];
 std::size_t bsj[N - 3];
 for(int j= 0; j <= N - 4; ++j) bsj[j]= std::size_t(N - j + 1) * 8, G[j]= huge_alloc((std::size_t(1) << (N - j - 4)) * bsj[j] * 4);
 {
  std::vector<std::vector<u32>> h(N + 1, std::vector<u32>(16, 0));
  for(int j= 0; j <= N; ++j) h[j][0]= e[j];
  for(int i= 0; i < 4; ++i)
   for(int j= 0; j + i < N; ++j)
    for(int s= 0; s < (1 << i); ++s) {
     u64 t= 0;
     for(int a= s;; a= (a - 1) & s) {
      t= (t + u64(h[j + 1][a]) * (g[(1 << i) | (s ^ a)] % k.P)) % k.P;
      if(!a) break;
     }
     h[j][(1 << i) | s]= u32(t);
    }
  for(int j= 0; j <= N - 4; ++j) init_block(G[j], h[j].data(), 0, 0, 5, k);
 }
 constexpr std::size_t BSF= N * 8;
 u32 *F= huge_alloc((std::size_t(1) << (N - 5)) * BSF * 4), *cF= small_alloc(col_blocks(N - 1) * BSF * 4);
 for(int i= 4; i < N; ++i) {
  const int K= std::min(i - 4, KC), J= N - i;
  const std::size_t cnt= std::size_t(1) << (i - 4);
  fwd_chunks(F, gs + (std::size_t(1) << i), i, 0, 0, i + 1, BSF, k);
  u32 *gp[N - 3], *gh[N - 3];
  auto step= [&](std::size_t q, const u32* Fb) {
   for(int j= 0; j <= J; ++j) gp[j]= G[j] + q * bsj[j];
   for(int j= 0; j < J; ++j) gh[j]= G[j] + (q + cnt) * bsj[j];
   comp_step(gp, gh, Fb, i, J, k);
  };
  if(i - 4 <= KC) {
   for(std::size_t q= 0; q < cnt; ++q) step(q, F + q * BSF);
  } else {
   const std::size_t hs= col_blocks(i);
   for(std::size_t j= 0; j < (std::size_t(1) << K); ++j) {
    fwd_column(cF, F, j, i, 0, 0, i + 1, BSF, k);
    for(std::size_t h= 0; h < hs; ++h) step((h << K) + j, cF + h * BSF);
   }
  }
 }
 constexpr int NS= N + 1;
 constexpr std::size_t BS= NS * 8;
 u32* G0= G[0];
 if constexpr(NB > KC) {
  const std::size_t hs= col_blocks(N);
  u32* cG= small_alloc(hs * BS * 4);
  for(std::size_t j= 0; j < (std::size_t(1) << KC); ++j) {
   for(std::size_t h= 0; h < hs; ++h) std::memcpy(cG + h * BS, G0 + ((h << KC) + j) * BS, BS * 4);
   inv_column(cG, G0, j, N, 0, 0, NS, BS, k);
  }
  std::free(cG);
 }
 const V Kc= cst(u32(u64(k.inv2pow(N - 1)) * k.R1 % k.P));
 inv_chunks(G0, N, 0, 0, NS, BS, k, [&](std::size_t q, const V* v) {
  for(int t= 0; t < 2; ++t) _mm256_storeu_si256((V*)(out + 16 * q + 8 * t), norm1(mont(v[t], Kc, k), k));
 });
 out[0]= e[0];
 for(int j= 0; j <= N - 4; ++j) std::free(G[j]);
 std::free(F), std::free(cF), std::free(gs);
}

// ---- 合成の転置 (power projection) ----
// F = ĝ'_i (R 倍、rank 0..i)。Zl[j], Zh[j] (j < J) は段 i+1 の双対の下と上 (rank 0..i+1)。Y[j] (j = 0..J) に段 i の双対 (rank 0..i) を書く。
SPS_AI void proj_step(u32* const* Zl, u32* const* Zh, u32* const* Y, const u32* F, int i, int J, const MK& k) {
 alignas(32) u32 Fe[(NMAX + 1) * 8], Fo[(NMAX + 1) * 8], De[2][(NMAX + 2) * 8], Do[2][(NMAX + 2) * 8], S[(NMAX + 2) * 8], Fr[(NMAX + 1) * 8], Fq[(NMAX + 1) * 8];
 split_lanes(F, i + 1, Fe, Fo, k, false);
 // 相関 corr_a = sum_{d=a+1}^{i+1} F_{d-1-a} D_d を dot で書くため、F を逆順に並べた Fr (Fr[x] = F_{i - x}) を作る
 for(int x= 0; x <= i; ++x) std::memcpy(Fr + 8 * x, Fe + 8 * (i - x), 32), std::memcpy(Fq + 8 * x, Fo + 8 * (i - x), 32);
 const V z= _mm256_setzero_si256();
 for(int j= 0; j <= J; ++j) {
  const int cur= j & 1, prv= cur ^ 1;
  if(j < J) {
   for(int r= 0; r <= i; ++r) st(S + 8 * r, add2(ld(Zl[j] + 8 * r), ld(Zh[j] + 8 * r), k));
   for(int d= 1; d <= i + 1; ++d) {
    V x= norm1(sub2(ld(Zl[j] + 8 * d), ld(Zh[j] + 8 * d), k), k);
    st(De[cur] + 8 * d, x), st(Do[cur] + 8 * d, _mm256_srli_epi64(x, 32));
   }
  }
  for(int a= 0; a <= i; ++a) {
   V y= j < J ? ld(S + 8 * a) : z;
   // corr_a = sum_{d=a+1}^{i+1} D_d Fr[i + 1 - d + a] ... F_{d-1-a} = Fr[i - (d - 1 - a)] = Fr[i + 1 + a - d]
   if(j >= 1) y= add2(y, dot(De[prv], Do[prv], Fr, Fq, a + 1, i + 1, i + 1 + a, z, z, k), k);
   st(Y[j] + 8 * a, y);
  }
 }
}
inline std::vector<std::vector<u32>> egf_powers(int n, const u32* g, u32 P) {
 const int L= 1 << n;
 std::vector<std::vector<u32>> c(n + 1, std::vector<u32>(L, 0));
 c[0][0]= 1;
 for(int m= 1; m <= n; ++m)
  for(int s= 1; s < L; ++s) {
   const int lo= s & -s, rest= s ^ lo;
   u64 t= 0;
   for(int a= rest;; a= (a - 1) & rest) {
    t= (t + u64(g[lo | a] % P) * c[m - 1][rest ^ a]) % P;
    if(!a) break;
   }
   c[m][s]= u32(t);
  }
 return c;
}
// V_k = sum_S w[S] (g^k / k!)[S] (k = 0..N、g[0] は 0)。
template <int N> std::vector<u32> egf_proj_n(const u32* g, const u32* w, const MK& k) {
 constexpr int NB= N - 4, NS= N + 1;
 u32* gs= huge_alloc(std::size_t(4) << N);
 scale_r(g, gs, std::size_t(1) << N, k);
 u32* D[N - 3];
 std::size_t bsj[N - 3];
 for(int j= 0; j <= N - 4; ++j) bsj[j]= std::size_t(N - j + 1) * 8, D[j]= huge_alloc((std::size_t(1) << (N - j - 4)) * bsj[j] * 4);
 {
  constexpr std::size_t BS= NS * 8;
  u32* X= D[0];
  fwd_chunks(X, w, N, 0, 0, NS, BS, k);
  if constexpr(NB > KC) {
   const std::size_t hs= col_blocks(N);
   u32* c= small_alloc(hs * BS * 4);
   for(std::size_t j= 0; j < (std::size_t(1) << KC); ++j) {
    fwd_column(c, X, j, N, 0, 0, NS, BS, k);
    for(std::size_t h= 0; h < hs; ++h) std::memcpy(X + ((h << KC) + j) * BS, c + h * BS, BS * 4);
   }
   std::free(c);
  }
 }
 constexpr std::size_t BSF= N * 8;
 u32 *F= huge_alloc((std::size_t(1) << (N - 5)) * BSF * 4), *cF= small_alloc(col_blocks(N - 1) * BSF * 4);
 for(int i= N - 1; i >= 4; --i) {
  const int K= std::min(i - 4, KC), J= N - i - 1;
  const std::size_t cnt= std::size_t(1) << (i - 4);
  fwd_chunks(F, gs + (std::size_t(1) << i), i, 0, 0, i + 1, BSF, k);
  u32 *zl[N - 3], *zh[N - 3], *yy[N - 3];
  auto step= [&](std::size_t q, const u32* Fb) {
   for(int j= 0; j <= J; ++j) zl[j]= D[j] + q * bsj[j], zh[j]= D[j] + (q + cnt) * bsj[j];
   for(int j= 0; j <= J + 1; ++j) yy[j]= D[j] + q * bsj[j];
   proj_step(zl, zh, yy, Fb, i, J + 1, k);
  };
  if(i - 4 <= KC) {
   for(std::size_t q= 0; q < cnt; ++q) step(q, F + q * BSF);
  } else {
   const std::size_t hs= col_blocks(i);
   for(std::size_t j= 0; j < (std::size_t(1) << K); ++j) {
    fwd_column(cF, F, j, i, 0, 0, i + 1, BSF, k);
    for(std::size_t h= 0; h < hs; ++h) step((h << K) + j, cF + h * BSF);
   }
  }
 }
 // 段 4 の双対 (各 j に塊 1 つ、rank 0..4) を 4 bit の素朴な形に戻し、g^m / m! と内積を取る
 std::vector<std::vector<u32>> c= egf_powers(4, g, k.P);
 std::vector<u32> vk(N + 1, 0);
 for(int j= 0; j <= N - 4; ++j) {
  alignas(32) u32 v[5][8], gt[16];
  for(int r= 0; r <= 4; ++r) st(v[r], norm1(lane_wht(ld(D[j] + 8 * r), k), k));
  for(int s= 0; s < 16; ++s) gt[s]= v[__builtin_popcount(s)][s >> 1];
  for(int m= 0; m <= 4 && j + m <= N; ++m) {
   u64 t= 0;
   for(int s= 0; s < 16; ++s) t= (t + u64(gt[s]) * c[m][s]) % k.P;
   vk[j + m]= u32((vk[j + m] + t) % k.P);
  }
 }
 const u64 sc= k.inv2pow(N - 1);
 for(auto& x: vk) x= u32(x * sc % k.P);
 for(int j= 0; j <= N - 4; ++j) std::free(D[j]);
 std::free(F), std::free(cF), std::free(gs);
 return vk;
}

// 実行時の n から N のテンプレートを選ぶ。lo..NMAX の外なら false を返す。
template <int LO, class Fn> inline bool dispatch(int n, Fn&& fn) {
 bool ok= false;
 [&]<int... I>(std::integer_sequence<int, I...>) { ((n == LO + I ? (fn.template operator()<LO + I>(), ok= true) : false), ...); }(std::make_integer_sequence<int, NMAX - LO + 1>());
 return ok;
}

// ---- 合成と power projection の前後の処理 ----
// e_k = P^{(k)}(c0) (k = 0..n)。P は M 個の係数 (2p 未満)。微分を n 回、それぞれを c0 で評価し、どちらも 8 個ずつ Montgomery で。
inline std::vector<u32> taylor(const u32* f, int M, u32 c0, int n, const MK& k) {
 std::vector<u32> e(n + 1, 0);
 if(M == 0) return e;
 const std::size_t L= (std::size_t(M) + 15) & ~std::size_t(15);
 u32 *d= small_alloc(L * 4 + 64), *pc= small_alloc(L * 4 + 64), *mr= small_alloc(L * 4 + 64);
 std::memset(d, 0, L * 4 + 64);
 for(int m= 0; m < M; ++m) d[m]= f[m] % k.P;
 {
  u64 x= k.R1;  // c0^m R
  for(std::size_t m= 0; m < L; ++m) pc[m]= u32(x), x= x * c0 % k.P;
  for(std::size_t m= 0; m < L; ++m) mr[m]= u32(u64(m + 1) % k.P * k.R1 % k.P);
 }
 for(int kk= 0, len= M; kk <= n && len > 0; ++kk, --len) {
  V acc= _mm256_setzero_si256();
  for(int m= 0; m < len; m+= 8) acc= add2(acc, red2(mont(_mm256_loadu_si256((const V*)(d + m)), ld(pc + m), k), k), k);
  alignas(32) u32 t[8];
  st(t, acc);
  u64 sum= 0;
  for(int r= 0; r < 8; ++r) sum+= t[r];
  e[kk]= u32(sum % k.P);
  for(int m= 0; m + 1 < len; m+= 8) _mm256_storeu_si256((V*)(d + m), norm1(mont(_mm256_loadu_si256((const V*)(d + m + 1)), ld(mr + m), k), k));
  for(int m= len - 1; m < len + 8 && m < int(L); ++m) d[m]= 0;
 }
 std::free(d), std::free(pc), std::free(mr);
 return e;
}
// ans_m = sum_{k=0}^{n} V_k u_k(m) (m = 0..M-1)、u_k(m) = (m)_k c0^{m-k} ((m)_k は下降階乗)。
// u_k(m+1) = c0 u_k(m) + k u_{k-1}(m) で m を 1 つずつ進め、k を 8 個ずつレーンに並べる (割り算を使わないので c0 は何でもよい)。
inline void proj_post(const std::vector<u32>& vk, u32 c0, int M, u32* ans, const MK& k) {
 const int n= int(vk.size()) - 1;
 assert(n < 24);
 alignas(32) u32 vr[24], kr[24], u0[24];
 for(int i= 0; i < 24; ++i) vr[i]= i <= n ? u32(u64(vk[i]) * k.R1 % k.P) : 0, kr[i]= i <= n ? u32(u64(i) * k.R1 % k.P) : 0, u0[i]= i == 0;
 V U[3]= {ld(u0), ld(u0 + 8), ld(u0 + 16)}, W[3]= {ld(vr), ld(vr + 8), ld(vr + 16)}, KR[3]= {ld(kr), ld(kr + 8), ld(kr + 16)};
 const V C= cst(u32(u64(c0 % k.P) * k.R1 % k.P)), sh= _mm256_setr_epi32(7, 0, 1, 2, 3, 4, 5, 6), z= _mm256_setzero_si256();
 for(int m= 0; m < M; ++m) {
  V s= add2(add2(red2(mont(U[0], W[0], k), k), red2(mont(U[1], W[1], k), k), k), red2(mont(U[2], W[2], k), k), k);
  alignas(32) u32 t[8];
  st(t, s);
  u64 sum= 0;
  for(int r= 0; r < 8; ++r) sum+= t[r];
  ans[m]= u32(sum % k.P);
  const V r0= _mm256_permutevar8x32_epi32(U[0], sh), r1= _mm256_permutevar8x32_epi32(U[1], sh), r2= _mm256_permutevar8x32_epi32(U[2], sh);
  const V p[3]= {_mm256_blend_epi32(r0, z, 0x01), _mm256_blend_epi32(r1, r0, 0x01), _mm256_blend_epi32(r2, r1, 0x01)};
  for(int i= 0; i < 3; ++i) U[i]= add2(red2(mont(U[i], C, k), k), red2(mont(p[i], KR[i], k), k), k);
 }
}

// ---- sps の各ヘッダから使う道具 ----
template <class T> struct is_zmod: std::false_type {};
template <unsigned M> struct is_zmod<ZMod<M>>: std::true_type {};
// T が法が奇数の ZMod で、n が速い核で受ける大きさ (lo <= n <= NMAX) か。
template <class T> inline bool fast(int n, int lo) {
 if constexpr(is_zmod<T>::value) return (T::mod() & 1) && T::mod() >= 3 && lo <= n && n <= NMAX;
 else return false;
}
// ZMod の列をそのままの値 (2M 未満) の u32 の列として読む。ZMod は u32 を 1 つだけ持つ標準レイアウトの型なので、先頭の要素と
// ポインタを相互に変換できる。
template <class T> inline const u32* raw(const std::vector<T>& v) {
 static_assert(std::is_standard_layout_v<T> && sizeof(T) == 4);
 return reinterpret_cast<const u32*>(v.data());
}
template <class T> inline u32* raw(std::vector<T>& v) {
 static_assert(std::is_standard_layout_v<T> && sizeof(T) == 4);
 return reinterpret_cast<u32*>(v.data());
}
template <class T> inline int log2_size(std::size_t N) {
 assert(N && !(N & (N - 1)));
 return __builtin_ctzll(N);
}
#undef SPS_AI
}  // namespace sps_xor
