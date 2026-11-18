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
// 範囲 (窓) だけを変換する。逆変換は上の bit からにして、最後の数段は答えに要る rank だけを持つ。変換の段は 2 段ずつ 1 回の読み書き
// で行い、区画の初めの段は入力を置く処理に、逆変換の最後の段は答えの取り出しに重ねる。塊を 2^9 個ずつの区画に分けて下の段を区画の
// 中で、上の段は列ごとに作業の配列で行い、列の初めと最後の段は大きい配列から窓の中の値だけを読み書きする。積は [0, p] に直してから
// 12 項までずつ u64 にためて Montgomery で還元し、隣り合う 2 つの出力をまとめて求める (片方の列の値をレジスタで次の項へ渡す)。
// exp、log、合成、power projection は、位置の更新が列の中で閉じるので、上の段を列ごとに全部続けて済ませる (段ごとに大きい配列を
// 読み書きしない)。区画をまたぐ段がある大きさでは、区画の段を済ませた大きい配列の塊を、持っている rank (KC + 5 = 14 個) だけに
// 詰めて持ち (合成を除く)、列ごとの段では今の列を計算している間に次の列で読む塊を L2 へ先に読む。procon-judge の Library Checker の 5 問で書き比べた xw1、xe1、xl1、xc1、xp1 の法を実行時の奇数に広げ、そこから
// 2 段ずつの変換、2 つずつの積和、列ごとの上の段を書き比べて入れたもの。
namespace sps_xor {
using u32= std::uint32_t;
using u64= std::uint64_t;
using V= __m256i;
#define SPS_AI [[gnu::always_inline]] inline
constexpr int KC= 9;     // 区画の塊の数は 2^KC
constexpr int NMAX= 22;  // 速い核で受ける集合の大きさの上限
constexpr int CW= KC + 5;                  // 区画の段を済ませた塊が持つ rank の数
constexpr std::size_t CWS= std::size_t(CW) * 8;  // 詰めた塊の幅 (u32 の数)
// 塊 (n 個の u32) を L2 へ先に読む
SPS_AI void pf_block(const u32* p, std::size_t n) {
 const std::uintptr_t e= std::uintptr_t(p + n);
 for(std::uintptr_t a= std::uintptr_t(p) & ~std::uintptr_t(63); a < e; a+= 64) __builtin_prefetch((const void*)a, 0, 2);
}
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
// 展開しない関数は、受けた MK を入口で手元の変数に写してから使う。__m256i の書き込みは何とでも別名になりうるので、参照のまま
// 使うと、gcc は書き込みのたびに定数をメモリから読み直す (変換の 1 段のループで 1 回に 2 回読んでいた)。
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
// 入力 16 個 (rank は lo + popcount) を置いてレーンの WHT まで済ませた 5 本 (rank lo + 0..4 の値)。lo は塊の上の bit の数 (+ 全体の
// rank のずれ) で、置き場所だけに効く。
SPS_AI void init_vec(V v[5], const u32* in, const MK& k) {
 V lo_= _mm256_loadu_si256((const V*)in), hi= _mm256_loadu_si256((const V*)(in + 8));
 __m256 lf= _mm256_castsi256_ps(lo_), hf= _mm256_castsi256_ps(hi);
 V E= _mm256_permute4x64_epi64(_mm256_castps_si256(_mm256_shuffle_ps(lf, hf, 0x88)), 0xD8);
 V O= _mm256_permute4x64_epi64(_mm256_castps_si256(_mm256_shuffle_ps(lf, hf, 0xDD)), 0xD8);
 const V M0= mpc(0), M1= mpc(1), M2= mpc(2), M3= mpc(3);
 v[0]= _mm256_broadcastd_epi32(_mm256_castsi256_si128(E));
 v[1]= lane_wht(_mm256_or_si256(_mm256_and_si256(E, M1), _mm256_and_si256(O, M0)), k);
 v[2]= lane_wht(_mm256_or_si256(_mm256_and_si256(E, M2), _mm256_and_si256(O, M1)), k);
 v[3]= lane_wht(_mm256_or_si256(_mm256_and_si256(E, M3), _mm256_and_si256(O, M2)), k);
 V o7= _mm256_permutevar8x32_epi32(O, cst(7));
 v[4]= _mm256_blend_epi32(o7, neg2(o7, k), 0x96);
}
SPS_AI void init_block(u32* X, const u32* in, int lo, int sb, int ns, const MK& k) {
 V v[5];
 init_vec(v, in, k);
 for(int j= 0; j < 5; ++j)
  if(int s= lo + j - sb; s >= 0 && s < ns) st(X + 8 * s, v[j]);
}
// 塊 2 つ (q が偶数の X と q + 1 の Y) の入力を置き、bit 0 の段まで済ませる。X は rank [lo, lo+4]、Y は [lo+1, lo+5] の値を持つ
// ので、どちらも [lo, lo+5] になる。lo は q の rank のずれ。
SPS_AI void init_pair(u32* X, u32* Y, const u32* in, int lo, int sb, int ns, const MK& k) {
 V a[5], b[5];
 init_vec(a, in, k), init_vec(b, in + 16, k);
 auto put= [&](int j, V x, V y) __attribute__((always_inline)) {
  if(int s= lo + j - sb; s >= 0 && s < ns) st(X + 8 * s, x), st(Y + 8 * s, y);
 };
 put(0, a[0], a[0]);
 for(int j= 1; j < 5; ++j) put(j, add2(a[j], b[j - 1], k), sub2(a[j], b[j - 1], k));
 put(5, b[4], neg2(b[4], k));
}
// 変換の段。読む塊 (S) と書く塊 (D) を分け、列の写しを初めと最後の段に重ねられるようにする (同じでもよい)。
// W は段の窓の幅で、区画の bit l の段なら W = 4 + l (列の段は l に区画の bit の数を足す)。
// 前向きの 1 段。S0 は rank [lo, lo+W]、S1 は [lo+1, lo+W+1] を持ち、D0 と D1 はどちらも [lo, lo+W+1] になる。
SPS_AI void fwd_pair(const u32* S0, const u32* S1, u32* D0, u32* D1, int lo, int W, int sb, int ns, const MK& k) {
 const int s0= lo - sb;
 if(s0 >= 0 && s0 < ns) {
  const V x= ld(S0 + 8 * s0);
  st(D0 + 8 * s0, x), st(D1 + 8 * s0, x);
 }
 const int a= std::max(s0 + 1, 0), b= std::min(s0 + W, ns - 1);
 for(int s= a; s <= b; ++s) {
  const V x= ld(S0 + 8 * s), y= ld(S1 + 8 * s);
  st(D0 + 8 * s, add2(x, y, k)), st(D1 + 8 * s, sub2(x, y, k));
 }
 if(int s= s0 + W + 1; s >= 0 && s < ns) {
  const V y= ld(S1 + 8 * s);
  st(D0 + 8 * s, y), st(D1 + 8 * s, neg2(y, k));
 }
}
// 逆向きの 1 段。S0 と S1 はどちらも rank [lo, lo+W+1] を持ち、D0 は [lo, lo+W]、D1 は [lo+1, lo+W+1] を残す。
SPS_AI void inv_pair(const u32* S0, const u32* S1, u32* D0, u32* D1, int lo, int W, int sb, int ns, const MK& k) {
 const int s0= lo - sb;
 if(s0 >= 0 && s0 < ns) st(D0 + 8 * s0, add2(ld(S0 + 8 * s0), ld(S1 + 8 * s0), k));
 const int a= std::max(s0 + 1, 0), b= std::min(s0 + W, ns - 1);
 for(int s= a; s <= b; ++s) {
  const V x= ld(S0 + 8 * s), y= ld(S1 + 8 * s);
  st(D0 + 8 * s, add2(x, y, k)), st(D1 + 8 * s, sub2(x, y, k));
 }
 if(int s= s0 + W + 1; s >= 0 && s < ns) st(D1 + 8 * s, sub2(ld(S0 + 8 * s), ld(S1 + 8 * s), k));
}
// 前向きの 2 段 (bit l と l+1、W = 4 + l)。Sm (m の下の bit が段 l、上の bit が段 l+1) は rank [lo+pc(m), lo+pc(m)+W] を持ち、
// Dm はどれも [lo, lo+W+2] になる。窓の端の 2 つずつの rank は、持たない塊の値を 0 として同じ式で求める。
SPS_AI void fwd_quad(const u32* S0, const u32* S1, const u32* S2, const u32* S3, u32* D0, u32* D1, u32* D2, u32* D3, int lo, int W, int sb, int ns, const MK& k) {
 const V z= _mm256_setzero_si256();
 auto bf= [&](int s, V x0, V x1, V x2, V x3) __attribute__((always_inline)) {
  const V a= add2(x0, x1, k), b= sub2(x0, x1, k), c= add2(x2, x3, k), d= sub2(x2, x3, k);
  st(D0 + 8 * s, add2(a, c, k)), st(D1 + 8 * s, add2(b, d, k)), st(D2 + 8 * s, sub2(a, c, k)), st(D3 + 8 * s, sub2(b, d, k));
 };
 auto in= [&](int s) __attribute__((always_inline)) { return s >= 0 && s < ns; };
 const int s0= lo - sb;
 if(in(s0)) {
  const V x= ld(S0 + 8 * s0);
  st(D0 + 8 * s0, x), st(D1 + 8 * s0, x), st(D2 + 8 * s0, x), st(D3 + 8 * s0, x);
 }
 if(const int s= s0 + 1; in(s)) bf(s, ld(S0 + 8 * s), ld(S1 + 8 * s), ld(S2 + 8 * s), z);
 const int a= std::max(s0 + 2, 0), b= std::min(s0 + W, ns - 1);
 for(int s= a; s <= b; ++s) bf(s, ld(S0 + 8 * s), ld(S1 + 8 * s), ld(S2 + 8 * s), ld(S3 + 8 * s));
 if(const int s= s0 + W + 1; in(s)) bf(s, z, ld(S1 + 8 * s), ld(S2 + 8 * s), ld(S3 + 8 * s));
 if(const int s= s0 + W + 2; in(s)) {
  const V x= ld(S3 + 8 * s), y= neg2(x, k);
  st(D0 + 8 * s, x), st(D1 + 8 * s, y), st(D2 + 8 * s, y), st(D3 + 8 * s, x);
 }
}
// 逆向きの 2 段 (bit l+1 のあと l、W = 4 + l)。Sm はどれも rank [lo, lo+W+2] を持ち、Dm は [lo+pc(m), lo+pc(m)+W] を残す。
SPS_AI void inv_quad(const u32* S0, const u32* S1, const u32* S2, const u32* S3, u32* D0, u32* D1, u32* D2, u32* D3, int lo, int W, int sb, int ns, const MK& k) {
 auto in= [&](int s) __attribute__((always_inline)) { return s >= 0 && s < ns; };
 const int s0= lo - sb;
 if(in(s0)) st(D0 + 8 * s0, add2(add2(ld(S0 + 8 * s0), ld(S2 + 8 * s0), k), add2(ld(S1 + 8 * s0), ld(S3 + 8 * s0), k), k));
 auto bf= [&](int s, bool k0, bool k3) __attribute__((always_inline)) {
  const V x0= ld(S0 + 8 * s), x1= ld(S1 + 8 * s), x2= ld(S2 + 8 * s), x3= ld(S3 + 8 * s);
  const V a= add2(x0, x2, k), c= sub2(x0, x2, k), b= add2(x1, x3, k), d= sub2(x1, x3, k);
  if(k0) st(D0 + 8 * s, add2(a, b, k));
  st(D1 + 8 * s, sub2(a, b, k)), st(D2 + 8 * s, add2(c, d, k));
  if(k3) st(D3 + 8 * s, sub2(c, d, k));
 };
 if(const int s= s0 + 1; in(s)) bf(s, true, false);
 const int a= std::max(s0 + 2, 0), b= std::min(s0 + W, ns - 1);
 for(int s= a; s <= b; ++s) bf(s, true, true);
 if(const int s= s0 + W + 1; in(s)) bf(s, false, true);
 if(const int s= s0 + W + 2; in(s)) st(D3 + 8 * s, sub2(sub2(ld(S0 + 8 * s), ld(S2 + 8 * s), k), sub2(ld(S1 + 8 * s), ld(S3 + 8 * s), k), k));
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
// レーンの WHT を済ませた rank c..c+4 の 5 本から、16 個の位置の rank |U| の値を取り出す (extract の後半)。
SPS_AI void pick(const V v[5], V q[2]) {
 const V M0= mpc(0), M1= mpc(1), M2= mpc(2), M3= mpc(3);
 V ev= _mm256_or_si256(_mm256_or_si256(_mm256_and_si256(v[0], M0), _mm256_and_si256(v[1], M1)), _mm256_or_si256(_mm256_and_si256(v[2], M2), _mm256_and_si256(v[3], M3)));
 V od= _mm256_or_si256(_mm256_or_si256(_mm256_and_si256(v[1], M0), _mm256_and_si256(v[2], M1)), _mm256_or_si256(_mm256_and_si256(v[3], M2), _mm256_and_si256(v[4], M3)));
 V ul= _mm256_unpacklo_epi32(ev, od), uh= _mm256_unpackhi_epi32(ev, od);
 q[0]= _mm256_permute2x128_si256(ul, uh, 0x20), q[1]= _mm256_permute2x128_si256(ul, uh, 0x31);
}
// 逆変換の bit 0 の段と取り出しを重ねる。X (q が偶数) と Y (q + 1) はどちらも rank [c, c+5] を持つ (c は q の rank のずれ)。
// X + Y の rank [c, c+4] から q の答えを、X - Y の rank [c+1, c+5] から q + 1 の答えを取り出す。
SPS_AI void extract_pair(const u32* X, const u32* Y, int c, int sb, int ns, V qx[2], V qy[2], const MK& k) {
 const V z= _mm256_setzero_si256();
 V vx[5], vy[5];
 for(int j= 0; j < 6; ++j) {
  const int s= c + j - sb;
  const bool h= s >= 0 && s < ns;
  const V x= h ? ld(X + 8 * s) : z, y= h ? ld(Y + 8 * s) : z;
  if(j < 5) vx[j]= h ? lane_wht(add2(x, y, k), k) : z;
  if(j > 0) vy[j - 1]= h ? lane_wht(sub2(x, y, k), k) : z;
 }
 pick(vx, qx), pick(vy, qy);
}
// 和 sum_{i} A_i B_{d-i} の u64 の和を、12 項までずつ還元して [0, 2p] で返す。Ae/Ao と Be/Bo は [0, p] の値の偶数と奇数のレーン。
// 項は a = a0..a1 (含む) で、B の添字は bo - a。e0 と o0 (p 未満) は初めの 12 項に足す。p < 2^30 なので、12 項 (p^2 以下) と
// p 未満の和の REDC は 4p 未満に収まる。項は 24 個まで (NMAX + 1 以下)。
// 内側のループは、読む場所を 2 本のポインタの増減だけで進める形にする (添字から場所を求め直すと、gcc 15 は 1 項に 4 命令多く出した)。
SPS_AI void dot_acc(V& e, V& o, const u32* xa, const u32* xo, const u32* ya, const u32* yo, int c) {
 for(int t= 0; t < c; ++t) {
  e= _mm256_add_epi64(e, mul64(ld(xa + 8 * t), ld(ya - 8 * t)));
  o= _mm256_add_epi64(o, mul64(ld(xo + 8 * t), ld(yo - 8 * t)));
 }
}
SPS_AI V dot(const u32* Ae, const u32* Ao, const u32* Be, const u32* Bo, int a0, int a1, int bo, V e0, V o0, const MK& k) {
 static_assert(NMAX + 1 <= 24);
 const u32 *xa= Ae + 8 * a0, *xo= Ao + 8 * a0, *ya= Be + 8 * (bo - a0), *yo= Bo + 8 * (bo - a0);
 const int len= a1 - a0 + 1;
 V e= e0, o= o0;
 dot_acc(e, o, xa, xo, ya, yo, len < 12 ? len : 12);
 V r= red2(redc2(e, o, k), k);
 if(len > 12) {
  V e2= _mm256_setzero_si256(), o2= e2;
  dot_acc(e2, o2, xa + 96, xo + 96, ya - 96, yo - 96, len - 12);
  r= add2(r, red2(redc2(e2, o2, k), k), k);
 }
 return r;
}
// 隣り合う 2 つの出力の積和で、出力 0 には X_u Y'_{u+1}、出力 1 には X_u Y'_u を足す (u = 0..c-1)。
// Y' は ya (出力 1 の u = 0 の値) から下へ進む列で、出力 0 の u 項の値は出力 1 の u+1 項の値と同じなので、レジスタで渡す。
// 1 項の読み込みは X の 2 本と Y の 2 本の 4 回で、出力 2 つに足す。
// Y' は ya から 8 c 個下まで読む。呼ぶ側はどれも ya - 8 c が配列の先頭以上になるように渡すが、gcc 15 は 12 項ずつに分けた枝
// (dot2 の c > 11) で c と ya の関係を追えず、exp_step の F の配列に負の添字で触ると誤って -Warray-bounds を出すので、ここだけ止める。
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#endif
SPS_AI void dot2_acc(V& e0, V& o0, V& e1, V& o1, const u32* xa, const u32* xo, const u32* ya, const u32* yo, int c) {
 V pe= ld(ya), po= ld(yo);
 for(int u= 0; u < c; ++u) {
  ya-= 8, yo-= 8;
  const V a= ld(xa), b= ld(xo), ye= ld(ya), yq= ld(yo);
  xa+= 8, xo+= 8;
  e0= _mm256_add_epi64(e0, mul64(a, ye)), o0= _mm256_add_epi64(o0, mul64(b, yq));
  e1= _mm256_add_epi64(e1, mul64(a, pe)), o1= _mm256_add_epi64(o1, mul64(b, po));
  pe= ye, po= yq;
 }
}
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
// 出力 0 = i0 + sum_u X_u Y'_{u+1}、出力 1 = i1 + sum_u X_u Y'_u (u = 0..c-1) に、掛けた 1 項 ex を片方 (EX1 なら出力 1) に足して
// [0, 2p] で返す。i0 と i1 (p 未満) は初めの組に、ex は最後の組に入れる。どの組も 12 項と p 未満の値までなので、c は 23 まで。
template <bool EX1> SPS_AI void dot2(const u32* xa, const u32* xo, const u32* ya, const u32* yo, int c, V i0e, V i0o, V i1e, V i1o, V exe, V exo, V& r0, V& r1, const MK& k) {
 if(c <= 11) {
  if constexpr(EX1) i1e= _mm256_add_epi64(i1e, exe), i1o= _mm256_add_epi64(i1o, exo);
  else i0e= _mm256_add_epi64(i0e, exe), i0o= _mm256_add_epi64(i0o, exo);
  dot2_acc(i0e, i0o, i1e, i1o, xa, xo, ya, yo, c);
  r0= red2(redc2(i0e, i0o, k), k), r1= red2(redc2(i1e, i1o, k), k);
 } else {
  dot2_acc(i0e, i0o, i1e, i1o, xa, xo, ya, yo, 12);
  const V z= _mm256_setzero_si256();
  V f0e= EX1 ? z : exe, f0o= EX1 ? z : exo, f1e= EX1 ? exe : z, f1o= EX1 ? exo : z;
  if(c > 12) dot2_acc(f0e, f0o, f1e, f1o, xa + 96, xo + 96, ya - 96, yo - 96, c - 12);
  r0= add2(red2(redc2(i0e, i0o, k), k), red2(redc2(f0e, f0o, k), k), k);
  r1= add2(red2(redc2(i1e, i1o, k), k), red2(redc2(f1e, f1o, k), k), k);
 }
}
SPS_AI void split_lanes(const u32* X, int cnt, u32* E, u32* O, const MK& k, bool zext) {
 const V lo32= _mm256_set1_epi64x(0xffffffffu);
 for(int i= 0; i < cnt; ++i) {
  V x= norm1(ld(X + 8 * i), k);
  st(E + 8 * i, zext ? _mm256_and_si256(x, lo32) : x), st(O + 8 * i, _mm256_srli_epi64(x, 32));
 }
}

// 塊の列を変換する。in は 2^m 個 (m >= 4、rank は popcount + roff)、X は 2^{m-4} 個の塊 (幅 bs)。区画の中の段だけ行う。
// cmp なら (nb > KC のとき)、区画 h の塊は rank popcount(h) + roff から CW 個を詰めて持つ (sb と ns は使わない)。
inline void fwd_chunks(u32* X, const u32* in, int m, int roff, int sb0, int ns0, std::size_t bs, const MK& k0, bool cmp= false) {
 const MK k= k0;
 const int nb= m - 4, K= std::min(nb, KC);
 const std::size_t cnt= std::size_t(1) << nb, csz= std::size_t(1) << K;
 for(std::size_t q0= 0; q0 < cnt; q0+= csz) {
  const int sb= cmp ? __builtin_popcountll(q0 >> K) + roff : sb0, ns= cmp ? CW : ns0;
  if(K == 0) {
   init_block(X + q0 * bs, in + 16 * q0, __builtin_popcountll(q0) + roff, sb, ns, k);
   continue;
  }
  // bit 0 の段は入力を置く処理に重ね、残りは 2 段ずつ (段の数が奇数なら最後の 1 段だけ 1 段で)
  for(std::size_t q= q0; q < q0 + csz; q+= 2) init_pair(X + q * bs, X + (q + 1) * bs, in + 16 * q, __builtin_popcountll(q) + roff, sb, ns, k);
  int kk= 1;
  for(; kk + 1 < K; kk+= 2) {
   const std::size_t d1= std::size_t(1) << kk, d2= d1 << 1;
   for(std::size_t qa= q0; qa < q0 + csz; qa+= d2 << 1) {
    const int lo= __builtin_popcountll(qa >> (kk + 2)) + roff;
    for(std::size_t q= qa; q < qa + d1; ++q) {
     u32 *x0= X + q * bs, *x1= X + (q + d1) * bs, *x2= X + (q + d2) * bs, *x3= X + (q + d1 + d2) * bs;
     fwd_quad(x0, x1, x2, x3, x0, x1, x2, x3, lo, 4 + kk, sb, ns, k);
    }
   }
  }
  if(kk < K)
   for(std::size_t q= q0; q < q0 + csz; ++q)
    if(!(q >> kk & 1)) {
     u32 *x= X + q * bs, *y= X + (q | std::size_t(1) << kk) * bs;
     fwd_pair(x, y, x, y, __builtin_popcountll(q >> (kk + 1)) + roff, 4 + kk, sb, ns, k);
    }
 }
}
// 列 j の塊 (h << K) + j を c に写して、上の段 (区画をまたぐ段) を変換する。
// 初めの 1 段か 2 段は、X から窓の中の値だけを読んで c に書く (写しの段を重ねる)。
// bsx は X の塊の幅 (0 なら bs)。cmp なら X の塊は詰めて持つので、場所 s = r - sb で詰めた場所 r - (popcount(h) + roff) を読むよう
// 塊の先頭をずらす。
inline void fwd_column(u32* c, const u32* X, std::size_t j, int m, int roff, int sb, int ns, std::size_t bs, const MK& k0, std::size_t bsx= 0, bool cmp= false) {
 const MK k= k0;
 const int nb= m - 4, K= std::min(nb, KC), H= nb - K;
 const std::size_t hs= std::size_t(1) << H;
 if(!bsx) bsx= bs;
 auto xb= [&](std::size_t h) __attribute__((always_inline)) { return X + ((h << K) + j) * bsx + (cmp ? 8 * (std::ptrdiff_t(sb) - __builtin_popcountll(h) - roff) : 0); };
 if(H == 1) {
  fwd_pair(xb(0), xb(1), c, c + bs, roff, 4 + K, sb, ns, k);
  return;
 }
 for(std::size_t h= 0; h < hs; h+= 4) fwd_quad(xb(h), xb(h + 1), xb(h + 2), xb(h + 3), c + h * bs, c + (h + 1) * bs, c + (h + 2) * bs, c + (h + 3) * bs, __builtin_popcountll(h >> 2) + roff, 4 + K, sb, ns, k);
 int kk= 2;
 for(; kk + 1 < H; kk+= 2) {
  const std::size_t d1= std::size_t(1) << kk, d2= d1 << 1;
  for(std::size_t ha= 0; ha < hs; ha+= d2 << 1) {
   const int lo= __builtin_popcountll(ha >> (kk + 2)) + roff;
   for(std::size_t h= ha; h < ha + d1; ++h) {
    u32 *x0= c + h * bs, *x1= c + (h + d1) * bs, *x2= c + (h + d2) * bs, *x3= c + (h + d1 + d2) * bs;
    fwd_quad(x0, x1, x2, x3, x0, x1, x2, x3, lo, 4 + K + kk, sb, ns, k);
   }
  }
 }
 if(kk < H)
  for(std::size_t h= 0; h < hs; ++h)
   if(!(h >> kk & 1)) {
    u32 *x= c + h * bs, *y= c + (h | std::size_t(1) << kk) * bs;
    fwd_pair(x, y, x, y, __builtin_popcountll(h >> (kk + 1)) + roff, 4 + K + kk, sb, ns, k);
   }
}
// 列の上の段の逆変換。最後の 1 段か 2 段は、c から読んで残す窓だけを X に書く (写しの段を重ねる)。
// bsx と cmp は fwd_column と同じ。cmp なら、詰めた塊のうち切り捨てる rank (sb より下と sb + ns 以上) の場所に 0 を書く。
inline void inv_column(u32* c, u32* X, std::size_t j, int m, int roff, int sb, int ns, std::size_t bs, const MK& k0, std::size_t bsx= 0, bool cmp= false) {
 const MK k= k0;
 const int nb= m - 4, K= std::min(nb, KC), H= nb - K;
 const std::size_t hs= std::size_t(1) << H;
 if(!bsx) bsx= bs;
 auto xb= [&](std::size_t h) __attribute__((always_inline)) { return X + ((h << K) + j) * bsx + (cmp ? 8 * (std::ptrdiff_t(sb) - __builtin_popcountll(h) - roff) : 0); };
 if(cmp) {
  const V z= _mm256_setzero_si256();
  for(std::size_t h= 0; h < hs; ++h) {
   u32* x= X + ((h << K) + j) * bsx;
   const int base= __builtin_popcountll(h) + roff;
   for(int s= 0; s < sb - base && s < CW; ++s) st(x + 8 * s, z);
   for(int s= std::max(sb + ns - base, 0); s < CW; ++s) st(x + 8 * s, z);
  }
 }
 if(H == 1) {
  inv_pair(c, c + bs, xb(0), xb(1), roff, 4 + K, sb, ns, k);
  return;
 }
 // 段 H-1..2 を c の中で (数が奇数なら初めの 1 段だけ 1 段で)、段 1 と 0 を X へ書きながら
 int kk= H - 1;
 if(H & 1) {
  for(std::size_t h= 0; h < hs; ++h)
   if(!(h >> kk & 1)) {
    u32 *x= c + h * bs, *y= c + (h | std::size_t(1) << kk) * bs;
    inv_pair(x, y, x, y, __builtin_popcountll(h >> (kk + 1)) + roff, 4 + K + kk, sb, ns, k);
   }
  --kk;
 }
 for(; kk >= 3; kk-= 2) {
  const int l= kk - 1;
  const std::size_t d1= std::size_t(1) << l, d2= d1 << 1;
  for(std::size_t ha= 0; ha < hs; ha+= d2 << 1) {
   const int lo= __builtin_popcountll(ha >> (l + 2)) + roff;
   for(std::size_t h= ha; h < ha + d1; ++h) {
    u32 *x0= c + h * bs, *x1= c + (h + d1) * bs, *x2= c + (h + d2) * bs, *x3= c + (h + d1 + d2) * bs;
    inv_quad(x0, x1, x2, x3, x0, x1, x2, x3, lo, 4 + K + l, sb, ns, k);
   }
  }
 }
 for(std::size_t h= 0; h < hs; h+= 4) inv_quad(c + h * bs, c + (h + 1) * bs, c + (h + 2) * bs, c + (h + 3) * bs, xb(h), xb(h + 1), xb(h + 2), xb(h + 3), __builtin_popcountll(h >> 2) + roff, 4 + K, sb, ns, k);
}
// 逆変換の区画の中の段と、塊ごとの取り出し。f(q, 取り出した 2 本) を呼ぶ。
// 段 K-1..1 を 2 段ずつ (数が奇数なら初めの 1 段だけ 1 段で) 行い、bit 0 の段は取り出しに重ねる。
// cmp なら、区画 h の塊は rank popcount(h) + roff から CW 個を詰めて持つ (切り捨てた rank の場所は 0)。
template <class Out> inline void inv_chunks(u32* X, int m, int roff, int sb0, int ns0, std::size_t bs, const MK& k0, const Out& out, bool cmp= false) {
 const MK k= k0;
 const int nb= m - 4, K= std::min(nb, KC);
 const std::size_t cnt= std::size_t(1) << nb, csz= std::size_t(1) << K;
 for(std::size_t q0= 0; q0 < cnt; q0+= csz) {
  const int sb= cmp ? __builtin_popcountll(q0 >> K) + roff : sb0, ns= cmp ? CW : ns0;
  if(K == 0) {
   V v[2];
   extract(X + q0 * bs, __builtin_popcountll(q0) + roff, sb, ns, v, k);
   out(q0, v);
   continue;
  }
  int kk= K - 1;
  if((K - 1) & 1) {
   for(std::size_t q= q0; q < q0 + csz; ++q)
    if(!(q >> kk & 1)) {
     u32 *x= X + q * bs, *y= X + (q | std::size_t(1) << kk) * bs;
     inv_pair(x, y, x, y, __builtin_popcountll(q >> (kk + 1)) + roff, 4 + kk, sb, ns, k);
    }
   --kk;
  }
  for(; kk >= 2; kk-= 2) {
   const int l= kk - 1;
   const std::size_t d1= std::size_t(1) << l, d2= d1 << 1;
   for(std::size_t qa= q0; qa < q0 + csz; qa+= d2 << 1) {
    const int lo= __builtin_popcountll(qa >> (l + 2)) + roff;
    for(std::size_t q= qa; q < qa + d1; ++q) {
     u32 *x0= X + q * bs, *x1= X + (q + d1) * bs, *x2= X + (q + d2) * bs, *x3= X + (q + d1 + d2) * bs;
     inv_quad(x0, x1, x2, x3, x0, x1, x2, x3, lo, 4 + l, sb, ns, k);
    }
   }
  }
  for(std::size_t q= q0; q < q0 + csz; q+= 2) {
   V v0[2], v1[2];
   extract_pair(X + q * bs, X + (q + 1) * bs, __builtin_popcountll(q) + roff, sb, ns, v0, v1, k);
   out(q, v0), out(q + 1, v1);
  }
 }
}
inline std::size_t col_blocks(int m) {
 const int nb= m - 4, K= std::min(nb, KC);
 return std::size_t(1) << (nb - K);
}
// f̂_i (rank 0..i、roff = sb = 0) の列 j の上の段 (i - 4 - KC 段) の前向きと逆向きの変換。0 段なら塊を写すだけ (fwd_column と
// inv_column は 1 段以上を受ける)。
// bsx は X の塊の幅。0 段の塊は rank 0..KC+4 の CW 個なので、詰めても詰めなくても同じ写しで済む。
inline void col_fwd(u32* c, const u32* X, std::size_t j, int i, std::size_t bs, const MK& k, std::size_t bsx, bool cmp) {
 if(i - 4 == KC) std::memcpy(c, X + j * bsx, std::size_t(i + 1) * 32);
 else fwd_column(c, X, j, i, 0, 0, i + 1, bs, k, bsx, cmp);
}
inline void col_inv(u32* c, u32* X, std::size_t j, int i, std::size_t bs, const MK& k, std::size_t bsx, bool cmp) {
 if(i - 4 == KC) std::memcpy(X + j * bsx, c, std::size_t(i + 1) * 32);
 else inv_column(c, X, j, i, 0, 0, i + 1, bs, k, bsx, cmp);
}
// X の塊 q (rank の起点 base) のうち、列の初めの段で読む窓 (rank [base, base + KC + 4] と [sb, sb + ns) の重なり) を先に読む。
SPS_AI void pf_win(const u32* X, std::size_t q, int base, int sb, int ns, std::size_t bsx, bool cmp) {
 if(cmp) pf_block(X + q * bsx, CWS);
 else {
  const int a= std::max(base - sb, 0), b= std::min(base - sb + CW, ns);
  if(a < b) pf_block(X + q * bsx + 8 * a, std::size_t(b - a) * 8);
 }
}
// 値 x (2p 未満) を R 倍して [0, p) に (n 個、8 の倍数)。
inline void scale_r(const u32* x, u32* y, std::size_t n, const MK& k0) {
 const MK k= k0;
 const V r2= cst(k.R2);
 for(std::size_t s= 0; s < n; s+= 8) _mm256_store_si256((V*)(y + s), norm1(mont(_mm256_loadu_si256((const V*)(x + s)), r2, k), k));
}

// ---- 部分集合畳み込み ----
// 塊ごとの rank の積 (長さ nm1 で切り詰める)。入力の rank r は場所 r - 1、積の rank r は場所 r - 2。結果は a に書く。
template <int NM1> SPS_AI void conv_block(u32* a, const u32* b, const MK& k) {
 alignas(32) u32 Ae[NM1 * 8], Ao[NM1 * 8], Be[NM1 * 8], Bo[NM1 * 8];
 split_lanes(a, NM1, Ae, Ao, k, false), split_lanes(b, NM1, Be, Bo, k, false);
 const V z= _mm256_setzero_si256();
 // 出力 s と s+1 をまとめる。t = 0..s の項を分け合い、出力 s+1 には t = s+1 の項 A_{s+1} B_0 を足す。
 int s= 0;
 for(; s + 1 < NM1; s+= 2) {
  V r0, r1;
  dot2<true>(Ae, Ao, Be + 8 * (s + 1), Bo + 8 * (s + 1), s + 1, z, z, z, z, mul64(ld(Ae + 8 * (s + 1)), ld(Be)), mul64(ld(Ao + 8 * (s + 1)), ld(Bo)), r0, r1, k);
  st(a + 8 * s, r0), st(a + 8 * (s + 1), r1);
 }
 if(s < NM1) st(a + 8 * s, dot(Ae, Ao, Be, Bo, 0, s, s, z, z, k));
}
template <int N> void convolve_n(const u32* f, const u32* g, u32* out, const MK& k0) {
 const MK k= k0;
 constexpr int NM1= N - 1, NB= N - 4, BS= NM1 * 8;
 constexpr bool CM= NB > KC;
 constexpr std::size_t XS= CM ? CWS : BS;
 const std::size_t nblk= std::size_t(1) << NB;
 u32 *A= huge_alloc(nblk * XS * 4), *B= huge_alloc(nblk * XS * 4);
 fwd_chunks(A, f, N, 0, 1, NM1, XS, k, CM), fwd_chunks(B, g, N, 0, 1, NM1, XS, k, CM);
 // 空集合は入力から外す (rank 0 は場所 -1 なので置かれない)
 if constexpr(NB <= KC) {
  for(std::size_t q= 0; q < nblk; ++q) conv_block<NM1>(A + q * BS, B + q * BS, k);
 } else {
  const std::size_t hs= col_blocks(N), K= KC;
  u32 *cA= small_alloc(hs * BS * 4), *cB= small_alloc(hs * BS * 4);
  for(std::size_t j= 0; j < (std::size_t(1) << K); ++j) {
   fwd_column(cA, A, j, N, 0, 1, NM1, BS, k, XS, CM), fwd_column(cB, B, j, N, 0, 1, NM1, BS, k, XS, CM);
   const bool nx= j + 1 < (std::size_t(1) << K);
   for(std::size_t h= 0; h < hs; ++h) {
    if(nx) {
     const std::size_t q= (h << K) + j + 1;
     const int base= __builtin_popcountll(h);
     pf_win(A, q, base, 1, NM1, XS, CM), pf_win(B, q, base, 1, NM1, XS, CM);
    }
    conv_block<NM1>(cA + h * BS, cB + h * BS, k);
   }
   inv_column(cA, A, j, N, 0, 2, NM1, BS, k, XS, CM);
  }
  std::free(cA), std::free(cB);
 }
 const V Kc= cst(u32(u64(k.inv2pow(N - 1)) * k.R2 % k.P)), F0= cst(u32(u64(f[0] % k.P) * k.R1 % k.P)), G0= cst(u32(u64(g[0] % k.P) * k.R1 % k.P));
 inv_chunks(A, N, 0, 2, NM1, XS, k, [&](std::size_t q, const V* v) {
  for(int t= 0; t < 2; ++t) {
   const u32 *fp= f + 16 * q + 8 * t, *gp= g + 16 * q + 8 * t;
   V fv= norm1(_mm256_loadu_si256((const V*)fp), k), gv= norm1(_mm256_loadu_si256((const V*)gp), k);
   V e= _mm256_add_epi64(_mm256_add_epi64(mul64(v[t], Kc), mul64(gv, F0)), mul64(fv, G0));
   V o= _mm256_add_epi64(_mm256_add_epi64(mul64(_mm256_srli_epi64(v[t], 32), Kc), mul64(_mm256_srli_epi64(gv, 32), F0)), mul64(_mm256_srli_epi64(fv, 32), G0));
   _mm256_storeu_si256((V*)(out + 16 * q + 8 * t), norm1(redc2(e, o, k), k));
  }
 }, CM);
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
 auto put= [&](int d, V D) __attribute__((always_inline)) {
  V g= d <= i ? ld(G + 8 * (d - 1)) : z;
  st(Ghi + 8 * (d - 1), sub2(g, D, k)), st(G + 8 * (d - 1), add2(g, D, k));
 };
 // D_d と D_{d+1} をまとめる。a = 1..d-1 の項を分け合い、D_{d+1} には a = d の項 G_d F_0 を足す (初めの値は F_{d-1} と F_d)。
 int d= 1;
 for(; d + 1 <= i + 1; d+= 2) {
  V r0, r1;
  dot2<true>(Ge + 8, Go + 8, Fe + 8 * (d - 1), Fo + 8 * (d - 1), d - 1, ld(Fe + 8 * (d - 1)), ld(Fo + 8 * (d - 1)), ld(Fe + 8 * d), ld(Fo + 8 * d), mul64(ld(Ge + 8 * d), ld(Fe)), mul64(ld(Go + 8 * d), ld(Fo)), r0, r1, k);
  put(d, r0), put(d + 1, r1);
 }
 if(d <= i + 1) put(d, dot(Ge, Go, Fe, Fo, 1, d - 1, d - 1, ld(Fe + 8 * (d - 1)), ld(Fo + 8 * (d - 1)), k));
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
// exp (b[0] は 0)。下の段 (i - 4 < KC) は G の 1 区画の中で行い、区画をまたぐ上の段は列 j ごとに全部続けて作業の配列の中で済ませる。
// G の列 j の塊 (h << KC) + j のうち、h = 2^H - 1 に下の段の G を、h ∈ [2^t - 1, 2^{t+1} - 1) に f̂_{KC+4+t} の区画の変換を置き
// (どちらも塊の幅は G と同じ)、列 j を済ませると逆変換の上の段の結果を列 j の全部の塊に書く。列 j の塊は列 j を済ませるときにだけ
// 読み書きするので、大きい配列は f̂_i の区画の変換と最後の答えの区画の変換を 1 回ずつ書いて読むだけになる。
template <int N> void exp_n(const u32* b, u32* out, const MK& k0) {
 const MK k= k0;
 constexpr int NS= N, BS= NS * 8, NB= N - 4, KL= NB < KC ? NB : KC, H= NB - KL;
 constexpr bool CM= H > 0;
 constexpr std::size_t CB= std::size_t(1) << KL, XS= CM ? CWS : BS;
 u32 *G= huge_alloc((std::size_t(1) << NB) * XS * 4), *bs= huge_alloc(std::size_t(4) << N);
 scale_r(b, bs, std::size_t(1) << N, k);
 u32* Glo= G + ((std::size_t(1) << H) - 1) * CB * XS;
 {
  alignas(32) u32 h[16], bb[16];
  for(int s= 0; s < 16; ++s) bb[s]= b[s] % k.P;
  naive_exp(4, bb, h, k.P);
  init_block(Glo, h, 0, 1, NS, k);
 }
 // 下の段 (区画の中で閉じる)
 u32* Fs= small_alloc((CB >> 1) * BS * 4);
 for(int i= 4; i < KL + 4; ++i) {
  const std::size_t cnt= std::size_t(1) << (i - 4);
  fwd_chunks(Fs, bs + (std::size_t(1) << i), i, 0, 0, i + 1, BS, k);
  for(std::size_t q= 0; q < cnt; ++q) exp_step(Glo + q * XS, Glo + (q + cnt) * XS, Fs + q * BS, i, k);
 }
 std::free(Fs);
 // 上の段。f̂_i の区画の変換を G の空いた塊に置き、列ごとに全部の段を済ませて逆変換の上の段の結果を G に書く
 if constexpr(H > 0) {
  auto Fi= [&](int i) { return G + ((std::size_t(1) << (i - 4 - KC)) - 1) * CB * XS; };
  for(int i= KC + 4; i < N; ++i) fwd_chunks(Fi(i), bs + (std::size_t(1) << i), i, 0, 0, i + 1, XS, k, CM);
  u32 *cG= small_alloc((std::size_t(1) << H) * BS * 4), *cF= small_alloc((std::size_t(1) << (H - 1)) * BS * 4);
  for(std::size_t j= 0; j < CB; ++j) {
   std::memcpy(cG, Glo + j * XS, std::min<std::size_t>(XS, BS) * 4);
   const bool nx= j + 1 < CB;
   if(nx) pf_block(Glo + (j + 1) * XS, XS);
   for(int i= KC + 4; i < N; ++i) {
    const std::size_t hs= std::size_t(1) << (i - 4 - KC);
    col_fwd(cF, Fi(i), j, i, BS, k, XS, CM);
    for(std::size_t h= 0; h < hs; ++h) {
     if(nx) pf_win(Fi(i), (h << KC) + j + 1, __builtin_popcountll(h), 0, i + 1, XS, CM);
     exp_step(cG + h * BS, cG + (h + hs) * BS, cF + h * BS, i, k);
    }
   }
   inv_column(cG, G, j, N, 0, 1, NS, BS, k, XS, CM);
  }
  std::free(cG), std::free(cF);
 }
 const V Kc= cst(u32(u64(k.inv2pow(N - 1)) * k.R1 % k.P));
 inv_chunks(G, N, 0, 1, NS, XS, k, [&](std::size_t q, const V* v) {
  for(int t= 0; t < 2; ++t) _mm256_storeu_si256((V*)(out + 16 * q + 8 * t), norm1(mont(v[t], Kc, k), k));
 }, CM);
 out[0]= 1;
 std::free(G), std::free(bs);
}

// ---- 割り算 (log と逆元と割り算の共通) ----
// Q_0 = A_0、Q_d = A_d - sum_{j=1}^{d} B_j Q_{d-j} (d = 1..D) を Q に書く (Q は A と同じ場所でもよい)。Be/Bo は B_1.. を 8 から置いた
// 偶数と奇数のレーン、Qe/Qo は Q の偶数と奇数のレーンの置き場。
SPS_AI void div_core(const u32* A, u32* Q, const u32* Be, const u32* Bo, u32* Qe, u32* Qo, int D, const MK& k) {
 const V z= _mm256_setzero_si256();
 st(Q, ld(A));
 split_lanes(Q, 1, Qe, Qo, k, false);
 auto fin= [&](int d, V S) __attribute__((always_inline)) {
  st(Q + 8 * d, sub2(ld(A + 8 * d), S, k));
  split_lanes(Q + 8 * d, 1, Qe + 8 * d, Qo + 8 * d, k, false);
 };
 // S_d = sum_{j<d} Q_j B_{d-j} と、S_{d+1} から Q_d の項を除いた sum_{j<d} Q_j B_{d+1-j} を Q_0..Q_{d-1} で分け合い、Q_d を求めてから
 // 後者の最後の組に B_1 Q_d を足して還元する。
 int d= 1;
 for(; d + 1 <= D; d+= 2) {
  const bool two= d >= 12;
  const u32 *ya= Be + 8 * (d + 1), *yo= Bo + 8 * (d + 1);
  V a0e= z, a0o= z, a1e= z, a1o= z, g0= z, g1= z;
  dot2_acc(a0e, a0o, a1e, a1o, Qe, Qo, ya, yo, two ? 12 : d);
  if(two) {
   g0= red2(redc2(a0e, a0o, k), k), g1= red2(redc2(a1e, a1o, k), k);
   a0e= a0o= a1e= a1o= z;
   if(d > 12) dot2_acc(a0e, a0o, a1e, a1o, Qe + 96, Qo + 96, ya - 96, yo - 96, d - 12);
  }
  V S0= red2(redc2(a0e, a0o, k), k);
  if(two) S0= add2(S0, g0, k);
  fin(d, S0);
  a1e= _mm256_add_epi64(a1e, mul64(ld(Be + 8), ld(Qe + 8 * d))), a1o= _mm256_add_epi64(a1o, mul64(ld(Bo + 8), ld(Qo + 8 * d)));
  V S1= red2(redc2(a1e, a1o, k), k);
  if(two) S1= add2(S1, g1, k);
  fin(d + 1, S1);
 }
 if(d <= D) fin(d, dot(Be, Bo, Qe, Qo, 1, d, d, z, z, k));
}
// A = Â' (R 倍、rank 0..i を場所 0..i)、B = B̂' (R 倍、rank 1..i を場所 0..i-1、rank 0 は 1)。Q' = A / B を rank i まで求めて A に書く。
// upd なら B を B̂'[y +] = B + t A に、上の半分 Bhi を B - t A にする (rank 1..i+1)。
SPS_AI void div_step(u32* A, u32* B, u32* Bhi, int i, bool upd, const MK& k) {
 alignas(32) u32 Be[(NMAX + 1) * 8], Bo[(NMAX + 1) * 8], Qe[(NMAX + 1) * 8], Qo[(NMAX + 1) * 8], Q[(NMAX + 1) * 8];
 split_lanes(B, i, Be + 8, Bo + 8, k, false);
 const V z= _mm256_setzero_si256();
 div_core(A, Q, Be, Bo, Qe, Qo, i, k);
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
// log (h[0] は 1)。大きい配列 Fb の列 j の塊のうち、h = 2^H - 1 に下の段の B (幅は Fb と同じ) を、h ∈ [2^t - 1, 2^{t+1} - 1) に段
// KC+4+t の f̂ の区画の変換を置く。上の段は列 j ごとに、B の列を作業の配列で伸ばしながら商を求め、商の列を逆変換の上の段まで
// 済ませて Fb の同じ場所に戻す。列を全部済ませてから、上の段ごとに区画の逆変換をして答えを取り出す。
template <int N> void log_n(const u32* h, u32* out, const MK& k0) {
 const MK k= k0;
 constexpr int NSB= N - 1, BSB= NSB * 8, NSF= N, BSF= NSF * 8, NB= N - 4, KL= NB < KC ? NB : KC, H= NB - KL;
 constexpr bool CM= H > 0;
 constexpr std::size_t CB= std::size_t(1) << KL, XS= CM ? CWS : BSF;
 u32 *Fb= huge_alloc((std::size_t(1) << NB) * XS * 4), *hs= huge_alloc(std::size_t(4) << N);
 scale_r(h, hs, std::size_t(1) << N, k);
 {
  u32 hh[16];
  for(int s= 0; s < 16; ++s) hh[s]= h[s] % k.P;
  naive_log(4, hh, out, k.P);
 }
 u32* Blo= Fb + ((std::size_t(1) << H) - 1) * CB * XS;
 init_block(Blo, hs, 0, 1, NSB, k);
 // 段 i の商の区画の逆変換をして、答えの [2^i, 2^{i+1}) を取り出す
 auto put= [&](int i, u32* F, std::size_t bs, bool cmp) {
  const V Kc= cst(k.inv2pow(i - 1));
  u32* o= out + (std::size_t(1) << i);
  inv_chunks(F, i, 0, 0, i + 1, bs, k, [&](std::size_t q, const V* v) {
   for(int t= 0; t < 2; ++t) _mm256_storeu_si256((V*)(o + 16 * q + 8 * t), norm1(mont(v[t], Kc, k), k));
  }, cmp);
 };
 // 下の段 (区画の中で閉じる)。B の塊の幅は XS
 u32* Fs= small_alloc((CB >> 1) * BSF * 4);
 for(int i= 4; i < KL + 4; ++i) {
  const std::size_t cnt= std::size_t(1) << (i - 4);
  fwd_chunks(Fs, hs + (std::size_t(1) << i), i, 0, 0, i + 1, BSF, k);
  for(std::size_t q= 0; q < cnt; ++q) div_step(Fs + q * BSF, Blo + q * XS, Blo + (q + cnt) * XS, i, i + 1 < N, k);
  put(i, Fs, BSF, false);
 }
 std::free(Fs);
 if constexpr(H > 0) {
  auto Fi= [&](int i) { return Fb + ((std::size_t(1) << (i - 4 - KC)) - 1) * CB * XS; };
  for(int i= KC + 4; i < N; ++i) fwd_chunks(Fi(i), hs + (std::size_t(1) << i), i, 0, 0, i + 1, XS, k, CM);
  u32 *cB= small_alloc((std::size_t(1) << (H - 1)) * BSB * 4), *cF= small_alloc((std::size_t(1) << (H - 1)) * BSF * 4);
  for(std::size_t j= 0; j < CB; ++j) {
   std::memcpy(cB, Blo + j * XS, std::min<std::size_t>(XS, BSB) * 4);
   const bool nx= j + 1 < CB;
   if(nx) pf_block(Blo + (j + 1) * XS, std::min<std::size_t>(XS, BSB));
   for(int i= KC + 4; i < N; ++i) {
    const std::size_t hc= std::size_t(1) << (i - 4 - KC);
    col_fwd(cF, Fi(i), j, i, BSF, k, XS, CM);
    for(std::size_t t= 0; t < hc; ++t) {
     if(nx) pf_win(Fi(i), (t << KC) + j + 1, __builtin_popcountll(t), 0, i + 1, XS, CM);
     div_step(cF + t * BSF, cB + t * BSB, cB + (t + hc) * BSB, i, i + 1 < N, k);
    }
    col_inv(cF, Fi(i), j, i, BSF, k, XS, CM);
   }
  }
  std::free(cB), std::free(cF);
  for(int i= KC + 4; i < N; ++i) put(i, Fi(i), XS, CM);
 }
 std::free(Fb), std::free(hs);
}
// 割り算 f / g (g[0] は 1)。変換の領域の各点で割って 1 回だけ逆変換する。
template <int N> void div_n(const u32* f, const u32* g, u32* out, const MK& k0) {
 const MK k= k0;
 constexpr int NS= N + 1, BS= NS * 8, NB= N - 4;
 constexpr bool CM= NB > KC;
 constexpr std::size_t XS= CM ? CWS : BS;
 const std::size_t nblk= std::size_t(1) << NB;
 u32 *A= huge_alloc(nblk * XS * 4), *B= huge_alloc(nblk * XS * 4), *gs= huge_alloc(std::size_t(4) << N), *fs= huge_alloc(std::size_t(4) << N);
 scale_r(g, gs, std::size_t(1) << N, k), scale_r(f, fs, std::size_t(1) << N, k);
 // A は f の R 倍 (rank 0..N)、B は g の R 倍 (rank 1..N を場所 0..N-1)
 fwd_chunks(A, fs, N, 0, 0, NS, XS, k, CM), fwd_chunks(B, gs, N, 0, 1, N, XS, k, CM);
 auto body= [&](u32* a, u32* b) __attribute__((always_inline)) {
  // Q' = A' / B' (B の rank 0 は R、R 倍のままの割り算は Q'_d = A'_d - REDC(sum B'_j Q'_{d-j}))
  alignas(32) u32 Be[(NMAX + 2) * 8], Bo[(NMAX + 2) * 8], Qe[(NMAX + 2) * 8], Qo[(NMAX + 2) * 8];
  split_lanes(b, N, Be + 8, Bo + 8, k, false);
  div_core(a, a, Be, Bo, Qe, Qo, N, k);
 };
 if constexpr(NB <= KC) {
  for(std::size_t q= 0; q < nblk; ++q) body(A + q * BS, B + q * BS);
 } else {
  const std::size_t hs= col_blocks(N);
  u32 *cA= small_alloc(hs * BS * 4), *cB= small_alloc(hs * BS * 4);
  for(std::size_t j= 0; j < (std::size_t(1) << KC); ++j) {
   fwd_column(cA, A, j, N, 0, 0, NS, BS, k, XS, CM), fwd_column(cB, B, j, N, 0, 1, N, BS, k, XS, CM);
   const bool nx= j + 1 < (std::size_t(1) << KC);
   for(std::size_t h= 0; h < hs; ++h) {
    if(nx) {
     const std::size_t q= (h << KC) + j + 1;
     const int base= __builtin_popcountll(h);
     pf_win(A, q, base, 0, NS, XS, CM), pf_win(B, q, base, 1, N, XS, CM);
    }
    body(cA + h * BS, cB + h * BS);
   }
   inv_column(cA, A, j, N, 0, 0, NS, BS, k, XS, CM);
  }
  std::free(cA), std::free(cB);
 }
 const V Kc= cst(k.inv2pow(N - 1));
 inv_chunks(A, N, 0, 0, NS, XS, k, [&](std::size_t q, const V* v) {
  for(int t= 0; t < 2; ++t) _mm256_storeu_si256((V*)(out + 16 * q + 8 * t), norm1(mont(v[t], Kc, k), k));
 }, CM);
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
  auto put= [&](int d, V D) __attribute__((always_inline)) {
   V g= d <= i ? ld(X + 8 * d) : z;
   st(Y + 8 * d, sub2(g, D, k)), st(X + 8 * d, add2(g, D, k));
  };
  // D_d = sum_{a<d} A_a F_{d-1-a} と D_{d+1} をまとめる。a = 0..d-1 の項を分け合い、D_{d+1} には a = d の項 A_d F_0 を足す。
  int d= 1;
  for(; d + 1 <= i + 1; d+= 2) {
   V r0, r1;
   dot2<true>(Ae, Ao, Fe + 8 * d, Fo + 8 * d, d, z, z, z, z, mul64(ld(Ae + 8 * d), ld(Fe)), mul64(ld(Ao + 8 * d), ld(Fo)), r0, r1, k);
   put(d, r0), put(d + 1, r1);
  }
  if(d <= i + 1) put(d, dot(Ae, Ao, Fe, Fo, 0, d - 1, d - 1, z, z, k));
 }
}
// sum_k e_k g^k / k! (g[0] は 0、e は n + 1 個)。G_0 の大きい配列の列 j の塊のうち、h = 2^H - 1 に下の段の G_0 を、
// h ∈ [2^t - 1, 2^{t+1} - 1) に ĝ_{KC+4+t} の区画の変換を置く。G_j (j >= 1) の下の段は別の配列 (塊 2^KC 個まで)。上の段は列 j ごとに、
// G_0..G_H の列を作業の配列で伸ばして済ませ、G_0 の列の逆変換の上の段の結果を G_0 に書く。
template <int N> void egf_comp_n(const u32* g, const u32* e, u32* out, const MK& k0) {
 const MK k= k0;
 constexpr int NB= N - 4, KL= NB < KC ? NB : KC, H= NB - KL, NS= N + 1;
 // 合成は詰めて持たない (procon-judge の run 44279978233 の Xeon 6975P-C で、詰めると 4 % から 6 % 遅かった。exp、log、power projection は
 // 縮むか変わらなかった)
 constexpr bool CM= false;
 constexpr std::size_t CB= std::size_t(1) << KL, BS= NS * 8, BSF= N * 8, XS= CM ? CWS : BS;
 u32 *gs= huge_alloc(std::size_t(4) << N), *G0= huge_alloc((std::size_t(1) << NB) * XS * 4);
 scale_r(g, gs, std::size_t(1) << N, k);
 u32* G[N - 3];
 std::size_t bsj[N - 3], lsj[N - 3];  // bsj は列の作業の配列の塊の幅、lsj は下の段の配列の塊の幅 (G_0 は XS)
 for(int j= 0; j <= N - 4; ++j) {
  bsj[j]= std::size_t(N - j + 1) * 8, lsj[j]= j == 0 ? XS : bsj[j];
  G[j]= j == 0 ? G0 + ((std::size_t(1) << H) - 1) * CB * XS : huge_alloc(std::min(std::size_t(1) << (N - j - 4), CB) * bsj[j] * 4);
 }
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
 u32 *gp[N - 3], *gh[N - 3];
 // 下の段 (区画の中で閉じる)
 u32* Fs= small_alloc((CB >> 1) * BSF * 4);
 for(int i= 4; i < KL + 4; ++i) {
  const int J= N - i;
  const std::size_t cnt= std::size_t(1) << (i - 4);
  fwd_chunks(Fs, gs + (std::size_t(1) << i), i, 0, 0, i + 1, BSF, k);
  for(std::size_t q= 0; q < cnt; ++q) {
   for(int j= 0; j <= J; ++j) gp[j]= G[j] + q * lsj[j];
   for(int j= 0; j < J; ++j) gh[j]= G[j] + (q + cnt) * lsj[j];
   comp_step(gp, gh, Fs + q * BSF, i, J, k);
  }
 }
 std::free(Fs);
 if constexpr(H > 0) {
  auto Fi= [&](int i) { return G0 + ((std::size_t(1) << (i - 4 - KC)) - 1) * CB * XS; };
  for(int i= KC + 4; i < N; ++i) fwd_chunks(Fi(i), gs + (std::size_t(1) << i), i, 0, 0, i + 1, XS, k, CM);
  u32* cG[H + 1];  // G_j の列は塊 2^{H-j} 個
  for(int j= 0; j <= H; ++j) cG[j]= small_alloc((std::size_t(1) << (H - j)) * bsj[j] * 4);
  u32* cF= small_alloc((std::size_t(1) << (H - 1)) * BS * 4);
  for(std::size_t c= 0; c < CB; ++c) {
   const bool nx= c + 1 < CB;
   for(int j= 0; j <= H; ++j) {
    std::memcpy(cG[j], G[j] + c * lsj[j], std::min(lsj[j], bsj[j]) * 4);
    if(nx) pf_block(G[j] + (c + 1) * lsj[j], std::min(lsj[j], bsj[j]));
   }
   for(int i= KC + 4; i < N; ++i) {
    const int J= N - i;
    const std::size_t hc= std::size_t(1) << (i - 4 - KC);
    col_fwd(cF, Fi(i), c, i, BS, k, XS, CM);
    for(std::size_t t= 0; t < hc; ++t) {
     if(nx) pf_win(Fi(i), (t << KC) + c + 1, __builtin_popcountll(t), 0, i + 1, XS, CM);
     for(int j= 0; j <= J; ++j) gp[j]= cG[j] + t * bsj[j];
     for(int j= 0; j < J; ++j) gh[j]= cG[j] + (t + hc) * bsj[j];
     comp_step(gp, gh, cF + t * BS, i, J, k);
    }
   }
   inv_column(cG[0], G0, c, N, 0, 0, NS, BS, k, XS, CM);
  }
  for(int j= 0; j <= H; ++j) std::free(cG[j]);
  std::free(cF);
 }
 const V Kc= cst(u32(u64(k.inv2pow(N - 1)) * k.R1 % k.P));
 inv_chunks(G0, N, 0, 0, NS, XS, k, [&](std::size_t q, const V* v) {
  for(int t= 0; t < 2; ++t) _mm256_storeu_si256((V*)(out + 16 * q + 8 * t), norm1(mont(v[t], Kc, k), k));
 }, CM);
 out[0]= e[0];
 for(int j= 1; j <= N - 4; ++j) std::free(G[j]);
 std::free(G0), std::free(gs);
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
  int a= 0;
  // corr_a と corr_{a+1} をまとめる。d = a+2..i+1 の項を分け合い、corr_a には d = a+1 の項 D_{a+1} Fr[i] を足す。
  if(j >= 1)
   for(; a + 1 <= i; a+= 2) {
    V r0, r1;
    dot2<false>(De[prv] + 8 * (a + 2), Do[prv] + 8 * (a + 2), Fr + 8 * i, Fq + 8 * i, i - a, z, z, z, z, mul64(ld(De[prv] + 8 * (a + 1)), ld(Fr + 8 * i)), mul64(ld(Do[prv] + 8 * (a + 1)), ld(Fq + 8 * i)), r0, r1, k);
    const V y0= j < J ? ld(S + 8 * a) : z, y1= j < J ? ld(S + 8 * (a + 1)) : z;
    st(Y[j] + 8 * a, add2(y0, r0, k)), st(Y[j] + 8 * (a + 1), add2(y1, r1, k));
   }
  for(; a <= i; ++a) {
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
// V_k = sum_S w[S] (g^k / k!)[S] (k = 0..N、g[0] は 0)。上の段 (区画をまたぐ段) があるときは、w の区画の変換 W と、上の段の ĝ_i の
// 区画の変換 (別の大きい配列 Fb の行 [2^t - 1, 2^{t+1} - 1)) を先に作る。列 j では W の列を前向きの上の段まで済ませて D_0 の列にし、
// 上の段 (i = N-1..KC+4) を作業の配列の中で済ませて、残った D_0..D_H の塊 1 つずつを下の段の配列の塊 j に書く。
template <int N> std::vector<u32> egf_proj_n(const u32* g, const u32* w, const MK& k0) {
 const MK k= k0;
 constexpr int NB= N - 4, KL= NB < KC ? NB : KC, H= NB - KL, NS= N + 1;
 constexpr std::size_t CB= std::size_t(1) << KL, BS= NS * 8, BSF= N * 8;
 u32* gs= huge_alloc(std::size_t(4) << N);
 scale_r(g, gs, std::size_t(1) << N, k);
 u32* D[N - 3];
 std::size_t bsj[N - 3];
 for(int j= 0; j <= N - 4; ++j) bsj[j]= std::size_t(N - j + 1) * 8, D[j]= huge_alloc(std::min(std::size_t(1) << (N - j - 4), CB) * bsj[j] * 4);
 u32 *zl[N - 3], *zh[N - 3], *yy[N - 3];
 if constexpr(H == 0) {
  fwd_chunks(D[0], w, N, 0, 0, NS, BS, k);
 } else {
  constexpr bool CM= true;
  constexpr std::size_t XS= CM ? CWS : BS, XF= CM ? CWS : BSF;
  u32 *W= huge_alloc((std::size_t(1) << NB) * XS * 4), *Fb= huge_alloc(((std::size_t(1) << H) - 1) * CB * XF * 4);
  auto Fi= [&](int i) { return Fb + ((std::size_t(1) << (i - 4 - KC)) - 1) * CB * XF; };
  fwd_chunks(W, w, N, 0, 0, NS, XS, k, CM);
  for(int i= KC + 4; i < N; ++i) fwd_chunks(Fi(i), gs + (std::size_t(1) << i), i, 0, 0, i + 1, XF, k, CM);
  u32* cD[H + 1];  // D_j の列は塊 2^{H-j} 個
  for(int j= 0; j <= H; ++j) cD[j]= small_alloc((std::size_t(1) << (H - j)) * bsj[j] * 4);
  u32* cF= small_alloc((std::size_t(1) << (H - 1)) * BSF * 4);
  for(std::size_t c= 0; c < CB; ++c) {
   fwd_column(cD[0], W, c, N, 0, 0, NS, BS, k, XS, CM);
   // 次の列で読む W の塊は、0 をここで、hc + t を段 i の t 番目で先に読む (hc = 2^{i-4-KC} で、全部で 2^H 個)
   const bool nx= c + 1 < CB;
   if(nx) pf_win(W, c + 1, 0, 0, NS, XS, CM);
   for(int i= N - 1; i >= KC + 4; --i) {
    const int J= N - i - 1;
    const std::size_t hc= std::size_t(1) << (i - 4 - KC);
    col_fwd(cF, Fi(i), c, i, BSF, k, XF, CM);
    for(std::size_t t= 0; t < hc; ++t) {
     if(nx) pf_win(Fi(i), (t << KC) + c + 1, __builtin_popcountll(t), 0, i + 1, XF, CM), pf_win(W, ((hc + t) << KC) + c + 1, __builtin_popcountll(hc + t), 0, NS, XS, CM);
     for(int j= 0; j <= J; ++j) zl[j]= cD[j] + t * bsj[j], zh[j]= cD[j] + (t + hc) * bsj[j];
     for(int j= 0; j <= J + 1; ++j) yy[j]= cD[j] + t * bsj[j];
     proj_step(zl, zh, yy, cF + t * BSF, i, J + 1, k);
    }
   }
   for(int j= 0; j <= H; ++j) std::memcpy(D[j] + c * bsj[j], cD[j], bsj[j] * 4);
  }
  for(int j= 0; j <= H; ++j) std::free(cD[j]);
  std::free(W), std::free(Fb), std::free(cF);
 }
 // 下の段 (区画の中で閉じる)
 u32* Fs= small_alloc((CB >> 1) * BSF * 4);
 for(int i= KL + 3; i >= 4; --i) {
  const int J= N - i - 1;
  const std::size_t cnt= std::size_t(1) << (i - 4);
  fwd_chunks(Fs, gs + (std::size_t(1) << i), i, 0, 0, i + 1, BSF, k);
  for(std::size_t q= 0; q < cnt; ++q) {
   for(int j= 0; j <= J; ++j) zl[j]= D[j] + q * bsj[j], zh[j]= D[j] + (q + cnt) * bsj[j];
   for(int j= 0; j <= J + 1; ++j) yy[j]= D[j] + q * bsj[j];
   proj_step(zl, zh, yy, Fs + q * BSF, i, J + 1, k);
  }
 }
 std::free(Fs);
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
 std::free(gs);
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
inline std::vector<u32> taylor(const u32* f, int M, u32 c0, int n, const MK& k0) {
 const MK k= k0;
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
inline void proj_post(const std::vector<u32>& vk, u32 c0, int M, u32* ans, const MK& k0) {
 const MK k= k0;
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
