#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <type_traits>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 整数の型 (法 2^64) の集合冪級数の exp と log の速い核。sps::exp と sps::log から呼ぶ内部の道具で、直接は使わない。
// 32 bit と 64 bit の整数 (符号付きも) を、法 2^64 の値として求める (32 bit は 64 bit で求めて切り詰める)。
// OR のランク付きのゼータ変換 F̂[y] = sum_{x ⊆ y} f[x] t^{|x|} を使う。部分集合の大きさより rank の大きい項 (重なる組から出る項) は、
// 足し算と掛け算を通しても大きさより rank の大きいままなので、途中で逆変換して rank を取り出さずに済み、rank n までで切れる。
//   exp: 変換の領域で Ĥ[y ∪ {i}] = Ĥ[y] (1 + t F̂_i[y]) を上の要素から順に求める (F̂_i は i を含む部分の i bit の変換)。
//   log: 各段の t F̂_i / F̂_lo (F̂_lo[∅] = 1) を、変換の領域の各点の冪級数の割り算で求める。
// どちらも割り算 (整数の) を使わない。逆変換は最後に 1 回だけで、位置 y の rank |y| の値が答え。
// 持ち方は位置 4 個を 1 本の __m256i (u64 × 4) に並べた塊で、塊ごとに rank の数だけ並べる。下の 2 bit はレーンの中で変換する。
// 位置の下の 11 bit (レーンの 2 bit と塊の 9 bit) を区画とし、区画の中の変換 (区画の段) は区画ごとに L2 の中で 2 段ずつ行う。上の bit
// は、下の 11 bit が同じ塊 (列) を作業の配列に集め、上の段の変換、各点の計算、上の段の逆変換を続けて行う。大きい配列 (huge page を
// 頼む) の列 c の塊は、区画の段の順変換の値 (区画の起点からの相対の rank 0 から popcount(c) + 2) と、区画の段の逆変換に要る値
// (相対の rank popcount(c) から 11) の多い方の幅だけを持つ。今の列を計算している間に、次の列の塊を L2 へ先に読む。
// u64 の積の和は、a_lo b_lo (vpmuludq) と (a_lo b_hi, a_hi b_lo) の下位 32 bit (vpmulld) を足して求める。
// n が 13 未満のときと区画 0 (下の 11 bit の部分) は、段ごとに i bit の変換をする形で求める。
namespace sps_or64 {
using V= __m256i;
#define SPS_OR64_AI [[gnu::always_inline]] inline
template <class T> constexpr bool fits= std::is_integral_v<T> && !std::is_same_v<T, bool> && (sizeof(T) == 4 || sizeof(T) == 8);
SPS_OR64_AI V add(V a, V b) { return _mm256_add_epi64(a, b); }
SPS_OR64_AI V sub(V a, V b) { return _mm256_sub_epi64(a, b); }
SPS_OR64_AI V hi32(V x) { return _mm256_srli_epi64(x, 32); }
SPS_OR64_AI V swap32(V x) { return _mm256_shuffle_epi32(x, 0xB1); }
// レーンの中の 2 bit のゼータ変換とメビウス変換
template <bool INV> SPS_OR64_AI V lane_zeta(V v) {
 const V m13= _mm256_setr_epi64x(0, -1, 0, -1), m23= _mm256_setr_epi64x(0, 0, -1, -1);
 V t= _mm256_and_si256(_mm256_permute4x64_epi64(v, 0xA0), m13);
 v= INV ? sub(v, t) : add(v, t);
 t= _mm256_and_si256(_mm256_permute4x64_epi64(v, 0x44), m23);
 return INV ? sub(v, t) : add(v, t);
}
// u64 の積の和 (段ごとの形): lo = sum a_lo b_lo、c = sum (a_lo b_hi + a_hi b_lo)。結果は lo + (c << 32)
SPS_OR64_AI void mac(V& lo, V& c, V a, V ah, V b, V bh) { lo= add(lo, _mm256_mul_epu32(a, b)), c= add(c, add(_mm256_mul_epu32(a, bh), _mm256_mul_epu32(ah, b))); }
SPS_OR64_AI V fin(V lo, V c) { return add(lo, _mm256_slli_epi64(c, 32)); }
// u64 の積の和 (区画と列の形): bs は b の 32 bit の半分を入れ替えたもの。c は 32 bit ごとの和。結果は lo + ((c_lo + c_hi) << 32)
SPS_OR64_AI void mac2(V& lo, V& c, V a, V b, V bs) { lo= add(lo, _mm256_mul_epu32(a, b)), c= _mm256_add_epi32(c, _mm256_mullo_epi32(a, bs)); }
SPS_OR64_AI V fin2(V lo, V c) { return add(lo, _mm256_slli_epi64(add(c, _mm256_srli_epi64(c, 32)), 32)); }
inline V* valloc(std::size_t n) { return static_cast<V*>(std::aligned_alloc(32, std::max<std::size_t>(n, 1) * sizeof(V))); }
inline V* huge_alloc(std::size_t n) {
 const std::size_t H= std::size_t(1) << 21;
 const std::size_t bytes= (n * sizeof(V) + H - 1) & ~(H - 1);
 void* p= std::aligned_alloc(H, bytes);
#ifdef __linux__
 madvise(p, bytes, MADV_HUGEPAGE);
#endif
 return static_cast<V*>(p);
}
SPS_OR64_AI void pf_block(const V* p, int n) {
 const std::uintptr_t e= std::uintptr_t(p + n);
 for(std::uintptr_t a= std::uintptr_t(p) & ~std::uintptr_t(63); a < e; a+= 64) __builtin_prefetch((const void*)a, 0, 2);
}
// ---- 段ごとの形 (n < 13 と区画 0) ----
// f の 2^m 個 (m >= 2) のランク付きのゼータ変換を F (塊 2^{m-2} 個、幅 S) に。塊 q は rank [0, popcount(q) + 2] を持つ。
template <class U> inline void rz_fwd(const U* f, V* F, int m, int S) {
 const std::size_t nb= std::size_t(1) << (m - 2);
 const V M0= _mm256_setr_epi64x(-1, 0, 0, 0), M1= _mm256_setr_epi64x(0, -1, -1, 0), M2= _mm256_setr_epi64x(0, 0, 0, -1);
 for(std::size_t q= 0; q < nb; ++q) {
  const V x= _mm256_loadu_si256((const V*)(f + 4 * q));
  V* b= F + q * S;
  const int c= __builtin_popcountll(q);
  b[c]= lane_zeta<false>(_mm256_and_si256(x, M0));
  b[c + 1]= lane_zeta<false>(_mm256_and_si256(x, M1));
  b[c + 2]= lane_zeta<false>(_mm256_and_si256(x, M2));
 }
 for(int kk= 0; kk < m - 2; ++kk) {
  const std::size_t d= std::size_t(1) << kk;
  for(std::size_t q0= 0; q0 < nb; q0+= 2 * d)
   for(std::size_t q= q0; q < q0 + d; ++q) {
    const V* s= F + q * S;
    V* t= F + (q + d) * S;
    const int lo= __builtin_popcountll(q >> kk), hi= __builtin_popcountll(q) + 2;
    t[lo]= s[lo];
    for(int r= lo + 1; r <= hi; ++r) t[r]= add(t[r], s[r]);
   }
 }
}
// 塊 q が rank [popcount(q), top(q)] を持つ F のメビウス変換をして、位置 y の rank popcount(y) を out に取り出す。
template <class U, class Top> inline void rz_inv(V* F, U* out, int m, int S, const Top& top) {
 const std::size_t nb= std::size_t(1) << (m - 2);
 for(int kk= 0; kk < m - 2; ++kk) {
  const std::size_t d= std::size_t(1) << kk;
  for(std::size_t q0= 0; q0 < nb; q0+= 2 * d)
   for(std::size_t q= q0; q < q0 + d; ++q) {
    const V* s= F + q * S;
    V* t= F + (q + d) * S;
    const int a= __builtin_popcountll(q) + 1, b= std::min(top(q), top(q + d));
    for(int r= a; r <= b; ++r) t[r]= sub(t[r], s[r]);
   }
 }
 const V M0= _mm256_setr_epi64x(-1, 0, 0, 0), M1= _mm256_setr_epi64x(0, -1, -1, 0), M2= _mm256_setr_epi64x(0, 0, 0, -1);
 for(std::size_t q= 0; q < nb; ++q) {
  const V* b= F + q * S;
  const int c= __builtin_popcountll(q), tp= top(q);
  V r= _mm256_and_si256(lane_zeta<true>(b[c]), M0);
  if(c + 1 <= tp) r= _mm256_or_si256(r, _mm256_and_si256(lane_zeta<true>(b[c + 1]), M1));
  if(c + 2 <= tp) r= _mm256_or_si256(r, _mm256_and_si256(lane_zeta<true>(b[c + 2]), M2));
  _mm256_storeu_si256((V*)(out + 4 * q), r);
 }
}
// h = exp(f) (f[0] = 0) の先頭 2^n 個
template <class U> inline void exp_small(const U* f, U* h, int n) {
 h[0]= 1;
 int i= 0;
 for(; i < n && i < 6; ++i) {  // h_hi[s] = sum_{x ⊆ s} f_hi[x] h[s \ x]
  const std::size_t L= std::size_t(1) << i;
  for(std::size_t s= 0; s < L; ++s) {
   U t= 0;
   for(std::size_t x= s;; x= (x - 1) & s) {
    t+= f[L | x] * h[s ^ x];
    if(!x) break;
   }
   h[L | s]= t;
  }
 }
 if(i == n) return;
 const int S= n;
 V *A= valloc((std::size_t(1) << (n - 3)) * S), *B= valloc((std::size_t(1) << (n - 3)) * S);
 alignas(32) V ah[32], bh[32];
 for(; i < n; ++i) {
  const std::size_t L= std::size_t(1) << i, nb= L >> 2;
  rz_fwd(h, A, i, S), rz_fwd(f + L, B, i, S);
  for(std::size_t q= 0; q < nb; ++q) {
   V *a= A + q * S, *b= B + q * S;
   const int pc= __builtin_popcountll(q), w= pc + 2, top= std::min(i, 2 * pc + 4);
   for(int e= 0; e <= w; ++e) ah[e]= hi32(a[e]), bh[e]= hi32(b[e]);
   for(int d= top; d >= pc; --d) {  // 上から書けば、まだ要る a[e] (e <= d) を壊さない
    V lo= _mm256_setzero_si256(), c= lo;
    const int e0= std::max(0, d - w), e1= std::min(d, w);
    for(int e= e0; e <= e1; ++e) mac(lo, c, a[e], ah[e], b[d - e], bh[d - e]);
    a[d]= fin(lo, c);
   }
  }
  rz_inv(A, h + L, i, S, [&](std::size_t q) { return std::min(i, 2 * __builtin_popcountll(q) + 4); });
 }
 std::free(A), std::free(B);
}
// g = log(f) (f[0] = 1) の先頭 2^n 個
template <class U> inline void log_small(const U* f, U* g, int n) {
 g[0]= 0;
 int i= 0;
 for(; i < n && i < 6; ++i) {  // g_hi[s] = f_hi[s] - sum_{∅≠x⊆s} f[x] g_hi[s \ x]
  const std::size_t L= std::size_t(1) << i;
  for(std::size_t s= 0; s < L; ++s) {
   U t= f[L | s];
   for(std::size_t x= s; x; x= (x - 1) & s) t-= f[x] * g[L | (s ^ x)];
   g[L | s]= t;
  }
 }
 if(i == n) return;
 const int S= n;
 // f の下の半分 (n - 1 bit) のゼータ変換を一度だけ (段 i ではその先頭の塊 2^{i-2} 個が i bit の変換)
 V *Bf= valloc((std::size_t(1) << (n - 3)) * S), *A= valloc((std::size_t(1) << (n - 3)) * S);
 rz_fwd(f, Bf, n - 1, S);
 alignas(32) V bh[32], q[32], qh[32];
 for(; i < n; ++i) {
  const std::size_t L= std::size_t(1) << i, nb= L >> 2;
  rz_fwd(f + L, A, i, S);
  for(std::size_t qq= 0; qq < nb; ++qq) {
   V* a= A + qq * S;
   const V* b= Bf + qq * S;
   const int w= __builtin_popcountll(qq) + 2;
   for(int e= 1; e <= w; ++e) bh[e]= hi32(b[e]);
   for(int d= 0; d <= i; ++d) {  // Q[d] = A[d] - sum_{e=1}^{min(d, w)} B[e] Q[d-e]
    V lo= _mm256_setzero_si256(), c= lo;
    for(int e= 1; e <= std::min(d, w); ++e) mac(lo, c, b[e], bh[e], q[d - e], qh[d - e]);
    const V ad= d <= w ? a[d] : _mm256_setzero_si256();
    q[d]= sub(ad, fin(lo, c)), qh[d]= hi32(q[d]);
   }
   for(int d= 0; d <= i; ++d) a[d]= q[d];
  }
  rz_inv(A, g + L, i, S, [&](std::size_t) { return i; });
 }
 std::free(A), std::free(Bf);
}
// ---- 区画と列の形 ----
constexpr int KC= 9, M0= KC + 2, NB= 1 << KC;
struct Layout {
 int off[NB + 1];  // 区画の中の塊の場所 (V の数)
 constexpr Layout(): off() {
  int o= 0;
  for(int c= 0; c < NB; ++c) {
   off[c]= o;
   const int p= __builtin_popcount(c);
   o+= std::max(p + 3, M0 + 1 - p);
  }
  off[NB]= o;
 }
};
inline constexpr Layout LY{};
constexpr int CS= LY.off[NB];  // 区画の幅 (V の数)
// 区画の段の順変換: 2^11 個の値 src から区画 dst へ。塊 q は相対の rank r を場所 r に置き、[0, popcount(q) + 2] を持つ。
template <class U> inline void chunk_fwd(const U* src, V* dst) {
 const V M0v= _mm256_setr_epi64x(-1, 0, 0, 0), M1v= _mm256_setr_epi64x(0, -1, -1, 0), M2v= _mm256_setr_epi64x(0, 0, 0, -1);
 for(int q= 0; q < NB; ++q) {
  const V x= _mm256_loadu_si256((const V*)(src + 4 * q));
  V* b= dst + LY.off[q];
  const int c= __builtin_popcount(q);
  b[c]= lane_zeta<false>(_mm256_and_si256(x, M0v));
  b[c + 1]= lane_zeta<false>(_mm256_and_si256(x, M1v));
  b[c + 2]= lane_zeta<false>(_mm256_and_si256(x, M2v));
 }
 int kk= 0;
 for(; kk + 2 <= KC; kk+= 2) {
  const int d= 1 << kk;
  for(int q0= 0; q0 < NB; q0+= 4 * d)
   for(int q= q0; q < q0 + d; ++q) {
    V *b0= dst + LY.off[q], *b1= dst + LY.off[q + d], *b2= dst + LY.off[q + 2 * d], *b3= dst + LY.off[q + 3 * d];
    // 前の窓: b0 [lo, hi]、b1 と b2 [lo + 1, hi + 1]、b3 [lo + 2, hi + 2]。後は b1 と b2 が [lo, hi + 1]、b3 が [lo, hi + 2]
    const int lo= __builtin_popcount(q >> kk), hi= __builtin_popcount(q) + 2;
    {
     const V v0= b0[lo];
     b1[lo]= v0, b2[lo]= v0, b3[lo]= v0;
    }
    {
     const V v0= b0[lo + 1], v1= b1[lo + 1], v2= b2[lo + 1];
     b1[lo + 1]= add(v1, v0), b2[lo + 1]= add(v2, v0), b3[lo + 1]= add(add(v1, v2), v0);
    }
    for(int r= lo + 2; r <= hi; ++r) {
     const V v0= b0[r], v1= b1[r], v2= b2[r], v3= b3[r];
     const V s01= add(v0, v1);
     b1[r]= s01, b2[r]= add(v2, v0), b3[r]= add(add(v3, v2), s01);
    }
    b3[hi + 1]= add(add(b3[hi + 1], b2[hi + 1]), b1[hi + 1]);
   }
 }
 for(; kk < KC; ++kk) {
  const int d= 1 << kk;
  for(int q0= 0; q0 < NB; q0+= 2 * d)
   for(int q= q0; q < q0 + d; ++q) {
    const V* s= dst + LY.off[q];
    V* t= dst + LY.off[q + d];
    const int lo= __builtin_popcount(q >> kk), hi= __builtin_popcount(q) + 2;
    t[lo]= s[lo];
    for(int r= lo + 1; r <= hi; ++r) t[r]= add(t[r], s[r]);
   }
 }
}
// 区画の段の逆変換: 塊 q は相対の rank [popcount(q), 11] を場所 rank - popcount(q) に持つ。位置 y の相対の rank popcount(y) を out に。
template <class U> inline void chunk_inv(V* src, U* out) {
 int kk= 0;
 for(; kk + 2 <= KC; kk+= 2) {
  const int d= 1 << kk;
  for(int q0= 0; q0 < NB; q0+= 4 * d)
   for(int q= q0; q < q0 + d; ++q) {
    // rank で引けるように、塊の先頭から起点をずらした場所
    const int c= __builtin_popcount(q);
    V *s0= src + LY.off[q] - c, *s1= src + LY.off[q + d] - (c + 1), *s2= src + LY.off[q + 2 * d] - (c + 1), *s3= src + LY.off[q + 3 * d] - (c + 2);
    // 要る範囲: s1 と s2 は [c + 1, T - 1]、s3 は [c + 2, T] (T は s0 に要る範囲の上端)
    const int T= c + KC - kk - __builtin_popcount(q >> kk) + 2;
    {
     const V v0= s0[c + 1];
     s1[c + 1]= sub(s1[c + 1], v0), s2[c + 1]= sub(s2[c + 1], v0);
    }
    for(int r= c + 2; r < T; ++r) {
     const V v0= s0[r], v1= s1[r], v2= s2[r], v3= s3[r];
     const V m1= sub(v1, v0);
     s1[r]= m1, s2[r]= sub(v2, v0), s3[r]= sub(sub(v3, v2), m1);
    }
    s3[T]= add(sub(sub(s3[T], s2[T]), s1[T]), s0[T]);
   }
 }
 for(; kk < KC; ++kk) {
  const int d= 1 << kk;
  for(int q0= 0; q0 < NB; q0+= 2 * d)
   for(int q= q0; q < q0 + d; ++q) {
    const int c= __builtin_popcount(q);
    const V* s= src + LY.off[q] - c;
    V* t= src + LY.off[q + d] - (c + 1);
    const int b= c + KC - kk - __builtin_popcount(q >> kk) + 2;
    for(int r= c + 1; r <= b; ++r) t[r]= sub(t[r], s[r]);
   }
 }
 const V M0v= _mm256_setr_epi64x(-1, 0, 0, 0), M1v= _mm256_setr_epi64x(0, -1, -1, 0), M2v= _mm256_setr_epi64x(0, 0, 0, -1);
 for(int q= 0; q < NB; ++q) {
  const V* b= src + LY.off[q];
  V r= _mm256_and_si256(lane_zeta<true>(b[0]), M0v);
  r= _mm256_or_si256(r, _mm256_and_si256(lane_zeta<true>(b[1]), M1v));
  r= _mm256_or_si256(r, _mm256_and_si256(lane_zeta<true>(b[2]), M2v));
  _mm256_storeu_si256((V*)(out + 4 * q), r);
 }
}
// 列の塊の窓を、塊の下の bit から上の段のゼータ変換
template <class Lo, class Hi> SPS_OR64_AI void col_zeta(V* Z, int S, int NJ, int nb, const Lo& lo0, const Hi& hi0) {
 for(int kk= 0; kk < nb; ++kk) {
  const int d= 1 << kk;
  for(int q0= 0; q0 < NJ; q0+= 2 * d)
   for(int q= q0; q < q0 + d; ++q) {
    const V* s= Z + std::size_t(q) * S;
    V* t= Z + std::size_t(q + d) * S;
    const int lo= lo0(q, kk), hi= hi0(q);
    t[lo]= s[lo];
    for(int r= lo + 1; r <= hi; ++r) t[r]= add(t[r], s[r]);
   }
 }
}
// 大きい配列の、区画 j と列 c の塊
SPS_OR64_AI V* blk(V* big, int j, int c) { return big + std::size_t(j) * CS + LY.off[c]; }
// 列 c の段の初め: 区画 j の塊の rank [popcount(j), popcount(j) + popcount(c) + 2] を作業の配列 Z の塊 j へ
SPS_OR64_AI void col_gather(V* big, V* Z, int S, int NJ, int c) {
 const int w0= __builtin_popcount(c) + 2;
 for(int j= 0; j < NJ; ++j) {
  const V* b= blk(big, j, c);
  V* z= Z + std::size_t(j) * S + __builtin_popcount(j);
  for(int s= 0; s <= w0; ++s) z[s]= b[s];
 }
}
// 列 c の段の終わり: 区画 j >= 1 の塊へ、rank [popcount(j) + popcount(c), popcount(j) + 11] を場所 0 から
SPS_OR64_AI void col_scatter(V* big, const V* Z, int S, int NJ, int c) {
 const int pcc= __builtin_popcount(c);
 for(int j= 1; j < NJ; ++j) {
  V* b= blk(big, j, c) - pcc;
  const V* z= Z + std::size_t(j) * S + __builtin_popcount(j);
  for(int s= pcc; s <= M0; ++s) b[s]= z[s];
 }
}
template <class U> inline void exp_big(const U* f, U* h, int n) {
 const int U_= n - M0, NJ= 1 << U_, S= n + 1;
 exp_small(f, h, M0);
 V* big= huge_alloc(std::size_t(NJ) * CS);
 chunk_fwd(h, big);
 for(int j= 1; j < NJ; ++j) chunk_fwd(f + (std::size_t(j) << M0), big + std::size_t(j) * CS);
 V* Z= valloc(std::size_t(NJ) * S);
 alignas(32) V Fl[32], Fs[32];
 for(int c= 0; c < NB; ++c) {
  const int pcc= __builtin_popcount(c), w0= pcc + 2;
  const int wn= c + 1 < NB ? LY.off[c + 2] - LY.off[c + 1] : 0;
  auto pf= [&](int j) {
   if(wn) pf_block(blk(big, j, c + 1), wn);
  };
  pf(0);
  col_gather(big, Z, S, NJ, c);
  for(int r= w0 + 1; r <= n; ++r) Z[r]= _mm256_setzero_si256();
  for(int i= M0; i < n; ++i) {
   const int L= 1 << (i - M0);
   V* ZL= Z + std::size_t(L) * S;
   // f_i の塊 (列の L..2L-1) を上の段で変換して t F̂_i に (rank は要素 i のぶん 1 ずれている)
   col_zeta(ZL, S, L, i - M0, [](int q, int kk) { return 1 + __builtin_popcount(q >> kk); }, [&](int q) { return 1 + __builtin_popcount(q) + w0; });
   for(int q= 0; q < L; ++q) {
    pf(L + q);
    const V* H= Z + std::size_t(q) * S;
    V* T= ZL + std::size_t(q) * S;
    const int wF= 1 + __builtin_popcount(q) + w0;
    for(int e= 1; e <= wF; ++e) Fl[e]= T[e], Fs[e]= swap32(T[e]);
    // T[r] = H[r] + sum_{e=1}^{min(r, wF)} H[r-e] F[e] を 2 つずつ (r と r + 1)。H の値はレジスタで次の項へ渡す
    T[0]= H[0];
    int r= 1;
    for(; r + 1 <= n; r+= 2) {
     V lo0= _mm256_setzero_si256(), c0= lo0, lo1= lo0, c1= lo0;
     V y= H[r];
     const int em= std::min(r, wF);
     for(int e= 1; e <= em; ++e) {
      const V x= H[r - e];
      mac2(lo1, c1, y, Fl[e], Fs[e]);
      mac2(lo0, c0, x, Fl[e], Fs[e]);
      y= x;
     }
     if(r + 1 <= wF) mac2(lo1, c1, y, Fl[r + 1], Fs[r + 1]);
     T[r]= add(H[r], fin2(lo0, c0)), T[r + 1]= add(H[r + 1], fin2(lo1, c1));
    }
    if(r == n) {
     V lo0= _mm256_setzero_si256(), c0= lo0;
     for(int e= 1; e <= std::min(r, wF); ++e) mac2(lo0, c0, H[r - e], Fl[e], Fs[e]);
     T[r]= add(H[r], fin2(lo0, c0));
    }
   }
  }
  // 上の段の逆変換。塊 j に要るのは rank [popcount(j) + popcount(c), popcount(j) + 11 + 残りの上の bit のうち j に無いものの数]
  for(int kk= 0; kk < U_; ++kk) {
   const int d= 1 << kk;
   for(int q0= 0; q0 < NJ; q0+= 2 * d)
    for(int q= q0; q < q0 + d; ++q) {
     const V* s= Z + std::size_t(q) * S;
     V* t= Z + std::size_t(q + d) * S;
     const int a= __builtin_popcount(q) + pcc + 1, b= __builtin_popcount(q) + M0 + U_ - kk - __builtin_popcount(q >> kk);
     for(int r= a; r <= b; ++r) t[r]= sub(t[r], s[r]);
    }
  }
  col_scatter(big, Z, S, NJ, c);
 }
 for(int j= 1; j < NJ; ++j) chunk_inv(big + std::size_t(j) * CS, h + (std::size_t(j) << M0));
 std::free(Z), std::free(big);
}
template <class U> inline void log_big(const U* f, U* g, int n) {
 const int U_= n - M0, NJ= 1 << U_, S= n + 1;
 log_small(f, g, M0);
 V* big= huge_alloc(std::size_t(NJ) * CS);
 for(int j= 0; j < NJ; ++j) chunk_fwd(f + (std::size_t(j) << M0), big + std::size_t(j) * CS);
 V* Z= valloc(std::size_t(NJ) * S);
 alignas(32) V Bl[32], Bs[32], Q[32];
 for(int c= 0; c < NB; ++c) {
  const int pcc= __builtin_popcount(c), w0= pcc + 2;
  const int wn= c + 1 < NB ? LY.off[c + 2] - LY.off[c + 1] : 0;
  auto pf= [&](int j) {
   if(wn) pf_block(blk(big, j, c + 1), wn);
  };
  pf(0);
  col_gather(big, Z, S, NJ, c);
  // 上の段を全部変換する。塊 j は rank [0, popcount(j) + popcount(c) + 2] を持つ
  col_zeta(Z, S, NJ, U_, [](int q, int kk) { return __builtin_popcount(q >> kk); }, [&](int q) { return __builtin_popcount(q) + w0; });
  // 段 i (上から): 列の 0..L-1 が F̂_lo、L..2L-1 が F̂_lo + t F̂_i。t Q_i = t F̂_i / F̂_lo を L..2L-1 に置いて、上の段を逆変換する
  for(int i= n - 1; i >= M0; --i) {
   const int L= 1 << (i - M0);
   V* ZL= Z + std::size_t(L) * S;
   for(int q= 0; q < L; ++q) {
    pf(L + q);
    const V* B= Z + std::size_t(q) * S;
    V* T= ZL + std::size_t(q) * S;
    const int wB= __builtin_popcount(q) + w0;
    for(int e= 1; e <= wB; ++e) Bl[e]= B[e], Bs[e]= swap32(B[e]);
    Q[0]= _mm256_setzero_si256();
    for(int r= 1; r <= i + 1; ++r) {
     V lo= _mm256_setzero_si256(), cc= lo;
     for(int e= std::min(r - 1, wB); e >= 1; --e) mac2(lo, cc, Q[r - e], Bl[e], Bs[e]);
     const V tf= r <= wB ? sub(T[r], B[r]) : r == wB + 1 ? T[r] : _mm256_setzero_si256();
     Q[r]= sub(tf, fin2(lo, cc));
    }
    for(int r= 0; r <= i + 1; ++r) T[r]= Q[r];
   }
   for(int kk= 0; kk < i - M0; ++kk) {
    const int d= 1 << kk;
    for(int q0= 0; q0 < L; q0+= 2 * d)
     for(int q= q0; q < q0 + d; ++q) {
      const V* s= ZL + std::size_t(q) * S;
      V* t= ZL + std::size_t(q + d) * S;
      const int a= __builtin_popcount(q) + 2 + pcc, b= 1 + __builtin_popcount(q) + i - kk - __builtin_popcount(q >> kk);
      for(int r= a; r <= b; ++r) t[r]= sub(t[r], s[r]);
     }
   }
  }
  col_scatter(big, Z, S, NJ, c);
 }
 for(int j= 1; j < NJ; ++j) chunk_inv(big + std::size_t(j) * CS, g + (std::size_t(j) << M0));
 std::free(Z), std::free(big);
}
// 整数の型 T の f から exp(f) か log(f) を法 2^64 で (32 bit は切り詰めて) 求める
template <bool LOG, class T> std::vector<T> run(const std::vector<T>& f, int n) {
 std::vector<T> h(f.size());
 auto go= [&](const auto* a, auto* b) {
  if(n < M0 + 2) LOG ? log_small(a, b, n) : exp_small(a, b, n);
  else LOG ? log_big(a, b, n) : exp_big(a, b, n);
 };
 if constexpr(sizeof(T) == 8) {
  using U= std::make_unsigned_t<T>;  // 符号付きと符号なしは同じ場所を読み書きしてよい
  go(reinterpret_cast<const U*>(f.data()), reinterpret_cast<U*>(h.data()));
 } else {
  std::vector<unsigned long long> a(f.size()), b(f.size());
  for(std::size_t i= 0; i < f.size(); ++i) a[i]= (unsigned long long)f[i];
  go(a.data(), b.data());
  for(std::size_t i= 0; i < f.size(); ++i) h[i]= T(b[i]);
 }
 return h;
}
#undef SPS_OR64_AI
}  // namespace sps_or64
