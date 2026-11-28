#pragma once
// AVX2 の実装 (procon-judge の self/ntt の 3 本目と 5 本目)。ntt.hpp が __AVX2__ か SIMDE_ENABLE_NATIVE_ALIASES の
// ときに include する。x86 では intrinsics、それ以外は SIMDe で、同じソースが通る。
// 値は Montgomery (R = 2^32) で持ち、dft の入口で直して idft の出口で戻す (dft と idft の入出力はふつうの表現)。
// dft は DIF で、評価点は 4 進の桁を逆にした順。idft は DIT。どちらも長さ nv の配列をベクトル (8 個の係数) で回す。
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>
#if defined(USE_SIMDE) || defined(SIMDE_ENABLE_NATIVE_ALIASES)
#include <simde/x86/avx2.h>
#if !defined(SIMDE_ENABLE_NATIVE_ALIASES)
#error "SIMDe のときは SIMDE_ENABLE_NATIVE_ALIASES も定義してください"
#endif
#define NEO_NTT_SIMDE 1
#elif defined(__AVX2__)
#include <immintrin.h>
#else
#error "AVX2 の intrinsics か SIMDe が要ります (-march=x86-64-v3 か -DSIMDE_ENABLE_NATIVE_ALIASES)"
#endif
#include <cstdlib>
#include "neo/base/core/algebra/ZMod.hpp"
#include "neo/base/core/convolve/naive.hpp"

// GCC の -O2 は段の内側の組み方が甘く、O3 と unroll-loops で self-ntt-convolution が 8% ほど速い
// (procon-judge の核の 12 本目で測った)。clang はこの pragma を知らないので、そのまま -O2 で組まれる。
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize("O3,unroll-loops")
#define NEO_NTT_AVX2_PUSHED 1
#endif

namespace neo_ntt_avx2 {
using u32= unsigned;
using u64= unsigned long long;
using V= __m256i;

// 小さい段をまとめて回すタイルの長さ (2 の冪)。4 KB から 8 KB くらい。
constexpr int dft_tile= 2048;
// 葉の積をまとめて回すタイルのベクトル数 (要素数では 8 倍)。
constexpr int tile_vec= 256;

// そろえた配列。表は 32 byte、作業の配列は 64 byte にそろえる (5 本目を測ったときと同じ条件)。
template <int ALIGN> struct abuf {
 u32* p= nullptr;
 int n= 0;
 abuf()= default;
 abuf(const abuf&)= delete;
 abuf& operator=(const abuf&)= delete;
 ~abuf() { std::free(p); }
 // 長さを n_ にして先頭を返す。すでにその長さなら何もしない (中身はそのまま)。
 u32* assign(int n_) {
  if (n == n_) return p;
  std::free(p);
  void* q= nullptr;
  if (posix_memalign(&q, ALIGN, size_t(n_) * 4)) std::abort();
  p= static_cast<u32*>(q), n= n_;
  return p;
 }
};

// n は奇数。n^{-1} mod 2^32。
constexpr u32 inv_mod2_32(u32 n) {
 u32 x= 1;
 for (int i= 0; i < 5; ++i) x*= 2 - n * x;
 return x;
}
template <unsigned P> constexpr u32 pow2_32_mod() {
 u64 r= 1;
 for (int i= 0; i < 32; ++i) r= r * 2 % P;
 return u32(r);
}

template <unsigned P> struct mc {
 static constexpr u32 p= P, p2= 2 * P;
 static constexpr u32 ni= u32(0) - inv_mod2_32(P);  // -P^{-1} mod 2^32
 static constexpr u32 r= pow2_32_mod<P>();          // 2^32 mod P
 static constexpr u32 r2= u32(u64(r) * r % P);      // 2^64 mod P
};

template <unsigned P> u32 spow(u32 a, u64 e) {  // ふつうの値のべき (a < P)
 u64 r= 1, b= a;
 for (; e; e>>= 1, b= b * b % P)
  if (e & 1) r= r * b % P;
 return u32(r);
}

// 原始根 (ふつうの値)。ntt.hpp の primitive_root を使い回す。
template <unsigned P> u32 groot() {
 static const u32 g= ntt_internal::primitive_root<P>().val();
 return g;
}

// Montgomery の掛け算 (スカラー)。a < 4P、b < P で、結果は P 未満。
template <unsigned P> constexpr u32 smul(u32 a, u32 b) {
 constexpr u32 ni= mc<P>::ni;
 u64 x= u64(a) * b;
 u64 z= (x + u64(u32(x) * ni) * P) >> 32;
 return u32(z >= P ? z - P : z);
}
// ふつうの値を Montgomery に直す (a < P でなくてもよい。a < 4P で結果 < P)。
template <unsigned P> u32 smont(u32 a) { return smul<P>(a, mc<P>::r2); }

template <unsigned P> u32 slow(u32 x) { return x >= mc<P>::p2 ? x - mc<P>::p2 : x; }
template <unsigned P> u32 sdiff(u32 x, u32 y) { return x + mc<P>::p2 - y; }

template <unsigned P> V vsplat(u32 x) { return _mm256_set1_epi32(int(x)); }
template <unsigned P> V vload(const u32* p) { return _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p)); }
template <unsigned P> void vstore(u32* p, V x) { _mm256_storeu_si256(reinterpret_cast<__m256i*>(p), x); }
template <unsigned P> __m128i vload128(const u32* p) { return _mm_loadu_si128(reinterpret_cast<const __m128i*>(p)); }
template <unsigned P> __m128i vsplat128(u32 x) { return _mm_set1_epi32(int(x)); }
// 128 bit を 2 本つないで 1 本の 8 レーンにする。
template <unsigned P> V vcombine128(__m128i lo, __m128i hi) { return _mm256_inserti128_si256(_mm256_castsi128_si256(lo), hi, 1); }
template <unsigned P> V vplus(V x, V y) { return _mm256_add_epi32(x, y); }
template <unsigned P> V vdiff(V x, V y) { return _mm256_sub_epi32(_mm256_add_epi32(x, vsplat<P>(mc<P>::p2)), y); }
// 4P 未満の値を 2P 未満に落とす。
template <unsigned P> V vlow(V x) { return _mm256_min_epu32(x, _mm256_sub_epi32(x, vsplat<P>(mc<P>::p2))); }
// 2P 未満の値を P 未満に落とす (表の値を正規化する)。
template <unsigned P> V vcanon(V x) { return _mm256_min_epu32(x, _mm256_sub_epi32(x, vsplat<P>(P))); }

// Montgomery の掛け算 (8 レーン)。x < 4P、w < P。結果は 2P 未満。
// 積 t = x w の下位 32 bit から q = t * ni mod 2^32 を作り、t + q P の上位 32 bit を取る。
// 偶数レーンは t を右に 32 ずらし、奇数レーンは o の上位 32 bit がもともと奇数レーンにあるので混ぜる。
// 奇数レーンは x と w を右に 32 ずらしてから掛ける (ずらさないと偶数レーンの値と混ざる)。
template <unsigned P> V vmul(V x, V w) {
 const V pv= vsplat<P>(P), niv= vsplat<P>(mc<P>::ni);
 const V xo= _mm256_srli_epi64(x, 32), wo= _mm256_srli_epi64(w, 32);
 const V t= _mm256_mul_epu32(x, w), to= _mm256_mul_epu32(xo, wo);
 const V qe= _mm256_mul_epu32(t, niv), qo= _mm256_mul_epu32(to, niv);
 const V e= _mm256_add_epi64(t, _mm256_mul_epu32(qe, pv));
 const V o= _mm256_add_epi64(to, _mm256_mul_epu32(qo, pv));
 return _mm256_blend_epi32(_mm256_srli_epi64(e, 32), o, 0xAA);
}
// w が全レーン同じとき (率や 4 乗根)。右にずらすのが 1 つ減る。
template <unsigned P> V vmuls(V x, V w) {
 const V pv= vsplat<P>(P), niv= vsplat<P>(mc<P>::ni);
 const V xo= _mm256_srli_epi64(x, 32);
 const V t= _mm256_mul_epu32(x, w), to= _mm256_mul_epu32(xo, w);
 const V qe= _mm256_mul_epu32(t, niv), qo= _mm256_mul_epu32(to, niv);
 const V e= _mm256_add_epi64(t, _mm256_mul_epu32(qe, pv));
 const V o= _mm256_add_epi64(to, _mm256_mul_epu32(qo, pv));
 return _mm256_blend_epi32(_mm256_srli_epi64(e, 32), o, 0xAA);
}

// ひねりを表に持ち、wi = w * (-P^{-1}) mod 2^32 を掛けて商を 1 命令 (vpmulld) で 8 レーン作る形。
// x < 4P、w < P、wi = w * ni。結果は 2P 未満。intrinsics の vpmulld は下位 32 bit だけなので符号は問わない。
// SIMDe の mullo は符号付きで掛けて UBSan に鳴るため、SIMDe のときは vmul / vmuls に回す。
#ifndef NEO_NTT_SIMDE
template <unsigned P> V vmult(V x, V w, V wi) {
 const V pv= vsplat<P>(P);
 const V xo= _mm256_srli_epi64(x, 32), wo= _mm256_srli_epi64(w, 32);
 const V q= _mm256_mullo_epi32(x, wi);
 const V t= _mm256_mul_epu32(x, w), to= _mm256_mul_epu32(xo, wo);
 const V e= _mm256_add_epi64(t, _mm256_mul_epu32(q, pv));
 const V o= _mm256_add_epi64(to, _mm256_mul_epu32(_mm256_srli_epi64(q, 32), pv));
 return _mm256_blend_epi32(_mm256_srli_epi64(e, 32), o, 0xAA);
}
// w が全レーン同じとき (率や 4 乗根、1/n)。
template <unsigned P> V vmults(V x, V w, V wi) {
 const V pv= vsplat<P>(P);
 const V xo= _mm256_srli_epi64(x, 32);
 const V q= _mm256_mullo_epi32(x, wi);
 const V t= _mm256_mul_epu32(x, w), to= _mm256_mul_epu32(xo, w);
 const V e= _mm256_add_epi64(t, _mm256_mul_epu32(q, pv));
 const V o= _mm256_add_epi64(to, _mm256_mul_epu32(_mm256_srli_epi64(q, 32), pv));
 return _mm256_blend_epi32(_mm256_srli_epi64(e, 32), o, 0xAA);
}
#else
template <unsigned P> V vmult(V x, V w, V) { return vmul<P>(x, w); }
template <unsigned P> V vmults(V x, V w, V) { return vmuls<P>(x, w); }
#endif

// 1 の原始 4 乗根 (順) とその逆。段ごとに作り直さない。
template <unsigned P> const V& im_vec() {
 static const V v= vsplat<P>(smont<P>(spow<P>(groot<P>(), u64(P - 1) / 4)));
 return v;
}
template <unsigned P> const V& ii_vec() {
 static const V v= vsplat<P>(smont<P>(spow<P>(groot<P>(), u64(P - 1) - u64(P - 1) / 4)));
 return v;
}
template <unsigned P> u32 im_scalar() {
 static const u32 v= smont<P>(spow<P>(groot<P>(), u64(P - 1) / 4));
 return v;
}
template <unsigned P> u32 ii_scalar() {
 static const u32 v= smont<P>(spow<P>(groot<P>(), u64(P - 1) - u64(P - 1) / 4));
 return v;
}
// n^{-1} mod P。n = 2^e ごとに残す。
template <unsigned P> u32 inv_n_scalar(int n) {
 static thread_local u32 tab[26]= {};
 const int e= __builtin_ctz(unsigned(n));
 if (!tab[e]) tab[e]= spow<P>(unsigned(n), u64(P - 2));
 return tab[e];
}

// 段のひねりの表。長さ len の段では h = len/4 として、tw[j] = ω_len^j、tw[h + j] = ω_len^{2j}、tw[2h + j] = ω_len^{3j}。
// 前半 3h が w (Montgomery)、後半 3h が wi = w * ni mod 2^32 (商を vpmulld で作るため)。
// 長さと向きごとに 1 回だけ作って残す (小さい n を何度も変換するときの作り直しを避ける)。
template <unsigned P> struct scratch {
 static inline thread_local abuf<32> tw[2][26];
};
// inv = false で順変換 (ω_len^j)、true で逆変換 (ω_len^{-j})。
template <unsigned P> const u32* build_tw(int len, int h, bool inv) {
 auto& tab= scratch<P>::tw[inv ? 1 : 0][__builtin_ctz(unsigned(len))];
 if (tab.n == 6 * h) return tab.p;
 tab.assign(6 * h);
 const u32 g= groot<P>();
 const u32 wl= smont<P>(spow<P>(g, inv ? u64(P - 1) - u64(P - 1) / unsigned(len) : u64(P - 1) / unsigned(len)));
 const u32 c2= smul<P>(wl, wl), c3= smul<P>(c2, wl);
 u32 a1[8], a2[8], a3[8];
 a1[0]= a2[0]= a3[0]= mc<P>::r;
 for (int j= 1; j < 8; ++j) a1[j]= smul<P>(a1[j - 1], wl), a2[j]= smul<P>(a2[j - 1], c2), a3[j]= smul<P>(a3[j - 1], c3);
 V w1= vload<P>(a1), w2= vload<P>(a2), w3= vload<P>(a3);
 const V r1= vsplat<P>(smul<P>(a1[7], wl));
 const V r2= vsplat<P>(smul<P>(a2[7], c2));
 const V r3= vsplat<P>(smul<P>(a3[7], c3));
 u32* p1= tab.p;
 u32* p2= p1 + h;
 u32* p3= p2 + h;
 // 表の値は P 未満に正規化する。vmuls の結果は 2P 未満なので、そのまま入れると
 // 段の掛け算の前提 (w < P) が崩れ、遅延還元の値が 2P を超えて 2^32 の桁あふれを起こす。
 for (int j= 0; j < h; j+= 8) {
  vstore<P>(p1 + j, w1), vstore<P>(p2 + j, w2), vstore<P>(p3 + j, w3);
  w1= vcanon<P>(vmuls<P>(w1, r1)), w2= vcanon<P>(vmuls<P>(w2, r2)), w3= vcanon<P>(vmuls<P>(w3, r3));
 }
 u32* q1= tab.p + 3 * h;
 u32* q2= q1 + h;
 u32* q3= q2 + h;
 const u32 niv= mc<P>::ni;
 for (int j= 0; j < h; ++j) q1[j]= u32(p1[j] * niv), q2[j]= u32(p2[j] * niv), q3[j]= u32(p3[j] * niv);
 return tab.p;
}

// 順変換の 1 段 (len の段、h = len/4 が 8 以上)。
template <unsigned P> void dft_stage_vec(u32* d, int n, int len) {
 const int h= len>> 2, hv= len>> 5;
 const u32* tab= build_tw<P>(len, h, false);
 const u32* t1= tab;
 const u32* t2= tab + h;
 const u32* t3= tab + 2 * h;
 const u32* u1w= tab + 3 * h;
 const u32* u2w= u1w + h;
 const u32* u3w= u2w + h;
 const V& im= im_vec<P>();
 const V imw= vsplat<P>(im_scalar<P>() * mc<P>::ni);
 for (int b= 0; b < n; b+= len) {
  u32* q= d + b;
  for (int jv= 0; jv < hv; ++jv) {
   u32* p= q + (jv << 3);
   const int j= jv << 3;
   const V w1= vload<P>(t1 + j), w2= vload<P>(t2 + j), w3= vload<P>(t3 + j);
   const V g1= vload<P>(u1w + j), g2= vload<P>(u2w + j), g3= vload<P>(u3w + j);
   const V u0= vload<P>(p), u1= vload<P>(p + h), u2= vload<P>(p + 2 * h), u3= vload<P>(p + 3 * h);
   const V t0= vplus<P>(u0, u2), t2v= vplus<P>(u1, u3);
   const V t1v= vdiff<P>(u0, u2), t3v= vdiff<P>(u1, u3);
   const V l0= vlow<P>(t0), l1= vlow<P>(t1v), l2= vlow<P>(t2v);
   const V pi= vmults<P>(t3v, im, imw);
   const V s1= vplus<P>(l1, pi), s3= vdiff<P>(l1, pi);
   vstore<P>(p, vlow<P>(vplus<P>(l0, l2)));
   vstore<P>(p + h, vmult<P>(s1, w1, g1));
   vstore<P>(p + 2 * h, vmult<P>(vdiff<P>(l0, l2), w2, g2));
   vstore<P>(p + 3 * h, vmult<P>(s3, w3, g3));
  }
 }
}

// 逆変換の 1 段 (len の段、h = len/4 が 8 以上)。
template <unsigned P> void idft_stage_vec(u32* d, int n, int len) {
 const int h= len>> 2, hv= len>> 5;
 const u32* tab= build_tw<P>(len, h, true);
 const u32* t1= tab;
 const u32* t2= tab + h;
 const u32* t3= tab + 2 * h;
 const u32* u1w= tab + 3 * h;
 const u32* u2w= u1w + h;
 const u32* u3w= u2w + h;
 const V& ii= ii_vec<P>();
 const V imw= vsplat<P>(ii_scalar<P>() * mc<P>::ni);
 for (int b= 0; b < n; b+= len) {
  u32* q= d + b;
  for (int jv= 0; jv < hv; ++jv) {
   u32* p= q + (jv << 3);
   const int j= jv << 3;
   const V w1= vload<P>(t1 + j), w2= vload<P>(t2 + j), w3= vload<P>(t3 + j);
   const V g1= vload<P>(u1w + j), g2= vload<P>(u2w + j), g3= vload<P>(u3w + j);
   const V u0= vload<P>(p);
   const V y1= vmult<P>(vload<P>(p + h), w1, g1);
   const V y2= vmult<P>(vload<P>(p + 2 * h), w2, g2);
   const V y3= vmult<P>(vload<P>(p + 3 * h), w3, g3);
   const V t0= vlow<P>(vplus<P>(u0, y2)), t2v= vlow<P>(vdiff<P>(u0, y2));
   const V t1v= vlow<P>(vplus<P>(y1, y3)), t3v= vmults<P>(vdiff<P>(y1, y3), ii, imw);
   vstore<P>(p, vlow<P>(vplus<P>(t0, t1v)));
   vstore<P>(p + h, vlow<P>(vplus<P>(t2v, t3v)));
   vstore<P>(p + 2 * h, vlow<P>(vdiff<P>(t0, t1v)));
   vstore<P>(p + 3 * h, vlow<P>(vdiff<P>(t2v, t3v)));
  }
 }
}

// 最後の基数 2 の段 (h = 1)。8 個 (4 対) ずつ。順変換にも逆変換にも使える (1/2 は掛けない)。
template <unsigned P> void dft_stage_r2(u32* d, int n) {
 int i= 0;
 for (; i + 8 <= n; i+= 8) {
  const V v= vload<P>(d + i);
  const V w= _mm256_shuffle_epi32(v, _MM_SHUFFLE(2, 3, 0, 1));  // [a1,a0,a3,a2]
  const V s= vlow<P>(vplus<P>(v, w));                          // 偶数レーン: a+b
  const V t= vlow<P>(vdiff<P>(w, v));                          // 奇数レーン: a-b
  vstore<P>(d + i, _mm256_blend_epi32(s, t, 0xAA));
 }
 for (; i < n; i+= 2) {
  const u32 u= d[i], v= d[i + 1];
  d[i]= slow<P>(u + v), d[i + 1]= slow<P>(sdiff<P>(u, v));
 }
}

// len = 4 の段 (h = 1)。2 ブロック (8 個) ずつ。128 bit レーンの中の 4 個で DFT4 をする。
// s = [u0+u2, u1+u3]、dd = [u0-u2, u1-u3] として、out0 = s0+s1、out2 = s0-s1、out1 = dd0+i dd1、out3 = dd0-i dd1 (ひねりは 1)。
// 値はレーン 0 か全レーンにできるので、最後にレーン 0〜3 へ並べ直す。順変換。
template <unsigned P> void dft_stage_h1(u32* d, int n) {
 const V& im= im_vec<P>();
 for (int b= 0; b < n; b+= 8) {
  const V v= vload<P>(d + b);
  const V sw= _mm256_shuffle_epi32(v, _MM_SHUFFLE(1, 0, 3, 2));  // [e2,e3,e0,e1]
  const V s= vlow<P>(vplus<P>(v, sw));                          // [t0,t1,t0,t1] (t0 = e0+e2、t1 = e1+e3)
  const V dd= vlow<P>(vdiff<P>(v, sw));                         // [t2,t3,-t2,-t3] (t2 = e0-e2、t3 = e1-e3)
  const V p= _mm256_shuffle_epi32(s, _MM_SHUFFLE(2, 3, 0, 1));   // [t1,t0,t1,t0]
  const V q= _mm256_shuffle_epi32(dd, _MM_SHUFFLE(2, 3, 0, 1));  // [t3,t2,t3,t2]
  const V iq= vmuls<P>(q, im);                                  // [i t3, i t2, ..]
  const V a0= vlow<P>(vplus<P>(s, p));                          // 全レーン out0 = t0+t1
  const V a1= _mm256_shuffle_epi32(vlow<P>(vplus<P>(dd, iq)), _MM_SHUFFLE(0, 0, 0, 0));   // out1 = t2+i t3
  const V a2= _mm256_shuffle_epi32(vlow<P>(vdiff<P>(s, p)), _MM_SHUFFLE(0, 0, 0, 0));      // out2 = t0-t1
  const V a3= _mm256_shuffle_epi32(vlow<P>(vdiff<P>(dd, iq)), _MM_SHUFFLE(0, 0, 0, 0));    // out3 = t2-i t3
  vstore<P>(d + b, _mm256_blend_epi32(_mm256_blend_epi32(a0, a1, 0x22), _mm256_blend_epi32(a2, a3, 0x88), 0xCC));
 }
}

// len = 4 の段 (h = 1)。逆変換。
template <unsigned P> void idft_stage_h1(u32* d, int n) {
 const V& ii= ii_vec<P>();
 for (int b= 0; b < n; b+= 8) {
  const V v= vload<P>(d + b);
  const V sw= _mm256_shuffle_epi32(v, _MM_SHUFFLE(1, 0, 3, 2));  // [e2,e3,e0,e1]
  const V s= vlow<P>(vplus<P>(v, sw));                          // [t0,t1,t0,t1]
  const V dd= vlow<P>(vdiff<P>(v, sw));                         // [t2,t3,-t2,-t3]
  const V p= _mm256_shuffle_epi32(s, _MM_SHUFFLE(2, 3, 0, 1));   // [t1,t0,t1,t0]
  const V q= _mm256_shuffle_epi32(dd, _MM_SHUFFLE(2, 3, 0, 1));  // [t3,t2,t3,t2]
  const V iq= vmuls<P>(q, ii);                                  // [ii t3, ii t2, ..]
  const V a0= vlow<P>(vplus<P>(s, p));                          // 全レーン x0 = t0+t1
  const V a1= _mm256_shuffle_epi32(vlow<P>(vplus<P>(dd, iq)), _MM_SHUFFLE(0, 0, 0, 0));   // x1 = t2+ii t3
  const V a2= _mm256_shuffle_epi32(vlow<P>(vdiff<P>(s, p)), _MM_SHUFFLE(0, 0, 0, 0));      // x2 = t0-t1
  const V a3= _mm256_shuffle_epi32(vlow<P>(vdiff<P>(dd, iq)), _MM_SHUFFLE(0, 0, 0, 0));    // x3 = t2-ii t3
  vstore<P>(d + b, _mm256_blend_epi32(_mm256_blend_epi32(a0, a1, 0x22), _mm256_blend_epi32(a2, a3, 0x88), 0xCC));
 }
}

// 小さい段 (h <= 4) のひねり (Montgomery) と、i / i^{-1} のスカラー版。長さと向きごとに 1 回だけ作る。
template <unsigned P> struct small_cache {
 static inline thread_local u32 a[2][6][3][4]= {};
 static inline thread_local bool built[2][6]= {};
 static inline thread_local u32 i4[2]= {0, 0};
 static inline thread_local bool built_i= false;
};
template <unsigned P> const u32* small_tw(int len, int h, bool inv) {
 using C= small_cache<P>;
 const int e= __builtin_ctz(unsigned(len));
 if (!C::built[inv ? 1 : 0][e]) {
  const u32 g= groot<P>();
  const u32 wl= smont<P>(spow<P>(g, inv ? u64(P - 1) - u64(P - 1) / unsigned(len) : u64(P - 1) / unsigned(len)));
  const u32 c2= smul<P>(wl, wl), c3= smul<P>(c2, wl);
  u32* t= &C::a[inv ? 1 : 0][e][0][0];
  u32* a1= t;
  u32* a2= t + 4;
  u32* a3= t + 8;
  a1[0]= a2[0]= a3[0]= mc<P>::r;
  for (int j= 1; j < 4; ++j) a1[j]= smul<P>(a1[j - 1], wl), a2[j]= smul<P>(a2[j - 1], c2), a3[j]= smul<P>(a3[j - 1], c3);
  C::built[inv ? 1 : 0][e]= true;
 }
 (void)h;
 return &C::a[inv ? 1 : 0][e][0][0];
}
template <unsigned P> u32 small_i(bool inv_i) {
 using C= small_cache<P>;
 if (!C::built_i) {
  const u32 g= groot<P>();
  C::i4[0]= smont<P>(spow<P>(g, u64(P - 1) / 4));
  C::i4[1]= smont<P>(spow<P>(g, u64(P - 1) - u64(P - 1) / 4));
  C::built_i= true;
 }
 return C::i4[inv_i ? 1 : 0];
}

// h = len/4 が 4 以下の段 (スカラー)。順変換。
template <unsigned P> void dft_stage_small(u32* d, int n, int len) {
 const int h= len>> 2;
 const u32* tw= small_tw<P>(len, h, false);
 const u32* a1= tw;
 const u32* a2= tw + 4;
 const u32* a3= tw + 8;
 const u32 im= small_i<P>(false);
 for (int b= 0; b < n; b+= len) {
  u32* q= d + b;
  for (int j= 0; j < h; ++j) {
   const u32 u0= q[j], u1= q[j + h], u2= q[j + 2 * h], u3= q[j + 3 * h];
   const u32 t0= u0 + u2, t2= u1 + u3;
   const u32 t1= sdiff<P>(u0, u2), t3= sdiff<P>(u1, u3);
   const u32 l0= slow<P>(t0), l1= slow<P>(t1), l2= slow<P>(t2);
   const u32 pi= smul<P>(t3, im);
   const u32 s1= l1 + pi, s3= sdiff<P>(l1, pi);
   q[j]= slow<P>(l0 + l2);
   q[j + h]= smul<P>(s1, a1[j]);
   q[j + 2 * h]= smul<P>(sdiff<P>(l0, l2), a2[j]);
   q[j + 3 * h]= smul<P>(s3, a3[j]);
  }
 }
}

// h = len/4 が 4 以下の段 (スカラー)。逆変換。
template <unsigned P> void idft_stage_small(u32* d, int n, int len) {
 const int h= len>> 2;
 const u32* tw= small_tw<P>(len, h, true);
 const u32* a1= tw;
 const u32* a2= tw + 4;
 const u32* a3= tw + 8;
 const u32 ii= small_i<P>(true);
 for (int b= 0; b < n; b+= len) {
  u32* q= d + b;
  for (int j= 0; j < h; ++j) {
   const u32 u0= q[j];
   const u32 y1= smul<P>(q[j + h], a1[j]), y2= smul<P>(q[j + 2 * h], a2[j]), y3= smul<P>(q[j + 3 * h], a3[j]);
   const u32 t0= slow<P>(u0 + y2), t2= slow<P>(sdiff<P>(u0, y2));
   const u32 t1= slow<P>(y1 + y3), t3= smul<P>(sdiff<P>(y1, y3), ii);
   q[j]= slow<P>(t0 + t1);
   q[j + h]= slow<P>(t2 + t3);
   q[j + 2 * h]= slow<P>(sdiff<P>(t0, t1));
   q[j + 3 * h]= slow<P>(sdiff<P>(t2, t3));
  }
 }
}

// len = 16 の段 (h = 4)。2 ブロック (32 個) ずつ。1 ブロックは 8 レーン 2 本 ([u0|u1] と [u2|u3]) なので、
// 2 ブロック分を 128 bit ごとに組み替えて、u0+u2 などを [b0|b1] の 8 レーンにする。順変換。
// 1 ブロックだけのタイル (n = 16) はスカラーに回す。
template <unsigned P> void dft_stage_h4(u32* d, int n) {
 if (n < 32) {
  dft_stage_small<P>(d, n, 16);
  return;
 }
 const u32* tw= small_tw<P>(16, 4, false);
 const V w1= vcombine128<P>(vload128<P>(tw), vload128<P>(tw));
 const V w2= vcombine128<P>(vload128<P>(tw + 4), vload128<P>(tw + 4));
 const V w3= vcombine128<P>(vload128<P>(tw + 8), vload128<P>(tw + 8));
 const V& im= im_vec<P>();
 for (int b= 0; b + 32 <= n; b+= 32) {
  u32* q= d + b;
  const V x0= vload<P>(q), x1= vload<P>(q + 8), x2= vload<P>(q + 16), x3= vload<P>(q + 24);
  const V s0= vplus<P>(x0, x1), s1= vplus<P>(x2, x3);
  const V d0= vdiff<P>(x0, x1), d1= vdiff<P>(x2, x3);
  const V t0= vlow<P>(_mm256_permute2x128_si256(s0, s1, 0x20));  // u0+u2 (b0 | b1)
  const V t2= vlow<P>(_mm256_permute2x128_si256(s0, s1, 0x31));  // u1+u3
  const V t1= vlow<P>(_mm256_permute2x128_si256(d0, d1, 0x20));  // u0-u2
  const V t3= _mm256_permute2x128_si256(d0, d1, 0x31);           // u1-u3
  const V pi= vmuls<P>(t3, im);
  const V a1= vplus<P>(t1, pi), a3= vdiff<P>(t1, pi);
  const V o0= vlow<P>(vplus<P>(t0, t2));
  const V o1= vmul<P>(a1, w1);
  const V o2= vmul<P>(vdiff<P>(t0, t2), w2);
  const V o3= vmul<P>(a3, w3);
  vstore<P>(q, _mm256_permute2x128_si256(o0, o1, 0x20));
  vstore<P>(q + 8, _mm256_permute2x128_si256(o2, o3, 0x20));
  vstore<P>(q + 16, _mm256_permute2x128_si256(o0, o1, 0x31));
  vstore<P>(q + 24, _mm256_permute2x128_si256(o2, o3, 0x31));
 }
}

// len = 16 の段 (h = 4)。逆変換。u0 には Montgomery の 1 を掛けて、u0 と y1 を 1 本にそろえる。
template <unsigned P> void idft_stage_h4(u32* d, int n) {
 if (n < 32) {
  idft_stage_small<P>(d, n, 16);
  return;
 }
 const u32* tw= small_tw<P>(16, 4, true);
 const V wa= vcombine128<P>(vsplat128<P>(mc<P>::r), vload128<P>(tw));
 const V wb= vcombine128<P>(vload128<P>(tw + 4), vload128<P>(tw + 8));
 const V& ii= ii_vec<P>();
 for (int b= 0; b + 32 <= n; b+= 32) {
  u32* q= d + b;
  const V x0= vload<P>(q), x1= vload<P>(q + 8), x2= vload<P>(q + 16), x3= vload<P>(q + 24);
  const V y0= vmul<P>(x0, wa), y1= vmul<P>(x1, wb);
  const V y2= vmul<P>(x2, wa), y3= vmul<P>(x3, wb);
  const V s0= vplus<P>(y0, y1), s1= vplus<P>(y2, y3);
  const V d0= vdiff<P>(y0, y1), d1= vdiff<P>(y2, y3);
  const V t0= vlow<P>(_mm256_permute2x128_si256(s0, s1, 0x20));  // u0+y2
  const V t1= vlow<P>(_mm256_permute2x128_si256(s0, s1, 0x31));  // y1+y3
  const V t2= vlow<P>(_mm256_permute2x128_si256(d0, d1, 0x20));  // u0-y2
  const V t3= vmuls<P>(_mm256_permute2x128_si256(d0, d1, 0x31), ii);  // (y1-y3) i^{-1}
  const V o0= vlow<P>(vplus<P>(t0, t1));
  const V o1= vlow<P>(vplus<P>(t2, t3));
  const V o2= vlow<P>(vdiff<P>(t0, t1));
  const V o3= vlow<P>(vdiff<P>(t2, t3));
  vstore<P>(q, _mm256_permute2x128_si256(o0, o1, 0x20));
  vstore<P>(q + 8, _mm256_permute2x128_si256(o2, o3, 0x20));
  vstore<P>(q + 16, _mm256_permute2x128_si256(o0, o1, 0x31));
  vstore<P>(q + 24, _mm256_permute2x128_si256(o2, o3, 0x31));
 }
}

// len = 8 の段 (h = 2)。2 ブロック (16 個) ずつ。1 ブロックが 8 レーン 1 本 ([u0|u1|u2|u3]、各 2 レーン) なので、
// 128 bit の入れ替えで u0+u2 と u1+u3 を同じ 128 bit にそろえ、2 ブロック分を 64 bit ずつ組にする。順変換。
// 1 ブロックだけのタイル (n = 8) はスカラーに回す。
template <unsigned P> void dft_stage_h2(u32* d, int n) {
 if (n < 16) {
  dft_stage_small<P>(d, n, 8);
  return;
 }
 const u32* tw= small_tw<P>(8, 2, false);
 const V w1= _mm256_set_epi32(tw[1], tw[0], tw[1], tw[0], tw[1], tw[0], tw[1], tw[0]);
 const V w2= _mm256_set_epi32(tw[5], tw[4], tw[5], tw[4], tw[5], tw[4], tw[5], tw[4]);
 const V w3= _mm256_set_epi32(tw[9], tw[8], tw[9], tw[8], tw[9], tw[8], tw[9], tw[8]);
 const V& im= im_vec<P>();
 for (int b= 0; b + 16 <= n; b+= 16) {
  u32* q= d + b;
  const V x0= vload<P>(q), x1= vload<P>(q + 8);
  const V z0= _mm256_permute2x128_si256(x0, x0, 0x01), z1= _mm256_permute2x128_si256(x1, x1, 0x01);
  const V s0= vplus<P>(x0, z0), s1= vplus<P>(x1, z1);
  const V d0= vdiff<P>(x0, z0), d1= vdiff<P>(x1, z1);
  const V t0= vlow<P>(_mm256_unpacklo_epi64(s0, s1));  // u0+u2 (b0 | b1)
  const V t2= vlow<P>(_mm256_unpackhi_epi64(s0, s1));  // u1+u3
  const V t1= vlow<P>(_mm256_unpacklo_epi64(d0, d1));  // u0-u2
  const V t3= _mm256_unpackhi_epi64(d0, d1);           // u1-u3
  const V pi= vmuls<P>(t3, im);
  const V a1= vplus<P>(t1, pi), a3= vdiff<P>(t1, pi);
  const V o0= vlow<P>(vplus<P>(t0, t2));
  const V o1= vmul<P>(a1, w1);
  const V o2= vmul<P>(vdiff<P>(t0, t2), w2);
  const V o3= vmul<P>(a3, w3);
  // o0 は [b0 | b1 | b0 | b1] (64 bit ごと) なので、b0 は unpacklo、b1 は unpackhi から取る。
  const V u01= _mm256_unpacklo_epi64(o0, o1), u23= _mm256_unpacklo_epi64(o2, o3);
  const V v01= _mm256_unpackhi_epi64(o0, o1), v23= _mm256_unpackhi_epi64(o2, o3);
  vstore<P>(q, _mm256_permute2x128_si256(u01, u23, 0x20));
  vstore<P>(q + 8, _mm256_permute2x128_si256(v01, v23, 0x20));
 }
}

// len = 8 の段 (h = 2)。逆変換。
template <unsigned P> void idft_stage_h2(u32* d, int n) {
 if (n < 16) {
  idft_stage_small<P>(d, n, 8);
  return;
 }
 const u32* tw= small_tw<P>(8, 2, true);
 const V mv= _mm256_set_epi32(tw[9], tw[8], tw[5], tw[4], tw[1], tw[0], int(mc<P>::r), int(mc<P>::r));
 const V& ii= ii_vec<P>();
 for (int b= 0; b + 16 <= n; b+= 16) {
  u32* q= d + b;
  const V x0= vload<P>(q), x1= vload<P>(q + 8);
  const V y0= vmul<P>(x0, mv), y1= vmul<P>(x1, mv);
  const V z0= _mm256_permute2x128_si256(y0, y0, 0x01), z1= _mm256_permute2x128_si256(y1, y1, 0x01);
  const V s0= vplus<P>(y0, z0), s1= vplus<P>(y1, z1);
  const V d0= vdiff<P>(y0, z0), d1= vdiff<P>(y1, z1);
  const V t0= vlow<P>(_mm256_unpacklo_epi64(s0, s1));
  const V t1= vlow<P>(_mm256_unpackhi_epi64(s0, s1));
  const V t2= vlow<P>(_mm256_unpacklo_epi64(d0, d1));
  const V t3= vmuls<P>(_mm256_unpackhi_epi64(d0, d1), ii);
  const V o0= vlow<P>(vplus<P>(t0, t1));
  const V o1= vlow<P>(vplus<P>(t2, t3));
  const V o2= vlow<P>(vdiff<P>(t0, t1));
  const V o3= vlow<P>(vdiff<P>(t2, t3));
  // o0 は [b0 | b1 | b0 | b1] (64 bit ごと) なので、b0 は unpacklo、b1 は unpackhi から取る。
  const V u01= _mm256_unpacklo_epi64(o0, o1), u23= _mm256_unpacklo_epi64(o2, o3);
  const V v01= _mm256_unpackhi_epi64(o0, o1), v23= _mm256_unpackhi_epi64(o2, o3);
  vstore<P>(q, _mm256_permute2x128_si256(u01, u23, 0x20));
  vstore<P>(q + 8, _mm256_permute2x128_si256(v01, v23, 0x20));
 }
}

// タイル (dft_tile 以下) の残りの段をまとめて回す。順変換。
template <unsigned P> void dft_stages_small(u32* d, int n) {
 int len= n;
 for (; len >= 32; len>>= 2) dft_stage_vec<P>(d, n, len);
 if (len == 16) dft_stage_h4<P>(d, n), len= 4;
 if (len == 8) dft_stage_h2<P>(d, n), len= 2;
 if (len == 4) {
  if (n >= 8) dft_stage_h1<P>(d, n);
  else dft_stage_small<P>(d, n, 4);
  len= 1;
 }
 if (len == 2) dft_stage_r2<P>(d, n);
}
// タイル (dft_tile 以下) の段を、小さい方から順に回す。逆変換。
template <unsigned P> void idft_stages_small(u32* d, int n) {
 const bool odd_e= __builtin_ctz(unsigned(n)) & 1;
 if (odd_e) {
  dft_stage_r2<P>(d, n);
  if (n >= 8) idft_stage_h2<P>(d, n);
  for (int len= 32; len <= n; len<<= 2) idft_stage_vec<P>(d, n, len);
 } else {
  if (n >= 8) idft_stage_h1<P>(d, n);
  else idft_stage_small<P>(d, n, 4);
  if (n >= 16) idft_stage_h4<P>(d, n);
  for (int len= 64; len <= n; len<<= 2) idft_stage_vec<P>(d, n, len);
 }
}

// DIF。a[i] に Σ_j a[j] ω^{j k} が入る。k は位置 i の 4 進の桁を逆にした数 (2 の指数が奇数なら最後の基数 2 の段で最上位の桁が最下位に来る)。
template <unsigned P> void to_mont(u32* d, int n) {
 const V r2v= vsplat<P>(mc<P>::r2);
 int i= 0;
 for (; i + 8 <= n; i+= 8) vstore<P>(d + i, vmuls<P>(vload<P>(d + i), r2v));
 for (; i < n; ++i) d[i]= smont<P>(d[i]);
}
template <unsigned P> void from_mont(u32* d, int n) {
 const V onev= vsplat<P>(1);
 int i= 0;
 for (; i + 8 <= n; i+= 8) vstore<P>(d + i, vmuls<P>(vload<P>(d + i), onev));
 for (; i < n; ++i) d[i]= smul<P>(d[i], 1);
}
// 段だけ。表現はそのまま (入力が Montgomery なら出力も Montgomery)。
template <unsigned P> void dft_stages(u32* d, int n) {
 if (n <= 1) return;
 // 大きい段 (タイルより長い段) は配列全体をなめる。タイルは段の長さに合わせて 4 で割って決める。
 int tile= n;
 while (tile > dft_tile) tile>>= 2;
 for (int len= n; len > tile; len>>= 2) dft_stage_vec<P>(d, n, len);
 // 残りの段はタイルごとに、タイルを L1 に載せたまま回す。
 for (int b= 0; b < n; b+= tile) dft_stages_small<P>(d + b, tile);
}
// 段だけ (逆変換)。段ごとの 1/4 は掛けない。
template <unsigned P> void idft_stages(u32* d, int n) {
 if (n <= 1) return;
 // 小さい段はタイルごとに、大きい段は配列全体をなめる。タイルは段の長さに合わせて 4 で割って決める。
 int tile= n;
 while (tile > dft_tile) tile>>= 2;
 for (int b= 0; b < n; b+= tile) idft_stages_small<P>(d + b, tile);
 for (int len= tile<< 2; len <= n; len<<= 2) idft_stage_vec<P>(d, n, len);
}

// Montgomery の値に c (ふつうの値) を掛ける。c = n^{-1} なら 1/n を掛けながらふつうの値に戻る。
template <unsigned P> void scale_by(u32* d, int n, u32 c) {
 const V cv= vsplat<P>(c);
 int i= 0;
 for (; i + 8 <= n; i+= 8) vstore<P>(d + i, vmuls<P>(vload<P>(d + i), cv));
 for (; i < n; ++i) d[i]= smul<P>(d[i], c);
}

// 以下、畳み込み (5 本目)。ひねりを全レーン同じにし、8 レーンを葉の 8 個の係数にして、段をベクトルの添字の上の
// DFT にする。葉の上の長さ 8 の DFT は、8 次の環 x^8 - w (w = η^{葉の位置の桁を逆にした数}) の積に置き換える。
// 値はふつうの表現のまま持ち、ひねりは Shoup の方法で掛ける。葉の積は Montgomery の MAC で、R^{-1} が 1 回付く。

// 32 bit の下位だけの積。intrinsics の vpmulld は下位 32 bit だけなので符号は問わない。
// SIMDe の mullo は符号付きで掛けて UBSan に鳴るので、SIMDe のときは 2 本の mul_epu32 で作る。
#ifndef NEO_NTT_SIMDE
template <unsigned P> V vpmulld(V a, V b) { return _mm256_mullo_epi32(a, b); }
#else
template <unsigned P> V vpmulld(V a, V b) {
 const V e= _mm256_mul_epu32(a, b);
 const V o= _mm256_mul_epu32(_mm256_srli_epi64(a, 32), _mm256_srli_epi64(b, 32));
 return _mm256_blend_epi32(e, _mm256_slli_epi64(o, 32), 0xAA);
}
#endif

// Shoup の掛け算 (w と wi は全レーン同じ)。x < 2^32、w < P、wi = floor(w 2^32 / P)。結果は 2P 未満。
template <unsigned P> V vshoups(V x, V w, V wi) {
 const V pv= vsplat<P>(P);
 const V qe= _mm256_mul_epu32(x, wi);
 const V qo= _mm256_mul_epu32(_mm256_srli_epi64(x, 32), wi);
 const V q= _mm256_blend_epi32(_mm256_srli_epi64(qe, 32), qo, 0xAA);
 return _mm256_sub_epi32(vpmulld<P>(x, w), vpmulld<P>(q, pv));
}

template <unsigned P> u32 shoup_q(u32 w) { return u32((u64(w) << 32) / P); }

template <unsigned P> u32 im_normal() {
 static const u32 v= spow<P>(groot<P>(), u64(P - 1) / 4);
 return v;
}
template <unsigned P> u32 ii_normal() {
 static const u32 v= spow<P>(groot<P>(), u64(P - 1) - u64(P - 1) / 4);
 return v;
}

// 段 (ベクトル数 L) のひねり。長さ m の DFT の、L の段の根は η_m^{m/L}。
// h = L/4 個ずつ、[w1|w2|w3|s1|s2|s3]。
template <unsigned P> struct stage_tables {
 static inline thread_local abuf<32> t[2][26];
};
template <unsigned P> const u32* stage_table(int m, int L, bool inv) {
 auto& tb= stage_tables<P>::t[inv ? 1 : 0][__builtin_ctz(unsigned(L))];
 const int h= L >> 2;
 if (tb.n != 6 * h) {
  tb.assign(6 * h);
  const u32 g= groot<P>();
  const u64 e= u64(P - 1) / unsigned(m) * unsigned(m / L);
  const u32 rl= spow<P>(g, inv ? u64(P - 1) - e : e);
  const u32 c2= u32(u64(rl) * rl % P), c3= u32(u64(c2) * rl % P);
  u32 a1= 1, a2= 1, a3= 1;
  u32* tp= tb.p;
  for (int j= 0; j < h; ++j) {
   tp[j]= a1, tp[h + j]= a2, tp[2 * h + j]= a3;
   tp[3 * h + j]= shoup_q<P>(a1), tp[4 * h + j]= shoup_q<P>(a2), tp[5 * h + j]= shoup_q<P>(a3);
   a1= u32(u64(a1) * rl % P), a2= u32(u64(a2) * c2 % P), a3= u32(u64(a3) * c3 % P);
  }
 }
 return tb.p;
}

// 葉の重み w = η^{rev(i)} の表 (m 個) と商。
template <unsigned P> struct leaf_tables {
 static inline thread_local abuf<32> t[26];
};
// i の桁 (4 進、奇数指数では最後が 2 進) を逆にした数。
template <unsigned P> int rev_index(int i, int e) {
 if (!(e & 1)) {
  int r= 0, x= i;
  for (int t= 0; t < e / 2; ++t) r= (r << 2) | (x & 3), x>>= 2;
  return r;
 }
 const int k= (e - 1) / 2;
 const int b= i & 1;
 int q[16], x= i >> 1;
 for (int t= 0; t < k; ++t) q[t]= x & 3, x>>= 2;
 int r= 0;
 for (int t= 0; t < k; ++t) r|= q[k - 1 - t] << (2 * t);
 return r | (b << (2 * k));
}
template <unsigned P> const u32* leaf_table(int m) {
 const int e= __builtin_ctz(unsigned(m));
 auto& tb= leaf_tables<P>::t[e];
 if (tb.n != 2 * m) {
  tb.assign(2 * m);
  const u32 g= groot<P>();
  const u32 eta= spow<P>(g, u64(P - 1) / unsigned(m));
  u32* w= tb.p;
  u32* q= w + m;
  // w[i] = eta^{rev(i)} を、i を 1 進めるごとの比で作る (1 つずつ pow を取らない)。
  // 比は、i の末尾の最大桁 (4 進の 3、2 進なら 1) の数 s で決まる。
  if (!(e & 1)) {
   const int k= e / 2;
   std::vector<u32> step(k);
   for (int s= 0; s < k; ++s) {
    const u64 d1= u64(1) << (2 * (k - 1 - s)), d2= u64(1) << (2 * (k - s));
    step[s]= spow<P>(eta, (d1 + d2) % unsigned(m));
   }
   u32 cur= 1;
   w[0]= 1;
   for (int i= 0; i < m - 1; ++i) {
    int s= __builtin_ctz(~unsigned(i)) >> 1;
    if (s >= k) s= k - 1;
    cur= u32(u64(cur) * step[s] % P);
    w[i + 1]= cur;
   }
  } else {
   const int k= (e - 1) / 2;
   const u32 half= spow<P>(eta, u64(m) / 2);  // eta^{4^k} = -1。b が 0 から 1 に進むときの比
   std::vector<u32> step(k);
   for (int s= 0; s < k; ++s) {
    const u64 d1= u64(1) << (2 * (k - 1 - s)), d2= u64(1) << (2 * (k - s));
    step[s]= spow<P>(eta, (d1 + d2) % unsigned(m));
   }
   u32 cur= 1;
   w[0]= 1;
   for (int i= 0; i < m - 1; ++i) {
    if (!(i & 1)) cur= u32(u64(cur) * half % P);
    else {
     int s= __builtin_ctz(~unsigned(i >> 1)) >> 1;
     if (s >= k) s= k - 1;
     cur= u32(u64(cur) * step[s] % P);
    }
    w[i + 1]= cur;
   }
  }
  for (int i= 0; i < m; ++i) q[i]= shoup_q<P>(w[i]);
 }
 return tb.p;
}

// 段 (ベクトル数 4) はひねりがすべて 1 なので、掛け算を省く。
template <unsigned P> void stage4v_h1(u32* d, int nv, bool inv) {
 const u32 imn= inv ? ii_normal<P>() : im_normal<P>();
 const V im= vsplat<P>(imn), imq= vsplat<P>(shoup_q<P>(imn));
 for (int b= 0; b < nv; b+= 4) {
  u32* p= d + 8 * b;
  const V u0= vload<P>(p), u1= vload<P>(p + 8), u2= vload<P>(p + 16), u3= vload<P>(p + 24);
  if (!inv) {
   const V t0= vplus<P>(u0, u2), t2= vplus<P>(u1, u3);
   const V t1= vdiff<P>(u0, u2), t3= vdiff<P>(u1, u3);
   const V l0= vlow<P>(t0), l1= vlow<P>(t1), l2= vlow<P>(t2);
   const V pi= vshoups<P>(vlow<P>(t3), im, imq);
   vstore<P>(p, vlow<P>(vplus<P>(l0, l2)));
   vstore<P>(p + 8, vlow<P>(vplus<P>(l1, pi)));
   vstore<P>(p + 16, vlow<P>(vdiff<P>(l0, l2)));
   vstore<P>(p + 24, vlow<P>(vdiff<P>(l1, pi)));
  } else {
   const V y2= vlow<P>(u2);
   const V y1= vlow<P>(u1), y3= vlow<P>(u3);
   const V t0= vlow<P>(vplus<P>(u0, y2)), t2= vlow<P>(vdiff<P>(u0, y2));
   const V t1= vlow<P>(vplus<P>(y1, y3)), t3= vshoups<P>(vlow<P>(vdiff<P>(y1, y3)), im, imq);
   vstore<P>(p, vlow<P>(vplus<P>(t0, t1)));
   vstore<P>(p + 8, vlow<P>(vplus<P>(t2, t3)));
   vstore<P>(p + 16, vlow<P>(vdiff<P>(t0, t1)));
   vstore<P>(p + 24, vlow<P>(vdiff<P>(t2, t3)));
  }
 }
}

// 段 (ベクトル数 L) の順変換 (DIF) と逆変換 (DIT)。L は 4 以上。
template <unsigned P> void stage4v(u32* d, int nv, int m, int L, bool inv) {
 const int h= L >> 2;
 if (h == 1) {
  stage4v_h1<P>(d, nv, inv);
  return;
 }
 const u32* t= stage_table<P>(m, L, inv);
 const u32 *w1= t, *w2= t + h, *w3= t + 2 * h;
 const u32 *s1= t + 3 * h, *s2= t + 4 * h, *s3= t + 5 * h;
 const u32 imn= inv ? ii_normal<P>() : im_normal<P>();
 const V im= vsplat<P>(imn), imq= vsplat<P>(shoup_q<P>(imn));
 for (int b= 0; b < nv; b+= L) {
  u32* p= d + 8 * b;
  for (int j= 0; j < h; ++j) {
   const V W1= vsplat<P>(w1[j]), W2= vsplat<P>(w2[j]), W3= vsplat<P>(w3[j]);
   const V S1= vsplat<P>(s1[j]), S2= vsplat<P>(s2[j]), S3= vsplat<P>(s3[j]);
   if (!inv) {
    const V u0= vload<P>(p + 8 * j), u1= vload<P>(p + 8 * (j + h)), u2= vload<P>(p + 8 * (j + 2 * h)), u3= vload<P>(p + 8 * (j + 3 * h));
    const V t0= vplus<P>(u0, u2), t2= vplus<P>(u1, u3);
    const V t1= vdiff<P>(u0, u2), t3= vdiff<P>(u1, u3);
    const V l0= vlow<P>(t0), l1= vlow<P>(t1), l2= vlow<P>(t2);
    const V pi= vshoups<P>(vlow<P>(t3), im, imq);
    const V a1= vplus<P>(l1, pi), a3= vdiff<P>(l1, pi);
    vstore<P>(p + 8 * j, vlow<P>(vplus<P>(l0, l2)));
    vstore<P>(p + 8 * (j + h), vshoups<P>(a1, W1, S1));
    vstore<P>(p + 8 * (j + 2 * h), vshoups<P>(vlow<P>(vdiff<P>(l0, l2)), W2, S2));
    vstore<P>(p + 8 * (j + 3 * h), vshoups<P>(a3, W3, S3));
   } else {
    const V y1= vshoups<P>(vload<P>(p + 8 * (j + h)), W1, S1);
    const V y2= vshoups<P>(vload<P>(p + 8 * (j + 2 * h)), W2, S2);
    const V y3= vshoups<P>(vload<P>(p + 8 * (j + 3 * h)), W3, S3);
    const V u0= vload<P>(p + 8 * j);
    const V t0= vlow<P>(vplus<P>(u0, y2)), t2= vlow<P>(vdiff<P>(u0, y2));
    const V t1= vlow<P>(vplus<P>(y1, y3)), t3= vshoups<P>(vlow<P>(vdiff<P>(y1, y3)), im, imq);
    vstore<P>(p + 8 * j, vlow<P>(vplus<P>(t0, t1)));
    vstore<P>(p + 8 * (j + h), vlow<P>(vplus<P>(t2, t3)));
    vstore<P>(p + 8 * (j + 2 * h), vlow<P>(vdiff<P>(t0, t1)));
    vstore<P>(p + 8 * (j + 3 * h), vlow<P>(vdiff<P>(t2, t3)));
   }
  }
 }
}

// 段 (ベクトル数 L) の、長さをコンパイル時に決めた版 (タイルの中の段)。
template <unsigned P, int L> void stage4v_c(u32* d, int nv, int m, bool inv) {
 constexpr int h= L >> 2;
 if constexpr (h == 1) {
  stage4v_h1<P>(d, nv, inv);
  return;
 }
 const u32* t= stage_table<P>(m, L, inv);
 const u32 *w1= t, *w2= t + h, *w3= t + 2 * h;
 const u32 *s1= t + 3 * h, *s2= t + 4 * h, *s3= t + 5 * h;
 const u32 imn= inv ? ii_normal<P>() : im_normal<P>();
 const V im= vsplat<P>(imn), imq= vsplat<P>(shoup_q<P>(imn));
 for (int b= 0; b < nv; b+= L) {
  u32* p= d + 8 * b;
#pragma GCC unroll 2
  for (int j= 0; j < h; ++j) {
   const V W1= vsplat<P>(w1[j]), W2= vsplat<P>(w2[j]), W3= vsplat<P>(w3[j]);
   const V S1= vsplat<P>(s1[j]), S2= vsplat<P>(s2[j]), S3= vsplat<P>(s3[j]);
   if (!inv) {
    const V u0= vload<P>(p + 8 * j), u1= vload<P>(p + 8 * (j + h)), u2= vload<P>(p + 8 * (j + 2 * h)), u3= vload<P>(p + 8 * (j + 3 * h));
    const V t0= vplus<P>(u0, u2), t2= vplus<P>(u1, u3);
    const V t1= vdiff<P>(u0, u2), t3= vdiff<P>(u1, u3);
    const V l0= vlow<P>(t0), l1= vlow<P>(t1), l2= vlow<P>(t2);
    const V pi= vshoups<P>(vlow<P>(t3), im, imq);
    const V a1= vplus<P>(l1, pi), a3= vdiff<P>(l1, pi);
    vstore<P>(p + 8 * j, vlow<P>(vplus<P>(l0, l2)));
    vstore<P>(p + 8 * (j + h), vshoups<P>(a1, W1, S1));
    vstore<P>(p + 8 * (j + 2 * h), vshoups<P>(vlow<P>(vdiff<P>(l0, l2)), W2, S2));
    vstore<P>(p + 8 * (j + 3 * h), vshoups<P>(a3, W3, S3));
   } else {
    const V y1= vshoups<P>(vload<P>(p + 8 * (j + h)), W1, S1);
    const V y2= vshoups<P>(vload<P>(p + 8 * (j + 2 * h)), W2, S2);
    const V y3= vshoups<P>(vload<P>(p + 8 * (j + 3 * h)), W3, S3);
    const V u0= vload<P>(p + 8 * j);
    const V t0= vlow<P>(vplus<P>(u0, y2)), t2= vlow<P>(vdiff<P>(u0, y2));
    const V t1= vlow<P>(vplus<P>(y1, y3)), t3= vshoups<P>(vlow<P>(vdiff<P>(y1, y3)), im, imq);
    vstore<P>(p + 8 * j, vlow<P>(vplus<P>(t0, t1)));
    vstore<P>(p + 8 * (j + h), vlow<P>(vplus<P>(t2, t3)));
    vstore<P>(p + 8 * (j + 2 * h), vlow<P>(vdiff<P>(t0, t1)));
    vstore<P>(p + 8 * (j + 3 * h), vlow<P>(vdiff<P>(t2, t3)));
   }
  }
 }
}

// 段 (ベクトル数 2) の順変換と逆変換。ひねりは 1。
template <unsigned P> void stage2v(u32* d, int nv) {
 for (int b= 0; b < nv; b+= 2) {
  const V u0= vload<P>(d + 8 * b), u1= vload<P>(d + 8 * (b + 1));
  vstore<P>(d + 8 * b, vlow<P>(vplus<P>(u0, u1)));
  vstore<P>(d + 8 * (b + 1), vlow<P>(vdiff<P>(u0, u1)));
 }
}

// タイル (NV ベクトル) の段の列。L = NV, NV/4, …, 4 (2 の指数が奇数なら 8 のあと 2)。
template <unsigned P, int NV, int L> struct tile_seq {
 static void fwd(u32* d, int m) { stage4v_c<P, L>(d, NV, m, false); tile_seq<P, NV, L / 4>::fwd(d, m); }
 static void inv(u32* d, int m) { tile_seq<P, NV, L / 4>::inv(d, m); stage4v_c<P, L>(d, NV, m, true); }
};
template <unsigned P, int NV> struct tile_seq<P, NV, 1> {
 static void fwd(u32*, int) {}
 static void inv(u32*, int) {}
};
template <unsigned P, int NV> struct tile_seq<P, NV, 2> {
 static void fwd(u32* d, int) { stage2v<P>(d, NV); }
 static void inv(u32* d, int) { stage2v<P>(d, NV); }
};
template <unsigned P, int NV> struct tile_seq<P, NV, 4> {
 static void fwd(u32* d, int) { stage4v_h1<P>(d, NV, false); }
 static void inv(u32* d, int) { stage4v_h1<P>(d, NV, true); }
};

// 葉の積 (mod x^8 - w)。窓 [w a, a] を 16 個の作業に置き、8 回の読み出しで回転を作る。
// 窓のストアと読み出しが近いと待つので、呼び出し側で 1 つ先の窓を作ってから今の積を取る。
// 積は Montgomery の MAC で、R^{-1} が 1 回付く。奇数レーンの積は 1 つ前の窓の読み出しの偶数レーンを使う。
template <unsigned P> void leaf_build(u32* buf, u32* cb, V ya, V yb, V wv, V wq) {
 vstore<P>(buf, vshoups<P>(ya, wv, wq));
 vstore<P>(buf + 8, ya);
 vstore<P>(cb, yb);
}
template <unsigned P> V leaf_mac(const u32* buf, const u32* cb) {
 const V pv= vsplat<P>(P), niv= vsplat<P>(mc<P>::ni);
 V prev= vload<P>(buf + 8);
 V e, o;
 {
  const V y0= vsplat<P>(cb[0]);
  e= _mm256_mul_epu32(prev, y0);
  o= _mm256_mul_epu32(_mm256_srli_epi64(prev, 32), y0);
 }
 for (int i= 1; i < 8; ++i) {
  const V x= vload<P>(buf + 8 - i);
  const V yi= vsplat<P>(cb[i]);
  e= _mm256_add_epi64(e, _mm256_mul_epu32(x, yi));
  o= _mm256_add_epi64(o, _mm256_mul_epu32(prev, yi));
  prev= x;
 }
 e= _mm256_add_epi64(e, _mm256_mul_epu32(_mm256_mul_epu32(e, niv), pv));
 o= _mm256_add_epi64(o, _mm256_mul_epu32(_mm256_mul_epu32(o, niv), pv));
 return _mm256_blend_epi32(_mm256_srli_epi64(e, 32), o, 0xAA);
}

// タイル (NV ベクトル) を、順変換の下の段、葉の積、逆変換の下の段まで回す。段の長さはコンパイル時に決まる。
template <unsigned P, int NV> void tile_work_c(u32* fa, u32* fb, int m, int B, bool square, V sc, V scq) {
 tile_seq<P, NV, NV>::fwd(fa, m);
 if (!square) tile_seq<P, NV, NV>::fwd(fb, m);
 const u32* lw= leaf_table<P>(m);
 const u32* lq= lw + m;
 // 窓は 1 つ先に作る (ストアした窓をすぐ読むと待つため)。
 alignas(32) u32 buf[2][16], cb[2][8];
 auto build= [&](int i, int slot) {
  const V ya= vcanon<P>(vload<P>(fa + 8 * i));
  const V yb= square ? ya : vcanon<P>(vload<P>(fb + 8 * i));
  leaf_build<P>(buf[slot], cb[slot], ya, yb, vsplat<P>(lw[B + i]), vsplat<P>(lq[B + i]));
 };
 build(0, 0);
 for (int i= 0; i < NV; ++i) {
  if (i + 1 < NV) build(i + 1, (i + 1) & 1);
  const V c= leaf_mac<P>(buf[i & 1], cb[i & 1]);
  vstore<P>(fa + 8 * i, vlow<P>(vshoups<P>(c, sc, scq)));
 }
 tile_seq<P, NV, NV>::inv(fa, m);
}

// タイル (nv ベクトル、全体では m ベクトル、先頭はベクトル B) を回す。nv は 16, 32, …, 256 のどれか。
template <unsigned P> void tile_work(u32* fa, u32* fb, int m, int nv, int B, bool square, V sc, V scq) {
 switch (nv) {
  case 16: tile_work_c<P, 16>(fa, fb, m, B, square, sc, scq); break;
  case 32: tile_work_c<P, 32>(fa, fb, m, B, square, sc, scq); break;
  case 64: tile_work_c<P, 64>(fa, fb, m, B, square, sc, scq); break;
  case 128: tile_work_c<P, 128>(fa, fb, m, B, square, sc, scq); break;
  case 256: tile_work_c<P, 256>(fa, fb, m, B, square, sc, scq); break;
  default: std::abort();
 }
}

// 深さ優先で区画を回す。
template <unsigned P> void rec(u32* fa, u32* fb, int nv, int m, int B, bool square, V sc, V scq) {
 if (nv <= tile_vec) {
  tile_work<P>(fa, fb, m, nv, B, square, sc, scq);
  return;
 }
 const int h= nv >> 2;
 if (square) stage4v<P>(fa, nv, m, nv, false);
 else stage4v<P>(fa, nv, m, nv, false), stage4v<P>(fb, nv, m, nv, false);
 for (int t= 0; t < 4; ++t) rec<P>(fa + 8 * t * h, fb + 8 * t * h, h, m, B + t * h, square, sc, scq);
 stage4v<P>(fa, nv, m, nv, true);
}

// 作業の配列。長さごとに残す (小さい畳み込みを何度も呼ぶ道具のため)。
template <unsigned P> struct conv_work {
 static inline thread_local abuf<64> fa, fb;
 static void reserve(int l) {
  if (fa.n < l) fa.assign(l), fb.assign(l);
 }
};

// 畳み込みの入口。
template <unsigned P> void convolve_entry(const ZMod<P>* a, int n, const ZMod<P>* b, int m_, int lo, int hi, ZMod<P>* out) {
 // 巡回の長さは L >= max(hi, n + m - 1 - lo) でよい。
 const int need= std::max(hi, n + m_ - 1 - lo);
 int l= 1;
 while (l < need) l<<= 1;
 if (std::min(n, m_) >= 60 && l >= 128 && (P - 1) % unsigned(l) == 0) {
  const int mv= l >> 3;
  const bool square= (a == b && n == m_);
  conv_work<P>::reserve(l);
  u32* fa= conv_work<P>::fa.p;
  u32* fb= conv_work<P>::fb.p;
  std::copy(reinterpret_cast<const u32*>(a), reinterpret_cast<const u32*>(a) + n, fa), std::fill(fa + n, fa + l, 0);
  if (!square) std::copy(reinterpret_cast<const u32*>(b), reinterpret_cast<const u32*>(b) + m_, fb), std::fill(fb + m_, fb + l, 0);
  // 葉の積が R^{-1} を 1 回付け、残りの逆の段が m 倍するので、R / m を掛ける (8 R / l)。
  const u32 c= u32(u64(mc<P>::r) * inv_n_scalar<P>(l) % P);
  const u32 c8= u32(u64(c) * 8 % P);
  rec<P>(fa, fb, mv, mv, 0, square, vsplat<P>(c8), vsplat<P>(shoup_q<P>(c8)));
  for (int k= lo; k < hi; ++k) out[k - lo]= ZMod<P>::raw(fa[k]);
  return;
 }
 Naive::convolve_range(a, n, b, m_, lo, hi, out);
}
}  // namespace neo_ntt_avx2

// ntt.hpp が使う口。
namespace ntt_internal {
static_assert(sizeof(ZMod<998244353>) == 4 && alignof(ZMod<998244353>) == 4, "AVX2 の実装は ZMod を u32 として読む");

template <unsigned P> void dft_impl(ZMod<P>* a, int n) {
 if (n <= 1) return;
 assert((n & (n - 1)) == 0 && (P - 1) % unsigned(n) == 0);
 u32* d= reinterpret_cast<u32*>(a);
 neo_ntt_avx2::to_mont<P>(d, n);
 neo_ntt_avx2::dft_stages<P>(d, n);
 neo_ntt_avx2::from_mont<P>(d, n);
}

// DIT。dft の逆。段ごとの 1/4 は掛けず、最後に 1/n を掛けながらふつうの値に戻す。
template <unsigned P> void idft_impl(ZMod<P>* a, int n) {
 if (n <= 1) return;
 assert((n & (n - 1)) == 0 && (P - 1) % unsigned(n) == 0);
 u32* d= reinterpret_cast<u32*>(a);
 neo_ntt_avx2::to_mont<P>(d, n);
 neo_ntt_avx2::idft_stages<P>(d, n);
 neo_ntt_avx2::scale_by<P>(d, n, neo_ntt_avx2::inv_n_scalar<P>(n));
}

// a[0..n) は長さ n の変換。係数に戻して 0 を詰め、長さ 2n の変換を a[0..2n) に作る (素直な形。変換の使い回しは、4 進の並びでは散らす手間が勝つので作っていない)。
// 途中は Montgomery のまま持つ: 逆変換の段のあとに Montgomery の 1/n を掛けて係数の表現にし、0 を詰めて順変換の段を回す。
template <unsigned P> void dft_doubling_impl(ZMod<P>* a, int n) {
 std::vector<ZMod<P>> c(a, a + n);
 u32* cd= reinterpret_cast<u32*>(c.data());
 neo_ntt_avx2::to_mont<P>(cd, n);
 neo_ntt_avx2::idft_stages<P>(cd, n);
 neo_ntt_avx2::scale_by<P>(cd, n, neo_ntt_avx2::smont<P>(neo_ntt_avx2::inv_n_scalar<P>(n)));
 std::fill(a, a + 2 * n, ZMod<P>());
 std::copy(c.begin(), c.end(), a);
 u32* d= reinterpret_cast<u32*>(a);
 neo_ntt_avx2::dft_stages<P>(d, 2 * n);
 neo_ntt_avx2::from_mont<P>(d, 2 * n);
}

// 畳み込み。Ntt::convolve_range が、最小の長さと法の条件を確かめてから呼ぶ。
template <unsigned P> void conv_kernel(const ZMod<P>* a, int n, const ZMod<P>* b, int m, int lo, int hi, ZMod<P>* out) {
 neo_ntt_avx2::convolve_entry<P>(a, n, b, m, lo, hi, out);
}
}  // namespace ntt_internal

#if defined(NEO_NTT_AVX2_PUSHED)
#pragma GCC pop_options
#undef NEO_NTT_AVX2_PUSHED
#endif
