#pragma once
// 単独の NTT の核 (スカラー、基数 4)。procon-judge の self/ntt で 2 本目として書いたものを、NeoLibrary の既定の NTT として移した。
// dft は DIF で bit 反転の並べ替えをせず、index の 4 進の桁を逆にした順に評価値を置く (2 の指数が奇数のときは最後に基数 2 の段を 1 つ)。
// idft は DIT でその逆。1/4 は段ごとに掛けず、最後に 1/n をまとめて掛ける。
// ひねりは段ごとの表に持ち、内側の繰り返しで ω を掛け直さない。
// 口は convolve_range<K>(a, n, b, m, lo, hi, out)、dft<P>、idft<P>、dft_doubling<P>、even_dft<P>、odd_dft<P>。
// 小さい積と、NTT に使えない係数の型 (ntt_friendly でない素数、実行時の法 ZMod<0>、ZMod 以外) は Naive に回す。
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <type_traits>
#include <vector>
#include "neo/base/core/algebra/ZMod.hpp"
#include "neo/base/core/convolve/naive.hpp"
#include "neo/base/core/nt/is_prime.hpp"

namespace ntt_internal {
using u32= unsigned;
using u64= unsigned long long;

// 単独の NTT の核を既定にしてよい素数か。2^20 までの変換を約束する。
// (名前を外に出さないのは、procon-judge の実験用の核 (_shared/ntt/_kernel_base.hpp) と同じ名前になるため。)
template <unsigned P> constexpr bool ntt_friendly= P > 2 && is_prime(P) && (P - 1) % (1u << 20) == 0;

template <unsigned P> ZMod<P> pow_mod(ZMod<P> a, u64 e) {
 ZMod<P> r= 1;
 while (e) {
  if (e & 1) r*= a;
  a*= a, e>>= 1;
 }
 return r;
}

// mod は素数。mod - 1 の素因数を集め、どの p でも g^((mod-1)/p) != 1 となる最小の g を返す。
template <unsigned P> const ZMod<P>& primitive_root() {
 static const ZMod<P> g= [] {
  u32 n= P - 1;
  std::vector<u32> ps;
  for (u32 p= 2; u64(p) * p <= n; ++p) {
   if (n % p) continue;
   ps.push_back(p);
   while (n % p == 0) n/= p;
  }
  if (n > 1) ps.push_back(n);
  for (u32 g= 2;; ++g) {
   bool ok= true;
   for (u32 p: ps)
    if (pow_mod(ZMod<P>(g), (P - 1) / p) == ZMod<P>(1)) {
     ok= false;
     break;
    }
   if (ok) return ZMod<P>(g);
  }
 }();
 return g;
}

// 1 の原始 4 乗根 ω^{(P-1)/4} とその逆。基数 4 の段の i の掛け算に使う。
template <unsigned P> const ZMod<P>& root4() {
 static const ZMod<P> v= pow_mod(primitive_root<P>(), u64(P - 1) / 4);
 return v;
}
template <unsigned P> const ZMod<P>& iroot4() {
 static const ZMod<P> v= root4<P>().inv();
 return v;
}

// 段のひねりの表。長さ len の段では h = len/4 として、tw[j] = ω_len^j、tw[h + j] = ω_len^{2j}、tw[2h + j] = ω_len^{3j}。
// 同じ段のすべてのブロックで使い回す。大きさは 3 * (n/4) まで。
template <unsigned P> struct scratch {
 static inline std::vector<ZMod<P>> tw;
};
template <unsigned P> ZMod<P>* get_tw(int h) {
 auto& v= scratch<P>::tw;
 if ((int)v.size() < 3 * h) v.resize(3 * h);
 return v.data();
}

// tw に 3 本の等比数列を詰める。wl は ω_len (逆変換では ω_len^{-1})。
template <unsigned P> void fill_tw(ZMod<P>* tw, int h, ZMod<P> wl) {
 ZMod<P>* w1= tw;
 ZMod<P>* w2= tw + h;
 ZMod<P>* w3= tw + 2 * h;
 const ZMod<P> wl2= wl * wl, wl3= wl2 * wl;
 w1[0]= w2[0]= w3[0]= ZMod<P>(1);
 for (int j= 1; j < h; ++j) w1[j]= w1[j - 1] * wl, w2[j]= w2[j - 1] * wl2, w3[j]= w3[j - 1] * wl3;
}

// DIF。a[i] に Σ_j a[j] ω^{j k} が入る。k は位置 i の 4 進の桁を逆にした数 (2 の指数が奇数なら、最後の基数 2 の段で最上位の桁が最下位に来る)。
template <unsigned P> void dft_impl(ZMod<P>* a, int n) {
 if (n <= 1) return;
 assert((n & (n - 1)) == 0 && (P - 1) % unsigned(n) == 0);
 const ZMod<P> g= primitive_root<P>();
 const auto fi= root4<P>().fixed();
 int len= n;
 for (; len >= 4; len>>= 2) {
  const int h= len>> 2;
  ZMod<P>* tw= get_tw<P>(h);
  fill_tw<P>(tw, h, pow_mod(g, u64(P - 1) / unsigned(len)));
  for (int i= 0; i < n; i+= len) {
   ZMod<P>* x0= a + i;
   ZMod<P>* x1= x0 + h;
   ZMod<P>* x2= x1 + h;
   ZMod<P>* x3= x2 + h;
   for (int j= 0; j < h; ++j) {
    const ZMod<P> u0= x0[j], u1= x1[j], u2= x2[j], u3= x3[j];
    const ZMod<P> t0= u0 + u2, t1= u0 - u2, t2= u1 + u3, t3= u1 - u3;
    const ZMod<P> y= t3 * fi;
    x0[j]= t0 + t2;
    x1[j]= (t1 + y) * tw[j];
    x2[j]= (t0 - t2) * tw[h + j];
    x3[j]= (t1 - y) * tw[2 * h + j];
   }
  }
 }
 if (len == 2) {
  // 長さ 2 の段 (最後の 1 段)。ひねりは 1。
  for (int i= 0; i < n; i+= 2) {
   const ZMod<P> u= a[i], v= a[i + 1];
   a[i]= u + v, a[i + 1]= u - v;
  }
 }
}

// DIT。dft の逆。段ごとの 1/4 (最後の基数 2 の段の 1/2) は掛けず、最後に 1/n を掛ける。
template <unsigned P> void idft_impl(ZMod<P>* a, int n) {
 if (n <= 1) return;
 assert((n & (n - 1)) == 0 && (P - 1) % unsigned(n) == 0);
 const ZMod<P> g= primitive_root<P>();
 const auto fii= iroot4<P>().fixed();
 if (__builtin_ctz(unsigned(n)) & 1) {
  // dft の最後の基数 2 の段を先に戻す。
  for (int i= 0; i < n; i+= 2) {
   const ZMod<P> u= a[i], v= a[i + 1];
   a[i]= u + v, a[i + 1]= u - v;
  }
 }
 // dft の基数 4 の段の長さは、2 の指数が偶数なら 4, 16, ..., n、奇数なら 8, 32, ..., n。逆はその昇順。
 for (int len= (__builtin_ctz(unsigned(n)) & 1) ? 8 : 4; len <= n; len<<= 2) {
  const int h= len>> 2;
  ZMod<P>* tw= get_tw<P>(h);
  fill_tw<P>(tw, h, pow_mod(g, u64(P - 1) - u64(P - 1) / unsigned(len)));
  for (int i= 0; i < n; i+= len) {
   ZMod<P>* x0= a + i;
   ZMod<P>* x1= x0 + h;
   ZMod<P>* x2= x1 + h;
   ZMod<P>* x3= x2 + h;
   for (int j= 0; j < h; ++j) {
    const ZMod<P> u0= x0[j];
    const ZMod<P> y1= x1[j] * tw[j];
    const ZMod<P> y2= x2[j] * tw[h + j];
    const ZMod<P> y3= x3[j] * tw[2 * h + j];
    const ZMod<P> t0= u0 + y2, t2= u0 - y2, t1= y1 + y3, t3= (y1 - y3) * fii;
    x0[j]= t0 + t1;
    x1[j]= t2 + t3;
    x2[j]= t0 - t1;
    x3[j]= t2 - t3;
   }
  }
 }
 const ZMod<P> iv= ZMod<P>(unsigned(n)).inv();
 for (int i= 0; i < n; ++i) a[i]*= iv;
}

// a[0..n) は長さ n の変換。係数に戻して 0 を詰め、長さ 2n の変換を a[0..2n) に作る (素直な形。変換の使い回しはまだ作っていない)。
template <unsigned P> void dft_doubling_impl(ZMod<P>* a, int n) {
 std::vector<ZMod<P>> c(a, a + n);
 idft_impl<P>(c.data(), n);
 std::fill(a, a + 2 * n, ZMod<P>());
 std::copy(c.begin(), c.end(), a);
 dft_impl<P>(a, 2 * n);
}

// 長さ n の変換から、偶数次の部分の長さ n/2 の変換を a[0..n/2) に作る。
// 並びが 4 進の桁の逆なので、対になる評価点は隣どうしではない (長さ n/2 が 4 の冪か、2 * 4 の冪かで変わる)。
template <unsigned P> void even_dft_impl(ZMod<P>* a, int n) {
 const int m= n>> 1;
 const ZMod<P> iv2= ZMod<P>(2).inv();
 if (__builtin_ctz(unsigned(n)) & 1) {
  // 長さ n/2 が 4 の冪。位置 2j と 2j + 1 が対 (上位の桁が 0 と 2)。
  for (int j= 0; j < m; ++j) a[j]= (a[2 * j] + a[2 * j + 1]) * iv2;
 } else {
  // 長さ n/2 が 2 * 4 の冪。位置 i と i + 2 が対 (i = 4 (j>>1) + (j&1) で、上位の桁が (j&1) と (j&1) + 2)。
  for (int j= 0; j < m; ++j) {
   const int i= ((j>> 1)<< 2) | (j & 1);
   a[j]= (a[i] + a[i + 2]) * iv2;
  }
 }
}

// 長さ n の変換から、奇数次の部分の長さ n/2 の変換を a[0..n/2) に作る。
// 評価点の番号 q は、出力の位置 j の 4 進の桁を逆にした数 (と、その最上位の桁)。ω_n^{-q} を段の差分で掛け進める。
template <unsigned P> void odd_dft_impl(ZMod<P>* a, int n) {
 const int e= __builtin_ctz(unsigned(n));
 const int m= n>> 1;
 const ZMod<P> g= primitive_root<P>();
 const ZMod<P> wi= pow_mod(g, u64(P - 1) - u64(P - 1) / unsigned(n));
 // 出力の位置を 1 進めるごとに掛ける因子。位取りの D 桁で、下から t 桁が 3 のとき δ_t = 5 * 4^{D-1-t} - 4^D だけ q が変わる。
 const int D= (e & 1) ? (e>> 1) : (e>> 1) - 1;
 ZMod<P> step[16];
 const u64 p4= D ? 1ull << (2 * D) : 1;
 const ZMod<P> w4p= pow_mod(g, u64(P - 1) / unsigned(n) * p4);  // ω_n^{+4^D}
 const ZMod<P> w4i= pow_mod(wi, p4);  // ω_n^{-4^D}
 for (int t= 0; t < D; ++t) step[t]= pow_mod(wi, 5ull << (2 * (D - 1 - t))) * w4p;
 if (e & 1) {
  // 長さ n/2 が 4 の冪。位置 2j と 2j + 1 が対で、q は j の 4 進の桁を逆にした数。
  ZMod<P> base= ZMod<P>(2).inv();
  for (int j= 0; j < m; ++j) {
   a[j]= (a[2 * j] - a[2 * j + 1]) * base;
   if (j + 1 < m) base*= step[__builtin_ctz(~unsigned(j))>> 1];
  }
 } else {
  // 長さ n/2 が 2 * 4 の冪。位置 4j と 4j + 2、4j + 1 と 4j + 3 が対で、q は 4 進の桁の逆に最上位の桁 (0 か 1) を足した数。
  const ZMod<P> iv2= ZMod<P>(2).inv();
  ZMod<P> base= iv2;
  const int mm= m>> 1;
  for (int j= 0; j < mm; ++j) {
   a[2 * j]= (a[4 * j] - a[4 * j + 2]) * base;
   a[2 * j + 1]= (a[4 * j + 1] - a[4 * j + 3]) * (base * w4i);
   if (j + 1 < mm) base*= step[__builtin_ctz(~unsigned(j))>> 1];
  }
 }
}

template <class K> struct static_ntt_mod: std::false_type {
 static constexpr unsigned mod= 0;
};
template <unsigned P> struct static_ntt_mod<ZMod<P>>: std::bool_constant<ntt_friendly<P>> {
 static constexpr unsigned mod= P;
};
}  // namespace ntt_internal

struct Ntt {
 // out[k - lo] = Σ_{i+j=k} a[i] b[j] (lo <= k < hi)。n >= 1、m >= 1、0 <= lo < hi <= n + m - 1。
 template <class K> static void convolve_range(const K* a, int n, const K* b, int m, int lo, int hi, K* out) {
  if constexpr (ntt_internal::static_ntt_mod<K>::value) {
   constexpr unsigned P= ntt_internal::static_ntt_mod<K>::mod;
   // 巡回の長さは L >= max(hi, n + m - 1 - lo) でよい (折り返した項が [lo, hi) に入らない)。
   const int need= std::max(hi, n + m - 1 - lo);
   int l= 1;
   while (l < need) l<<= 1;
   if (std::min(n, m) >= 60 && (P - 1) % unsigned(l) == 0) {
    std::vector<K> fa(a, a + n), fb;
    fa.resize(l);
    if (a == b && n == m) {
     ntt_internal::dft_impl<P>(fa.data(), l);
     for (int i= 0; i < l; ++i) fa[i]*= fa[i];
    } else {
     fb.assign(b, b + m), fb.resize(l);
     ntt_internal::dft_impl<P>(fa.data(), l), ntt_internal::dft_impl<P>(fb.data(), l);
     for (int i= 0; i < l; ++i) fa[i]*= fb[i];
    }
    ntt_internal::idft_impl<P>(fa.data(), l);
    for (int k= lo; k < hi; ++k) out[k - lo]= fa[k];
    return;
   }
  }
  Naive::convolve_range(a, n, b, m, lo, hi, out);
 }

 template <unsigned P> requires ntt_internal::ntt_friendly<P> static void dft(ZMod<P>* a, int n) { ntt_internal::dft_impl<P>(a, n); }
 template <unsigned P> requires ntt_internal::ntt_friendly<P> static void idft(ZMod<P>* a, int n) { ntt_internal::idft_impl<P>(a, n); }
 template <unsigned P> requires ntt_internal::ntt_friendly<P> static void dft_doubling(ZMod<P>* a, int n) { ntt_internal::dft_doubling_impl<P>(a, n); }
 template <unsigned P> requires ntt_internal::ntt_friendly<P> static void even_dft(ZMod<P>* a, int n) { ntt_internal::even_dft_impl<P>(a, n); }
 template <unsigned P> requires ntt_internal::ntt_friendly<P> static void odd_dft(ZMod<P>* a, int n) { ntt_internal::odd_dft_impl<P>(a, n); }
};

// 既定の核。ntt_friendly な素数だけを受け持つ (この宣言が無いと、核を include し忘れたときに黙って Naive で動く)。
template <unsigned P> requires ntt_internal::ntt_friendly<P> Ntt default_kernel(const ZMod<P>*);
