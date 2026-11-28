#pragma once
// 任意の法 (素数でも合成数でもよい) の畳み込みの核。単独の NTT の核 Ntt をいくつかの素数で回し、Garner の方法で法 M の値に戻す。
// 素数の本数は値の上界から選ぶ。上界が素数の積で覆えない大きさは Naive に回す。
// 位の高い使い方として、素数ごとの畳み込みの前後に掛け算を挟む道 (each と garner) を持つ。Hurwitz 級数環の積
// (neo/base/core/poly/hurwitz.hpp) が、二項係数を素因数に分けて戻すのに使う。
// 口はほかの核と同じ convolve_range<ZMod<M>>(a, n, b, m, lo, hi, out)。K は ZMod (M は素数でも合成数でもよい)。
// 素数は 2^20 以上の 2 冪を p - 1 が持つ、2^30 未満のもの。2 の指数の大きい順に並べる (長さ len の畳み込みに使える
// 素数が先頭からの本数になる)。2^20 までの長さを約束する。
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <type_traits>
#include <vector>
#include "neo/base/core/algebra/ZMod.hpp"
#include "neo/base/core/convolve/convolve.hpp"
#include "neo/base/core/convolve/naive.hpp"
#include "neo/base/core/convolve/ntt.hpp"

namespace crt_internal {
using u32= unsigned;
using u64= unsigned long long;
using u128= unsigned __int128;

// 2^23, 2^22, 2^21, 2^20 の 2 冪を (p - 1) が持つ素数。同じ 2 の指数の中では大きい順。
constexpr int MAX_PRIMES= 12;
constexpr u32 PRIMES[MAX_PRIMES]= {
 998244353,  897581057,  880803841,   // 2^23
 985661441,  943718401,  935329793,   // 2^22
 1012924417, 1004535809, 975175681,   // 2^21
 1053818881, 1051721729, 1045430273,  // 2^20
};
constexpr int CAP_BITS[MAX_PRIMES]= {23, 23, 23, 22, 22, 22, 21, 21, 21, 20, 20, 20};

constexpr u64 pow_mod_c(u64 a, u64 e, u64 p) {
 u64 r= 1;
 while (e) {
  if (e & 1) r= r * a % p;
  a= a * a % p, e>>= 1;
 }
 return r;
}

// 6 本の u64 の桁で 384 bit までの非負整数を持つ。値の上界の比較にだけ使う。
// 桁あふれしたら飽和させ (すべての桁を ~0 にする)、どの素数の積よりも大きい扱いにする。
struct Big {
 u64 w[6];
 constexpr Big(u64 v= 0): w{v, 0, 0, 0, 0, 0} {}
 constexpr void mul_small(u64 m) {
  u64 carry= 0;
  for (int i= 0; i < 6; ++i) {
   const u128 t= u128(w[i]) * m + carry;
   w[i]= u64(t), carry= u64(t >> 64);
  }
  if (carry) w[0]= w[1]= w[2]= w[3]= w[4]= w[5]= ~0ull;
 }
 constexpr int cmp(const Big& o) const {
  for (int i= 5; i >= 0; --i)
   if (w[i] != o.w[i]) return w[i] < o.w[i] ? -1 : 1;
  return 0;
 }
};

// 先頭 t 本の素数の積の表。値は 2^360 未満なので 6 桁に収まる。
struct PrefixTable {
 Big pre[MAX_PRIMES + 1];
 constexpr PrefixTable(): pre{} {
  pre[0]= Big(1);
  for (int t= 0; t < MAX_PRIMES; ++t) pre[t + 1]= pre[t], pre[t + 1].mul_small(PRIMES[t]);
 }
};
inline constexpr PrefixTable PT= PrefixTable();

// Garner の定数。W_t = q_0 ... q_{t-1} として、WM[s][t] = W_t mod q_s (t < s)、IW[t] = W_t^{-1} mod q_t。
struct GarnerTable {
 u32 wm[MAX_PRIMES][MAX_PRIMES];
 u32 iw[MAX_PRIMES];
 constexpr GarnerTable(): wm{}, iw{} {
  for (int s= 0; s < MAX_PRIMES; ++s) {
   u64 w= 1;
   for (int t= 0; t < s; ++t) {
    wm[s][t]= u32(w % PRIMES[s]);
    w= w * PRIMES[t] % PRIMES[s];
   }
  }
  iw[0]= 1;
  for (int t= 1; t < MAX_PRIMES; ++t) {
   u64 w= 1;
   for (int s= 0; s < t; ++s) w= w * PRIMES[s] % PRIMES[t];
   iw[t]= u32(pow_mod_c(w, u64(PRIMES[t]) - 2, PRIMES[t]));
  }
 }
};
inline constexpr GarnerTable GT= GarnerTable();

// 長さ len の畳み込み (NTT の長さは len 以上の 2 冪) に使える、先頭からの素数の本数。
constexpr int usable(int len) {
 int l= 1;
 while (l < len) l<<= 1;
 int t= 0;
 while (t < MAX_PRIMES && (1 << CAP_BITS[t]) >= l) ++t;
 return t;
}
}  // namespace crt_internal

// Ntt は単独の NTT の核 (convolve_range<K> と、P を受ける形の conv_at を持つ)。
template <class Ntt> struct crt_engine {
 using u32= crt_internal::u32;
 static constexpr int max_primes= crt_internal::MAX_PRIMES;
 static constexpr u32 prime(int t) { return crt_internal::PRIMES[t]; }

 // 値の上界 bound_len^(1 + k) * M^2 (< 2^384) を覆うのに要る、長さ sz の畳み込みに使える素数の本数。
 // 覆えなければ 0。k = 0 は全体の畳み込み (上界 min(n, m) M^2)、k = |S| は Hurwitz 級数環の積。
 static int count_for(int bound_len, int sz, int k, u32 M) {
  using namespace crt_internal;
  Big b(1);
  for (int i= 0; i <= k; ++i) b.mul_small(u64(bound_len));
  b.mul_small(M), b.mul_small(M);
  const int use= usable(sz);
  for (int t= 1; t <= use; ++t)
   if (b.cmp(PT.pre[t]) < 0) return t;
  return 0;
 }

 // 素数 P の下での全体の積 (長さ sz)。a と b の値は P 未満。
 template <u32 P> static void conv_at(const u32* a, int n, const u32* b, int m, int sz, u32* out) {
  using Z= ZMod<P>;
  std::vector<Z> fa(n), fb(m), fc(sz);
  for (int i= 0; i < n; ++i) fa[i]= Z::raw(a[i]);
  for (int i= 0; i < m; ++i) fb[i]= Z::raw(b[i]);
  Ntt::template convolve_range<Z>(fa.data(), n, fb.data(), m, 0, sz, fc.data());
  for (int i= 0; i < sz; ++i) out[i]= fc[i].val();
 }
 // 素数を実行時の番号で選ぶ形。
 static void conv_at(int t, const u32* a, int n, const u32* b, int m, int sz, u32* out) {
  switch (t) {
   case 0: conv_at<prime(0)>(a, n, b, m, sz, out); break;
   case 1: conv_at<prime(1)>(a, n, b, m, sz, out); break;
   case 2: conv_at<prime(2)>(a, n, b, m, sz, out); break;
   case 3: conv_at<prime(3)>(a, n, b, m, sz, out); break;
   case 4: conv_at<prime(4)>(a, n, b, m, sz, out); break;
   case 5: conv_at<prime(5)>(a, n, b, m, sz, out); break;
   case 6: conv_at<prime(6)>(a, n, b, m, sz, out); break;
   case 7: conv_at<prime(7)>(a, n, b, m, sz, out); break;
   case 8: conv_at<prime(8)>(a, n, b, m, sz, out); break;
   case 9: conv_at<prime(9)>(a, n, b, m, sz, out); break;
   case 10: conv_at<prime(10)>(a, n, b, m, sz, out); break;
   case 11: conv_at<prime(11)>(a, n, b, m, sz, out); break;
   default: assert(false);
  }
 }
 // 先頭 cnt 本の素数それぞれで f<T, P>() を呼ぶ (T は番号、P は素数)。前後に掛け算を挟む道が使う。
 template <class F> static void each(int cnt, F&& f) {
  switch (cnt) {
   case 12: f.template operator()<11, prime(11)>(); [[fallthrough]];
   case 11: f.template operator()<10, prime(10)>(); [[fallthrough]];
   case 10: f.template operator()<9, prime(9)>(); [[fallthrough]];
   case 9: f.template operator()<8, prime(8)>(); [[fallthrough]];
   case 8: f.template operator()<7, prime(7)>(); [[fallthrough]];
   case 7: f.template operator()<6, prime(6)>(); [[fallthrough]];
   case 6: f.template operator()<5, prime(5)>(); [[fallthrough]];
   case 5: f.template operator()<4, prime(4)>(); [[fallthrough]];
   case 4: f.template operator()<3, prime(3)>(); [[fallthrough]];
   case 3: f.template operator()<2, prime(2)>(); [[fallthrough]];
   case 2: f.template operator()<1, prime(1)>(); [[fallthrough]];
   case 1: f.template operator()<0, prime(0)>(); [[fallthrough]];
   default: assert(cnt >= 0 && cnt <= max_primes);
  }
 }
 // r[t] = x mod prime(t) (t < cnt) から、x mod M を求める。x は先頭 cnt 本の素数の積未満の一意な値。
 static u32 garner(const u32* r, int cnt, u32 M) {
  using namespace crt_internal;
  if (M == 1) return 0;
  u32 c[MAX_PRIMES];
  c[0]= r[0];
  u64 acc= u64(r[0]) % M, w= 1;  // w = W_t (t = 0 のとき 1)
  for (int t= 1; t < cnt; ++t) {
   const u32 q= PRIMES[t];
   u64 xq= 0;
   for (int s= 0; s < t; ++s) xq= (xq + u64(c[s]) * GT.wm[t][s]) % q;
   c[t]= u32((u64(r[t]) + q - xq) % q * GT.iw[t] % q);
   w= w * PRIMES[t - 1] % M;  // W_t = W_{t-1} q_{t-1}
   acc= (acc + u64(c[t]) * w) % M;
  }
  return u32(acc);
 }

 // out[k - lo] = Σ_{i+j=k} a[i] b[j] (lo <= k < hi)。n >= 1、m >= 1、0 <= lo < hi <= n + m - 1。
 template <class K> static void convolve_range(const K* a, int n, const K* b, int m, int lo, int hi, K* out) {
  static_assert(is_zmod_v<K>, "Crt は ZMod の係数だけ");
  using namespace crt_internal;
  const u32 M= K::mod();
  if (M == 1) {
   std::fill(out, out + (hi - lo), K());
   return;
  }
  if (std::min(n, m) < 60) {
   Naive::convolve_range(a, n, b, m, lo, hi, out);
   return;
  }
  const int sz= n + m - 1;
  const int cnt= count_for(std::min(n, m), sz, 0, M);
  if (cnt == 0) {
   Naive::convolve_range(a, n, b, m, lo, hi, out);
   return;
  }
  std::vector<u32> xa(n), xb(m), rs(size_t(cnt) * sz);
  for (int t= 0; t < cnt; ++t) {
   const u32 q= PRIMES[t];
   for (int i= 0; i < n; ++i) {
    const u32 x= a[i].val();
    xa[i]= x >= q ? x - q : x;  // M < 2q (q >= 880803841、M < 2^30)
   }
   for (int i= 0; i < m; ++i) {
    const u32 x= b[i].val();
    xb[i]= x >= q ? x - q : x;
   }
   conv_at(t, xa.data(), n, xb.data(), m, sz, rs.data() + size_t(t) * sz);
  }
  u32 r[MAX_PRIMES];
  for (int k= lo; k < hi; ++k) {
   for (int t= 0; t < cnt; ++t) r[t]= rs[size_t(t) * sz + k];
   out[k - lo]= K::raw(garner(r, cnt, M));
  }
 }
};

using Crt= crt_engine<Ntt>;

// 既定の核。Ntt に使える素数は ntt_friendly の付いた宣言 (ntt.hpp) が優先される。
template <unsigned M> Crt default_kernel(const ZMod<M>*);
// 既定の多素数の核。Hurwitz 級数環の道具が、合成数の法の道で使うエンジンをこの宣言から探す。
template <unsigned M> Crt default_crt_engine(const ZMod<M>*);
