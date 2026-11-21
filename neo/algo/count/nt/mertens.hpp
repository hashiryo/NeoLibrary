#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#include <simde/x86/fma.h>
#endif
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstring>
#include <vector>
// Mertens 関数 M(x) = Σ_{d≤x} μ(d) を、N の床の値 ⌊N/q⌋ (q ≤ R = ⌊√N⌋) すべてで求める (Θ(N^{2/3}))。
// ⌊N/q⌋ ≤ L (L ≈ 1.6 N^{2/3}) の q は、30 と互いに素な数だけを詰めた区画ごとの μ の篩で、30 と互いに素な数についての μ の和
// C を床の値で書き留め、M(x) = Σ_{a|30} μ(a) C(⌊x/a⌋) で戻す (平方因子の無い n は n = a m、a = gcd(n, 30) と 1 通りに書ける)。
// 素数 p ≥ 7 の倍数 p m は m mod 30 の 8 通りごとに、詰めた添字で 8p 飛びの等差数列になる。√L より大きい素因数は、足した
// log2 p の sc 倍の和が log2 x より小さいことで見分ける (sc は L < 2^25 なら 5、それ以上は 4 で、sc log2 L < 128。閾値は区画の
// 左端で決めるので、区画は左端の 2 倍まで)。最後の走査は
// AVX2 で 32 個ずつ。残りの q ≤ K = ⌊N/(L+1)⌋ は v = ⌊N/q⌋ について q の大きい順に、分け目 u = max(⌊R/q⌋, 0.4 √v) で
//   M(v) = 1 - Σ_{j=2}^{u} M(⌊v/j⌋) - Σ_{t≤t1} μ(t) ⌊v/t⌋ + u M(t1)  (t1 = ⌊v/(u+1)⌋ ≤ R + 2)
// とする。⌊v/j⌋ = ⌊N/(qj)⌋ なので、qj ≤ R なら求めた表を、そうでなければ R 以下の M の表を引く (AVX2 で 4 個ずつ添字を作り
// gather で引く)。後ろの和は AVX2 の double 4 本で、floor と FMA で double のまま足す (和は 2^53 未満)。⌊v/t⌋ は
// (v + 1/2) (1/t) の切り捨てで、v < 2^51 なら誤差が 1/(2t) 未満なので割り算と一致する。procon-judge の
// yosupo-sum-of-totient-function の 7 回目の tot7_l16c40g を写したもの。
namespace mertens_internal {
using u8= unsigned char;
using i8= signed char;
using u32= unsigned;
using i64= long long;
using u64= unsigned long long;
inline u64 isqrt(u64 n) {
 u64 r= (u64)std::sqrt((double)n);
 while(r * r > n) --r;
 while((r + 1) * (r + 1) <= n) ++r;
 return r;
}
// 30 と互いに素な 8 個の剰余と、その番号 (無ければ 255)、r 以下の個数
inline constexpr std::array<u8, 8> RES30= {1, 7, 11, 13, 17, 19, 23, 29};
inline constexpr auto IDX30= [] {
 std::array<u8, 30> a{};
 for(auto& x: a) x= 255;
 for(int i= 0; i < 8; ++i) a[RES30[i]]= (u8)i;
 return a;
}();
inline constexpr auto CNT30= [] {
 std::array<u8, 30> a{};
 int c= 0;
 for(int r= 0; r < 30; ++r) c+= IDX30[r] != 255, a[r]= (u8)c;
 return a;
}();
inline u64 cand_le(u64 x) { return 8 * (x / 30) + CNT30[x % 30]; }  // x 以下の 30 と互いに素な正の数の個数
struct Table {
 u64 N, R;
 std::vector<i8> mu;      // μ(x)、x ≤ R + 3 (+ α)
 std::vector<int> M;      // M(x)、x ≤ R + 3 (+ α)
 std::vector<i64> W;      // W[q] = M(⌊N/q⌋)、1 ≤ q ≤ R + 1
 std::vector<double> rt;  // 1/t
};
// N ≥ 1
inline Table floor_table(u64 N) {
 Table T;
 T.N= N;
 const u64 R= T.R= isqrt(N);
 u64 L= (u64)(1.6 * std::cbrt((double)N) * std::cbrt((double)N));
 L= std::max<u64>(L, R + 1);
 L= std::min<u64>((L / 30 + 1) * 30 - 1, (1ull << 32) - 31);
 const double sc= L < (1u << 25) ? 5 : 4;  // log の倍率 (sc log2 L < 128)
 const u64 K= N / (L + 1);                 // q ≤ K は ⌊N/q⌋ > L で再帰、q > K は篩で書き留める
 // 篩は 30 の倍数 nb (> R + 3) から始める。nb 未満の μ、M、30 と互いに素な数についての和 Cs は小さい篩で
 const u64 nb= (R + 4 + 29) / 30 * 30, RS= nb - 1;
 auto& mus= T.mu;
 auto& Ms= T.M;
 mus.assign(RS + 1, 1), Ms.assign(RS + 1, 0);
 std::vector<int> Cs(RS + 1, 0);
 {
  std::vector<u32> pr, spf(RS + 1, 0);
  mus[0]= 0;
  for(u64 i= 2; i <= RS; ++i) {
   if(!spf[i]) spf[i]= (u32)i, pr.push_back((u32)i), mus[i]= -1;
   for(u32 p: pr) {
    if((u64)p * i > RS) break;
    spf[p * i]= p;
    if(p == spf[i]) {  // p | i
     mus[p * i]= 0;
     break;
    }
    mus[p * i]= (i8)-mus[i];
   }
  }
  for(u64 i= 1; i <= RS; ++i) Ms[i]= Ms[i - 1] + mus[i], Cs[i]= Cs[i - 1] + (IDX30[i % 30] != 255 ? mus[i] : 0);
 }
 auto& rt= T.rt;
 rt.assign(RS + 1, 0);
 for(u64 t= 1; t <= RS; ++t) rt[t]= 1.0 / (double)t;
 const double Nh= (double)N + 0.5;
 auto fl= [&](u64 q) { return (u64)(Nh * rt[q]); };  // ⌊N/q⌋ (q ≤ RS)
 auto& W= T.W;
 W.assign(R + 2, 0);
 // 篩: 30 と互いに素な n ∈ [nb, L] の μ を区画ごとに。床の値 ⌊N/q⌋ (q は R から K+1 へ、値は小さい順) を通るところで C を
 // Cf[q] に書き留める
 std::vector<i64> Cf(R + 2, 0);
 u64 q= R;
 while(q > K && fl(q) < nb) --q;
 if(q > K && nb <= L) {
  const u32 sq= (u32)isqrt(L);
  struct Prog {
   u32 c, i;  // 数列の始めのブロックと、ブロックの中の番号
  };
  std::vector<u32> pr;
  std::vector<u8> lgp;
  std::vector<std::array<Prog, 8>> pg, pg2;
  {
   std::vector<u8> c(sq + 1, 0);
   for(u32 i= 2; i <= sq; ++i)
    if(!c[i]) {
     if(i >= 7) {
      pr.push_back(i), lgp.push_back((u8)(0x80 + std::lround(std::log2((double)i) * sc)));
      std::array<Prog, 8> a, b;
      const u64 i2= (u64)i * i;
      for(int j= 0; j < 8; ++j) a[j]= {(u32)((u64)i * RES30[j] / 30), IDX30[(u64)i * RES30[j] % 30]}, b[j]= {(u32)(i2 * RES30[j] / 30), IDX30[i2 * RES30[j] % 30]};
      pg.push_back(a), pg2.push_back(b);
     }
     for(u32 j= i * i; j <= sq; j+= i) c[j]= 1;
    }
  }
  constexpr u64 B= 1u << 13;  // 区画のブロック数 (8B 個の数)
  std::vector<u8> sv(8 * B + 32), zf(8 * B + 32);
  i64 C= Cs[nb - 1];
  u64 nx= fl(q);
  auto rec= [&]() {
   Cf[q]= C, --q;
   nx= q > K ? fl(q) : ~0ull;
  };
  for(u64 lb= nb / 30; lb * 30 <= L;) {
   const u64 nblk= std::min({B, (L + 1) / 30 - lb, lb});  // 区画は左端の 2 倍まで
   const u64 lo= lb * 30, hi= lo + 30 * nblk, len= 8 * nblk;
   std::memset(sv.data(), 0, len), std::memset(zf.data(), 0, len);
   for(size_t k= 0; k < pr.size(); ++k) {
    const u64 p= pr[k];
    if(p >= hi) break;
    const u8 l= lgp[k];
    for(int j= 0; j < 8; ++j) {
     const u64 c= pg[k][j].c, t0= lb > c ? (lb - c + p - 1) / p : 0;
     for(u64 ci= 8 * (p * t0 + c - lb) + pg[k][j].i; ci < len; ci+= 8 * p) sv[ci]= (u8)(sv[ci] + l);
    }
    const u64 p2= p * p;
    if(p2 < hi)
     for(int j= 0; j < 8; ++j) {
      const u64 c= pg2[k][j].c, t0= lb > c ? (lb - c + p2 - 1) / p2 : 0;
      for(u64 ci= 8 * (p2 * t0 + c - lb) + pg2[k][j].i; ci < len; ci+= 8 * p2) zf[ci]= 1;
     }
   }
   // 素数 p ≤ √L も上の輪で -1 になる (log が log p なので反転しない)
   const int th= (int)std::lround((std::log2((double)lo) - 0.5 * std::log2((double)L) * 0.5) * sc);
   const __m256i zero= _mm256_setzero_si256(), m7f= _mm256_set1_epi8(0x7F), one= _mm256_set1_epi8(1), thv= _mm256_set1_epi8((char)th);
   alignas(32) i8 mu8[32];
   auto target= [&]() -> i64 { return nx == ~0ull ? (i64)len : (i64)cand_le(nx) - (i64)(8 * lb) - 1; };
   i64 tg= target();
   while(q > K && tg < 0) rec(), tg= target();
   for(u64 i= 0; i < len; i+= 32) {
    const __m256i v= _mm256_loadu_si256((const __m256i*)&sv[i]), z= _mm256_loadu_si256((const __m256i*)&zf[i]);
    const __m256i sm= _mm256_cmpgt_epi8(zero, v), fm= _mm256_cmpgt_epi8(thv, _mm256_and_si256(v, m7f));
    const __m256i mu= _mm256_andnot_si256(_mm256_cmpeq_epi8(z, one), _mm256_or_si256(_mm256_xor_si256(sm, fm), one));
    const u64 e= std::min<u64>(32, len - i);
    if((i64)(i + e) <= tg && e == 32) {  // 書き留める番号を含まない 32 個は (μ + 1) の和から 32 を引く
     const __m256i sad= _mm256_sad_epu8(_mm256_add_epi8(mu, one), zero);
     const __m128i h= _mm_add_epi64(_mm256_castsi256_si128(sad), _mm256_extracti128_si256(sad, 1));
     C+= (i64)(_mm_cvtsi128_si64(h) + _mm_extract_epi64(h, 1)) - 32;
     continue;
    }
    _mm256_store_si256((__m256i*)mu8, mu);
    for(u64 d= 0; d < e; ++d) {
     C+= mu8[d];
     while(q > K && tg == (i64)(i + d)) rec(), tg= target();
    }
   }
   lb+= nblk;
  }
 }
 // q > K の W[q] = Σ_{a | 30} μ(a) C(⌊N/(aq)⌋)。aq ≤ R で ⌊N/(aq)⌋ ≥ nb なら Cf、そうでなければ Cs
 for(u64 qq= K + 1; qq <= R; ++qq) {
  const u64 x= fl(qq);
  if(x <= RS) {
   W[qq]= Ms[x];
   continue;
  }
  auto cv= [&](u64 a, u64 xa) -> i64 { return a * qq <= R && xa >= nb ? Cf[a * qq] : Cs[xa]; };
  W[qq]= cv(1, x) - cv(2, x / 2) - cv(3, x / 3) - cv(5, x / 5) + cv(6, x / 6) + cv(10, x / 10) + cv(15, x / 15) - cv(30, x / 30);
 }
 // q ≤ K は q の大きい順に再帰で
 {
  // μ(t) ≠ 0 の t の列 (μ(t)、1/t を double で) と、x 以下の個数 cnt[x]
  std::vector<double> sfm, sfr;
  std::vector<u32> cnt(RS + 1, 0);
  for(u64 t= 1; t <= RS; ++t) {
   if(mus[t]) sfm.push_back(mus[t]), sfr.push_back(rt[t]);
   cnt[t]= (u32)sfm.size();
  }
  for(int pad= 0; pad < 4; ++pad) sfm.push_back(0), sfr.push_back(0);
  for(u64 qi= K; qi >= 1; --qi) {
   const u64 v= fl(qi), s= isqrt(v);
   const u64 u= std::min(s, std::max<u64>({R / qi, (u64)(0.4 * (double)s), 1})), t1= v / (u + 1);
   const double vh= (double)v + 0.5;
   i64 r= 1 + (i64)u * Ms[t1];
   // Σ_{j=2}^{u} M(⌊v/j⌋)
   u64 j= 2;
   for(; j <= u && qi * j <= R; ++j) r-= W[qi * j];
   {
    const __m256d vv= _mm256_set1_pd(vh);
    __m128i acc= _mm_setzero_si128();
    for(; j + 3 <= u; j+= 4) acc= _mm_add_epi32(acc, _mm_i32gather_epi32(Ms.data(), _mm256_cvttpd_epi32(_mm256_mul_pd(vv, _mm256_loadu_pd(&rt[j]))), 4));
    alignas(16) int h[4];
    _mm_store_si128((__m128i*)h, acc);
    r-= (i64)h[0] + h[1] + h[2] + h[3];
   }
   for(; j <= u; ++j) r-= Ms[(u64)(vh * rt[j])];
   // Σ_{t≤t1} μ(t) ⌊v/t⌋ を 4 本ずつ (列の後ろは μ = 0 で埋めてある)
   {
    const size_t na= cnt[t1];
    const __m256d vv= _mm256_set1_pd(vh);
    __m256d acc= _mm256_setzero_pd();
    for(size_t a= 0; a < na; a+= 4) acc= _mm256_fmadd_pd(_mm256_loadu_pd(&sfm[a]), _mm256_floor_pd(_mm256_mul_pd(vv, _mm256_loadu_pd(&sfr[a]))), acc);
    // na を 4 の倍数に切り上げた分は t > t1 を含むので引き戻す
    for(size_t a= na; a < (na + 3) / 4 * 4; ++a) r+= (i64)(sfm[a] * std::floor(vh * sfr[a]));
    alignas(32) double h[4];
    _mm256_store_pd(h, acc);
    r-= (i64)h[0] + (i64)h[1] + (i64)h[2] + (i64)h[3];  // 本ごとの和は 2^53 未満
   }
   W[qi]= r;
  }
 }
 return T;
}
}
// Mertens 関数 M(N) = Σ_{d≤N} μ(d)。N < 2^50。Θ(N^{2/3}) 時間、Θ(N^{1/2}) 空間
inline long long mertens(unsigned long long N) {
 assert(N < (1ull << 50));
 if(N == 0) return 0;
 return mertens_internal::floor_table(N).W[1];
}
