#pragma once
#if defined(__AVX2__) || defined(__BMI2__)
#include <immintrin.h>
#endif
#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <vector>
// Deléglise–Rivat の方法で n 以下の素数の個数を数える。y = α n^{1/3}、z = n / y、a = π(y)、c = 8 (p_c = 19) とし、
//   π(n) = S1 + S2 + a - 1 - P2
//   S1 = Σ_{m ≤ y, 平方因子なし, lpf(m) > p_c} μ(m) φ(n / m, c)                          (ordinary leaves)
//   S2 = Σ_{c < b < a} Σ_{y / p_b < m ≤ y, 平方因子なし, lpf(m) > p_b} -μ(m) φ(n / (p_b m), b - 1)  (special leaves)
//   P2 = Σ_{y < p_b ≤ √n} (π(n / p_b) - b + 1)
// を足す。φ(t, b) は t 以下で最初の b 個の素数のどれでも割り切れない数の個数。special leaves は t = n / (p_b m) で分け、
// t ≥ p_b^2 (hard) は [1, z] を区間ごとに篩って p_1, ..., p_{b-1} を除いた残りを数え、p_b ≤ t < p_b^2 (easy) は
// φ = π(t) - b + 2 なので π の表を引き、t < p_b (trivial) は φ = 1 なので個数だけ足す。p_b > max(n^{1/3}, √y) の b は
// trivial しかなく、和を式で足す。m が素数 q の easy leaves のうち q > √(n / p) の部分は、Gourdon の反転
// Σ_{α<q≤β} π(N/q) = π(β)π(N/β) - π(α)π(N/α) + Σ_{N/β<r≤N/α} π(N/r) で √(n / p) 以下の r の和に直し、sparse の部分と
// 同じ表引きを重みつきで使い回す。合成数の m の easy leaves は素数 2 つの積なので、素数の 2 重の走査でたどる。
//
// 篩は 30 の車輪のビットで持つ (1 バイトに 30k + {1, 7, 11, 13, 17, 19, 23, 29})。7、11、13 の倍数は 1001 バイト周期の型を
// 写し、17 と 19 の倍数は 17 (19) バイト周期の型を 8 バイトずつ AND し、23 以上は倍数の車輪の 1 周 8 個を並べて消す。
// hard leaves の合成数の m は、b ごとに lpf(m) > p_b のものだけを t の昇順に並べておく。π の表は同じ篩の走査で作り、
// 奇数だけで持つ (64 個の数ごとに u64 で、下位 32 bit が奇素数の印、上位 32 bit が手前の奇素数の個数)。
// 割り算 N / d は、N d < 2^64 のとき floor(2^64 / d) + 1 との掛け算の上位で正確に求まる。この条件と、π の表の添字を
// 32 bit に収める条件から、n ≤ 10^14 に限る。
//
// AVX2 があれば、easy leaves の素数の q の表引きと ordinary leaves を 4 つずつ回す。商は 1 / d を 1 ulp 上げた
// double との掛け算の floor で求め (N < 2^50 で正確)、π の表は gather で引く。π(t) は表の値 e の上位 32 bit と、e を
// 63 - ((t - 1) / 2 mod 32) だけ左にずらした値の popcount (pshufb の 4 bit の表と psadbw) の和に 1 を足したもの。
namespace prime_pi_internal {
using u8= unsigned char;
using u32= unsigned;
using u64= unsigned long long;
using i64= long long;
using u128= unsigned __int128;
inline u64 isqrt(u64 n) {
 u64 r= (u64)std::sqrt((double)n);
 while(r * r > n) --r;
 while((r + 1) * (r + 1) <= n) ++r;
 return r;
}
inline u64 icbrt(u64 n) {
 u64 r= (u64)std::cbrt((double)n);
 while(r * r * r > n) --r;
 while((r + 1) * (r + 1) * (r + 1) <= n) ++r;
 return r;
}
inline u64 mulhi(u64 a, u64 b) { return (u64)(((u128)a * b) >> 64); }
// N d < 2^64 のとき mulhi(N, magic(d)) = N / d
inline u64 magic(u64 d) { return ~u64(0) / d + 1; }
// n < 2^53 なら double で割る。正しく丸めた商は floor(n / d) か 1 大きいだけなので、1 回直せばよい。
inline u64 fdiv(u64 n, u64 d) {
 if(n < (u64(1) << 53)) {
  const u64 q= (u64)(i64)((double)n / (double)d);
  return q - (q * d > n);
 }
 return n / d;
}
// 1 / d を丸めてから 1 ulp 上げた値。n < 2^50 なら floor(n inv_up(d)) = floor(n / d)。正の有限の double は、ビット列に 1 を
// 足すと次に大きい値になる (std::nextafter は関数の呼び出しになる)。
inline double inv_up(double d) { return std::bit_cast<double>(std::bit_cast<u64>(1.0 / d) + 1); }
inline constexpr u32 WR[8]= {1, 7, 11, 13, 17, 19, 23, 29};
inline constexpr std::array<u8, 30> BI= [] {
 std::array<u8, 30> t{};
 for(auto& v: t) v= 255;
 for(u32 i= 0; i < 8; ++i) t[WR[i]]= (u8)i;
 return t;
}();
// MASK[j]: 240 個の塊の中で j 以下の数のビット
inline constexpr std::array<u64, 240> MASK= [] {
 std::array<u64, 240> t{};
 for(u32 j= 0; j < 240; ++j)
  for(u32 k= 0; k < 8; ++k)
   for(u32 i= 0; i < 8; ++i)
    if(30 * k + WR[i] <= j) t[j]|= u64(1) << (8 * k + i);
 return t;
}();
// p ≡ WR[c] (mod 30) の倍数 p m (m ≡ WR[i]) のビットの位置と、m が 30 進む 1 周の先頭からのバイトの差の端数
inline constexpr std::array<std::array<u8, 8>, 8> WBIT= [] {
 std::array<std::array<u8, 8>, 8> t{};
 for(u32 c= 0; c < 8; ++c)
  for(u32 i= 0; i < 8; ++i) t[c][i]= BI[WR[c] * WR[i] % 30];
 return t;
}();
inline constexpr std::array<std::array<u8, 8>, 8> WOFF= [] {
 std::array<std::array<u8, 8>, 8> t{};
 for(u32 c= 0; c < 8; ++c)
  for(u32 i= 0; i < 8; ++i) t[c][i]= (u8)(WR[c] * WR[i] / 30);
 return t;
}();
// 奇数だけの π の表。t[k] の下位 32 bit の bit j は 64 k + 2 j + 1 が素数か、上位 32 bit は 64 k 未満の奇素数の個数。
// 引くときに 240 で割らずに済み、位置は shift だけで決まる。
struct PiOdd {
 std::unique_ptr<u64[]> t;  // 篩の走査ですべて書くので 0 で埋めない
 // 2 ≤ n < 2^32
 u64 operator()(u64 n) const {
  const u32 m= (u32)n - 1;
  const u64 e= t[m >> 6];
  const u32 pc= (u32)std::popcount((u32)e << (~m >> 1 & 31)) + 1;  // u32 で足すと符号の拡張が入らない
  return (e >> 32) + pc;
 }
};
// 30 の車輪の 1 バイトを、奇数の位置 (30 k + r の r = 1, 7, ..., 29 は (r - 1) / 2 番目) の 15 bit に広げる
inline constexpr std::array<uint16_t, 256> EXP15= [] {
 std::array<uint16_t, 256> t{};
 for(u32 v= 0; v < 256; ++v)
  for(u32 i= 0; i < 8; ++i)
   if(v >> i & 1) t[v]|= (uint16_t)(1u << (WR[i] - 1) / 2);
 return t;
}();
// 区間の車輪のバイト列 seg を奇数だけの表 out[0, ne) にする。running は手前の奇素数の個数で、書いた分だけ進める。
// first なら車輪の外の 3 と 5 を足す。BMI2 があれば 4 バイトを pdep で 60 bit に広げ、32 バイト (480 bit) から 15 語を作る。
inline void to_odd(const u8* seg, u32 ne, u64* out, u64& running, bool first) {
#ifdef __BMI2__
 constexpr u64 M60= 0x4B69ull | 0x4B69ull << 15 | 0x4B69ull << 30 | 0x4B69ull << 45;  // 0x4B69 は位置 {0,3,5,6,8,9,11,14}
 for(u32 g= 0, o= 0; o < ne; ++g) {
  u64 q[8]= {};
  for(u32 k= 0; k < 8; ++k) {
   u32 b;
   std::memcpy(&b, seg + 32 * g + 4 * k, 4);
   const u64 v= _pdep_u64(b, M60);
   const u32 off= 60 * k;
   q[off / 64]|= v << (off % 64);
   if(off % 64 > 4) q[off / 64 + 1]|= v >> (64 - off % 64);
  }
  if(first && g == 0) q[0]|= 6;
  for(u32 j= 0; j < 15 && o < ne; ++j, ++o) {
   const u32 v= (u32)(q[j / 2] >> (32 * (j % 2)));
   out[o]= running << 32 | v, running+= (u64)std::popcount(v);
  }
 }
#else
 // 奇数の位置の列は 1 バイトにつき 15 bit 進む。32 bit たまるごとに書く。
 u64 acc= first ? 6 : 0;
 u32 nacc= 0;
 for(u32 k= 0, o= 0; o < ne; ++k) {
  acc|= (u64)EXP15[seg[k]] << nacc, nacc+= 15;
  if(nacc >= 32) {
   const u32 v= (u32)acc;
   out[o++]= running << 32 | v, running+= (u64)std::popcount(v), acc>>= 32, nacc-= 32;
  }
 }
#endif
}
// φ(t, 6) を 30030 周期の表で引き、φ(t, 8) = φ(t, 6) - φ(t / 17, 6) - φ(t / 19, 6) + φ(t / 323, 6)
struct Phi8 {
 std::vector<uint16_t> tab;
 Phi8(): tab(30032) {  // 32 bit の gather で 2 バイト先まで読むので余分に持つ
  std::vector<u8> co(30030, 1);
  co[0]= 0;
  for(u32 p: {2u, 3u, 5u, 7u, 11u, 13u})
   for(u32 j= p; j < 30030; j+= p) co[j]= 0;
  u32 c= 0;
  for(u32 i= 0; i < 30030; ++i) tab[i]= (uint16_t)(c+= co[i]);
 }
 u64 phi6(u64 t) const { return t / 30030 * 5760 + tab[t % 30030]; }
 u64 operator()(u64 t) const { return phi6(t) - phi6(t / 17) - phi6(t / 19) + phi6(t / 323); }
};
// 篩う素数の状態。pos は今の区間の先頭からのバイトの位置、i は倍数 p m の m の車輪の位置。
struct SievingPrime {
 u32 pos, a;
 u8 cls, i;
};
// 区間の [0, nbytes) バイトにある p の倍数を消し、COUNT なら消えた数を返す。消えた数は局所変数に数える。
template <u32 C, bool COUNT> inline u32 cross_cls(u8* seg, u32 nbytes, SievingPrime& s) {
 const u32 a= s.a, p= 30 * a + WR[C];
 const u32 off[8]= {0, a * 6 + WOFF[C][1], a * 10 + WOFF[C][2], a * 12 + WOFF[C][3], a * 16 + WOFF[C][4], a * 18 + WOFF[C][5], a * 22 + WOFF[C][6], a * 28 + WOFF[C][7]};
 u32 gone= 0;
 auto clr= [&](u32 q, u32 bit) {
  if constexpr(COUNT) {
   const u8 v= seg[q];
   gone+= (v >> bit) & 1;
   seg[q]= (u8)(v & ~(1u << bit));
  } else {
   seg[q]&= (u8) ~(1u << bit);
  }
 };
 u32 i= s.i, base= s.pos - off[i];  // 周の先頭 (区間より前なら 2^32 で回る)
 for(; i < 8; ++i) {
  const u32 q= base + off[i];
  if(q >= nbytes) return s.pos= q - nbytes, s.i= (u8)i, gone;
  clr(q, WBIT[C][i]);
 }
 for(base+= p; base + off[7] < nbytes; base+= p) {
  clr(base, WBIT[C][0]), clr(base + off[1], WBIT[C][1]), clr(base + off[2], WBIT[C][2]), clr(base + off[3], WBIT[C][3]);
  clr(base + off[4], WBIT[C][4]), clr(base + off[5], WBIT[C][5]), clr(base + off[6], WBIT[C][6]), clr(base + off[7], WBIT[C][7]);
 }
 for(i= 0;; ++i) {
  const u32 q= base + off[i];
  if(q >= nbytes) return s.pos= q - nbytes, s.i= (u8)i, gone;
  clr(q, WBIT[C][i]);
 }
}
template <bool COUNT> inline u32 cross_off(u8* seg, u32 nbytes, SievingPrime& s) {
 switch(s.cls) {
  case 0: return cross_cls<0, COUNT>(seg, nbytes, s);
  case 1: return cross_cls<1, COUNT>(seg, nbytes, s);
  case 2: return cross_cls<2, COUNT>(seg, nbytes, s);
  case 3: return cross_cls<3, COUNT>(seg, nbytes, s);
  case 4: return cross_cls<4, COUNT>(seg, nbytes, s);
  case 5: return cross_cls<5, COUNT>(seg, nbytes, s);
  case 6: return cross_cls<6, COUNT>(seg, nbytes, s);
  default: return cross_cls<7, COUNT>(seg, nbytes, s);
 }
}
#if defined(__AVX2__)
// Σ_{l < i ≤ r} π(floor(n inv[i]))。4 つずつ AVX2 で回し、端は 1 つずつ。
inline i64 sum_pi4(const PiOdd& pt, const double* inv, u64 l, u64 r, double n) {
 u64 i= l + 1;
 i64 s= 0;
 if(i + 3 <= r) {
  const __m256d vn= _mm256_set1_pd(n);
  const __m256i one= _mm256_set1_epi64x(1), c31= _mm256_set1_epi64x(31), c63= _mm256_set1_epi64x(63), m4= _mm256_set1_epi8(0x0f);
  const __m256i lut= _mm256_setr_epi8(0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4);
  const long long* tab= (const long long*)pt.t.get();
  __m256i acc= _mm256_setzero_si256();
  for(; i + 3 <= r; i+= 4) {
   const __m256i t= _mm256_cvtepu32_epi64(_mm256_cvttpd_epi32(_mm256_mul_pd(vn, _mm256_loadu_pd(inv + i))));
   const __m256i m= _mm256_sub_epi64(t, one);
   const __m256i e= _mm256_i64gather_epi64(tab, _mm256_srli_epi64(m, 6), 8);
   // 下位 32 bit の bit j (j ≤ (m / 2) mod 32) だけが残り、上位 32 bit の個数は押し出される
   const __m256i v= _mm256_sllv_epi64(e, _mm256_sub_epi64(c63, _mm256_and_si256(_mm256_srli_epi64(m, 1), c31)));
   const __m256i pc= _mm256_add_epi8(_mm256_shuffle_epi8(lut, _mm256_and_si256(v, m4)), _mm256_shuffle_epi8(lut, _mm256_and_si256(_mm256_srli_epi16(v, 4), m4)));
   acc= _mm256_add_epi64(acc, _mm256_add_epi64(_mm256_sad_epu8(pc, _mm256_setzero_si256()), _mm256_srli_epi64(e, 32)));
  }
  alignas(32) u64 a4[4];
  _mm256_store_si256((__m256i*)a4, acc);
  s= (i64)(a4[0] + a4[1] + a4[2] + a4[3] + (i - l - 1));  // 1 つにつき 2 の分の 1
 }
 for(; i <= r; ++i) s+= (i64)pt((u64)(i64)(n * inv[i]));
 return s;
}
// ordinary leaves の Σ s_i φ(floor(x inv[i]), 8) (0 ≤ i < n)。s_i は sg が空なら 1、そうでなければ sg[i] (±1)。
// 商はどれも、整数の値を持つ double と 1 ulp 上げた逆数の掛け算の floor で求める (どの値も 2^50 未満なので正確)。
// φ(u, 6) = (u / 30030) 5760 + tab[u mod 30030] の表は 32 bit の gather で引いて下位 16 bit を取る。
inline i64 s1_sum4(const Phi8& ph, const double* inv, const int8_t* sg, size_t n, u64 x) {
 size_t i= 0;
 i64 s= 0;
 if(n >= 4) {
  const __m256d vx= _mm256_set1_pd((double)x), i17= _mm256_set1_pd(inv_up(17)), i19= _mm256_set1_pd(inv_up(19));
  const __m256d i323= _mm256_set1_pd(inv_up(323)), i30030= _mm256_set1_pd(inv_up(30030)), c30030= _mm256_set1_pd(30030);
  const __m256i c5760= _mm256_set1_epi64x(5760);
  const __m128i m16= _mm_set1_epi32(0xffff);
  const int* tab= (const int*)ph.tab.data();
  auto phi6= [&](__m256d u) {
   const __m256d q= _mm256_floor_pd(_mm256_mul_pd(u, i30030));
   const __m128i r= _mm256_cvttpd_epi32(_mm256_sub_pd(u, _mm256_mul_pd(q, c30030)));  // q 30030 < 2^53 なので正確
   const __m128i tv= _mm_and_si128(_mm_i32gather_epi32(tab, r, 2), m16);
   return _mm256_add_epi64(_mm256_mul_epi32(_mm256_cvtepi32_epi64(_mm256_cvttpd_epi32(q)), c5760), _mm256_cvtepu32_epi64(tv));
  };
  __m256i acc= _mm256_setzero_si256();
  for(; i + 4 <= n; i+= 4) {
   const __m256d t= _mm256_floor_pd(_mm256_mul_pd(vx, _mm256_loadu_pd(inv + i)));
   __m256i f= _mm256_sub_epi64(_mm256_add_epi64(phi6(t), phi6(_mm256_floor_pd(_mm256_mul_pd(t, i323)))), _mm256_add_epi64(phi6(_mm256_floor_pd(_mm256_mul_pd(t, i17))), phi6(_mm256_floor_pd(_mm256_mul_pd(t, i19)))));
   if(sg) {
    int32_t s4;
    std::memcpy(&s4, sg + i, 4);
    const __m256i ng= _mm256_cmpgt_epi64(_mm256_setzero_si256(), _mm256_cvtepi8_epi64(_mm_cvtsi32_si128(s4)));
    f= _mm256_sub_epi64(_mm256_xor_si256(f, ng), ng);
   }
   acc= _mm256_add_epi64(acc, f);
  }
  alignas(32) i64 a4[4];
  _mm256_store_si256((__m256i*)a4, acc);
  s= a4[0] + a4[1] + a4[2] + a4[3];
 }
 for(; i < n; ++i) {
  const i64 f= (i64)ph((u64)(i64)((double)x * inv[i]));
  s+= sg && sg[i] < 0 ? -f : f;
 }
 return s;
}
#endif
// 奇数だけの篩で数える。小さい x 用。
inline u64 pi_small(u64 n) {
 if(n < 2) return 0;
 u64 h= (n - 1) / 2, cnt= 1;
 std::vector<u8> comp(h + 1, 0);
 for(u64 i= 1; i <= h; ++i)
  if(!comp[i]) {
   ++cnt;
   for(u64 p= 2 * i + 1, j= (p * p) / 2; j <= h; j+= p) comp[j]= 1;
  }
 return cnt;
}
inline u64 prime_pi(u64 x) {
 if(x < 100000) return pi_small(x);
 constexpr u32 c= 8;
 const double L= std::log((double)x);
 const double alpha= std::max(1.0, 1.5 * (((0.00148918 * L - 0.0691909) * L + 1.00165) * L + 0.372253));
 const u64 x13= icbrt(x), sx= isqrt(x);
 const u64 y= std::min(sx, std::max(x13 + 1, (u64)(alpha * (double)x13))), z= x / y, sz= isqrt(z);
 // y + 2 までの素数 (1 始まり)。反転で引く r は √(x / p) + 2 以下で、反転があるのは √(x / p) < y のときだけ。
 std::vector<u32> primes;
 {
  const u64 L2= y + 2, h= (L2 - 1) / 2;
  std::vector<u8> comp(h + 1, 0);
  for(u64 i= 1; (2 * i + 1) * (2 * i + 1) <= L2; ++i)
   if(!comp[i])
    for(u64 p= 2 * i + 1, j= (p * p) / 2; j <= h; j+= p) comp[j]= 1;
  primes.resize(h + 2), primes[0]= 0, primes[1]= 2;
  size_t k= 2;
  for(u64 i= 1; i <= h; ++i) primes[k]= (u32)(2 * i + 1), k+= !comp[i];
  primes.resize(k);
 }
 auto npr= [&](u64 n) -> u64 { return (u64)(std::upper_bound(primes.begin() + 1, primes.end(), (u32)n) - primes.begin()) - 1; };
 const u64 a= npr(y), nsz= npr(sz);
 // y 以下の平方因子のない合成数で lpf > p_c のもの (昇順)。素数の積を深さ優先でたどって表に lpf と μ を書き、小さい順に詰める。
 std::vector<u32> cm, clp;
 std::vector<int8_t> cmu;
 {
  std::vector<uint16_t> tag(y + 1, 0);  // lpf | (μ > 0 ? 0x8000 : 0)。合成数の lpf は √y 以下
  size_t cnt= 0;
  auto dfs= [&](auto& self, u64 m, size_t i, u32 lp, int mu) -> void {
   for(; i <= a; ++i) {
    const u64 mm= m * primes[i];
    if(mm > y) break;
    if(lp) tag[mm]= (uint16_t)(lp | (mu > 0 ? 0x8000 : 0)), ++cnt;
    self(self, mm, i + 1, lp ? lp : primes[i], -mu);
   }
  };
  dfs(dfs, 1, c + 1, 0, 1);
  cm.resize(cnt + 1);
  size_t k= 0;
  for(u64 m= 1; m <= y; ++m) cm[k]= (u32)m, k+= tag[m] != 0;
  cm.resize(cnt), clp.resize(cnt), cmu.resize(cnt);
  for(size_t i= 0; i < cnt; ++i) clp[i]= tag[cm[i]] & 0x7fff, cmu[i]= (tag[cm[i]] & 0x8000) ? -1 : 1;
 }
 std::vector<u64> pd(primes.size());
 for(size_t i= 1; i < primes.size(); ++i) pd[i]= magic(primes[i]);
 std::vector<double> cinv(cm.size());  // 1 / m を 1 ulp 上げた値。N < 2^50 なら floor(N cinv) = floor(N / m)
 for(size_t i= 0; i < cm.size(); ++i) cinv[i]= inv_up(cm[i]);
#if defined(__AVX2__)
 std::vector<double> pinv(primes.size());  // 1 / p を 1 ulp 上げた値。N < 2^50 なら floor(N pinv) = floor(N / p)
 for(size_t i= 1; i < primes.size(); ++i) pinv[i]= inv_up(primes[i]);
#endif
 // hard leaves: b ごとに m ∈ (max(y/p, p), min(y, x/p^3)]。素数の m は primes の添字 (plo, pcur] で持ち、上から下る
 // (t が増える向き)。範囲が空でないなら p^2 ≤ z。
 u64 bmax= c;
 std::vector<u64> N(a + 1, 0);
 std::vector<u32> pcur(a + 1, 0), plo(a + 1, 0);
 for(u64 b= c + 1; b < a; ++b) {
  const u64 p= primes[b];
  if(p * p > z) break;
  const u64 lo= std::max(y / p, p), hi= std::min(y, x / (p * p * p));
  if(lo >= hi) continue;
  bmax= b, N[b]= x / p, plo[b]= (u32)npr(lo), pcur[b]= (u32)npr(hi);
 }
 // 合成数の m の hard leaves は、b の分を ht[hcb[b], hcb[b + 1]) に t の昇順で並べ、μ(m) > 0 なら (-μ φ を足すので
 // 引く) bit 31 を立てる。合成数の lpf は √y 以下なので葉のある b は p_b^2 < y に限られ、そこでは lo_b = y / p_b が b に
 // ついて減る。m を大きい順にたどると、m が葉になる b の範囲 (lo_b < m ≤ hi_b、p_b < lpf(m)) の両端は増えるだけなので、
 // 2 本の指で求まる。この範囲の差分から b ごとの個数を数えて場所を取り、書くのは lpf > p_b の合成数の列 alive から
 // lpf = p_b のものを抜きながら b を進め、(lo_b, hi_b] を後ろから写す。
 std::vector<u32> hcb(bmax + 2, 0), ht;
 {
  u64 bsq= c;  // p_b^2 < y となる最後の b
  while(bsq < bmax && (u64)primes[bsq + 1] * primes[bsq + 1] < y) ++bsq;
  std::vector<u32> pis(isqrt(y) + 1);  // v ≤ √y の π(v)
  for(size_t v= 0, i= 1; v < pis.size(); ++v) {
   while(i < primes.size() && primes[i] <= v) ++i;
   pis[v]= (u32)(i - 1);
  }
  std::vector<u64> lo(bsq + 1), hi(bsq + 1);
  for(u64 b= c + 1; b <= bsq; ++b) lo[b]= y / primes[b], hi[b]= std::min(y, x / ((u64)primes[b] * primes[b] * primes[b]));
  std::vector<i64> d(bsq + 2, 0);
  for(u64 i= cm.size(), bl= c + 1, bh= c; i-- > 0;) {
   const u64 m= cm[i];
   while(bh < bsq && hi[bh + 1] >= m) ++bh;
   while(bl <= bsq && lo[bl] >= m) ++bl;
   const u64 bt= std::min<u64>(bh, pis[clp[i]] - 1);
   if(bl <= bt) ++d[bl], --d[bt + 1];
  }
  for(u64 b= c + 1, run= 0; b <= bmax; ++b) hcb[b + 1]= hcb[b] + (u32)(b <= bsq ? run+= d[b] : 0);
  ht.resize(hcb[bmax + 1]);
  std::vector<u32> alive(cm.size());
  for(size_t i= 0; i < cm.size(); ++i) alive[i]= (u32)i;
  for(u64 b= c + 1; b <= bsq; ++b) {
   size_t k= 0;
   for(size_t j= 0; j < alive.size(); ++j) alive[k]= alive[j], k+= clp[alive[j]] > primes[b];
   alive.resize(k);
   auto gt= [&](u32 v, u32 i) { return v < cm[i]; };
   const auto i1= std::upper_bound(alive.begin(), alive.end(), (u32)lo[b], gt), i2= std::upper_bound(i1, alive.end(), (u32)hi[b], gt);
   u32* out= ht.data() + hcb[b];
   for(auto it= i2; it != i1;) --it, *out++= (u32)(i64)((double)N[b] * cinv[*it]) | (cmu[*it] > 0 ? 1u << 31 : 0);
  }
 }
 std::vector<u32> hcur(hcb.begin(), hcb.end());
 // 区間の篩。区間の先頭 low は 240 の倍数で、バイト k は 30 (low / 30 + k) + WR[i]。
 constexpr u32 SW= 2048;
 constexpr u64 SEGN= 240 * (u64)SW;
 std::vector<u64> seg64(SW + 1, 0);
 u8* seg= reinterpret_cast<u8*>(seg64.data());
 std::vector<u8> pat(1001);  // 7, 11, 13 の倍数を除く型
 for(u32 k= 0; k < 1001; ++k)
  for(u32 i= 0; i < 8; ++i) {
   const u64 n= 30 * (u64)k + WR[i];
   if(n % 7 && n % 11 && n % 13) pat[k]|= (u8)(1u << i);
  }
 // 17 と 19 の倍数を除く型。バイト g の型は g mod 17 (19) で決まるので、8 バイトずつ AND する。先頭から 8 バイト読める
 // ように 1 周より 7 バイト長く持つ。
 u8 p17[24], p19[26];
 for(u32 g= 0; g < 26; ++g) {
  u8 v17= 0, v19= 0;
  for(u32 i= 0; i < 8; ++i) {
   const u64 n= 30 * (u64)g + WR[i];
   v17|= (u8)((n % 17 != 0) << i), v19|= (u8)((n % 19 != 0) << i);
  }
  if(g < 24) p17[g]= v17;
  p19[g]= v19;
 }
 // 17 以上の素数で篩う。p_c までと hard leaves のある素数は φ のために p 自身から、残りは π の表のためだけなので p^2 から消す。
 const u64 nsp= std::max<u64>(nsz, c);
 std::vector<SievingPrime> sp(std::max(nsp, bmax) + 1);
 for(u64 b= c + 1; b <= nsp; ++b) {
  const u64 p= primes[b], m0= b <= std::max<u64>(bmax, c) ? 1 : p;
  sp[b]= SievingPrime{(u32)(p * m0 / 30), (u32)(p / 30), BI[p % 30], BI[m0 % 30]};
 }
 std::vector<i64> phi(bmax + 1, 0);
 std::vector<u32> pre(SW + 1);
 PiOdd pt;
 pt.t.reset(new u64[z / 64 + 1]);
 u64 running= 0;
 i64 s2h= 0;
 for(u64 low= 0; low <= z; low+= SEGN) {
  const u64 high= std::min(low + SEGN, z + 1);
  const u32 nbytes= (u32)((high - low + 29) / 30), nw= (nbytes + 7) / 8;
  for(u32 k= 0, r= (u32)(low / 30 % 1001); k < nbytes; r= 0) {
   const u32 len= std::min<u32>(1001 - r, nbytes - k);
   std::memcpy(seg + k, pat.data() + r, len), k+= len;
  }
  std::memset(seg + nbytes, 0, 8 * (SW + 1) - nbytes);
  seg64[(high - 1 - low) / 240]&= MASK[(high - 1 - low) % 240];  // high 以上の数のビットを落とす
  u64 cnt= 0;
  for(u32 k= 0, o17= (u32)(low / 30 % 17), o19= (u32)(low / 30 % 19); k < nw; ++k) {  // 17 と 19 は数えずに消す
   u64 a17, a19;
   std::memcpy(&a17, p17 + o17, 8), std::memcpy(&a19, p19 + o19, 8);
   cnt+= (u64)std::popcount(seg64[k]&= a17 & a19);
   o17+= 8, o17-= o17 >= 17 ? 17 : 0, o19+= 8, o19-= o19 >= 19 ? 19 : 0;
  }
  for(u64 b= c + 1; b <= bmax; ++b) {
   // この区間に入る葉は t < high のもの。素数の m は添字 (pj, pih]、合成数は ht の [hc0, hk)。
   const u32 pih= pcur[b], hc0= hcur[b], hc1= hcb[b + 1];
   if(pih > plo[b] || hc0 < hc1) {
    const u64 Nb= N[b];
    u32 pj= pih;
    if(pih > plo[b] && mulhi(Nb, pd[pih]) < high) {
     const u32 mb= (u32)(N[b] / high);  // t < high ⇔ m > mb
     pj= (u32)(std::upper_bound(primes.begin() + plo[b] + 1, primes.begin() + pih + 1, mb) - primes.begin()) - 1;
    }
    const u32 hk= (u32)(std::partition_point(ht.begin() + hc0, ht.begin() + hc1, [&](u32 v) { return (v & 0x7fffffff) < high; }) - ht.begin());
    if(pj < pih || hc0 < hk) {
     // 最も大きい t の語まで、語ごとの累積を作ってから葉を数える
     u64 tmax= 0;
     if(pj < pih) tmax= mulhi(Nb, pd[pj + 1]);
     if(hc0 < hk) tmax= std::max<u64>(tmax, ht[hk - 1] & 0x7fffffff);
     const u32 wl= (u32)((tmax - low) / 240);
     u32 run= 0;
     for(u32 w= 0; w <= wl; ++w) pre[w]= run, run+= (u32)std::popcount(seg64[w]);
     i64 sum= 0;
     for(u32 i= pih; i > pj; --i) {
      const u32 u= (u32)(mulhi(Nb, pd[i]) - low), w= u / 240;
      sum+= (i64)pre[w] + std::popcount(seg64[w] & MASK[u - 240 * w]);
     }
     s2h+= sum + (i64)(pih - pj) * phi[b], pcur[b]= pj;
     for(u32 k= hc0; k < hk; ++k) {
      const u32 u= (ht[k] & 0x7fffffff) - (u32)low, w= u / 240;
      const i64 ph= phi[b] + (i64)pre[w] + std::popcount(seg64[w] & MASK[u - 240 * w]);
      const i64 ng= -(i64)(ht[k] >> 31);  // 0 か -1。条件式で書くと gcc は分岐にする
      s2h+= (ph ^ ng) - ng;
     }
     hcur[b]= hk;
    }
   }
   phi[b]+= (i64)cnt;
   cnt-= cross_off<true>(seg, nbytes, sp[b]);  // p_b の倍数を消し、消えた数を cnt から引く
  }
  for(u64 b= bmax + 1; b <= nsz; ++b) cross_off<false>(seg, nbytes, sp[b]);
  if(low == 0) {
   seg[0]&= (u8)~1u;                                                                                      // 1 は素数でない
   for(u64 b= 4; b <= std::max<u64>(bmax, c); ++b) seg[primes[b] / 30]|= (u8)(1u << BI[primes[b] % 30]);  // 消した素数を戻す
  }
  to_odd(seg, (u32)((high - low + 63) / 64), pt.t.get() + low / 64, running, low == 0);
 }
 // easy leaves と trivial leaves。p_b > max(x^{1/3}, √y) の b は葉が q ∈ (p_b, y] の trivial だけで、b の分は a - b。
 i64 s2e= 0;
 u64 b0= c + 1;
 while(b0 < a && !(primes[b0] * primes[b0] * primes[b0] > x && y < primes[b0] * (primes[b0] + 1))) ++b0;
 if(b0 < a) s2e+= (i64)((a - b0) * (a - b0 + 1) / 2);
 for(u64 b= c + 1; b < b0; ++b) {
  const u64 p= primes[b];
  const u64 dv= pd[b];
  const u64 Nb= mulhi(x, dv), Nbp= mulhi(Nb, dv), lo= std::max({mulhi(y, dv), p, mulhi(Nbp, dv)});  // 素数の q も合成数の m も lo より大きい
  if(lo >= y) continue;
  // 合成数の m = q1 q2 (p < q1 < q2、lo < m ≤ y)。t ≥ p (p < √y ≤ √z) なので trivial はなく、μ(m) = 1 なので引く。
  for(u64 j= b + 1; primes[j] * primes[j + 1] <= y; ++j) {
   const u64 q1= primes[j], k1= npr(y / q1);
   const u64 Nq= mulhi(Nb, pd[j]);
   for(u64 k= std::max<u64>(j, npr(lo / q1)) + 1; k <= k1; ++k) s2e-= (i64)pt(mulhi(Nq, pd[k])) - (i64)b + 2;
  }
  // 素数の q ∈ (lo, y]: q ≤ he なら easy、q > he なら trivial
  const u64 he= std::min(y, Nbp);
  if(he < y) s2e+= (i64)(a - pt(std::max(lo, he)));
  if(lo >= he) continue;
  const u64 s= isqrt(Nb), ilo= pt(lo), ihe= pt(he), imid= std::min(ihe, pt(std::max(lo, std::min(he, s))));
  auto acc= [&](u64 l, u64 r) -> i64 {  // π(Nb / q_i) の (l, r] の和
#if defined(__AVX2__)
   return sum_pi4(pt, pinv.data(), l, r, (double)Nb);
#else
   i64 v= 0;
   for(u64 i= l + 1; i <= r; ++i) v+= (i64)pt(mulhi(Nb, pd[i]));
   return v;
#endif
  };
  // sparse: q ∈ (lo, min(he, s)]。clustered: q ∈ (al, he] を反転して r ∈ (Nb / he, Nb / al] の和にし、
  // r の添字 (rlo, rhi] のうち sparse の範囲に入る部分は、同じ表引きを 2 回足す。
  i64 sum= (i64)(imid - ilo) * (2 - (i64)b);
  const u64 al= std::max(lo, s);
  u64 rlo= imid, rhi= imid;
  if(al < he) {
   const u64 pal= pt(al), nb= Nb / he, na= Nb / al;
   sum+= (i64)(ihe * pt(nb)) - (i64)(pal * pt(na)) + (i64)(ihe - pal) * (2 - (i64)b), rlo= pt(nb), rhi= pt(na);
  }
  const u64 o1= std::clamp(rlo, ilo, imid), o2= std::clamp(rhi, o1, imid);
  sum+= acc(ilo, o1) + 2 * acc(o1, o2) + acc(o2, imid);
  if(rlo < ilo) sum+= acc(rlo, std::min(rhi, ilo));
  if(rhi > imid) sum+= acc(std::max(rlo, imid), rhi);
  s2e+= sum;
 }
 // ordinary leaves: m = 1、p_c < q ≤ y の素数、合成数
 Phi8 phic;
 i64 s1= (i64)phic(x);
#if defined(__AVX2__)
 s1+= s1_sum4(phic, cinv.data(), cmu.data(), cm.size(), x) - s1_sum4(phic, pinv.data() + c + 1, nullptr, a - c, x);
#else
 for(u64 i= c + 1; i <= a; ++i) s1-= (i64)phic(fdiv(x, primes[i]));
 for(size_t i= 0; i < cm.size(); ++i) {
  const i64 ph= (i64)phic((u64)(i64)((double)x * cinv[i]));
  s1+= cmu[i] > 0 ? ph : -ph;
 }
#endif
 // P2: y < p ≤ √x の素数を π の表のビットから列挙する
 i64 p2= 0;
 for(u64 w= (y + 1) / 64, b= a; w <= sx / 64; ++w)
  for(u32 bits= (u32)pt.t[w]; bits; bits&= bits - 1) {
   const u64 p= 64 * w + 2 * (u32)std::countr_zero(bits) + 1;
   if(p <= y) continue;
   if(p > sx) break;
   p2+= (i64)pt(fdiv(x, p)) - (i64)(++b) + 1;
  }
 return (u64)(s1 + s2h + s2e + (i64)a - 1 - p2);
}
}
// n 以下の素数の個数 π(n)。0 ≤ n ≤ 10^14。
inline unsigned long long prime_pi(unsigned long long n) {
 assert(n <= 100000000000000ull);
 return prime_pi_internal::prime_pi(n);
}
