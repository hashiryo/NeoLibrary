#pragma once
#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cmath>
#include <cstring>
#include <memory>
#include <utility>
#include <vector>
#if defined(__x86_64__) && defined(__BMI2__)
#include <immintrin.h>
#endif
// エラトステネスの篩を 30 の車輪のビット (1 バイトに 30k + {1, 7, 11, 13, 17, 19, 23, 29}) で 32 KB の区画ごとに回す。
//
// 7 から 163 までの素数は、2 つか 3 つずつまとめた型 (周期は積のバイト数で 64 KB 以下) を 16 本作っておき、8 本ずつの AND を
// 32 バイトずつ書く (GCC のベクトル拡張で、x64 は AVX2、arm は NEON になる)。
//
// √N 以下の残りの素数 p は、倍数 p m の m を車輪で進めて消す。素数ごとに次に消す周の先頭だけを持ち、区画の中に先頭がある周を
// 丸ごと消す。最後の周は区画の後ろへはみ出すので、区画の後ろに余白を持って全部 1 にしておき、はみ出しを次の区画の先頭に
// AND する (余白が次の区画より長くても、余白の後ろ半分は次の区画の余白に重なるので、そのまま先へ送られる)。周の途中から
// 始めないので、素数ごとの出入りで輪の出口の分岐が 1 回しか外れない。p ≤ 4000 は m を 210 の車輪で進めて 7 の倍数を
// 飛ばす。m が 210 進む 7 周 (48 個) をまとめて消し、周の先頭を m ≡ 1 (mod 210) にそろえると、飛ばす位置は p によらず、
// 周 r の中で (2r + w_k) mod 7 = 0 となる k になる。p > 4000 は 30 の車輪の 1 周 8 個ずつ消す。素数は p mod 30 の類ごとに
// 並べ、類ごとに展開した関数を呼ぶ。
//
// A 個おきの取り出しは区画ごとに 2 段で行う。1 段目は語の popcount の累積 c から、番号が c 未満の取り出しの数
// floor((c + A - 1 - t0) / A) を ceil(2^64 / A) との掛け算の上位で求め、取り出す素数がある語の番号を列に書く (分岐なし)。
// 2 段目は列の各項で語の中の k 番目のビットを引く (BMI2 があれば pdep、無ければバイトごとの累積と 256 通りの表)。
namespace enumerate_primes_internal {
using u8= unsigned char;
using u32= unsigned;
using u64= unsigned long long;
using u128= unsigned __int128;
inline constexpr u32 WR[8]= {1, 7, 11, 13, 17, 19, 23, 29};
inline constexpr std::array<u8, 30> BI= [] {
 std::array<u8, 30> t{};
 for(auto& v: t) v= 255;
 for(u32 i= 0; i < 8; ++i) t[WR[i]]= (u8)i;
 return t;
}();
// p ≡ WR[c] (mod 30) の倍数 p m (m ≡ WR[i]) のビットの位置と、周の先頭からのバイトの差の端数
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
// 語のビット b が表す数の、語の先頭からの差。64 は捨てる値。
inline constexpr std::array<u8, 65> BITVAL= [] {
 std::array<u8, 65> t{};
 for(u32 b= 0; b < 64; ++b) t[b]= (u8)(30 * (b >> 3) + WR[b & 7]);
 return t;
}();
#if !(defined(__x86_64__) && defined(__BMI2__))
// バイト b の中で r 番目に立っているビットの位置
inline constexpr std::array<std::array<u8, 8>, 256> SEL8= [] {
 std::array<std::array<u8, 8>, 256> t{};
 for(u32 b= 0; b < 256; ++b)
  for(u32 i= 0, r= 0; i < 8; ++i)
   if(b >> i & 1) t[b][r++]= (u8)i;
 return t;
}();
#endif
// 語 x の中で k 番目 (0 始まり、k < 64) に立っているビットの位置。立っているビットが k 個以下なら 0 から 64 の何か。
inline u32 kth_bit(u64 x, u32 k) {
#if defined(__x86_64__) && defined(__BMI2__)
 return (u32)std::countr_zero(_pdep_u64(u64(1) << k, x));
#else
 constexpr u64 L8= 0x0101010101010101ull, H8= 0x8080808080808080ull;
 u64 s= x - ((x >> 1) & 0x5555555555555555ull);
 s= (s & 0x3333333333333333ull) + ((s >> 2) & 0x3333333333333333ull);
 s= (s + (s >> 4)) & 0x0f0f0f0f0f0f0f0full;
 const u64 pre= s * L8;                                 // バイト i は、バイト 0 から i までの popcount
 const u64 le= ((k * L8 | H8) - pre) & H8;              // 累積が k 以下のバイトに最上位ビットを立てる
 const u32 sh= (u32)(((le >> 7) * L8 >> 56) * 8) & 63;  // k 番目のビットがあるバイトの位置
 const u32 before= (u32)((pre << 8) >> sh) & 0xff;      // それより前のバイトの popcount
 return sh + SEL8[(x >> sh) & 0xff][(k - before) & 7];
#endif
}
// 篩う素数の状態。pos は次に消す周 (210 の車輪なら 7 周のまとまり) の先頭の、今の区画の先頭からのバイトの位置。a = p / 30。
struct SievingPrime {
 u32 pos, a;
};
// 区画の [0, nbytes) に先頭がある周を 8 個ずつ丸ごと消す (p ≡ WR[C] mod 30)。最後の周は p バイト未満はみ出す。
template <u32 C> inline void cross30(u8* seg, u32 nbytes, SievingPrime& s) {
 const u32 a= s.a, p= 30 * a + WR[C];
 const u32 o1= a * 6 + WOFF[C][1], o2= a * 10 + WOFF[C][2], o3= a * 12 + WOFF[C][3], o4= a * 16 + WOFF[C][4];
 const u32 o5= a * 18 + WOFF[C][5], o6= a * 22 + WOFF[C][6], o7= a * 28 + WOFF[C][7];
 constexpr u8 m0= (u8) ~(1u << WBIT[C][0]), m1= (u8) ~(1u << WBIT[C][1]), m2= (u8) ~(1u << WBIT[C][2]), m3= (u8) ~(1u << WBIT[C][3]);
 constexpr u8 m4= (u8) ~(1u << WBIT[C][4]), m5= (u8) ~(1u << WBIT[C][5]), m6= (u8) ~(1u << WBIT[C][6]), m7= (u8) ~(1u << WBIT[C][7]);
 u32 base= s.pos;
 for(; base < nbytes; base+= p) {
  u8* q= seg + base;
  q[0]&= m0, q[o1]&= m1, q[o2]&= m2, q[o3]&= m3, q[o4]&= m4, q[o5]&= m5, q[o6]&= m6, q[o7]&= m7;
 }
 s.pos= base - nbytes;
}
// 7 周のまとまりの周 R で、周の先頭 q から差 o の位置を消す。m が 7 の倍数になる k ((2R + WR[k]) mod 7 = 0) は飛ばす。
template <u32 C, u32 R, u32... K> inline void row210(u8* q, const u32* o, std::integer_sequence<u32, K...>) { ((((2 * R + WR[K]) % 7 != 0) ? (void)(q[o[K]]&= (u8) ~(1u << WBIT[C][K])) : (void)0), ...); }
// 区画の [0, nbytes) に先頭がある 7 周のまとまり (m ≡ 1 mod 210 から) を丸ごと消す。最後のまとまりは 7p バイト未満はみ出す。
template <u32 C> inline void cross210(u8* seg, u32 nbytes, SievingPrime& s) {
 const u32 a= s.a, p= 30 * a + WR[C];
 const u32 o[8]= {0, a * 6 + WOFF[C][1], a * 10 + WOFF[C][2], a * 12 + WOFF[C][3], a * 16 + WOFF[C][4], a * 18 + WOFF[C][5], a * 22 + WOFF[C][6], a * 28 + WOFF[C][7]};
 constexpr auto K8= std::make_integer_sequence<u32, 8>{};
 u32 base= s.pos;
 for(; base < nbytes; base+= 7 * p) {
  u8* q= seg + base;
  row210<C, 0>(q, o, K8), q+= p;
  row210<C, 1>(q, o, K8), q+= p;
  row210<C, 2>(q, o, K8), q+= p;
  row210<C, 3>(q, o, K8), q+= p;
  row210<C, 4>(q, o, K8), q+= p;
  row210<C, 5>(q, o, K8), q+= p;
  row210<C, 6>(q, o, K8);
 }
 s.pos= base - nbytes;
}
// st[beg[c], beg[c + 1]) が類 c の素数。W210 なら 210 の車輪で消す。
template <bool W210, u32 C= 0> inline void cross_all(u8* seg, u32 nbytes, SievingPrime* st, const u32* beg) {
 for(u32 k= beg[C]; k < beg[C + 1]; ++k) {
  if constexpr(W210) cross210<C>(seg, nbytes, st[k]);
  else cross30<C>(seg, nbytes, st[k]);
 }
 if constexpr(C < 7) cross_all<W210, C + 1>(seg, nbytes, st, beg);
}
// 7 から 163 までの素数の倍数を除く型。バイト g の型は g mod (まとめた素数の積) で決まる。1 周の先に CH + 32 バイト足して
// 持ち、途中から CH バイトを続けて読めるようにする。
struct Presieve {
 static constexpr u32 PMAX= 163, CH= 2048, G= 8;  // CH は 1 回に書くバイト数、G は 1 回の読み書きで AND する型の本数
 using V= u8 __attribute__((vector_size(32)));
 struct Pat {
  std::vector<u8> b;
  u32 period, step;  // step = CH mod period
 };
 std::vector<Pat> pats;  // G の倍数本 (足りない分は全部 1 の型)
 Presieve() {
  std::vector<u32> ps;
  for(u32 p= 7; p <= PMAX; p+= 2) {
   bool pr= true;
   for(u32 d= 3; d * d <= p; d+= 2) pr&= p % d != 0;
   if(pr) ps.push_back(p);
  }
  for(std::size_t i= 0; i < ps.size();) {
   std::size_t j= i + 1;
   u64 prod= ps[i];
   while(j < ps.size() && j < i + 3 && prod * ps[j] <= 65536) prod*= ps[j++];
   add(std::vector<u32>(ps.begin() + (std::ptrdiff_t)i, ps.begin() + (std::ptrdiff_t)j));
   i= j;
  }
  while(pats.size() % G) add({});
 }
 void add(const std::vector<u32>& g) {
  Pat P;
  P.period= 1;
  for(u32 p: g) P.period*= p;
  const u32 len= P.period + CH + 32;
  P.b.assign(len, 0xff);
  for(u32 p: g) {
   std::vector<u8> one(p, 0xff);  // p の倍数を除く 1 周 p バイトの型
   for(u32 i= 0; i < 8; ++i) {
    u32 h= 0;
    while((30 * h + WR[i]) % p) ++h;
    one[h]&= (u8) ~(1u << i);
   }
   for(u32 o= 0; o < len; o+= p) {
    const u32 n= std::min(p, len - o);
    for(u32 k= 0; k < n; ++k) P.b[o + k]&= one[k];
   }
  }
  P.step= CH % P.period;
  pats.push_back(std::move(P));
 }
 // seg[0, len) に、先頭のバイトが g0 の型を書く。len を 32 の倍数に切り上げた所まで書く。
 void fill(u8* seg, u32 len, u64 g0) const {
  const std::size_t np= pats.size();
  u32 off[32];
  for(std::size_t k= 0; k < np; ++k) off[k]= (u32)(g0 % pats[k].period);
  auto ld= [](const u8* p) {
   V v;
   std::memcpy(&v, p, 32);
   return v;
  };
  for(u32 c= 0; c < len; c+= CH) {
   const u32 nv= (std::min(CH, len - c) + 31) / 32;
   u8* out= seg + c;
   for(std::size_t k= 0; k < np; k+= G) {
    const u8* p[G];
    for(u32 t= 0; t < G; ++t) p[t]= pats[k + t].b.data() + off[k + t];
    if(k == 0)
     for(u32 j= 0; j < nv; ++j) {
      V v= ld(p[0] + 32 * j);
#pragma GCC unroll 16
      for(u32 t= 1; t < G; ++t) v&= ld(p[t] + 32 * j);
      std::memcpy(out + 32 * j, &v, 32);
     }
    else
     for(u32 j= 0; j < nv; ++j) {
      V v= ld(out + 32 * j);
#pragma GCC unroll 16
      for(u32 t= 0; t < G; ++t) v&= ld(p[t] + 32 * j);
      std::memcpy(out + 32 * j, &v, 32);
     }
   }
   for(std::size_t k= 0; k < np; ++k) {
    off[k]+= pats[k].step;
    if(off[k] >= pats[k].period) off[k]-= pats[k].period;
   }
  }
 }
};
inline const Presieve& presieve() {
 static const Presieve pre;
 return pre;
}
// [7, N] を SW 語 (8 SW バイト) の区画ごとに篩い、区画ごとに f(語の列, 語の数, 区画の先頭の数) を呼ぶ。語の列は N より
// 大きい数と 1 のビットを落としてある。PSPLIT 以下の素数は 210 の車輪で消す。7 ≤ N。
template <u32 SW= 4096, u32 PSPLIT= 4000, class F> inline void sieve(u32 N, F&& f) {
 static_assert(SW % 4 == 0 && 7 * PSPLIT <= 8 * SW);
 constexpr u32 PMAX= Presieve::PMAX;
 u32 sq= (u32)std::sqrt((double)N);
 while((u64)sq * sq > N) --sq;
 while((u64)(sq + 1) * (sq + 1) <= N) ++sq;
 std::vector<u32> sp;  // PMAX より大きく √N 以下の素数 (奇数だけの篩)
 {
  std::vector<u8> comp(sq / 2 + 1, 0);
  for(u32 i= 1; (2 * i + 1) * (2 * i + 1) <= sq; ++i)
   if(!comp[i])
    for(u32 j= (2 * i + 1) * (2 * i + 1) / 2; j <= sq / 2; j+= 2 * i + 1) comp[j]= 1;
  for(u32 i= 3; i <= sq / 2; ++i)
   if(!comp[i] && 2 * i + 1 > PMAX) sp.push_back(2 * i + 1);
 }
 // 小さい素数と大きい素数を、それぞれ類ごとに並べる。小さい素数は m = 210 floor(p / 210) + 1 のまとまりから、大きい素数は
 // p^2 を含む周から消す (その手前の倍数も合成数)。p < 210 は m = 1 で p 自身も消すので、最初の区画で戻す。
 std::vector<SievingPrime> st(sp.size());
 u32 bs[9]= {}, bl[9]= {};
 {
  const u32 nsmall= (u32)(std::upper_bound(sp.begin(), sp.end(), PSPLIT) - sp.begin());
  u32 cs[8]= {}, cl[8]= {};
  for(u32 k= 0; k < sp.size(); ++k) ++(k < nsmall ? cs : cl)[BI[sp[k] % 30]];
  for(u32 c= 0; c < 8; ++c) bs[c + 1]= bs[c] + cs[c];
  bl[0]= nsmall;
  for(u32 c= 0; c < 8; ++c) bl[c + 1]= bl[c] + cl[c];
  u32 ps[8], pl[8];
  for(u32 c= 0; c < 8; ++c) ps[c]= bs[c], pl[c]= bl[c];
  for(u32 k= 0; k < sp.size(); ++k) {
   const u32 p= sp[k], c= BI[p % 30];
   if(k < nsmall) st[ps[c]++]= SievingPrime{7 * p * (p / 210) + p / 30, p / 30};
   else st[pl[c]++]= SievingPrime{p / 30 * (p + 1), p / 30};
  }
 }
 const Presieve& pre= presieve();
 // 区画の後ろの余白 (はみ出しは max(7 PSPLIT, √N) バイト未満) と、次の区画の先頭に AND するはみ出し
 const u32 SL= (std::max(sq, 7 * PSPLIT) + 64) / 64 * 64;
 std::vector<u64> seg64(SW + SL / 8 + 8), carry(SL / 8, ~u64(0));
 u8* const seg= reinterpret_cast<u8*>(seg64.data());
 constexpr u64 SEGN= 240 * (u64)SW;
 for(u64 low= 0; low <= N; low+= SEGN) {
  const u64 high= std::min<u64>(low + SEGN, (u64)N + 1);
  const u32 nbytes= (u32)((high - low + 29) / 30), nw= (nbytes + 7) / 8;
  pre.fill(seg, nbytes, low / 30);
  std::memset(seg + nbytes, 0xff, SL);
  cross_all<true>(seg, nbytes, st.data(), bs);
  for(u32 i= 0; i < SL / 8; ++i) seg64[i]&= carry[i];
  cross_all<false>(seg, nbytes, st.data(), bl);
  std::memcpy(carry.data(), seg + nbytes, SL);
  if(low == 0) {
   seg[0]&= (u8)~1u;                                            // 1 は素数でない
   for(u32 p= 7; p <= std::max(PMAX, 210u) && p <= N; p+= 2) {  // 型と 210 の車輪の m = 1 で消した素数を戻す
    bool pr= p % 3 && p % 5;
    for(u32 d= 7; d * d <= p; d+= 2) pr&= p % d != 0;
    if(pr) seg[p / 30]|= (u8)(1u << BI[p % 30]);
   }
  }
  // high 以上の数のビットを落とす
  const u64 last= high - 1 - low;
  u64 mlast= 0;
  for(u32 k= 0; k < 8; ++k)
   for(u32 i= 0; i < 8; ++i)
    if(30 * k + WR[i] <= last % 240) mlast|= u64(1) << (8 * k + i);
  seg64[last / 240]&= mlast;
  f((const u64*)seg64.data(), nw, (u32)low);
 }
}
// π(N) の上からの見積もり (N > 1 で π(N) < 1.26 N / ln N)
inline std::size_t pi_bound(u32 N) { return N < 17 ? 8 : (std::size_t)(1.26 * N / std::log((double)N)) + 8; }
inline std::vector<u32> all(u32 N) {
 std::vector<u32> res;
 res.reserve(pi_bound(N));
 for(u32 p: {2u, 3u, 5u})
  if(p <= N) res.push_back(p);
 if(N < 7) return res;
 sieve(N, [&](const u64* s, u32 nw, u32 low) {
  u32 c= 0;
  for(u32 w= 0; w < nw; ++w) c+= (u32)std::popcount(s[w]);
  const std::size_t z= res.size();
  res.resize(z + c);
  u32* out= res.data() + z;
  for(u32 w= 0; w < nw; ++w)
   for(u64 x= s[w]; x; x&= x - 1) *out++= low + 240 * w + BITVAL[std::countr_zero(x)];
 });
 return res;
}
// 語 W[0, n) の各項 (語の番号と、その語より前の個数) で、語の中の k 番目のビットを引いて sel に足す。項 i の番号は t0 + i A。
inline void select_words(const u64* W, u32 n, u32 t0, u32 A, u32 low, const u64* s, std::vector<u32>& sel) {
 const std::size_t z= sel.size();
 sel.resize(z + n);
 u32* out= sel.data() + z;
 for(u32 i= 0; i < n; ++i) {
  const u32 w= (u32)W[i], k= t0 + i * A - (u32)(W[i] >> 32);
  out[i]= low + 240 * w + BITVAL[kth_bit(s[w], k)];
 }
}
inline std::pair<u32, std::vector<u32>> select(u32 N, u32 A, u32 B) {
 if(A == 1) {
  std::vector<u32> v= all(N);
  return {(u32)v.size(), std::move(v)};
 }
 std::vector<u32> sel;
 sel.reserve(pi_bound(N) / A + 8);
 u64 cnt= 0, next= B;  // cnt は数えた素数の個数、next は次に取り出す番号 (0 始まり)
 for(u32 p: {2u, 3u, 5u})
  if(p <= N) {
   if(cnt == next) sel.push_back(p), next+= A;
   ++cnt;
  }
 if(N < 7) return {(u32)cnt, sel};
 constexpr u32 SW= 4096;
 const std::unique_ptr<u64[]> wc(new u64[4 * SW + 64]);  // 取り出す素数がある語の番号と、その語より前の個数
 const u64 MA= ~u64(0) / A + 1;                          // ceil(2^64 / A)
 sieve<SW>(N, [&](const u64* s, u32 nw, u32 low) {
  if(A >= 16 && A < (1u << 24)) {
   // 区画の中で取り出すのは 4 nw + 1 個以下なので、wc は溢れない。c1 + K < 2^32 なので掛け算の割り算は正確
   const u32 t0= (u32)(next - cnt), K= A - 1 - t0;  // t0 < A は区画の中の最初の番号
   u32 c= 0, n= 0;
   u64* const W= wc.get();
   for(u32 w= 0; w < nw; ++w) {
    const u32 c1= c + (u32)std::popcount(s[w]);
    const u64 e= (u64)c << 32 | w;
    W[n]= e;
    const u32 n1= (u32)((u128)(c1 + K) * MA >> 64);  // 番号が c1 未満の取り出しの数
    if(n1 > n + 1) [[unlikely]]
     for(u32 j= n + 1; j < n1; ++j) W[j]= e;  // 同じ語の 2 個目から
    n= n1, c= c1;
   }
   select_words(W, n, t0, A, low, s, sel);
   next= cnt + t0 + (u64)n * A, cnt+= c;
  } else {
   u64 cn= cnt, nx= next;
   for(u32 w= 0; w < nw; ++w) {
    const u64 x= s[w];
    const u32 pc= (u32)std::popcount(x);
    for(; nx < cn + pc; nx+= A) sel.push_back(low + 240 * w + BITVAL[kth_bit(x, (u32)(nx - cn))]);
    cn+= pc;
   }
   cnt= cn, next= nx;
  }
 });
 return {(u32)cnt, sel};
}
}
// N 以下の素数を小さい順にすべて返す。0 ≤ N < 2^32。
inline std::vector<unsigned> enumerate_primes(unsigned N) { return enumerate_primes_internal::all(N); }
// N 以下の素数の個数 π(N) と、小さい方から数えて B 番目 (0 始まり) から A 個おきの素数を返す。0 ≤ N < 2^32、0 ≤ B < A。
inline std::pair<unsigned, std::vector<unsigned>> enumerate_primes(unsigned N, unsigned A, unsigned B) {
 assert(B < A);
 return enumerate_primes_internal::select(N, A, B);
}
