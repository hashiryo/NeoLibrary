#pragma once
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <numeric>
#include <utility>
#include <vector>
#include "neo/number_theory/factorize.hpp"
#include "neo/number_theory/inv_gcd.hpp"
// x^k ≡ y (mod m) となる最小の k を求める。
//
// x と m が互いに素でない部分: k ≤ 32 は直に調べる。そこで見つからなければ、g = gcd(x, m) > 1 の間、y が g で割り切れる
// ことを確かめて m、y を g で割り、d に x / g を掛けて、段の数 c を数える (c ≤ 31)。残りは d x^{k'} ≡ y (mod m')、
// gcd(x, m') = 1 で、答えは k = c + k' (k' は最小の解)。
//
// 互いに素な部分は Pohlig-Hellman。m' を factorize_pairs で分解して φ(m') の素因数分解を作り、x の位数 n を φ(m') から
// 素因数を外していって求める。t = y d^{-1} について t^n ≠ 1 なら解は無い。n の素因数 r^e ごとに、γ = x^{n / r^e}、
// τ = t^{n / r^e} として γ^u = τ の u を r 進の桁ごとに求め、CRT でまとめる (まとめた値が [0, n) の最小の解)。桁は位数 r の
// ζ = γ^{r^{e-1}} の離散対数で、r ≤ 64 なら表を順に見る。それより大きければ baby-step giant-step で引く。表は 4 バイトの項
// (鍵のハッシュの 16 bit の指紋と添字) で B = √r / 2 の 8 倍以上の大きさにとり、ほとんど空きにしておく (表に無い鍵を引く
// たびの分岐が外れない)。r > 2^22 で (Z/m')^* が巡回群 (m' が素数冪かその 2 倍) のときは、Pollard の rho 法で引く
// (32 個の乗数の r-adding walk を 16 本並べ、区別点に着いた値を小さい表に入れて、別の歩みとぶつかったら解く)。大きい表を
// 読み書きしないので速い。巡回群でない m' では引く元が部分群に無いことがあり、rho 法が終わらないので baby-step giant-step
// を使う (見つからなければ解は無い)。剰余の掛け算は Barrett (m' は偶数でもよい)。
namespace discrete_log_internal {
using u32= unsigned;
using u64= unsigned long long;
using u128= unsigned __int128;
struct Barrett {
 u64 m, im;  // im = floor(2^64 / m)。m ≥ 2、m < 2^32
 explicit Barrett(u64 m_): m(m_), im((u64)(((u128)1 << 64) / m_)) {}
 u64 mul(u64 a, u64 b) const {
  const u64 t= a * b, q= (u64)((u128)t * im >> 64);
  const u64 r= t - q * m;
  return r >= m ? r - m : r;
 }
 // 掛ける値はマスクで選ぶ (gcc は三項演算子を分岐にする)
 u64 pow(u64 b, u64 e) const {
  u64 r= 1;
  for(; e; e>>= 1, b= mul(b, b)) r= mul(r, 1 ^ ((1 ^ b) & (0 - (e & 1))));
  return r;
 }
};
inline std::vector<u32>& table_buf() {
 static std::vector<u32> buf;
 return buf;
}
// 位数 r の ζ の冪 ζ^j (j < B) の表。c = ζ^l の l を引く (無ければ -1)。項は (指紋 << 16) | (j + 1)、0 は空き。
struct OrderRLog {
 const Barrett& M;
 u64 zeta, r, B= 0, H= 0, step= 0;  // step = ζ^{-B}
 int lg= 0;
 u32* tab= nullptr;
 OrderRLog(const Barrett& M_, u64 zeta_, u64 r_): M(M_), zeta(zeta_), r(r_) {
  if(r <= 64) return;
  B= (((u64)std::sqrt((double)r) >> 1) + 8) & ~7ull;
  lg= 4;
  while((1ull << lg) < 8 * B) ++lg;
  H= 1ull << lg;
  auto& buf= table_buf();
  if(buf.size() < H) buf.resize(H);
  tab= buf.data();
  std::memset(tab, 0, H * sizeof(u32));
  u64 zp[8];
  zp[0]= 1;
  for(int d= 1; d < 8; ++d) zp[d]= M.mul(zp[d - 1], zeta);
  const u64 z8= M.mul(zp[7], zeta);
  for(u64 j= 0; j < B; j+= 8)
   for(int d= 0; d < 8; ++d) ins(zp[d], j + d), zp[d]= M.mul(zp[d], z8);
  step= M.pow(zp[0], r - 1);  // zp[0] = ζ^B
 }
 static u64 hash(u64 k) { return k * 0x9e3779b97f4a7c15ull; }
 void ins(u64 k, u64 j) {
  const u64 h= hash(k);
  u64 i= h >> (64 - lg);
  const u32 e= (u32)(h >> 16 & 0xffff0000u) | (u32)(j + 1);
  while(tab[i]) i= (i + 1) & (H - 1);
  tab[i]= e;
 }
 long long look(u64 k) const {
  const u64 h= hash(k);
  const u32 fp= (u32)(h >> 16 & 0xffff0000u);
  for(u64 i= h >> (64 - lg); tab[i]; i= (i + 1) & (H - 1))
   if((tab[i] & 0xffff0000u) == fp) {
    const u64 j= (tab[i] & 0xffffu) - 1;
    if(M.pow(zeta, j) == k) return (long long)j;
   }
  return -1;
 }
 long long operator()(u64 c) const {
  if(r <= 64) {
   u64 v= 1;
   for(u64 l= 0; l < r; ++l, v= M.mul(v, zeta))
    if(v == c) return (long long)l;
   return -1;
  }
  u64 sp[8], u[8];
  sp[0]= 1;
  for(int d= 1; d < 8; ++d) sp[d]= M.mul(sp[d - 1], step);
  const u64 s8= M.mul(sp[7], step);
  for(int d= 0; d < 8; ++d) u[d]= M.mul(c, sp[d]);
  for(u64 i= 0; i <= r / B + 8; i+= 8) {
   for(int d= 0; d < 8; ++d)
    if(const long long j= look(u[d]); j >= 0) return (long long)(((i + d) * B + (u64)j) % r);
   for(int d= 0; d < 8; ++d) u[d]= M.mul(u[d], s8);
  }
  return -1;
 }
};
// 素数 r の位数の ζ について、c = ζ^l の l を Pollard の rho 法で求める。c は ζ の冪とする。
inline u64 rho_log(const Barrett& M, u64 zeta, u64 c, u64 r) {
 if(c == 1) return 0;
 constexpr int K= 32, W= 16;
 u64 seed= r * 0x9e3779b97f4a7c15ull ^ c ^ (zeta << 7);
 auto rnd= [&] {
  u64 z= (seed+= 0x9e3779b97f4a7c15ull);
  z= (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
  z= (z ^ (z >> 27)) * 0x94d049bb133111ebull;
  return z ^ (z >> 31);
 };
 u64 mx[K], ma[K], mb[K];
 for(int k= 0; k < K; ++k) ma[k]= rnd() % r, mb[k]= rnd() % r, mx[k]= M.mul(M.pow(zeta, ma[k]), M.pow(c, mb[k]));
 // 区別点の割合 2^{-t}: 1 本の歩みが区別点に着くまで √r / (8 W) 歩ほど
 int t= 4;
 while((1ull << (2 * t + 2)) * 64 * W * W <= r) ++t;
 const u64 dmask= (1ull << t) - 1, maxlen= 20ull << t;
 constexpr int LG= 12;
 constexpr u32 EMPTY= ~0u;
 std::vector<u32> key(1u << LG, EMPTY), va(1u << LG), vb(1u << LG);  // 区別点の表 (x → a, b)
 int stored= 0;
 u64 x[W], a[W], b[W], len[W];
 auto start= [&](int w) {
  a[w]= rnd() % r, b[w]= rnd() % r, len[w]= 0;
  x[w]= M.mul(M.pow(zeta, a[w]), M.pow(c, b[w]));
 };
 for(int w= 0; w < W; ++w) start(w);
 for(;;)
  for(int w= 0; w < W; ++w) {
   const int k= (int)((x[w] * 0xc2b2ae3d27d4eb4full) >> 59);
   x[w]= M.mul(x[w], mx[k]);
   a[w]+= ma[k], b[w]+= mb[k];  // 法 r では割らない (区別点で割る)
   ++len[w];
   if(((x[w] * 0x9e3779b97f4a7c15ull) >> 40 & dmask) == 0) [[unlikely]] {
    a[w]%= r, b[w]%= r;
    u32 h= (u32)((x[w] * 0xd6e8feb86659fd93ull) >> (64 - LG));
    while(key[h] != EMPTY && key[h] != x[w]) h= (h + 1) & ((1u << LG) - 1);
    if(key[h] == EMPTY) {
     if(stored * 2 < (1 << LG)) key[h]= (u32)x[w], va[h]= (u32)a[w], vb[h]= (u32)b[w], ++stored;
     len[w]= 0;
    } else if(vb[h] != b[w]) {
     // ζ^{a'} c^{b'} = ζ^{a} c^{b} → c^{b - b'} = ζ^{a' - a}
     const u64 l= (u64)((u128)((va[h] + r - a[w]) % r) * ::inv_mod<u64>((b[w] + r - vb[h]) % r, r) % r);
     if(M.pow(zeta, l) == c) return l;
     start(w);
    } else start(w);
   } else if(len[w] > maxlen) [[unlikely]]
    start(w);
  }
}
// gcd(x, m) = 1、m ≥ 2 で x^k ≡ t (mod m) となる最小の k (無ければ -1)
inline long long coprime_log(u64 x, u64 t, u64 m) {
 const Barrett M(m);
 std::vector<std::pair<u64, int>> phi;  // φ(m) の素因数分解
 u64 ph= 1;
 auto add= [&](u64 q, int e) {
  for(auto& [r, c]: phi)
   if(r == q) return void(c+= e);
  phi.emplace_back(q, e);
 };
 const auto fm= factorize_pairs(m);
 const bool cyclic= fm.size() == 1 ? fm[0].first != 2 || fm[0].second <= 2 : fm.size() == 2 && fm[0].first == 2 && fm[0].second == 1;
 for(auto [p, e]: fm) {
  for(int i= 1; i < e; ++i) ph*= p;
  if(e > 1) add(p, e - 1);
  ph*= p - 1;
  for(auto [q, f]: factorize_pairs(p - 1)) add(q, f);
 }
 u64 n= ph;  // x の位数
 for(auto& [r, c]: phi)
  while(c > 0 && M.pow(x, n / r) == 1) n/= r, --c;
 if(M.pow(t, n) != 1) return -1;
 u64 X= 0, mod= 1;
 auto crt= [&](u64 xr, u64 re) {  // X ≡ X (mod mod)、≡ xr (mod re) をまとめる
  const u64 k= (u64)((u128)((xr + re - X % re) % re) * ::inv_mod<u64>(mod % re, re) % re);
  X+= mod * k, mod*= re;
 };
 for(auto [r, e]: phi) {
  if(e == 0) continue;
  u64 re= 1;
  for(int i= 0; i < e; ++i) re*= r;
  const u64 gam= M.pow(x, n / re), tau= M.pow(t, n / re), zeta= M.pow(gam, re / r);
  if(cyclic && r > (1u << 22)) {
   crt(rho_log(M, zeta, tau, r), re);  // r^2 > 2^44 > n なので e = 1
   continue;
  }
  const OrderRLog lg(M, zeta, r);
  u64 xr= 0, ri= 1, ginv_ri= M.pow(gam, re - 1);  // ginv_ri = γ^{-r^i}
  u64 cur= tau;                                   // τ γ^{-xr}
  for(int i= 0; i < e; ++i) {
   const long long l= lg(M.pow(cur, re / ri / r));
   if(l < 0) return -1;
   if(l) cur= M.mul(cur, M.pow(ginv_ri, (u64)l)), xr+= (u64)l * ri;
   ginv_ri= M.pow(ginv_ri, r), ri*= r;
  }
  crt(xr, re);
 }
 return (long long)X;
}
}
// x^k ≡ y (mod m) となる最小の k ≥ 0 を返す。無ければ -1。0^0 = 1 とする (m = 1 なら 0)。1 ≤ m < 2^32。
inline long long discrete_log(unsigned long long x, unsigned long long y, unsigned long long m) {
 using namespace discrete_log_internal;
 assert(1 <= m && m < (1ull << 32));
 if(m == 1) return 0;
 x%= m, y%= m;
 {  // k ≤ 32 は直に調べる
  u64 v= 1;
  for(int k= 0; k <= 32; ++k, v= v * x % m)
   if(v == y) return k;
 }
 u64 mm= m, b= y, d= 1;
 int c= 0;
 for(u64 g; (g= std::gcd(x, mm)) > 1;) {
  if(b % g) return -1;
  mm/= g, b/= g, ++c;
  d= d * (x / g) % mm;
 }
 if(mm == 1) return -1;  // k ≥ c の解は k ≤ 32 で見つかっているはず
 const u64 t= (u64)((u128)b * ::inv_mod<u64>(d % mm, mm) % mm);
 const long long k= coprime_log(x % mm, t, mm);
 return k < 0 ? -1 : k + c;
}
