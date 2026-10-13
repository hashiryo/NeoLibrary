#pragma once
#include <array>
#include <cassert>
#include <cmath>
#include <numeric>
#include <tuple>
#include <utility>
#include <vector>
#include "neo/number_theory/factorize.hpp"
#include "neo/number_theory/inv_gcd.hpp"
// 素数 p を法として X^k ≡ y となる X を 1 つ求める。m = p - 1、g = gcd(k, m) とすると、解があるのは y^(m / g) = 1 のとき
// (y = 0 と k = 0 は別に扱う)。k = g u とすると gcd(u, m / g) = 1 なので、v = u^{-1} mod (m / g) として、y^v の任意の
// g 乗根 W が答えになる (W^k = (y^v)^u = y)。y^(m / g) と y^v は 2 乗の列を共有して求める。
//
// g の素因数 r ごとに、r^e (r^e ‖ g) 乗根を取っていく。m = r^s t (gcd(r, t) = 1) として、群を r の Sylow 部分群 S (位数
// r^s) と位数 t の部分 T に分ける。a の T の成分 a_T = a^{e_T} (e_T ≡ 1 mod t、≡ 0 mod r^s) の q = r^e 乗根は、q が t と
// 互いに素なので x_T = a_T^{q^{-1} mod t} と 1 つに決まる。S の成分 a_S = a^{e_S} は、r 乗非剰余 ρ から z = ρ^t (位数 r^s)
// を作り、a_S = z^L の L を Pohlig-Hellman で r 進の桁ごとに求める。桁は位数 r の ζ = z^{r^{s-1}} = ρ^{m/r} の離散対数で、
// r ≤ 64 は表を順に見て、それより大きい r は baby-step giant-step で引く (表に c ζ^j を入れ、ζ^{iB} を引く。掛け算の鎖は
// 4 本に分けて並べる)。a_S は q 乗なので L は r^e で割り切れ、x_S = z^{L / r^e}。s = e なら x_S = 1。T の成分を a_T の冪で
// 選ぶので、x = x_T x_S は残りの素因数の冪乗根も持つ。
//
// 剰余の掛け算は Montgomery で、p < 2^30 は R = 2^32 で値を 2p 未満のまま持ち、p < 2^62 は R = 2^64 で値を 2p 未満のまま、
// それ以上は p 未満に保つ。累乗は右から左の 2 進で、掛ける値はマスクの算術で選ぶ (分岐なし)。
namespace kth_root_mod_internal {
using u32= unsigned;
using u64= unsigned long long;
using u128= unsigned __int128;
struct Mont32 {
 using T= u32;
 u32 p, pinv, r1, r2;  // pinv = p^{-1} mod 2^32、r1 = 2^32 mod p、r2 = 2^64 mod p
 explicit Mont32(u64 p_): p((u32)p_) {
  u32 x= p;
  for(int i= 0; i < 5; ++i) x*= 2 - p * x;
  pinv= x, r1= (u32)((1ull << 32) % p), r2= (u32)((u64)r1 * r1 % p);
 }
 u32 mul(u32 a, u32 b) const {
  const u64 t= (u64)a * b;
  return (u32)(t >> 32) + p - (u32)((u64)((u32)t * pinv) * p >> 32);
 }
 u32 norm(u32 a) const { return a >= p ? a - p : a; }
};
template <bool LAZY> struct Mont64 {
 using T= u64;
 u64 p, pinv, r1, r2;  // pinv = p^{-1} mod 2^64、r1 = 2^64 mod p、r2 = 2^128 mod p
 explicit Mont64(u64 p_): p(p_) {
  u64 x= p;
  for(int i= 0; i < 6; ++i) x*= 2 - p * x;
  pinv= x, r1= (u64)(((u128)1 << 64) % p), r2= (u64)((u128)r1 * r1 % p);
 }
 u64 mul(u64 a, u64 b) const {
  const u128 t= (u128)a * b;
  const u64 hi= (u64)(t >> 64), m= (u64)((u128)((u64)t * pinv) * p >> 64);
  if constexpr(LAZY) return hi + p - m;
  else return hi >= m ? hi - m : hi - m + p;
 }
 u64 norm(u64 a) const {
  if constexpr(LAZY) return a >= p ? a - p : a;
  else return a;
 }
};
template <class M> struct Ops: M {
 using T= typename M::T;
 using M::M, M::mul, M::norm, M::r1, M::r2, M::p;
 T to(u64 a) const { return mul((T)a, r2); }
 u64 from(T a) const { return norm(mul(a, 1)); }
 bool eq(T a, T b) const { return norm(a) == norm(b); }
 T sel(T b, u64 bit) const { return r1 ^ ((r1 ^ b) & (T(0) - (T)(bit & 1))); }
 T pow(T b, u64 e) const {
  T r= r1;
  for(; e; e>>= 1, b= mul(b, b)) r= mul(r, sel(b, e));
  return r;
 }
 // (b^e1, b^e2)。2 乗の列を共有する。
 std::pair<T, T> pow2(T b, u64 e1, u64 e2) const {
  T x= r1, y= r1;
  for(; e1 | e2; e1>>= 1, e2>>= 1, b= mul(b, b)) x= mul(x, sel(b, e1)), y= mul(y, sel(b, e2));
  return {x, y};
 }
};
// a^{-1} mod m (gcd(a, m) = 1、m ≥ 1)
inline u64 inv_m(u64 a, u64 m) { return m == 1 ? 0 : (m >> 32 ? ::inv_mod<u64>(a % m, m) : (u64)::inv_mod<u32>((u32)(a % m), (u32)m)); }
// 256 未満の奇素数と、その逆元 (mod 2^32) と floor((2^32 - 1) / p)
struct SmallPrime {
 u32 p, inv, lim;
};
inline constexpr auto SMALLP= [] {
 std::array<SmallPrime, 53> t{};
 int k= 0;
 for(u32 p= 3; p < 256; p+= 2) {
  bool pr= true;
  for(u32 d= 3; d * d <= p; d+= 2) pr&= p % d != 0;
  if(!pr) continue;
  u32 x= p;
  for(int i= 0; i < 5; ++i) x*= 2 - p * x;
  t[k++]= {p, x, ~0u / p};
 }
 return t;
}();
// 位数 r の元 zeta について、c = zeta^l の l (0 ≤ l < r) を求める
template <class O> u64 dlog_order_r(const O& M, typename O::T zeta, typename O::T c, u64 r) {
 using T= typename O::T;
 if(r <= 64) {
  T v= M.r1;
  for(u64 l= 0; l < r; ++l, v= M.mul(v, zeta))
   if(M.eq(v, c)) return l;
  return 0;
 }
 const u64 B= ((u64)std::sqrt(r * 0.5) + 4) & ~3ull;  // 4 の倍数
 // 表は B の 4 倍以上の 2 の冪 (B ≤ 128 ならスタックの 512 個)
 int lg= 9;
 while((1ull << lg) < 4 * B) ++lg;
 const u64 H= 1ull << lg;
 constexpr T EMPTY= ~T(0);
 T skey[512];
 u32 sval[512];
 std::vector<T> hkey;
 std::vector<u32> hval;
 T* key= skey;
 u32* val= sval;
 if(H > 512) hkey.assign(H, EMPTY), hval.resize(H), key= hkey.data(), val= hval.data();
 else
  for(int i= 0; i < 512; ++i) skey[i]= EMPTY;
 auto hash= [&](T k) { return (u64)(((u64)k * 0x9e3779b97f4a7c15ull) >> (64 - lg)); };
 auto ins= [&](T v, u64 j) {
  const T k= M.norm(v);
  u64 h= hash(k);
  while(key[h] != EMPTY && key[h] != k) h= (h + 1) & (H - 1);
  if(key[h] == EMPTY) key[h]= k, val[h]= (u32)j;
 };
 // baby-step: c ζ^j (j < B)
 const T z2= M.mul(zeta, zeta), z3= M.mul(z2, zeta), z4= M.mul(z3, zeta);
 T v0= c, v1= M.mul(c, zeta), v2= M.mul(c, z2), v3= M.mul(c, z3);
 for(u64 j= 0; j < B; j+= 4) {
  ins(v0, j), ins(v1, j + 1), ins(v2, j + 2), ins(v3, j + 3);
  v0= M.mul(v0, z4), v1= M.mul(v1, z4), v2= M.mul(v2, z4), v3= M.mul(v3, z4);
 }
 // giant-step: ζ^{iB} (i ≥ 1) が c ζ^j と一致すれば c = ζ^{iB - j}
 const T g1= M.pow(zeta, B), g2= M.mul(g1, g1), g3= M.mul(g2, g1), g4= M.mul(g3, g1);
 T u[4]= {g1, g2, g3, g4};
 auto look= [&](T v) -> long long {
  const T k= M.norm(v);
  for(u64 h= hash(k); key[h] != EMPTY; h= (h + 1) & (H - 1))
   if(key[h] == k) return val[h];
  return -1;
 };
 for(u64 i= 1; i <= r / B + 4; i+= 4) {
  for(int d= 0; d < 4; ++d)
   if(const long long j= look(u[d]); j >= 0) return (u64)(((u128)(i + d) * B + (u128)r * B - (u64)j) % r);
  for(int d= 0; d < 4; ++d) u[d]= M.mul(u[d], g4);
 }
 return 0;
}
// a (q = r^e 乗剰余、Montgomery 形) の r^e 乗根 (Montgomery 形)。m = p - 1。
template <class O> typename O::T root_prime_power(const O& M, typename O::T a, u64 r, int e, u64 m) {
 using T= typename O::T;
 u64 t= m;
 int s= 0;
 while(t % r == 0) t/= r, ++s;
 u64 rs= 1, q= 1;
 for(int i= 0; i < s; ++i) rs*= r;
 for(int i= 0; i < e; ++i) q*= r;
 // T の成分の q 乗根: a^{e_T w}、e_T = r^s (r^s)^{-1} mod t、w = q^{-1} mod t
 const u64 v1= inv_m(rs % t, t), eT= (u64)((u128)rs * v1 % m), w= s == e ? v1 : inv_m(q % t, t);
 const u64 eTw= (u64)((u128)eT * w % m);
 if(s == e) return M.pow(a, eTw);
 // S の成分 a_S = a^{e_S}、e_S = t (t^{-1} mod r^s)
 const u64 eS= (u64)((u128)t * inv_m(t % rs, rs) % m);
 const auto [xT, aS]= M.pow2(a, eTw, eS);
 // r 乗非剰余 ρ から z = ρ^t (位数 r^s)、ζ = z^{r^{s-1}} = ρ^{m / r} (位数 r)
 u64 rho= 1;
 T zeta, z;
 do std::tie(zeta, z)= M.pow2(M.to(++rho), m / r, t);
 while(M.eq(zeta, M.r1));
 // a_S = z^L の L を r 進で求める (下の e 桁は 0)。最後の桁のあとは a_S z^{-L} を更新しない
 u64 L= 0, re= q;        // re = r^i
 T cur= aS, zinv_re= 0;  // a_S z^{-L} と z^{-r^i}
 for(int i= e; i < s; ++i) {
  const T c= M.pow(cur, rs / re / r);  // (a_S z^{-L})^{r^{s-1-i}}
  const u64 l= dlog_order_r(M, zeta, c, r);
  L+= l * re;
  if(i + 1 < s) {
   if(i == e) zinv_re= M.pow(M.pow(z, rs - 1), q);
   if(l) cur= M.mul(cur, M.pow(zinv_re, l));
   zinv_re= M.pow(zinv_re, r);
  }
  re*= r;
 }
 return M.mul(xT, M.pow(z, L / q));
}
template <class O> u64 solve(u64 k, u64 y, u64 p) {
 const O M(p);
 const u64 m= p - 1, g= std::gcd(k, m);
 const u64 v= inv_m(k / g % (m / g), m / g);
 const auto [chk, yv]= M.pow2(M.to(y), m / g, v);
 if(!M.eq(chk, M.r1)) return ~0ull;
 auto z= yv;
 if(g < (1u << 16)) {
  // 256 未満の素数で試し割り (g < 2^16 なので、残りは 1 か素数)
  u32 h= (u32)g;
  if(const int c= __builtin_ctz(h)) h>>= c, z= root_prime_power(M, z, 2, c, m);
  for(const auto& [r, inv, lim]: SMALLP) {
   if(r * r > h) break;
   if(h * inv > lim) continue;
   int e= 0;
   do h*= inv, ++e;
   while(h * inv <= lim);
   z= root_prime_power(M, z, r, e, m);
  }
  if(h > 1) z= root_prime_power(M, z, h, 1, m);
 } else
  for(auto [r, e]: factorize_pairs(g)) z= root_prime_power(M, z, r, e, m);
 return M.from(z);
}
}
// 素数 p を法として X^k ≡ y となる 0 ≤ X < p を 1 つ返す。無ければ -1。0^0 = 1 とする。2 ≤ p < 2^63 で、p が素数で
// あることは確かめない。
inline long long kth_root_mod(unsigned long long k, unsigned long long y, unsigned long long p) {
 using namespace kth_root_mod_internal;
 assert(2 <= p && p < (1ull << 63));
 y%= p;
 if(y == 0) return k > 0 ? 0 : -1;
 if(k == 0) return y == 1 ? 1 : -1;
 if(p == 2) return 1;
 u64 x;
 if(p < (1u << 30)) x= solve<Ops<Mont32>>(k, y, p);
 else if(p < (1ull << 62)) x= solve<Ops<Mont64<true>>>(k, y, p);
 else x= solve<Ops<Mont64<false>>>(k, y, p);
 return x == ~0ull ? -1 : (long long)x;
}
