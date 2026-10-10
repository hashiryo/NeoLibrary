#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <type_traits>
#include <utility>
#include <vector>
#include "neo/number_theory/gcd.hpp"
#include "neo/number_theory/inv_gcd.hpp"
#include "neo/number_theory/is_prime.hpp"
// 97 以下の素数で試し割りしてから、残りの合成数の因数を 1 つずつ見つけて再帰する。平方数は平方根で割り、
// 2^48 未満は Pollard rho (Brent、c の違う 2 本の列を同時に進める)、それより大きければ ECM で探す。
//
// ECM は Montgomery 曲線 B y^2 = x^3 + A x^2 + x の x 座標だけで計算する (X:Z の射影座標)。曲線は Suyama の形で、
// 群の位数が 12 で割り切れる。σ = 6, 7, ... に u = σ^2 - 5、v = 4σ、点 (u^3 : v^3)、
// (A + 2) / 4 = (v - u)^3 (3u + v) / (16 u^3 v) で、分母の逆元は inv_gcd で求める。
// stage 1 は B1 = 150 以下の素数冪の積 (2 と 3 は 3 個と 2 個余分に掛ける) を 64 bit ずつに束ねて Montgomery ladder で掛け、
// stage 2 は D = 210 の baby-step giant-step で、m = 1..16 と D と互いに素な b < D / 2 の組のうち、mD ± b のどちらかが
// 素数になる組の X_m Z_b - X_b Z_m を 4 本の積に掛けていく (B2 ≈ 3400)。
// n < 2^60 では σ と σ + 1 の 2 本の曲線を同じ手順で同時に回す (1 本だと鎖の長さで決まり、掛け算器が空く)。
// 2 本分の逆元と gcd は、積を取って 1 回にまとめる。
//
// 剰余の掛け算は R = 2^64 の Montgomery。n < 2^60 なら値を 2n 未満のまま持ち、掛け算に渡すだけの和と差は
// 4n 未満のまま補正しない (4n 未満どうしの積の還元は 2n 未満に収まる)。n >= 2^60 は値を n 未満に保つ。
namespace factorize_internal {
using u64= unsigned long long;
using u128= unsigned __int128;
constexpr u64 inv64(u64 n) {
 u64 x= (3 * n) ^ 2, y= 1 - n * x;
 x*= 1 + y, y*= y;
 x*= 1 + y, y*= y;
 x*= 1 + y, y*= y;
 return x * (1 + y);
}
// n < 2^60。値は [0, 2n)、add と sub の結果は [0, 4n) で掛け算にだけ渡す。
struct Lazy {
 u64 n, ninv, n2;
 constexpr explicit Lazy(u64 n_): n(n_), ninv(inv64(n_)), n2(2 * n_) {}
 constexpr u64 mul(u64 a, u64 b) const {
  u128 t= (u128)a * b;
  return u64(t >> 64) + n - u64((u128)(u64(t) * ninv) * n >> 64);
 }
 constexpr u64 add(u64 a, u64 b) const { return a + b; }
 constexpr u64 sub(u64 a, u64 b) const { return a - b + n2; }
 // a + b - c (a, b, c < 2n) を、2n 未満の値との積に渡せる形で返す ((0, 6n)、6n 2n < n 2^64)。
 constexpr u64 add_sub(u64 a, u64 b, u64 c) const { return a + b - c + n2; }
 constexpr u64 norm(u64 x) const {
  x= x >= n2 ? x - n2 : x;
  return x >= n ? x - n : x;
 }
};
// 奇数の n < 2^64。値は [0, n)。
struct Strict {
 u64 n, ninv;
 constexpr explicit Strict(u64 n_): n(n_), ninv(inv64(n_)) {}
 constexpr u64 mul(u64 a, u64 b) const {
  u128 t= (u128)a * b;
  u64 hi= u64(t >> 64), m= u64((u128)(u64(t) * ninv) * n >> 64);
  return hi >= m ? hi - m : hi - m + n;
 }
 constexpr u64 add(u64 a, u64 b) const { return a >= n - b ? a - (n - b) : a + b; }
 constexpr u64 sub(u64 a, u64 b) const { return a >= b ? a - b : a - b + n; }
 constexpr u64 add_sub(u64 a, u64 b, u64 c) const { return sub(add(a, b), c); }
 constexpr u64 norm(u64 x) const { return x; }
};
template <class M> constexpr u64 from_mont(const M& m, u64 x) { return m.norm(m.mul(x, 1)); }
struct Pt {
 u64 X, Z;
};
template <class M> constexpr Pt xdbl(const M& m, Pt P, u64 a24) {
 const u64 ss= m.mul(m.add(P.X, P.Z), m.add(P.X, P.Z)), dd= m.mul(m.sub(P.X, P.Z), m.sub(P.X, P.Z)), t= m.sub(ss, dd);
 return {m.mul(ss, dd), m.mul(t, m.add(dd, m.mul(a24, t)))};
}
// P + Q。D = P - Q。
template <class M> constexpr Pt xadd(const M& m, Pt P, Pt Q, Pt D) {
 const u64 u= m.mul(m.sub(P.X, P.Z), m.add(Q.X, Q.Z)), v= m.mul(m.add(P.X, P.Z), m.sub(Q.X, Q.Z)), w= m.add(u, v), y= m.sub(u, v);
 return {m.mul(D.Z, m.mul(w, w)), m.mul(D.X, m.mul(y, y))};
}
constexpr int B1= 150, MH= 16, D= 210;
// stage 1 で掛ける数: B1 以下の素数冪の積 (2 と 3 は 3 個と 2 個余分に掛ける。micro-ecm にならった) を 64 bit ずつに束ねたもの。
struct Stage1 {
 std::array<u64, 8> c{};
 int len= 0;
};
inline constexpr Stage1 STAGE1= [] {
 Stage1 s;
 u64 cur= 1;
 for(int p= 2; p <= B1; ++p) {
  bool pr= true;
  for(int d= 2; d * d <= p; ++d)
   if(p % d == 0) pr= false;
  if(!pr) continue;
  int e= (p == 2 ? 3 : p == 3 ? 2 : 0);
  for(int q= p; q <= B1; q*= p) ++e;
  for(int i= 0; i < e; ++i) {
   if(cur > ~0ull / p) s.c[s.len++]= cur, cur= 1;
   cur*= p;
  }
 }
 s.c[s.len++]= cur;
 return s;
}();
// [k]P (k >= 2)。入れ替えはビットで選ぶ cmov で書き、分岐にしない。
template <class M> constexpr Pt ladder(const M& m, Pt P, u64 k, u64 a24) {
 Pt R0= P, R1= xdbl(m, P, a24);
 u64 sw= 0;
 for(int i= 62 - __builtin_clzll(k); i >= 0; --i) {
  const u64 bit= (k >> i) & 1, c= bit ^ sw;
  sw= bit;
  const Pt A{c ? R1.X : R0.X, c ? R1.Z : R0.Z}, B{c ? R0.X : R1.X, c ? R0.Z : R1.Z};
  R1= xadd(m, A, B, P), R0= xdbl(m, A, a24);
 }
 return sw ? R1 : R0;
}
// 2 本の曲線に同じ手順で同時に [k] を掛ける。1 本の ladder は 1 段が剰余乗算 3 回分の鎖で決まって掛け算器が空くので、
// それを 2 本目で埋める。ビットが同じなので入れ替えは共通で、xor とマスクで書く (gcc 15 は三項演算子だと分岐にした)。
template <class M> constexpr void ladder2(const M& m, Pt& Pa, Pt& Pb, u64 k, u64 a24a, u64 a24b) {
 Pt R0a= Pa, R1a= xdbl(m, Pa, a24a), R0b= Pb, R1b= xdbl(m, Pb, a24b);
 u64 sw= 0;
 for(int i= 62 - __builtin_clzll(k); i >= 0; --i) {
  const u64 bit= (k >> i) & 1, mk= 0 - (bit ^ sw);
  sw= bit;
  const u64 dax= (R0a.X ^ R1a.X) & mk, daz= (R0a.Z ^ R1a.Z) & mk, dbx= (R0b.X ^ R1b.X) & mk, dbz= (R0b.Z ^ R1b.Z) & mk;
  const Pt Aa{R0a.X ^ dax, R0a.Z ^ daz}, Ba{R1a.X ^ dax, R1a.Z ^ daz}, Ab{R0b.X ^ dbx, R0b.Z ^ dbz}, Bb{R1b.X ^ dbx, R1b.Z ^ dbz};
  R1a= xadd(m, Aa, Ba, Pa), R0a= xdbl(m, Aa, a24a);
  R1b= xadd(m, Ab, Bb, Pb), R0b= xdbl(m, Ab, a24b);
 }
 Pa= sw ? R1a : R0a, Pb= sw ? R1b : R0b;
}
inline constexpr std::array<int, 24> BABY= {1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103};
// m ごとに、mD ± b のどちらかが (B1, MH D + 103] の素数になる BABY の添字の集合 (384 組のうち 318 組)。
inline constexpr auto MASK= [] {
 std::array<unsigned, MH + 1> mask{};
 auto isp= [](int x) {
  for(int d= 2; d * d <= x; ++d)
   if(x % d == 0) return false;
  return x > 1;
 };
 for(int k= 1; k <= MH; ++k)
  for(int j= 0; j < 24; ++j)
   if(int lo= k * D - BABY[j], hi= k * D + BABY[j]; (lo > B1 && isp(lo)) || (hi > B1 && isp(hi))) mask[k]|= 1u << j;
 return mask;
}();
// D = 210 の baby-step giant-step で、使う組の X_m Z_b - X_b Z_m の積を返す。
template <class M> constexpr u64 stage2(const M& m, Pt Q, u64 a24) {
 // odd[b >> 1] = [b]Q。b ≡ 1, 5 (mod 6) の 2 本の鎖を [b] = [b - 6] + [6] で伸ばす。差は鎖の 2 つ前で、最初は [5] と [1]。
 std::array<Pt, 53> odd{};
 const Pt Q2= xdbl(m, Q, a24), Q3= xadd(m, Q2, Q, Q), Q5= xadd(m, Q3, Q2, Q), Q6= xdbl(m, Q3, a24);
 odd[0]= Q, odd[2]= Q5;
 Pt pa= Q5, ca= Q, pb= Q, cb= Q5;
 for(int b= 7; b <= 103; b+= 6) {
  const Pt na= xadd(m, ca, Q6, pa);
  odd[b >> 1]= na, pa= ca, ca= na;
  if(b + 4 <= 101) {
   const Pt nb= xadd(m, cb, Q6, pb);
   odd[(b + 4) >> 1]= nb, pb= cb, cb= nb;
  }
 }
 std::array<u64, 24> bx{}, bz{}, bxz{};
 for(int j= 0; j < 24; ++j) {
  const Pt B= odd[BABY[j] >> 1];
  bx[j]= B.X, bz[j]= B.Z, bxz[j]= m.mul(B.X, B.Z);
 }
 const Pt G1= xdbl(m, xadd(m, odd[51], Q2, odd[50]), a24);  // [210] = 2 ([103] + [2])、差 [101]
 Pt Gp= G1, Gc= xdbl(m, G1, a24);
 u64 acc[4]= {1, 1, 1, 1};
 unsigned c= 0;
 // (X_m - X_b)(Z_m + Z_b) + X_b Z_b - X_m Z_m = X_m Z_b - X_b Z_m。積は 4 本に分けて鎖をつながない。
 auto pairs= [&](Pt G, unsigned mk) {
  const u64 gxz= m.mul(G.X, G.Z);
  for(; mk; mk&= mk - 1, ++c) {
   const int j= __builtin_ctz(mk);
   acc[c & 3]= m.mul(acc[c & 3], m.add_sub(m.mul(m.sub(G.X, bx[j]), m.add(G.Z, bz[j])), bxz[j], gxz));
  }
 };
 pairs(Gp, MASK[1]), pairs(Gc, MASK[2]);
 for(int k= 3; k <= MH; ++k) {
  const Pt Gn= xadd(m, Gc, G1, Gp);
  Gp= Gc, Gc= Gn, pairs(Gc, MASK[k]);
 }
 return m.mul(m.mul(acc[0], acc[1]), m.mul(acc[2], acc[3]));
}
// σ の曲線の値。逆元はまだ取らない。
struct Raw {
 u64 U3, V3, num, den;
};
template <class M> constexpr Raw curve(const M& m, u64 n, u64 sigma, u64 r2) {
 auto to= [&](u64 x) { return m.mul(x, r2); };
 const u64 U= m.norm(to((sigma * sigma - 5) % n)), V= m.norm(to(4 * sigma % n));
 const u64 U3= m.mul(m.mul(U, U), U), V3= m.mul(m.mul(V, V), V), vu= m.sub(V, U), w= m.add(m.add(U, U), m.add(U, V));
 return {U3, V3, m.mul(m.mul(m.mul(vu, vu), vu), w), m.mul(m.mul(U3, V), to(16))};
}
// 奇数の合成数 n (平方数でない) の自明でない因数。見つかるまで曲線を替える。DUAL なら 2 本ずつ回し、
// (A + 2)/4 = num / den の逆元は 2 本分の den の積の 1 回の inv_gcd でまとめ、gcd も 2 本の積でまとめて取る。
template <class M, bool DUAL> constexpr u64 ecm(u64 n) {
 const M m(n);
 const u64 r2= u64(((u128)((0 - n) % n) << 64) % n), one= m.mul(1, r2);
 // za zb の gcd。n ごと割れたときだけ 1 本ずつ見る。因数が無ければ 0、どれも n ごとなら n。
 auto check= [&](u64 za, u64 zb) -> u64 {
  const u64 g= gcd(m.norm(m.mul(za, zb)), n);
  if(g != n) return g == 1 ? 0 : g;
  for(u64 z: {za, zb})
   if(const u64 h= gcd(m.norm(z), n); h != 1 && h != n) return h;
  return n;
 };
 for(u64 sigma= 6;; sigma+= 1 + DUAL) {
  const Raw ca= curve(m, n, sigma, r2), cb= DUAL ? curve(m, n, sigma + 1, r2) : Raw{0, 0, 0, one};
  const auto [g, inv]= inv_gcd(from_mont(m, m.mul(ca.den, cb.den)), n);
  if(g != 1) {
   if(g != n) return g;
   continue;
  }
  const u64 it= m.mul(inv, r2), a24a= m.mul(ca.num, m.mul(it, cb.den));
  Pt Pa{ca.U3, ca.V3};
  if constexpr(DUAL) {
   const u64 a24b= m.mul(cb.num, m.mul(it, ca.den));
   Pt Pb{cb.U3, cb.V3};
   for(int i= 0; i < STAGE1.len; ++i) ladder2(m, Pa, Pb, STAGE1.c[i], a24a, a24b);
   if(const u64 g1= check(Pa.Z, Pb.Z); g1) {
    if(g1 != n) return g1;
    continue;
   }
   if(const u64 g2= check(stage2(m, Pa, a24a), stage2(m, Pb, a24b)); g2 && g2 != n) return g2;
  } else {
   for(int i= 0; i < STAGE1.len; ++i) Pa= ladder(m, Pa, STAGE1.c[i], a24a);
   if(const u64 g1= gcd(m.norm(Pa.Z), n); g1 != 1) {
    if(g1 != n) return g1;
    continue;
   }
   if(const u64 g2= gcd(m.norm(stage2(m, Pa, a24a)), n); g2 != 1 && g2 != n) return g2;
  }
 }
}
// Pollard rho (Brent)。奇数の合成数 n < 2^60。c の違う 2 本の列を同時に進め、128 段ごとに差の積の gcd を取る。
constexpr u64 rho(u64 n) {
 const Lazy m(n);
 constexpr u64 S= 128;
 for(u64 c1= 1, c2= 2;; c1+= 2, c2+= 2) {
  u64 z1= c1, z2= c2, y1= 0, y2= 0, x1= 0, x2= 0, g= 1;
  for(u64 k= S;; k<<= 1) {
   // z は (c, 2n + c)、取っておいた z は c より大きいので、x - z は (0, 4n) に収まる
   x1= z1 + 2 * n, x2= z2 + 2 * n;
   for(u64 j= 0; j < k && g == 1; j+= S) {
    y1= z1, y2= z2;
    u64 q1= 1, q2= 1;
    for(u64 i= 0; i < S; ++i) {
     z1= m.mul(z1, z1) + c1, z2= m.mul(z2, z2) + c2;
     q1= m.mul(q1, x1 - z1), q2= m.mul(q2, x2 - z2);
    }
    g= gcd(m.norm(m.mul(q1, q2)), n);
   }
   if(g == 1) continue;
   if(g != n) return g;
   // 1 段ずつやり直す
   for(int s= 0; s < 2; ++s) {
    u64 y= s ? y2 : y1, x= s ? x2 : x1, c= s ? c2 : c1;
    for(u64 i= 0; i < S; ++i) {
     y= m.mul(y, y) + c;
     if(u64 gg= gcd(m.norm(x - y), n); gg != 1) {
      if(gg != n) return gg;
      break;
     }
    }
   }
   break;  // どちらの列も n ごと割れたので c を替える
  }
 }
}
constexpr u64 isqrt(u64 n) {
 u64 r;
 if(std::is_constant_evaluated()) {
  r= n < 2 ? n : 1ull << ((64 - __builtin_clzll(n)) / 2 + 1);
  for(u64 s; (s= (r + n / r) / 2) < r;) r= s;
 } else r= std::sqrt((double)n);
 if(r > 0xffffffffull) r= 0xffffffffull;
 while(r * r > n) --r;
 while(r < 0xffffffffull && (r + 1) * (r + 1) <= n) ++r;
 return r;
}
// n は奇数で 97 以下の素因数を持たない。
constexpr void rec(u64 n, std::vector<u64>& out) {
 if(n == 1) return;
 if(n < 101 * 101 || is_prime(n)) return out.push_back(n);
 if(const u64 r= isqrt(n); r * r == n) return rec(r, out), rec(r, out);
 const u64 d= n < (1ull << 48) ? rho(n) : n < (1ull << 60) ? ecm<Lazy, true>(n) : ecm<Strict, false>(n);
 rec(d, out), rec(n / d, out);
}
struct SmallPrime {
 u64 p, inv, lim;
};
inline constexpr auto SMALL= [] {
 std::array<SmallPrime, 24> t{};
 int k= 0;
 for(u64 p= 3; p < 100; p+= 2) {
  bool pr= true;
  for(u64 d= 3; d * d <= p; d+= 2)
   if(p % d == 0) pr= false;
  if(pr) t[k++]= {p, inv64(p), ~0ull / p};
 }
 return t;
}();
}
// n の素因数を重複込みで昇順に並べて返す。n <= 1 なら空。
constexpr std::vector<unsigned long long> factorize(unsigned long long n) {
 using namespace factorize_internal;
 std::vector<u64> out;
 if(n <= 1) return out;
 const int z= __builtin_ctzll(n);
 out.insert(out.end(), z, 2);
 n>>= z;
 for(const auto& [p, inv, lim]: SMALL) {
  if(p * p > n) break;
  while(n * inv <= lim) out.push_back(p), n*= inv;
 }
 rec(n, out);
 std::sort(out.begin(), out.end());
 return out;
}
// n の素因数分解を (素数, 指数) の組で、素数の昇順に返す。n <= 1 なら空。
constexpr std::vector<std::pair<unsigned long long, int>> factorize_pairs(unsigned long long n) {
 std::vector<std::pair<unsigned long long, int>> out;
 for(auto p: factorize(n)) {
  if(out.empty() || out.back().first != p) out.emplace_back(p, 1);
  else ++out.back().second;
 }
 return out;
}
