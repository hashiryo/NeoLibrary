#pragma once
#include <array>
#include <cassert>
#include <utility>
// 素数 p を法として X^2 ≡ y となる X を 1 つ求める。p の剰余で道を分ける。
//
// p ≡ 3 (mod 4): X = y^{(p+1)/4} とし、X^2 = y を確かめる (違えば平方非剰余)。
// p ≡ 5 (mod 8): Atkin の方法。v = (2y)^{(p-5)/8}、i = 2y v^2 (i^2 = -1)、X = y v (i - 1) とし、X^2 = y を確かめる。
// p ≡ 1 (mod 8): p - 1 = 2^s Q (Q は奇数)。
//  s ≤ 8 なら Tonelli-Shanks。平方非剰余 z は、相互法則で (q / p) = (p mod q / q) となるので、64 未満の奇素数 q ごとの
//  平方剰余の表から累乗を使わずに探す。t = y^Q、R = y^{(Q+1)/2}、c = z^Q は 1 つの輪で求める。
//  s > 8 なら Lucas 列。Jacobi 記号で平方非剰余を除いてから、a = 2^j で a^2 - y が平方非剰余になる j を探す (候補は
//  2 つずつ並べて求める)。α = 1 + ω / a (ω^2 = a^2 - y) は和が 2、積が y / a^2 の F_{p^2} の元で、α^{(p+1)/2} は F_p に
//  入り 2 乗すると y / a^2 になるので、X = a V_{(p+1)/2} / 2 (V は P' = 2、Q' = y / a^2 の Lucas 列)。V は
//  (V_k, V_{k+1}, Q'^k, Q'^{k+1}) の梯子で求め、1 段の掛け算は 4 回で、どれも前の段から 1 回ぶんの深さしかない。選択は
//  ビットの分岐にする (同じ p が続くと当たる)。a = 2^j (j < 8) で見つからなければ、a を 1 から順に探して F_{p^2} の
//  累乗を取る (Cipolla の方法)。
// Jacobi 記号は二進の形で、差 a - n の ctz を |a - n| と min(a, n) の選択と並べて取る (x86 では選択を cmov で書く。gcc は
// 三項演算子を分岐にする)。
// 剰余の掛け算は Montgomery で、p < 2^30 は R = 2^32 で値を 2p 未満のまま持ち、p < 2^62 は R = 2^64 で値を 2p 未満の
// まま、それ以上は p 未満に保つ。累乗は右から左の 2 進で、掛ける値はマスクの算術で選ぶ (分岐なし)。
namespace sqrt_mod_internal {
using u32= unsigned;
using u64= unsigned long long;
using u128= unsigned __int128;
struct Mont32 {
 using T= u32;
 using S= int;
 static constexpr bool LAZY= true;
 u32 p, pinv, r1;  // pinv = p^{-1} mod 2^32、r1 = 2^32 mod p
 explicit Mont32(u64 p_): p((u32)p_) {
  u32 x= (3 * p) ^ 2;  // 下位 5 ビットが正しい
  for(int i= 0; i < 3; ++i) x*= 2 - p * x;
  pinv= x, r1= (0u - p) % p;
 }
 u32 mul(u32 a, u32 b) const {
  const u64 t= (u64)a * b;
  return (u32)(t >> 32) + p - (u32)((u64)((u32)t * pinv) * p >> 32);
 }
 // a b R^{-1} - c (c < 2p) を [0, 2p) に。上位の項に p - c を先に足し、和は (-2p, 2p) に入るので負なら 2p を足す
 u32 mulsub(u32 a, u32 b, u32 c) const {
  const u64 t= (u64)a * b;
  const u32 r= (u32)(t >> 32) + (p - c) - (u32)((u64)((u32)t * pinv) * p >> 32);
  return r + (2 * p & (u32)((int)r >> 31));
 }
 u32 to(u64 a) const { return (u32)((a << 32) % p); }
};
template <bool L> struct Mont64 {
 using T= u64;
 using S= long long;
 static constexpr bool LAZY= L;
 u64 p, pinv, r1, r2;  // pinv = p^{-1} mod 2^64、r1 = 2^64 mod p、r2 = 2^128 mod p
 explicit Mont64(u64 p_): p(p_) {
  u64 x= (3 * p) ^ 2;
  for(int i= 0; i < 4; ++i) x*= 2 - p * x;
  pinv= x, r1= (u64)(((u128)1 << 64) % p), r2= (u64)((u128)r1 * r1 % p);
 }
 u64 mul(u64 a, u64 b) const {
  const u128 t= (u128)a * b;
  const u64 hi= (u64)(t >> 64), m= (u64)((u128)((u64)t * pinv) * p >> 64);
  if constexpr(L) return hi + p - m;
  else return hi >= m ? hi - m : hi - m + p;
 }
 u64 mulsub(u64 a, u64 b, u64 c) const {
  if constexpr(L) {
   const u128 t= (u128)a * b;
   const u64 r= (u64)(t >> 64) + (p - c) - (u64)((u128)((u64)t * pinv) * p >> 64);
   return r + (2 * p & (u64)((long long)r >> 63));
  } else {
   const u64 x= mul(a, b);
   return x - c + (p & (0ull - (u64)(x < c)));
  }
 }
 u64 to(u64 a) const { return mul(a, r2); }
};
template <class M> struct Ops: M {
 using T= typename M::T;
 using M::M, M::mul, M::p, M::r1;
 static constexpr T W= M::LAZY ? 2 : 1;  // 値は [0, W p)
 T norm(T a) const {
  if constexpr(M::LAZY) return a >= p ? a - p : a;
  else return a;
 }
 u64 from(T a) const { return norm(mul(a, 1)); }
 bool eq(T a, T b) const { return norm(a) == norm(b); }
 // 和と差と半分。補正はマスクで足す
 T add(T a, T b) const {
  const T t= a + b - W * p;
  return t + (W * p & (T)((typename M::S)t >> (8 * sizeof(T) - 1)));
 }
 T sub(T a, T b) const { return a - b + (W * p & (T(0) - (T)(a < b))); }
 T half(T a) const { return (a + (p & (T(0) - (a & 1)))) >> 1; }
 T sel(T b, u64 bit) const { return r1 ^ ((r1 ^ b) & (T(0) - (T)(bit & 1))); }
 T pow(T b, u64 e) const {
  T r= r1;
  for(; e; e>>= 1, b= mul(b, b)) r= mul(r, sel(b, e));
  return r;
 }
 // (b^e1, b^e2, c^e1)、e2 ≤ e1
 std::array<T, 3> pow3(T b, T c, u64 e1, u64 e2) const {
  T x= r1, y= r1, w= r1;
  for(; e1; e1>>= 1, e2>>= 1, b= mul(b, b), c= mul(c, c)) x= mul(x, sel(b, e1)), y= mul(y, sel(b, e2)), w= mul(w, sel(c, e1));
  return {x, y, w};
 }
};
// 二進の Jacobi 記号の 1 段。a と n は奇数で a ≠ n。sg のビット 0 が符号。
// (a / n) は、a > n なら (a - n / n)、a < n なら相互法則で (n - a / a) に符号 (a ≡ n ≡ 3 (mod 4) なら反転) を掛けたもの。
// |a - n| の 2 の冪を外すと、z が奇数で新しい n が ≡ 3, 5 (mod 8) のとき反転する。
inline void jacobi_step(u64& a, u64& n, u64& sg) {
 u64 d= a - n;
 const u64 z= __builtin_ctzll(d);
 sg^= (u64)(a < n) & (a & n) >> 1;
 u64 nn= n;
#if defined(__x86_64__)
 const u64 e= n - a;
 asm("cmp %[n], %[a]\n\tcmovb %[e], %[d]\n\tcmovb %[a], %[nn]" : [d] "+&r"(d), [nn] "+&r"(nn) : [a] "r"(a), [n] "r"(n), [e] "r"(e) : "cc");
#else
 const u64 m= 0ull - (u64)(a < n);
 nn^= (a ^ n) & m, d= (d ^ m) - m;
#endif
 sg^= z & ((nn >> 1) ^ (nn >> 2));
 a= d >> z, n= nn;
}
inline int jacobi_sign(u64 n, u64 sg) { return n == 1 ? 1 - 2 * (int)(sg & 1) : 0; }
// a の 2 の冪を外す (a は 0 でない)
inline u64 jacobi_init(u64& a, u64 n) {
 const u64 z= __builtin_ctzll(a);
 a>>= z;
 return z & ((n >> 1) ^ (n >> 2));
}
// Jacobi 記号 (a / n)、n は奇数、a < n
inline int jacobi(u64 a, u64 n) {
 if(!a) return n == 1;
 u64 sg= jacobi_init(a, n);
 while(a != n) jacobi_step(a, n, sg);
 return jacobi_sign(n, sg);
}
// (a1 / n) と (a2 / n) を 2 本並べて求める
inline std::pair<int, int> jacobi2(u64 a1, u64 a2, u64 n) {
 if(!a1 || !a2) return {jacobi(a1, n), jacobi(a2, n)};
 u64 n1= n, n2= n, s1= jacobi_init(a1, n), s2= jacobi_init(a2, n);
 while(a1 != n1 && a2 != n2) jacobi_step(a1, n1, s1), jacobi_step(a2, n2, s2);
 while(a1 != n1) jacobi_step(a1, n1, s1);
 while(a2 != n2) jacobi_step(a2, n2, s2);
 return {jacobi_sign(n1, s1), jacobi_sign(n2, s2)};
}
// 64 未満の奇素数 q と、q の平方剰余の集合のビット (bit r が立てば r は平方剰余)
struct QrTable {
 u32 q;
 u64 qr;
};
inline constexpr auto QRT= [] {
 std::array<QrTable, 17> t{};
 int k= 0;
 for(u32 q= 3; q < 64; q+= 2) {
  bool pr= true;
  for(u32 d= 3; d * d <= q; d+= 2) pr&= q % d != 0;
  if(!pr) continue;
  u64 m= 0;
  for(u32 x= 1; x < q; ++x) m|= 1ull << (x * x % q);
  t[k++]= {q, m};
 }
 return t;
}();
// P' = 2、Q' = q の Lucas 列の V_e (e は 3 以上の奇数)
template <class O> typename O::T lucas_v(const O& M, typename O::T q, u64 e) {
 using T= typename O::T;
 const T two= M.add(M.r1, M.r1);
 T v0= two, v1= M.sub(M.add(two, two), M.add(q, q)), q0= q, q1= M.mul(q, q);  // k = 1
 for(int i= 62 - __builtin_clzll(e); i > 0; --i) {
  const T mx= M.mulsub(v0, v1, M.add(q0, q0)), qmx= M.mul(q0, q1);
  if(e >> i & 1) v1= M.mulsub(v1, v1, M.add(q1, q1)), q1= M.mul(q1, q1), v0= mx, q0= qmx;
  else v0= M.mulsub(v0, v0, M.add(q0, q0)), q0= M.mul(q0, q0), v1= mx, q1= qmx;
 }
 return M.mulsub(v0, v1, M.add(q0, q0));  // e の最下位ビットは 1
}
// a を 1 から順に探す Cipolla の方法 (y は平方剰余)
template <class O> u64 cipolla(const O& M, typename O::T y, u64 p) {
 using T= typename O::T;
 T a= M.r1, w;
 for(;; a= M.add(a, M.r1))
  if(w= M.sub(M.mul(a, a), y); jacobi(M.from(w), p) == -1) break;
 T x0= M.r1, x1= 0, b0= a, b1= M.r1;  // 結果と底 (x0 + x1 ω、ω^2 = w)
 for(u64 e= (p + 1) / 2; e; e>>= 1) {
  const T m= T(0) - (T)(e & 1), c0= M.r1 ^ ((M.r1 ^ b0) & m), c1= b1 & m;
  const T n0= M.add(M.mul(x0, c0), M.mul(M.mul(x1, c1), w)), n1= M.add(M.mul(x0, c1), M.mul(x1, c0));
  x0= n0, x1= n1;
  const T s0= M.add(M.mul(b0, b0), M.mul(M.mul(b1, b1), w)), s1= M.mul(M.add(b0, b0), b1);
  b0= s0, b1= s1;
 }
 return M.from(x0);
}
// 0 < Y < P、P は奇素数。無ければ ~0
template <class O> u64 solve(u64 Y, u64 P) {
 using T= typename O::T;
 const O M(P);
 const T y= M.to(Y);
 if((P & 3) == 3) {
  const T x= M.pow(y, (P + 1) / 4);
  return M.eq(M.mul(x, x), y) ? M.from(x) : ~0ull;
 }
 if((P & 7) == 5) {
  const T y2= M.add(y, y), v= M.pow(y2, (P - 5) / 8);
  const T i= M.mul(y2, M.mul(v, v));
  const T x= M.mul(M.mul(y, v), M.sub(i, M.r1));
  return M.eq(M.mul(x, x), y) ? M.from(x) : ~0ull;
 }
 const int s= __builtin_ctzll(P - 1);
 const u64 Q= (P - 1) >> s;
 if(s <= 8) {
  u64 z= 0;
  for(const auto& [q, qr]: QRT)
   if(const u64 r= P % q; r && !(qr >> r & 1)) {
    z= q;
    break;
   }
  if(!z)
   for(z= 67; jacobi(z, P) != -1;) ++z;
  auto [t, R, c]= M.pow3(y, M.to(z), Q, (Q + 1) / 2);
  {  // t^{2^{s-1}} ≠ 1 なら平方非剰余
   T u= t;
   for(int i= 0; i < s - 1; ++i) u= M.mul(u, u);
   if(!M.eq(u, M.r1)) return ~0ull;
  }
  for(int m= s; !M.eq(t, M.r1);) {
   int i= 0;
   T u= t;
   while(!M.eq(u, M.r1)) u= M.mul(u, u), ++i;
   T b= c;
   for(int j= 0; j < m - i - 1; ++j) b= M.mul(b, b);
   m= i, c= M.mul(b, b), t= M.mul(t, c), R= M.mul(R, b);
  }
  return M.from(R);
 }
 if(jacobi(Y, P) != 1) return ~0ull;
 auto cand= [&](int j) {  // 4^j - Y mod P
  const u64 c= (1ull << (2 * j)) % P;
  return c >= Y ? c - Y : c + P - Y;
 };
 int j= -1;
 for(int k= 0; j < 0 && k < 8; k+= 2) {
  const auto [e1, e2]= jacobi2(cand(k), cand(k + 1), P);
  j= e1 == -1 ? k : e2 == -1 ? k + 1 : -1;
 }
 if(j < 0) return cipolla(M, y, P);
 T q= y;
 for(int i= 0; i < 2 * j; ++i) q= M.half(q);  // Y / a^2
 const T v= lucas_v(M, q, (P + 1) / 2);
 T x= j == 0 ? M.half(v) : v;  // a V / 2 = 2^{j-1} V
 for(int i= 1; i < j; ++i) x= M.add(x, x);
 return M.from(x);
}
}
// 素数 p を法として X^2 ≡ y となる 0 ≤ X < p を 1 つ返す (2 つの根のどちらを返すかは決めない)。無ければ -1。
// 2 ≤ p < 2^63 で、p が素数であることは確かめない。
inline long long sqrt_mod(unsigned long long y, unsigned long long p) {
 using namespace sqrt_mod_internal;
 assert(2 <= p && p < (1ull << 63));
 if(y >= p) y%= p;
 if(y == 0 || p == 2) return (long long)y;
 u64 x;
 if(p < (1u << 30)) x= solve<Ops<Mont32>>(y, p);
 else if(p < (1ull << 62)) x= solve<Ops<Mont64<true>>>(y, p);
 else x= solve<Ops<Mont64<false>>>(y, p);
 return x == ~0ull ? -1 : (long long)x;
}
