#pragma once
#include <cassert>
#include <iostream>
#include <type_traits>
#include <utility>
// Z/MZ の元。法 M は 1 <= M < 2^30 で、偶奇も素数かどうかも問わない。ZMod<0> は、法を実行時に set_mod で決める型になる。
// 値はそのままの値を 2M 未満で持ち (Barrett の冗長な表現)、掛け算は Barrett で還元する。前計算 X = floor((2^64 - 1) / M) なら、
// 2M 未満どうしの積 (4M^2 < 2^62) の剰余は 2M 未満に収まり、M = 1 でも X が u64 に入る。
// 片方が決まった積は fixed() で前計算してから掛ける。法が奇数の定数なら Plantard の前計算を、そのままの値に掛ける形で持つ
// (w を Plantard の表現 -w R^2 mod M に直してから rho = M^{-1} mod 2^64 を掛けておくと、還元の結果がそのままの値になる)。
// 偶数の定数なら Shoup の方法で、w' = floor(w 2^32 / M) を持つ。ZMod<0> は fixed() を作るときに法の偶奇を見て、奇数なら Plantard、
// 偶数なら Shoup にし、掛け算のたびに作ったときの印で分岐する (Shoup は掛け算 3 回で、Plantard の 2 回より throughput が落ちる CPU がある)。
// 方式は procon-judge の self/modulo-* の問題で選んだ。記録は algo-notes の notes/modular_arithmetic/modint_problems.md。
namespace zmod_internal {
using u32= unsigned;
using u64= unsigned long long;
using i64= long long;
using u128= unsigned __int128;
using i128= __int128;
// inv64、redc、odd、inv_gcd32 は neo/number_theory/inv_gcd.hpp の 32 bit の道の写し (W = 32 に決めた形)。直すときは両方を直す。
// 数論のヘッダを読まないのは、数論のヘッダを直すたびに ZMod を使う提出が全部測り直しになるのを避けるため。
// n は奇数。n^{-1} mod 2^64。
constexpr u64 inv64(u64 n) {
 u64 x= (3 * n) ^ 2, y= 1 - n * x;
 x*= 1 + y, y*= y;
 x*= 1 + y, y*= y;
 x*= 1 + y, y*= y;
 return x * (1 + y);
}
constexpr int GAP= 8;
// v 2^{-c} mod M を [0, M] で返す。0 <= c <= 64、v <= M、M は奇数、Minv = M^{-1} mod 2^64。
constexpr u64 redc(u64 v, int c, u64 M, u64 Minv) {
 u128 t= (u128)v << (64 - c);
 u64 hi= u64(t >> 64), mh= u64((u128)(u64(t) * Minv) * M >> 64);
 return hi >= mh ? hi - mh : hi - mh + M;
}
struct Res {
 u64 g, x, M;  // A x ≡ g (mod B)、0 <= x < M = B / g
};
// B は奇数、Binv = B^{-1} mod 2^64。A, B < 2^32 で、k < 64 なので還元は 1 回で済む。
[[gnu::always_inline]] constexpr Res odd(u64 A, u64 B, u64 Binv) {
 if(A > B && __builtin_clzll(B) - __builtin_clzll(A) > GAP) A-= u32(A) / u32(B) * B;
 if(A == 0) return {B, 0, 1};
 int k= __builtin_ctzll(A);
 u64 a= B, b= A >> k, ca= 0, cb= 1, neg= 0;
 if(__builtin_clzll(b) - __builtin_clzll(a) > GAP) {
  u64 q= u32(a) / u32(b), r= a - q * b;
  if(r == 0) a= b, ca= q - 1;
  else {
   int c= __builtin_ctzll(r);
   a= r >> c, ca= q, cb<<= c, k+= c;
  }
 }
 while(a != b) {
  u64 d= a - b, e= b - a;
  int c= __builtin_ctzll(d);
  u64 mn= a < b ? a : b, t= a > b ? d : e, cm= a > b ? cb : ca;
  neg+= a < b;
  ca+= cb, cb= cm << c, k+= c;
  a= t >> c, b= mn;
 }
 const u64 M= ca + cb, Minv= Binv * b;
 u64 x= redc(cb, k, M, Minv);
 x= x >= M ? x - M : x;
 return {b, (neg & 1) && x ? M - x : x, M};
}
// g = gcd(a, b) と、a x ≡ g (mod b)、0 <= x < b / g の x を返す。a, b < 2^32、b >= 1。
constexpr std::pair<u64, u64> inv_gcd32(u64 a, u64 b) {
 const int z= __builtin_ctzll(a | b);
 const u64 a1= a >> z, b1= b >> z;
 const bool sw= !(b1 & 1);
 const u64 B= sw ? a1 : b1, Binv= inv64(B);
 const Res r= odd(sw ? b1 : a1, B, Binv);
 if(!sw) return {r.g << z, r.x};
 const u64 t= (b1 * r.x - r.g) * Binv, m= b1 * (r.M * Binv);
 const u64 x= m - t;
 return {r.g << z, x >= m ? x - m : x};
}
// ZMod<0> の法。u32 の static 変数で持つと、u32 の配列への書き込みのたびに GCC が法を読み直すので、u64 で持つ。
struct Dyn {
 static inline u64 m= 1, x= ~0ull, k= 2;  // k = 2m - (2^64 mod m)
 static inline u64 rho= 1, r4= 0;  // m が奇数のとき m^{-1} mod 2^64 と 2^128 mod m (Plantard の前計算)
};
}
template <unsigned MOD> class ZMod {
 using u32= zmod_internal::u32;
 using u64= zmod_internal::u64;
 using i64= zmod_internal::i64;
 using u128= zmod_internal::u128;
 using i128= zmod_internal::i128;
 using Dyn= zmod_internal::Dyn;
 static_assert(MOD < (1u << 30), "法は 1 以上 2^30 未満 (0 なら実行時の法)");
 static constexpr bool DYN= MOD == 0;
 static constexpr bool PLANTARD= !DYN && (MOD & 1);
 static constexpr u64 RHO= zmod_internal::inv64(MOD | 1);  // 法が奇数の定数のときだけ使う
 u32 v;  // 2M 未満
 static constexpr u64 X() {
  if constexpr(DYN) return Dyn::x;
  else return u64(-1) / MOD;
 }
 static constexpr u32 M2() { return mod() << 1; }
 // 2^62 未満の n の剰余を 2M 未満で返す。
 static constexpr u32 rem62(u64 n) { return u32(n - u64((u128(n) * X()) >> 64) * mod()); }
 // u64 全体の n では q の誤差が 2 まであって n - q M は 3M 未満なので、2M 以上なら 2M を引く。
 static constexpr u32 rem64(u64 n) {
  u32 r= rem62(n);
  return r >= M2() ? r - M2() : r;
 }
 template <class T> static constexpr u32 conv(T n) {
  if constexpr(sizeof(T) > 8) {  // 128 bit は遅い割り算のままにする
   if constexpr(std::is_signed_v<T>) {
    i128 r= i128(n) % mod();
    return u32(r < 0 ? r + mod() : r);
   } else return u32(u128(n) % mod());
  } else if constexpr(!DYN) {  // コンパイラが定数の割り算を掛け算とシフトに直す
   if constexpr(std::is_signed_v<T>) {
    i64 r= i64(n) % i64(MOD);
    return u32(r < 0 ? r + MOD : r);
   } else return u32(u64(n) % MOD);
  } else if constexpr(sizeof(T) <= 4) {
   if constexpr(std::is_signed_v<T>) {
    u32 r= rem62(n < 0 ? u64(-i64(n)) : u64(n));
    return n < 0 && r ? M2() - r : r;
   } else return rem62(u64(n));
  } else if constexpr(std::is_signed_v<T>) {  // u64(n) = n + 2^64 (n < 0) の 2^64 の分を k で戻す
   u32 r= rem64(u64(n)) + (u32(i64(n) >> 63) & u32(Dyn::k));
   return r >= M2() ? r - M2() : r;
  } else return rem64(u64(n));
 }
 // Plantard の還元。2^64 / phi 未満の w に -w R^-2 mod M を M 未満で返す (R = 2^32、法は奇数)。
 static constexpr u32 plantard(u64 w, u64 rho, u64 m) { return u32((u128((w * rho) | u32(-1)) * m) >> 64); }
public:
 static constexpr u32 mod() {
  if constexpr(DYN) return u32(Dyn::m);
  else return MOD;
 }
 // ZMod<0> の法を決める。1 <= m < 2^30。それまでに作った値と fixed() の前計算は使えなくなる。
 static void set_mod(u32 m) {
  static_assert(DYN, "set_mod は ZMod<0> だけ");
  assert(1 <= m && m < (1u << 30));
  Dyn::m= m, Dyn::x= u64(-1) / m, Dyn::k= 2ull * m - (u64(-1) % m + 1) % m;
  Dyn::rho= m & 1 ? zmod_internal::inv64(m) : 0, Dyn::r4= m & 1 ? u64(-u128(m) % m) : 0;
 }
 constexpr ZMod(): v(0) {}
 template <class T, std::enable_if_t<std::is_integral_v<T>, int> = 0> constexpr ZMod(T n): v(conv(n)) {}
 // 2M 未満の x をそのまま値にする (還元しない)。
 static constexpr ZMod raw(u32 x) {
  ZMod r;
  r.v= x;
  return r;
 }
 constexpr u32 val() const { return v >= mod() ? v - mod() : v; }
 constexpr ZMod& operator+=(ZMod r) { return v+= r.v, v= v >= M2() ? v - M2() : v, *this; }
 constexpr ZMod& operator-=(ZMod r) { return v-= r.v, v= int(v) < 0 ? v + M2() : v, *this; }
 constexpr ZMod& operator*=(ZMod r) { return v= rem62(u64(v) * r.v), *this; }
 constexpr ZMod& operator/=(ZMod r) { return *this*= r.inv(); }
 constexpr ZMod operator-() const { return ZMod() - *this; }
 friend constexpr ZMod operator+(ZMod l, ZMod r) { return l+= r; }
 friend constexpr ZMod operator-(ZMod l, ZMod r) { return l-= r; }
 friend constexpr ZMod operator*(ZMod l, ZMod r) { return l*= r; }
 friend constexpr ZMod operator/(ZMod l, ZMod r) { return l/= r; }
 friend constexpr bool operator==(ZMod l, ZMod r) { return l.val() == r.val(); }
 friend constexpr bool operator!=(ZMod l, ZMod r) { return l.val() != r.val(); }
 constexpr ZMod pow(u64 e) const {
  ZMod r= raw(1), b= *this;
  for(; e; e>>= 1, b*= b)
   if(e & 1) r*= b;
  return r;
 }
 // 割り算を使わない拡張 gcd (neo/number_theory/inv_gcd.hpp の写し)。gcd(val(), M) = 1 でなければ assert で止める。
 constexpr ZMod inv() const {
  const auto [g, x]= zmod_internal::inv_gcd32(val(), mod());
  assert(g == 1);
  return raw(u32(x));
 }
 // 片方が決まった積の前計算。p.fixed() を一度作り、x * P や x *= P で何度も掛ける。Shoup の前計算には割り算が要るので、
 // ループの外で作る。ZMod<0> では、作ったあとに set_mod で法を変えると使えなくなる。
 struct Fixed {
  u32 w;  // M 未満
  u64 b;  // Plantard なら (-w R^2 mod M) rho、Shoup なら floor(w 2^32 / M)
  bool pl;  // Plantard で前計算したか (ZMod<0> でだけ見る)
 };
 constexpr Fixed fixed() const {
  const u32 w= val();
  if constexpr(PLANTARD) {
   constexpr u32 R4= u32(-u128(MOD) % MOD);  // R^4 mod M。plantard(w R^4) = -w R^2 mod M
   u64 b= u64(plantard(u64(w) * R4, RHO, MOD)) * RHO;
   // 法が定数だと、GCC は b を w と定数の積と見て a b を (a w) rho に組み替え、掛け算を 3 回に戻すので、b を隠す。
   if(!std::is_constant_evaluated()) asm("" : "+r"(b));
   return {w, b, true};
  } else if constexpr(DYN) {
   if(Dyn::m & 1) {
    u64 b= u64(plantard(u64(w) * Dyn::r4, Dyn::rho, Dyn::m)) * Dyn::rho;
    if(!std::is_constant_evaluated()) asm("" : "+r"(b));
    return {w, b, true};
   }
  }
  return {w, (u64(w) << 32) / mod(), false};
 }
 friend constexpr ZMod operator*(ZMod a, Fixed f) {
  if constexpr(PLANTARD) return raw(u32((u128((u64(a.v) * f.b) | u32(-1)) * MOD) >> 64));
  else if constexpr(DYN) {
   if(f.pl) return raw(u32((u128((u64(a.v) * f.b) | u32(-1)) * Dyn::m) >> 64));
  }
  return raw(a.v * f.w - u32((u64(a.v) * f.b) >> 32) * mod());
 }
 friend constexpr ZMod operator*(Fixed f, ZMod a) { return a * f; }
 friend constexpr ZMod& operator*=(ZMod& a, Fixed f) { return a= a * f; }
 friend std::ostream& operator<<(std::ostream& os, ZMod r) { return os << r.val(); }
 friend std::istream& operator>>(std::istream& is, ZMod& r) {
  i64 n;
  is >> n, r= ZMod(n);
  return is;
 }
};
