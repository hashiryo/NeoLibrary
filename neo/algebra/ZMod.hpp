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
// それ以外 (偶数の定数と ZMod<0>) は Shoup の方法で、w' = floor(w 2^32 / M) を持つ。
// 方式は procon-judge の self/modulo-* の問題で選んだ。記録は algo-notes の notes/modular_arithmetic/modint_problems.md。
namespace zmod_internal {
using u32= unsigned;
using u64= unsigned long long;
using i64= long long;
using u128= unsigned __int128;
using i128= __int128;
// n は奇数。n^{-1} mod 2^64 (x = n は 3 bit 正しく、1 回で正しい桁が倍になる)。
constexpr u64 inv64(u64 n) {
 u64 x= n;
 for(int i= 0; i < 5; ++i) x*= 2 - n * x;
 return x;
}
// ZMod<0> の法。u32 の static 変数で持つと、u32 の配列への書き込みのたびに GCC が法を読み直すので、u64 で持つ。
struct Dyn {
 static inline u64 m= 1, x= ~0ull, k= 2;  // k = 2m - (2^64 mod m)
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
 // Plantard の還元。2^64 / phi 未満の w に -w R^-2 mod M を M 未満で返す (R = 2^32、法は奇数の定数)。
 static constexpr u32 plantard(u64 w) { return u32((u128((w * RHO) | u32(-1)) * MOD) >> 64); }
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
 // 拡張 Euclid の互除法。gcd(val(), M) = 1 でなければ assert で止める。
 constexpr ZMod inv() const {
  u32 a= val(), b= mod();
  int x= 1, y= 0;
  while(b) {
   u32 q= a / b;
   a-= q * b, x-= int(q) * y;
   std::swap(a, b), std::swap(x, y);
  }
  assert(a == 1);
  return raw(u32(x < 0 ? x + int(mod()) : x));
 }
 // 片方が決まった積の前計算。p.fixed() を一度作り、x * P や x *= P で何度も掛ける。Shoup の前計算には割り算が要るので、
 // ループの外で作る。ZMod<0> では、作ったあとに set_mod で法を変えると使えなくなる。
 struct Fixed {
  u32 w;  // M 未満
  u64 b;  // Plantard なら (-w R^2 mod M) rho、Shoup なら floor(w 2^32 / M)
 };
 constexpr Fixed fixed() const {
  const u32 w= val();
  if constexpr(PLANTARD) {
   constexpr u32 R4= u32(-u128(MOD) % MOD);  // R^4 mod M。plantard(w R^4) = -w R^2 mod M
   u64 b= u64(plantard(u64(w) * R4)) * RHO;
   // 法が定数だと、GCC は b を w と定数の積と見て a b を (a w) rho に組み替え、掛け算を 3 回に戻すので、b を隠す。
   if(!std::is_constant_evaluated()) asm("" : "+r"(b));
   return {w, b};
  } else return {w, (u64(w) << 32) / mod()};
 }
 friend constexpr ZMod operator*(ZMod a, Fixed f) {
  if constexpr(PLANTARD) return raw(u32((u128((u64(a.v) * f.b) | u32(-1)) * MOD) >> 64));
  else return raw(a.v * f.w - u32((u64(a.v) * f.b) >> 32) * mod());
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
