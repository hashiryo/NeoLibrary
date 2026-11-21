#pragma once
#include <cassert>
#include <type_traits>
#include <utility>
// Stein の方法の拡張版を Kaliski の almost inverse と同じ形で書く。値を 2^c で割る代わりにもう一方の係数を 2^c 倍し、
// 最後に 2^{-k} を Montgomery の還元で掛けて戻す。値は neo の gcd と同じく
// (大きいほう, 小さいほう) → (差 >> c, 小さいほう) の形で持つので、鎖は引き算、tzcnt、シフトの 3 段のまま。
//
// B を奇数として、2 つの枠 (値 a, 係数 ca) と (値 b, 係数 cb) は、符号 σ_a = -σ_b と k について
//   A ca ≡ σ_a a 2^k,  A cb ≡ σ_b b 2^k  (mod B),  B = a cb + b ca
// を保つ。最後の式の項がすべて非負なので係数はずっと B 以下に収まり、u64 全体の入力を 128 bit の値なしで扱える。
// 終わりは a = b = g で、B / g = ca + cb、(A / g) cb ≡ σ_b 2^k (mod B / g) となる。
// 桁数の差が大きいときは、A %= B か互除法の 1 段 (B = q A + r、ca = q) を先に済ませる。どちらも上の式を保つ。
namespace inv_gcd_internal {
using u32= unsigned;
using u64= unsigned long long;
using u128= unsigned __int128;
constexpr int GAP= 8;
constexpr u64 inv64(u64 n) {
 u64 x= (3 * n) ^ 2, y= 1 - n * x;
 x*= 1 + y, y*= y;
 x*= 1 + y, y*= y;
 x*= 1 + y, y*= y;
 return x * (1 + y);
}
// v 2^{-c} mod M を [0, M] で返す。0 <= c <= 64、v <= M、M は奇数、Minv = M^{-1} mod 2^64。
constexpr u64 redc(u64 v, int c, u64 M, u64 Minv) {
 u128 t= (u128)v << (64 - c);
 u64 hi= u64(t >> 64), mh= u64((u128)(u64(t) * Minv) * M >> 64);
 return hi >= mh ? hi - mh : hi - mh + M;
}
// 割り算は幅に合わせる (32 bit の割り算のほうが速い CPU がある)。
template <int W> constexpr u64 quo(u64 a, u64 b) {
 if constexpr(W == 32) return u32(a) / u32(b);
 else return a / b;
}
struct Res {
 u64 g, x, M;  // A x ≡ g (mod B)、0 <= x < M = B / g
};
// B は奇数、Binv = B^{-1} mod 2^64。W は入力の幅 (32 か 64) で、k < 2W なので W = 32 なら還元は 1 回で済む。
// gcc は大きさを見て展開しないことがあり、そのときは結果をメモリ経由で返して 1 回 3 ns ほど遅くなるので、展開を強いる。
template <int W> [[gnu::always_inline]] constexpr Res odd(u64 A, u64 B, u64 Binv) {
 if(A > B && __builtin_clzll(B) - __builtin_clzll(A) > GAP) A-= quo<W>(A, B) * B;
 if(A == 0) return {B, 0, 1};
 int k= __builtin_ctzll(A);
 u64 a= B, b= A >> k, ca= 0, cb= 1, neg= 0;
 if(__builtin_clzll(b) - __builtin_clzll(a) > GAP) {
  u64 q= quo<W>(a, b), r= a - q * b;
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
  neg+= a < b;  // 入れ替えの回数。足し算で書くと adc 1 つになる
  ca+= cb, cb= cm << c, k+= c;
  a= t >> c, b= mn;
 }
 const u64 M= ca + cb, Minv= Binv * b;
 u64 x;
 if constexpr(W == 32) x= redc(cb, k, M, Minv);
 else x= redc(redc(cb, k >> 1, M, Minv), k - (k >> 1), M, Minv);
 x= x >= M ? x - M : x;
 return {b, (neg & 1) && x ? M - x : x, M};
}
// odd を呼ぶところを 1 か所にする。3 か所に分けると、gcc は odd をインライン展開せず、結果をメモリ経由で返していた。
template <int W> constexpr std::pair<u64, u64> inv_gcd(u64 a, u64 b) {
 // 共通の 2 の冪を除く (b が奇数なら z = 0)。a = 0 なら a1 = 0 で、odd が (b1, 0) を返す。
 const int z= __builtin_ctzll(a | b);
 const u64 a1= a >> z, b1= b >> z;
 // b1 が偶数なら a1 は奇数なので、役を入れ替えて b1 y ≡ g1 (mod a1) を解き、a1 x + b1 y = g1 から x を出す。
 const bool sw= !(b1 & 1);
 const u64 B= sw ? a1 : b1, Binv= inv64(B);
 const Res r= odd<W>(sw ? b1 : a1, B, Binv);
 if(!sw) return {r.g << z, r.x};
 // t = (b1 y - g1) / a1 は割り切れて -1 <= t < b1 / g1 なので、mod 2^64 で a1^{-1} を掛けて求まる。
 const u64 t= (b1 * r.x - r.g) * Binv, m= b1 * (r.M * Binv);  // m = b1 / g1 (g1^{-1} = (a1 / g1) a1^{-1})
 const u64 x= m - t;
 return {r.g << z, x >= m ? x - m : x};
}
}
// g = gcd(a, b) と、a x ≡ g (mod b)、0 <= x < b / g の x を返す。b >= 1。
template <class T> constexpr std::pair<T, T> inv_gcd(T a, T b) {
 static_assert(std::is_integral_v<T> && sizeof(T) <= 8);
 assert(b >= 1);
 using U= std::make_unsigned_t<T>;
 if constexpr(std::is_signed_v<T>) {
  if(a < 0 && (a%= b) < 0) a+= b;
 }
 auto [g, x]= inv_gcd_internal::inv_gcd<sizeof(T) <= 4 ? 32 : 64>(U(a), U(b));
 return {T(g), T(x)};
}
// a x ≡ 1 (mod m)、0 <= x < m の x を返す。m >= 1 で、gcd(a, m) = 1 でなければ assert で止める。
template <class T> constexpr T inv_mod(T a, T m) {
 auto [g, x]= inv_gcd(a, m);
 assert(g == 1);
 return x;
}
