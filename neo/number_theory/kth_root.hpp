#pragma once
#include <bit>
#include <cassert>
#include <cmath>
#include <vector>
// ⌊a^{1/k}⌋ (0 ≤ a < 2^64、1 ≤ k ≤ 64)。浮動小数点で見積もった r を、整数の r^k と (r + 1)^k で 1 だけ直す。
// 直しは r - [r^k > a] + [(r + 1)^k ≤ a] の算術で、分岐にしない (a が r^k のすぐそばだと、double に丸めた a からの見積もりが
// 1 大きく出るかが半々に変わり、分岐が外れる)。
// - k = 2 は sqrt、k = 3 は x = 2^e m (e = 3q + s) の (2^s m)^{1/3} を、m の上位 7 bit の表と 3 次の多項式で出す。
// - k ≥ 4 は log2 と exp2 を表と 2 次の多項式で出す (答えは 65536 未満なので相対誤差 2^-20 ほどで足りる)。
// - k ≥ 5 は答えが 7132 以下なので、r^k - 1 (r^k ≥ 2^64 なら 2^64 - 1) の表を引く。r^k - 1 で持つのは、2^64 - 1 で頭打ちに
//   すると a = 2^64 - 1 のとき (r + 1)^k ≤ a を取り違えるため。k = 3、4 は掛け算で作る。
// 表は初めて呼んだときに作る (1 万項ほど)。procon-judge の yosupo-kth-root-integer の 6 回目の kr6 の形。
namespace kth_root_internal {
using u32= unsigned;
using u64= unsigned long long;
using i64= long long;
using u128= unsigned __int128;
struct Tables {
 double lg[128], inv[128], ex[64], invk[65], cb[3][128];
 u32 off[65];
 std::vector<u64> pw;
 Tables() {
  for(int j= 0; j < 128; ++j) {
   const double c= 1.0 + (j + 0.5) / 128;
   lg[j]= std::log2(c), inv[j]= 1.0 / c;
   for(int s= 0; s < 3; ++s) cb[s][j]= std::cbrt(std::ldexp(c, s));
  }
  for(int j= 0; j < 64; ++j) ex[j]= std::exp2(j / 64.0);
  for(int k= 1; k <= 64; ++k) invk[k]= 1.0 / k;
  u32 o= 0;
  for(int k= 5; k <= 64; ++k) {
   off[k]= o;
   for(u64 r= 0;; ++r) {
    u128 p= 1;
    bool sat= false;
    for(int i= 0; i < k && !sat; ++i) p*= r, sat= (p >> 64) != 0;
    pw.push_back(sat ? ~0ull : u64(p) - 1), ++o;
    if(sat) break;
   }
   pw.push_back(~0ull), ++o;
  }
 }
};
inline const Tables& tables() {
 static const Tables t;
 return t;
}
// a ≥ 2^63 なら 1 bit ずらしてから符号付きで変換して 2 倍する (相対誤差 2^-53)
inline double to_double(u64 a) {
 const u64 s= a >> 63;
 return double(i64(a >> s)) * double(1 + s);
}
// x^{1/3} (1 ≤ x < 2^64)。誤差 2^-36 ほど
inline double cbrt3(double x, const Tables& T) {
 const u64 b= std::bit_cast<u64>(x);
 const u32 e= u32(b >> 52) - 1023, q= e / 3, s= e - 3 * q, j= u32(b >> 45) & 127;
 const double m= std::bit_cast<double>((b & 0x000fffffffffffffull) | 0x3ff0000000000000ull);
 const double r= m * T.inv[j] - 1.0;
 const double v= T.cb[s][j] * (1.0 + r * (1.0 / 3 + r * (-1.0 / 9 + r * (5.0 / 81))));
 return std::bit_cast<double>(std::bit_cast<u64>(v) + (u64(q) << 52));
}
// x^{1/k} (1 ≤ x < 2^64、k ≥ 4)。相対誤差 2^-20 ほど
inline double root_lo(double x, u32 k, const Tables& T) {
 const u64 b= std::bit_cast<u64>(x);
 const int e= int(b >> 52) - 1023;
 const u32 j= u32(b >> 45) & 127;
 const double m= std::bit_cast<double>((b & 0x000fffffffffffffull) | 0x3ff0000000000000ull);
 const double r= m * T.inv[j] - 1.0;
 const double y= (double(e) + T.lg[j] + r * (1.4426950408889634 + r * -0.7213475204444817)) * T.invk[k];
 const i64 n= i64(y * 64.0);
 const double f= (y - double(n) * (1.0 / 64)) * 0.6931471805599453;
 const double v= T.ex[n & 63] * (1.0 + f * (1.0 + f * 0.5));
 return std::bit_cast<double>(std::bit_cast<u64>(v) + (u64(n >> 6) << 52));
}
}
inline unsigned long long kth_root(unsigned long long a, unsigned k) {
 using namespace kth_root_internal;
 assert(1 <= k && k <= 64);
 if(k == 1 || a <= 1) return a;
 const Tables& T= tables();
 const double x= to_double(a);
 if(k == 2) {
  u64 r= u64(i64(std::sqrt(x)));
  r= r < 0xffffffffull ? r : 0xfffffffeull;  // (r + 1)^2 を u64 に収める
  const u64 s= r + 1;
  return r - u64(r * r > a) + u64(s * s <= a);
 }
 if(k == 3) {
  const u64 r= u64(i64(cbrt3(x, T)));
  const u128 p= u128(r * r) * r, s= r + 1, q= u128(s * s) * s;
  return r - u64(p > a) + u64(q <= a);
 }
 const u64 r= u64(i64(root_lo(x, k, T)));
 if(k == 4) {
  const u64 r2= r * r, s= r + 1, s2= s * s;
  const u128 p= u128(r2) * r2, q= u128(s2) * s2;
  return r - u64(p > a) + u64(q <= a);
 }
 const u64* p= T.pw.data() + T.off[k];
 return r - u64(p[r] >= a) + u64(p[r + 1] < a);
}
