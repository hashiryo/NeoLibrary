#pragma once
#include <array>
#include <bit>
#include <cmath>
#include <type_traits>
// Baillie-PSW: 底 2 の強擬素数判定と extra strong Lucas 判定 (Q = 1、P は 3, 4, 5, ... で (P^2-4 / n) = -1 となる最初のもの)。
// 2^64 未満に反例が無いことは、Feitsma の底 2 の擬素数の一覧で確かめられている (Math::Prime::Util 0.31 の Changes)。
// 2 つの判定を 1 つのループで同時に回す。Lucas 側は V 列だけを梯子で求め、1 段は独立な積 2 回になる。
namespace is_prime_internal {
using u64= unsigned long long;
using u32= unsigned;
using u128= unsigned __int128;
using i64= long long;
template <u64 p> constexpr u32 divisible(u64 n) {
 constexpr u64 inv= [] {
  u64 x= p;
  for(int i= 0; i < 6; ++i) x*= 2 - p * x;
  return x;
 }();
 return n * inv <= ~0ull / p;
}
constexpr u64 inv64(u64 n) {
 u64 x= (3 * n) ^ 2, y= 1 - n * x;
 x*= 1 + y, y*= y;
 x*= 1 + y, y*= y;
 x*= 1 + y, y*= y;
 return x * (1 + y);
}
constexpr u32 SMALL_P[8]= {2, 3, 5, 7, 11, 13, 17, 19};
// P^2-4 の平方因子を除いた部分に SMALL_P[i] が現れるかの bit。P <= 20 では素因数がすべて 19 以下。
constexpr u32 d_mask(u32 P) {
 u32 D= P * P - 4, m= 0;
 for(int i= 0; i < 8; ++i) {
  int e= 0;
  while(D % SMALL_P[i] == 0) D/= SMALL_P[i], ++e;
  if(e & 1) m|= 1u << i;
 }
 return m;
}
// (p / n) = -1 となる p の組から、(P^2-4 / n) = -1 となる最初の P (3..20)。無ければ 0。
inline constexpr auto FIRST_P= [] {
 std::array<unsigned char, 256> t{};
 for(u32 j= 0; j < 256; ++j)
  for(u32 P= 3; P <= 20; ++P)
   if(std::popcount(d_mask(P) & j) & 1) {
    t[j]= P;
    break;
   }
 return t;
}();
// (p / n) = -1 なら 1。n は奇数で p で割り切れない。
template <u32 p> constexpr u32 jneg(u64 n) {
 constexpr u64 qr= [] {
  u64 m= 0;
  for(u32 x= 1; x < p; ++x) m|= 1ull << (x * x % p);
  return m;
 }();
 return u32(~qr >> (n % p) & 1) ^ u32((p & 3) == 3 && (n & 3) == 3);
}
constexpr int jacobi(u64 a, u64 n) {
 int t= 1;
 for(a%= n; a;) {
  int z= std::countr_zero(a);
  a>>= z;
  if((z & 1) && ((n & 7) == 3 || (n & 7) == 5)) t= -t;
  if((a & n & 3) == 3) t= -t;
  u64 r= n % a;
  n= a, a= r;
 }
 return n == 1 ? t : 0;
}
constexpr bool is_square(u64 n) {
 u64 r;
 if(std::is_constant_evaluated()) {
  r= n < 2 ? n : 1ull << ((64 - std::countl_zero(n)) / 2 + 1);
  for(u64 s; (s= (r + n / r) / 2) < r;) r= s;
 } else r= std::sqrt((double)n);
 if(r > 0xffffffffull) r= 0xffffffffull;
 while(r * r > n) --r;
 while(r < 0xffffffffull && (r + 1) * (r + 1) <= n) ++r;
 return r * r == n;
}
// Lucas 判定の P。n が平方数か、P^2-4 と共通の素因数を持てば 0 (合成数)。
constexpr u64 select_P(u64 n) {
 u32 j= u32(0x28 >> (n & 7) & 1) | jneg<3>(n) << 1 | jneg<5>(n) << 2 | jneg<7>(n) << 3 | jneg<11>(n) << 4 | jneg<13>(n) << 5 | jneg<17>(n) << 6 | jneg<19>(n) << 7;
 if(u64 P= FIRST_P[j]; P) return P;
 if(is_square(n)) return 0;
 for(u64 P= 21;; ++P) {
  int jj= jacobi(P * P - 4, n);
  if(jj == -1) return P;
  if(jj == 0 && (P * P - 4) % n) return 0;
 }
}
// 指数の上 T 桁は、2^e と V_k (k <= 2^T) を小さい整数のまま求めてから Montgomery 表現に直す。
constexpr int T= 3;
inline constexpr auto LUCAS_V= [] {
 std::array<std::array<u64, (1 << T) + 1>, 21> t{};
 for(u64 P= 3; P <= 20; ++P) {
  t[P][0]= 2, t[P][1]= P;
  for(int k= 2; k <= (1 << T); ++k) t[P][k]= P * t[P][k - 1] - t[P][k - 2];
 }
 return t;
}();
// b R mod n を [0, 2n) で返す (R = 2^64、one = R mod n、b < 2^40)。商 b one / n を浮動小数点で見積もる (±1 しかずれない)。
constexpr u64 to_mont(u64 b, u64 one, u64 n, double inv_n) {
 u64 q= u64(i64(double(i64(one)) * double(i64(b)) * inv_n)), r= b * one - q * n;
 return i64(r) < 0 ? r + n : r;
}
// n < 2^60。値を 3n 未満のまま持ち、a b < 9 n^2 < n R を保つ。
constexpr bool bpsw_fast(u64 n, u64 P) {
 const u64 ninv= inv64(n), one= (0 - n) % n, mone= n - one, two= 2 * one >= n ? 2 * one - n : 2 * one;
 // a b R^{-1} + c mod n を (0, 3n) で返す。nc = n + c (0 <= c <= n)。
 auto mul= [n, ninv](u64 a, u64 b, u64 nc) {
  u128 t= (u128)a * b;
  return u64(t >> 64) + nc - u64((u128)(u64(t) * ninv) * n >> 64);
 };
 auto norm= [n](u64 v) {
  v= v >= n ? v - n : v;
  return v >= n ? v - n : v;
 };
 const double inv_n= 1.0 / double(i64(n));
 u64 Pm= to_mont(P, one, n, inv_n);
 Pm= Pm >= n ? Pm - n : Pm;
 const u64 nc2= 2 * n - two, ncP= 2 * n - Pm;
 // 底 2 は n-1 = dA 2^sA の dA、Lucas は n+1 = dL 2^sL の dL を、先頭に 0 を詰めて同じ桁数に揃えて上から読む。
 // 0 のビットでは、底 2 の x = 1 も Lucas の (V_0, V_1) = (2, P) も変わらない。
 const int sA= std::countr_zero(n - 1), sL= std::countr_zero(n + 1);
 const u64 dA= (n - 1) >> sA, dL= (n + 1) >> sL;
 const int L= 64 - std::countl_zero(dA > dL ? dA : dL), top= L > T ? L - T : 0;
 const u64 eA= dA >> top, eL= dL >> top;
 u64 pa= top ? dA << (64 - top) : 0, pl= top ? dL << (64 - top) : 0;
 // 底 2 の段は、片方の因数を 2 倍してから掛けて x^2 2^bit を積 1 回で求める (x < 2n、2x x < 8 n^2)。
 // Lucas は (V_k, V_{k+1}) = prev ? (a, s) : (s, a) で持ち、V_{2k+1} = V_k V_{k+1} - P と、V_k か V_{k+1} の 2 乗 - 2 を求める。
 u64 x= to_mont(1ull << eA, one, n, inv_n), a, s, prev= 1;
 auto stepL= [&](u64 b) {
  const u64 sel= b == prev ? s : a, a2= mul(a, s, ncP);
  s= mul(sel, sel, nc2), a= a2, prev= b;
 };
 if(P <= 20) a= to_mont(LUCAS_V[P][eL], one, n, inv_n), s= to_mont(LUCAS_V[P][eL + 1], one, n, inv_n);
 else {
  a= two, s= Pm;
  for(int i= T; i--;) stepL(eL >> i & 1);
 }
 for(int i= 0; i < top; ++i, pa<<= 1, pl<<= 1) x= mul(x, x << (pa >> 63), n), stepL(pl >> 63);
 // 底 2: x ≡ ±1 か、2 乗を sA-1 回までして -1 が出れば通る。
 // Lucas: U_d ≡ 0 かつ V_d ≡ ±2 (Q = 1 では V_d ≡ ±2 のとき V_{d+1} ≡ ±P と同じ) か、0 <= r < sL-1 で V_{d 2^r} ≡ 0 なら通る。
 x= norm(x);
 u64 V= norm(prev ? a : s), W= norm(prev ? s : a);
 bool okA= x == one || x == mone, okL= (V == two && W == Pm) || (V == n - two && W == n - Pm) || (sL >= 2 && V == 0);
 for(int ra= sA - 1, rl= sL - 2;;) {
  if(!okA && (ra <= 0 || x == one)) return false;
  if(!okL && (rl <= 0 || V == two)) return false;
  if(okA && okL) return true;
  if(!okA) x= norm(mul(x, x, n)), okA= x == mone, --ra;
  if(!okL) V= norm(mul(V, V, nc2)), okL= V == 0, --rl;
 }
}
// n >= 2^60。値を n 未満に保つ。
constexpr bool bpsw_strict(u64 n, u64 P) {
 const u64 ninv= inv64(n), one= (0 - n) % n, mone= n - one;
 auto mul= [n, ninv](u64 a, u64 b) {
  u128 t= (u128)a * b;
  u64 m= (u128)(u64(t) * ninv) * n >> 64, hi= t >> 64;
  return hi >= m ? hi - m : hi - m + n;
 };
 auto add= [n](u64 a, u64 b) { return a >= n - b ? a - (n - b) : a + b; };
 auto sub= [n](u64 a, u64 b) { return a >= b ? a - b : a + (n - b); };
 const u64 two= add(one, one), Pm= (u128)P * one % n;
 const int sA= std::countr_zero(n - 1), sL= std::countr_zero(n + 1);
 const u64 dA= (n - 1) >> sA, dL= (n + 1) >> sL;
 const int L= 64 - std::countl_zero(dA > dL ? dA : dL);
 u64 pa= dA << (64 - L), pl= dL << (64 - L), x= one, a= two, s= Pm, prev= 1;
 for(int i= 0; i < L; ++i, pa<<= 1, pl<<= 1) {
  x= mul(x, x), x= add(x, x & (0 - (pa >> 63)));
  const u64 b= pl >> 63, sel= b == prev ? s : a, a2= sub(mul(a, s), Pm);
  s= sub(mul(sel, sel), two), a= a2, prev= b;
 }
 u64 V= prev ? a : s, W= prev ? s : a;
 bool okA= x == one || x == mone, okL= (V == two && W == Pm) || (V == n - two && W == n - Pm) || (sL >= 2 && V == 0);
 for(int ra= sA - 1, rl= sL - 2;;) {
  if(!okA && (ra <= 0 || x == one)) return false;
  if(!okL && (rl <= 0 || V == two)) return false;
  if(okA && okL) return true;
  if(!okA) x= mul(x, x), okA= x == mone, --ra;
  if(!okL) V= sub(mul(V, V), two), okL= V == 0, --rl;
 }
}
}
constexpr bool is_prime(unsigned long long n) {
 using namespace is_prime_internal;
 if(n < 64) return 0x28208a20a08a28acull >> n & 1;
 if(u32(~n & 1) | divisible<3>(n) | divisible<5>(n) | divisible<7>(n) | divisible<11>(n) | divisible<13>(n) | divisible<17>(n) | divisible<19>(n) | divisible<23>(n) | divisible<29>(n) | divisible<31>(n) | divisible<37>(n)) return false;
 if(n < 41 * 41) return true;
 const u64 P= select_P(n);
 if(!P) return false;
 return n < (1ull << 60) ? bpsw_fast(n, P) : bpsw_strict(n, P);
}
