#pragma once
#include <type_traits>
// Stein の方法 (binary GCD)。大きさの差が大きいときだけ、先に大きいほうを小さいほうで 1 回割る。
// Stein は大きい数と小さい数の組で 1 段に 2 bit ほどしか縮まないので、割り算 1 回のほうが速い。
namespace gcd_internal {
template <class U> constexpr int ctz(U x) {
 if constexpr(sizeof(U) == 8) return __builtin_ctzll(x);
 else return __builtin_ctz(x);
}
template <class U> constexpr int clz(U x) {
 if constexpr(sizeof(U) == 8) return __builtin_clzll(x);
 else return __builtin_clz(x);
}
template <class U> constexpr U gcd(U a, U b) {
 U mx= a > b ? a : b, mn= a > b ? b : a;
 if(mn == 0) return mx;
 if(clz(mn) - clz(mx) > 8) {
  mx%= mn;
  if(mx == 0) return mn;
 }
 const int z= ctz(mx | mn);
 mx>>= ctz(mx), mn>>= ctz(mn);
 // 差の末尾の 0 を大小の入れ替えより先に数えるので、鎖は引き算、tzcnt、シフトの 3 段になる。
 while(mx != mn) {
  U d= mx - mn;
  int s= ctz(d);
  bool f= mx > mn;
  mn= f ? mn : mx, mx= (f ? d : -d) >> s;
 }
 return mx << z;
}
}
// |a| と |b| の最大公約数。gcd(0, 0) = 0。
template <class T> constexpr T gcd(T a, T b) {
 static_assert(std::is_integral_v<T> && sizeof(T) <= 8);
 using U= std::conditional_t<sizeof(T) <= 4, unsigned, unsigned long long>;
 if constexpr(std::is_signed_v<T>) return T(gcd_internal::gcd<U>(a < 0 ? U(0) - U(a) : U(a), b < 0 ? U(0) - U(b) : U(b)));
 else return T(gcd_internal::gcd<U>(a, b));
}
