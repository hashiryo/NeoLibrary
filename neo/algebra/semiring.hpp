#pragma once
#include <limits>
#include <type_traits>
// 半環の値の型。和と積と 0 と 1 (静的メンバ関数 zero() と one()) を持ち、行列 (neo/linalg/Matrix.hpp、SquareMatrix.hpp) の
// 要素やセグメント木の値にそのまま使える。
//   MinPlus<T>: 和は min、積は + (最短路)。0 は INF、1 は 0。
//   MaxPlus<T>: 和は max、積は + (最長路)。0 は -INF、1 は 0。
//   MaxMin<T>: 和は max、積は min (道の幅の最大)。0 は T の最小、1 は T の最大。
//   OrAnd: 和は or、積は and (到達)。
// 整数の INF は T の最大の半分で、積は片方が INF (MaxPlus は -INF) ならそのまま INF にし、和が INF を超えれば INF に丸める。
// 負の重みがあっても INF は INF のまま。浮動小数点数の INF は infinity。
namespace semiring_internal {
template <class T> constexpr T inf() {
 if constexpr(std::is_floating_point_v<T>) return std::numeric_limits<T>::infinity();
 else return std::numeric_limits<T>::max() / 2;
}
}
template <class T> struct MinPlus {
 static constexpr T INF= semiring_internal::inf<T>();
 T x;
 constexpr MinPlus(): x(INF) {}
 constexpr MinPlus(T x): x(x) {}
 static constexpr MinPlus zero() { return MinPlus(INF); }
 static constexpr MinPlus one() { return MinPlus(T(0)); }
 constexpr MinPlus& operator+=(MinPlus r) { return x= r.x < x ? r.x : x, *this; }
 constexpr MinPlus& operator*=(MinPlus r) {
  if constexpr(std::is_floating_point_v<T>) x+= r.x;
  else {
   // max(x, r.x, x + r.x) >= INF なら INF。GCC は条件式の選択を分岐にするので、和は符号なしで求め、マスクで選ぶ。
   using U= std::make_unsigned_t<T>;
   const T s= T(U(x) + U(r.x)), t= x > r.x ? x : r.x;
   const U m= U(0) - U((t > s ? t : s) >= INF);
   x= T((U(s) & ~m) | (U(INF) & m));
  }
  return *this;
 }
 friend constexpr MinPlus operator+(MinPlus l, MinPlus r) { return l+= r; }
 friend constexpr MinPlus operator*(MinPlus l, MinPlus r) { return l*= r; }
 friend constexpr bool operator==(MinPlus l, MinPlus r) { return l.x == r.x; }
 friend constexpr bool operator!=(MinPlus l, MinPlus r) { return l.x != r.x; }
};
template <class T> struct MaxPlus {
 static constexpr T INF= semiring_internal::inf<T>();
 T x;
 constexpr MaxPlus(): x(-INF) {}
 constexpr MaxPlus(T x): x(x) {}
 static constexpr MaxPlus zero() { return MaxPlus(-INF); }
 static constexpr MaxPlus one() { return MaxPlus(T(0)); }
 constexpr MaxPlus& operator+=(MaxPlus r) { return x= r.x > x ? r.x : x, *this; }
 constexpr MaxPlus& operator*=(MaxPlus r) {
  if constexpr(std::is_floating_point_v<T>) x+= r.x;
  else {
   using U= std::make_unsigned_t<T>;
   const T s= T(U(x) + U(r.x)), t= x < r.x ? x : r.x;
   const U m= U(0) - U((t < s ? t : s) <= -INF);
   x= T((U(s) & ~m) | (U(-INF) & m));
  }
  return *this;
 }
 friend constexpr MaxPlus operator+(MaxPlus l, MaxPlus r) { return l+= r; }
 friend constexpr MaxPlus operator*(MaxPlus l, MaxPlus r) { return l*= r; }
 friend constexpr bool operator==(MaxPlus l, MaxPlus r) { return l.x == r.x; }
 friend constexpr bool operator!=(MaxPlus l, MaxPlus r) { return l.x != r.x; }
};
template <class T> struct MaxMin {
 T x;
 constexpr MaxMin(): x(std::numeric_limits<T>::lowest()) {}
 constexpr MaxMin(T x): x(x) {}
 static constexpr MaxMin zero() { return MaxMin(std::numeric_limits<T>::lowest()); }
 static constexpr MaxMin one() { return MaxMin(std::numeric_limits<T>::max()); }
 constexpr MaxMin& operator+=(MaxMin r) { return x= r.x > x ? r.x : x, *this; }
 constexpr MaxMin& operator*=(MaxMin r) { return x= r.x < x ? r.x : x, *this; }
 friend constexpr MaxMin operator+(MaxMin l, MaxMin r) { return l+= r; }
 friend constexpr MaxMin operator*(MaxMin l, MaxMin r) { return l*= r; }
 friend constexpr bool operator==(MaxMin l, MaxMin r) { return l.x == r.x; }
 friend constexpr bool operator!=(MaxMin l, MaxMin r) { return l.x != r.x; }
};
struct OrAnd {
 bool x;
 constexpr OrAnd(): x(false) {}
 constexpr OrAnd(bool x): x(x) {}
 static constexpr OrAnd zero() { return OrAnd(false); }
 static constexpr OrAnd one() { return OrAnd(true); }
 constexpr OrAnd& operator+=(OrAnd r) { return x|= r.x, *this; }
 constexpr OrAnd& operator*=(OrAnd r) { return x&= r.x, *this; }
 friend constexpr OrAnd operator+(OrAnd l, OrAnd r) { return l+= r; }
 friend constexpr OrAnd operator*(OrAnd l, OrAnd r) { return l*= r; }
 friend constexpr bool operator==(OrAnd l, OrAnd r) { return l.x == r.x; }
 friend constexpr bool operator!=(OrAnd l, OrAnd r) { return l.x != r.x; }
};
