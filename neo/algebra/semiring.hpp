#pragma once
#include <cstddef>
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
 // 行列の積の速い道。p と q の値がすべて 0 以上 INF 以下なら、積は min(x + y, INF) と同じになる (INF に 0 以上を足しても INF
 // 以上)。INF は 2^k - 1 の形なので、範囲に入るかは値の or が INF の外の bit を持たないかで分かる。
 static constexpr bool fast_range(const MinPlus* p, const MinPlus* q, std::size_t n) {
  if constexpr(std::is_floating_point_v<T>) return false;
  else {
   using U= std::make_unsigned_t<T>;
   U o= 0;
   for(std::size_t i= 0; i < n; ++i) o|= U(p[i].x) | U(q[i].x);
   return !(o & ~U(INF));
  }
 }
 static constexpr MinPlus fast_mul(MinPlus l, MinPlus r) {
  const T s= l.x + r.x;
  return MinPlus(s < INF ? s : INF);
 }
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
 // 行列の積の速い道。p と q の値がすべて -INF 以上 0 以下なら、積は max(x + y, -INF) と同じになる。
 static constexpr bool fast_range(const MaxPlus* p, const MaxPlus* q, std::size_t n) {
  if constexpr(std::is_floating_point_v<T>) return false;
  else {
   using U= std::make_unsigned_t<T>;
   U o= 0;
   for(std::size_t i= 0; i < n; ++i) o|= (U(0) - U(p[i].x)) | (U(0) - U(q[i].x));
   return !(o & ~U(INF));
  }
 }
 static constexpr MaxPlus fast_mul(MaxPlus l, MaxPlus r) {
  const T s= l.x + r.x;
  return MaxPlus(s > -INF ? s : -INF);
 }
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
