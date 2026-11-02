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
// 整数の MinPlus と MaxPlus は、行列の積の速い道を持つ。積を INF の比べのない和 x + y (fast_mul) にして和 (min か max) で足し、
// 最後に成分ごとに fast_fix で INF に直す。掛ける 2 つの行列の値がどれも fast_one_sided の範囲 (MinPlus は 0 以上 INF 以下、
// MaxPlus は -INF 以上 0 以下) か、どれも fast_bounded の範囲 (±INF か、絶対値が INF / 4 以下) にあれば、厳密な積と同じ値になる。
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
 // 行列の積の速い道。fast_one_sided は p の値がすべて 0 以上 INF 以下か。このとき積は min(x + y, INF) と同じで (INF に 0 以上を
 // 足しても INF 以上)、min は単調なので、和 x + y の min を取ってから INF 以上を INF にしてよい。INF は 2^k - 1 の形なので、
 // 範囲に入るかは値の or が INF の外の bit を持たないかで分かる。
 static constexpr bool fast_one_sided(const MinPlus* p, std::size_t n) { return fast_one_sided(p, p, n); }
 // p と q の n 個ずつを 1 本のループで調べる (大きさを型で決める行列の積で、判定を 1 回の or の畳み込みにする)。
 static constexpr bool fast_one_sided(const MinPlus* p, const MinPlus* q, std::size_t n) {
  if constexpr(std::is_floating_point_v<T>) return false;
  else {
   using U= std::make_unsigned_t<T>;
   U o= 0;
   for(std::size_t i= 0; i < n; ++i) o|= U(p[i].x) | U(q[i].x);
   return !(o & ~U(INF));
  }
 }
 // fast_bounded は p の値がすべて INF か、絶対値が L = INF / 4 以下か。このとき有限の値どうしの和は絶対値が INF / 2 未満で、
 // INF を含む和は INF - L 以上 (2 INF 以下なので溢れない) なので、和の min を取ってから INF / 2 以上を INF にすれば、厳密な積の
 // min と同じになる。
 static constexpr bool fast_bounded(const MinPlus* p, std::size_t n) {
  if constexpr(std::is_floating_point_v<T> || !std::is_signed_v<T>) return false;
  else {
   using U= std::make_unsigned_t<T>;
   constexpr U L= U(INF) >> 2;
   bool ok= true;
   for(std::size_t i= 0; i < n; ++i) ok&= (U(p[i].x) == U(INF)) | (U(U(p[i].x) + L) <= 2 * L);
   return ok;
  }
 }
 static constexpr MinPlus fast_mul(MinPlus l, MinPlus r) { return MinPlus(T(l.x + r.x)); }
 // 速い道で足した和を直す。bounded は fast_bounded の範囲で掛けたか。
 static constexpr MinPlus fast_fix(MinPlus v, bool bounded) {
  const T t= bounded ? T(INF >> 1) : INF;
  return MinPlus(v.x >= t ? INF : v.x);
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
 // 行列の積の速い道。MinPlus の符号を逆にしたもので、fast_one_sided は p の値がすべて -INF 以上 0 以下か、fast_bounded は
 // すべて -INF か、絶対値が INF / 4 以下か。
 static constexpr bool fast_one_sided(const MaxPlus* p, std::size_t n) { return fast_one_sided(p, p, n); }
 static constexpr bool fast_one_sided(const MaxPlus* p, const MaxPlus* q, std::size_t n) {
  if constexpr(std::is_floating_point_v<T> || !std::is_signed_v<T>) return false;
  else {
   using U= std::make_unsigned_t<T>;
   U o= 0;
   for(std::size_t i= 0; i < n; ++i) o|= (U(0) - U(p[i].x)) | (U(0) - U(q[i].x));
   return !(o & ~U(INF));
  }
 }
 static constexpr bool fast_bounded(const MaxPlus* p, std::size_t n) {
  if constexpr(std::is_floating_point_v<T> || !std::is_signed_v<T>) return false;
  else {
   using U= std::make_unsigned_t<T>;
   constexpr U L= U(INF) >> 2;
   bool ok= true;
   for(std::size_t i= 0; i < n; ++i) ok&= (U(p[i].x) == U(-INF)) | (U(U(p[i].x) + L) <= 2 * L);
   return ok;
  }
 }
 static constexpr MaxPlus fast_mul(MaxPlus l, MaxPlus r) { return MaxPlus(T(l.x + r.x)); }
 static constexpr MaxPlus fast_fix(MaxPlus v, bool bounded) {
  const T t= bounded ? T(INF >> 1) : INF;
  return MaxPlus(v.x <= -t ? -INF : v.x);
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
