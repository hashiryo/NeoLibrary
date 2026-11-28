#pragma once
#include <cmath>
#include <istream>
#include <ostream>
#include <type_traits>
// 2 次元の点。K は整数でも浮動小数でもよい。整数なら dot と cross を 1 段広い型 (int -> long long -> __int128)
// で計算するので、判定は exact になる。浮動小数なら sgn を EPS で判定する。
// 比較 (operator< と operator==) は座標の比較そのもので、EPS を入れない。並べ替えの順序が推移的であることと、
// 同じ座標の点を重複として扱うことを優先する。位置の判定は sgn(cross(...)) のように sgn を通して書く。
namespace geo_internal {
template <class K> struct wide {
 using type= K;
};
template <> struct wide<char> {
 using type= short;
};
template <> struct wide<unsigned char> {
 using type= unsigned short;
};
template <> struct wide<short> {
 using type= int;
};
template <> struct wide<unsigned short> {
 using type= unsigned;
};
template <> struct wide<int> {
 using type= long long;
};
template <> struct wide<unsigned> {
 using type= unsigned long long;
};
template <> struct wide<long long> {
 using type= __int128;
};
template <> struct wide<unsigned long long> {
 using type= unsigned __int128;
};
template <> struct wide<float> {
 using type= double;
};
template <> struct wide<double> {
 using type= long double;
};
}
template <class K> using wide_t= typename geo_internal::wide<K>::type;
// x の符号。浮動小数は絶対値が EPS 以下なら 0 とみなす。
template <class K> constexpr int sgn(K x) {
 if constexpr(std::is_floating_point_v<K>) {
  constexpr K eps= K(1e-9);
  return x < -eps ? -1 : x > eps;
 } else return x < 0 ? -1 : x > 0;
}
template <class K> struct Point {
 K x, y;
 constexpr Point(K x= K(), K y= K()): x(x), y(y) {}
 constexpr Point& operator+=(const Point& p) { return x+= p.x, y+= p.y, *this; }
 constexpr Point& operator-=(const Point& p) { return x-= p.x, y-= p.y, *this; }
 constexpr Point& operator*=(K a) { return x*= a, y*= a, *this; }
 constexpr Point& operator/=(K a) { return x/= a, y/= a, *this; }
 constexpr Point operator+(const Point& p) const { return {x + p.x, y + p.y}; }
 constexpr Point operator-(const Point& p) const { return {x - p.x, y - p.y}; }
 constexpr Point operator*(K a) const { return {x * a, y * a}; }
 constexpr Point operator/(K a) const { return {x / a, y / a}; }
 friend constexpr Point operator*(K a, const Point& p) { return {a * p.x, a * p.y}; }
 constexpr Point operator-() const { return {-x, -y}; }
 // 90 度回す。
 constexpr Point operator!() const { return {-y, x}; }
 constexpr bool operator<(const Point& p) const { return x != p.x ? x < p.x : y < p.y; }
 constexpr bool operator>(const Point& p) const { return p < *this; }
 constexpr bool operator<=(const Point& p) const { return !(p < *this); }
 constexpr bool operator>=(const Point& p) const { return !(*this < p); }
 constexpr bool operator==(const Point& p) const { return x == p.x && y == p.y; }
 constexpr bool operator!=(const Point& p) const { return x != p.x || y != p.y; }
 friend std::istream& operator>>(std::istream& is, Point& p) { return is >> p.x >> p.y; }
 friend std::ostream& operator<<(std::ostream& os, const Point& p) { return os << "(" << p.x << ", " << p.y << ")"; }
};
template <class K> constexpr wide_t<K> dot(const Point<K>& p, const Point<K>& q) {
 return wide_t<K>(p.x) * q.x + wide_t<K>(p.y) * q.y;
}
// 左に曲がる (反時計回り) なら正、右なら負。
template <class K> constexpr wide_t<K> cross(const Point<K>& p, const Point<K>& q) {
 return wide_t<K>(p.x) * q.y - wide_t<K>(p.y) * q.x;
}
template <class K> constexpr wide_t<K> norm2(const Point<K>& p) { return dot(p, p); }
template <class K> constexpr wide_t<K> dist2(const Point<K>& p, const Point<K>& q) { return norm2(p - q); }
template <class K> long double dist(const Point<K>& p, const Point<K>& q) { return std::sqrt((long double)dist2(p, q)); }
// アフィン変換。p -> A p + b。
template <class K> struct Affine {
 K a00= 1, a01= 0, a10= 0, a11= 1;
 Point<K> b;
 constexpr Point<K> operator()(const Point<K>& p) const { return {a00 * p.x + a01 * p.y + b.x, a10 * p.x + a11 * p.y + b.y}; }
 constexpr Affine operator*(const Affine& r) const {
  return {a00 * r.a00 + a01 * r.a10, a00 * r.a01 + a01 * r.a11, a10 * r.a00 + a11 * r.a10, a10 * r.a01 + a11 * r.a11, (*this)(r.b)};
 }
};
enum class Ccw : int { counter_clockwise= 0, clockwise= 1, online_back= 2, online_front= 3, on_segment= 4 };
constexpr const char* ccw_name(Ccw c) {
 switch(c) {
  case Ccw::counter_clockwise: return "COUNTER_CLOCKWISE";
  case Ccw::clockwise: return "CLOCKWISE";
  case Ccw::online_back: return "ONLINE_BACK";
  case Ccw::online_front: return "ONLINE_FRONT";
  default: return "ON_SEGMENT";
 }
}
inline std::ostream& operator<<(std::ostream& os, Ccw c) { return os << ccw_name(c); }
// p0 -> p1 -> p2 の向き。浮動小数は外積を正規化してから sgn を通す。
template <class K> Ccw ccw(const Point<K>& p0, const Point<K>& p1, const Point<K>& p2) {
 Point<K> a= p1 - p0, b= p2 - p0;
 int s;
 if constexpr(std::is_floating_point_v<K>) s= sgn(sgn(cross(a, b) / std::sqrt((long double)norm2(a) * (long double)norm2(b))));
 else s= sgn(cross(a, b));
 if(s) return s > 0 ? Ccw::counter_clockwise : Ccw::clockwise;
 auto d= dot(a, b);
 if(sgn(d) < 0) return Ccw::online_back;
 return sgn(d - norm2(a)) > 0 ? Ccw::online_front : Ccw::on_segment;
}
