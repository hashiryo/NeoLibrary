#pragma once
#include <cmath>
#include <numbers>
#include "neo/base/core/geo/Point.hpp"
// 角度の変換と回転、偏角の比較。AngleComp は整数の座標で exact に偏角を比べる (-pi, pi] の順。
inline long double radian_to_degree(long double r) { return r * 180.0 / std::numbers::pi_v<long double>; }
inline long double degree_to_radian(long double d) { return d * std::numbers::pi_v<long double> / 180.0; }
inline long double normalize_radian(long double r) {
 constexpr long double pi= std::numbers::pi_v<long double>;
 return r= std::fmod(r + pi, 2 * pi), r > 0 ? r - pi : r + pi;
}
template <class K> long double angle(const Point<K>& p) { return std::atan2((long double)p.y, (long double)p.x); }
template <class K> long double angle(const Point<K>& p, const Point<K>& q) { return std::atan2((long double)cross(p, q), (long double)dot(p, q)); }
template <class K> Affine<K> rotate(long double theta) {
 long double c= std::cos(theta), s= std::sin(theta);
 return {K(c), K(-s), K(s), K(c), Point<K>{0, 0}};
}
template <class K> Affine<K> rotate(const Point<K>& p, long double theta) {
 long double c= std::cos(theta), s= std::sin(theta);
 return {K(c), K(-s), K(s), K(c), Point<K>{K(p.x - c * p.x + s * p.y), K(p.y - s * p.x - c * p.y)}};
}
template <class K> Affine<K> rotate90(const Point<K>& p) { return {0, -1, 1, 0, p - !p}; }
template <class K> class AngleComp {
 using P= Point<K>;
 static int quad(const P& p) {
  if(int s= sgn(p.y); s) return s;
  return sgn(p.x) < 0 ? 2 : 0;
 }
public:
 bool operator()(const P& p, const P& q) const {
  if(int a= quad(p), b= quad(q); a != b) return a < b;
  return cross(p, q) > 0;
 }
};
