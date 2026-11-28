#pragma once
#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>
#include "neo/base/core/geo/Circle.hpp"
#include "neo/base/core/geo/Polygon.hpp"
// 円と円、円と多角形の共通部分の面積。
namespace intersection_area_internal {
template <class R> R angle(const Point<R>& p, const Point<R>& q) { return std::atan2((long double)cross(p, q), (long double)dot(p, q)); }
}
template <class R> R intersection_area(const Circle<R>& c, const Circle<R>& d) {
 R l= dist2(d.o, c.o), g= c.r + d.r;
 if(sgn(l - g * g) >= 0) return 0;
 R h= std::abs(c.r - d.r), a= c.r * c.r, b= d.r * d.r;
 if(sgn(l - h * h) <= 0) return std::min(a, b) * std::numbers::pi_v<R>;
 R m= std::sqrt(l);
 // 三角形の面積は Heron の式で作る。接する近くで acos の引数が 1 に近づくと精度を落とし、
 // 半径の 2 乗 (1e8 など) を掛けたときに効くので、角度は atan2(高さ, 射影) で作る。
 R tri4= (g + m) * (m + h) * (m - h) * (g - m);
 R K= std::sqrt(std::max(R(0), tri4)) / 4;
 R ht= 2 * K / m;
 R z= (l + a - b) / (m * 2);
 R alpha= std::atan2(ht, z), beta= std::atan2(ht, m - z);
 return a * alpha + b * beta - 2 * K;
}
template <class R> R intersection_area(const Circle<R>& c, const Polygon<R>& g) {
 using P= Point<R>;
 auto sub= [&](const P& p, const P& q) {
  P u= q - p;
  R a= norm2(u), b= dot(u, p) / a, e= c.r * c.r, d= b * b - (norm2(p) - e) / a;
  if(sgn(d) <= 0) return e * intersection_area_internal::angle(p, q);
  d= std::sqrt(d);
  P x= p + std::clamp(-b - d, R(0), R(1)) * u, y= p + std::clamp(-b + d, R(0), R(1)) * u;
  return cross(x, y) + e * (intersection_area_internal::angle(p, x) + intersection_area_internal::angle(y, q));
 };
 std::vector<P> ps;
 for(const auto& p: g) ps.emplace_back(p - c.o);
 R ret= 0;
 for(int i= ps.size(); i--;) ret+= sub(ps[i], ps[g.next(i)]);
 return ret / 2;
}
template <class R> R intersection_area(const Polygon<R>& g, const Circle<R>& c) { return intersection_area(c, g); }
