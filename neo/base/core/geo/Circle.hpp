#pragma once
#include <cmath>
#include <numbers>
#include <vector>
#include "neo/base/core/geo/Segment.hpp"
// 円。中心 o と半径 r。浮動小数を受ける。
template <class R> struct Circle {
 using P= Point<R>;
 P o;
 R r;
 Circle() {}
 Circle(const P& o, R r): o(o), r(r) {}
 R area() const { return r * r * std::numbers::pi_v<R>; }
 // +1: 内, 0: 上, -1: 外
 int where(const P& p) const { return sgn(r * r - dist2(p, o)); }
 // +1: 2 点で交わる, 0: 接する, -1: 交わらない
 int where(const Line<R>& l) const { return sgn(r * r - dist2(l, o)); }
 // c がこの円の中に収まるか。
 bool in(const Circle& c) const {
  R a= c.r - r;
  return sgn(a) <= 0 && sgn(dist2(o, c.o) - a * a) <= 0;
 }
 std::vector<Line<R>> tangent(const P& p) const {
  P d= p - o, e= !d;
  R b= norm2(d), a= b - r * r;
  if(int s= sgn(a); s < 0) return {};
  else if(s == 0) return {{p, e}};
  d*= r, e*= std::sqrt(a);
  return {Line<R>(p, !(d + e)), Line<R>(p, !(d - e))};
 }
 friend std::ostream& operator<<(std::ostream& os, const Circle& c) { return os << c.o << " " << c.r; }
};
// 2 点で交わるなら 2 点、接するなら 1 点、交わらないなら空、同じ円なら円周上の 3 点。
// 円 c の周上を p -> q と反時計回りに結んだ弧が円 d の内部にあるようにしてある。
template <class R> std::vector<Point<R>> cross_points(const Circle<R>& c, const Circle<R>& d) {
 Point<R> v= d.o - c.o;
 R g= norm2(v), a= c.r - d.r, b= c.r + d.r;
 if(!sgn(g)) {
  if(sgn(a)) return {};
  return {{c.o.x + c.r, c.o.y}, {c.o.x - c.r, c.o.y}, {c.o.x, c.o.y + c.r}};
 }
 int in= sgn(g - a * a), out= sgn(g - b * b);
 if(in < 0 || out > 0) return {};
 if(!in) return {(c.r * d.o - d.r * c.o) / a};
 if(!out) return {(c.r * d.o + d.r * c.o) / b};
 R e= (a * b + g) / (g * 2);
 Point<R> q= c.o + e * v, n= !v * std::sqrt(c.r * c.r / g - e * e);
 return {q - n, q + n};
}
// 円と直線の交点。接するなら 1 点、交わらないなら空。
template <class R> std::vector<Point<R>> cross_points(const Circle<R>& c, const Line<R>& l) {
 Point<R> v= l.p - c.o;
 R a= norm2(l.d), b= dot(l.d, v) / a, d= b * b - (norm2(v) - c.r * c.r) / a;
 int s= sgn(d);
 if(s < 0) return {};
 if(!s) return {l.p - b * l.d};
 d= std::sqrt(d);
 return {l.p - (b + d) * l.d, l.p - (b - d) * l.d};
}
template <class R> std::vector<Point<R>> cross_points(const Line<R>& l, const Circle<R>& c) { return cross_points(c, l); }
template <class R> std::vector<Point<R>> cross_points(const Circle<R>& c, const Segment<R>& s) {
 Point<R> u= s.q - s.p, v= s.p - c.o;
 R a= norm2(u), b= dot(u, v) / a, d= b * b - (norm2(v) - c.r * c.r) / a;
 int t= sgn(d);
 if(t < 0) return {};
 if(!t && sgn(b) <= 0 && sgn(1 + b) >= 0) return {s.p - b * u};
 d= std::sqrt(d), a= -b - d, b= -b + d;
 std::vector<Point<R>> ps;
 if(0 <= sgn(a) && sgn(a - 1) <= 0) ps.emplace_back(s.p + a * u);
 if(0 <= sgn(b) && sgn(b - 1) <= 0) ps.emplace_back(s.p + b * u);
 return ps;
}
template <class R> std::vector<Point<R>> cross_points(const Segment<R>& s, const Circle<R>& c) { return cross_points(c, s); }
// 三角形 ABC の外接円。
template <class R> Circle<R> circumscribed_circle(const Point<R>& A, const Point<R>& B, const Point<R>& C) {
 Point<R> u= !(B - A), v= C - A, o= (A + B + dot(C - B, v) / dot(u, v) * u) / 2;
 return {o, std::sqrt((long double)dist2(A, o))};
}
// 三角形 ABC の内接円。半径は 2 * 面積 / 周長で作る。細い三角形で (s - a)(s - b)(s - c) を計算すると桁落ちするため。
template <class R> Circle<R> inscribed_circle(const Point<R>& A, const Point<R>& B, const Point<R>& C) {
 R a= std::sqrt((long double)dist2(B, C)), b= std::sqrt((long double)dist2(C, A)), c= std::sqrt((long double)dist2(A, B)), s= (a + b + c) / 2;
 R r= std::abs((long double)cross(B - A, C - A)) / (a + b + c);
 return {(a * A + b * B + c * C) / (s * 2), r};
}
// 2 円の共通接線。
template <class R> std::vector<Line<R>> common_tangent(const Circle<R>& c, const Circle<R>& d) {
 Point<R> u= d.o - c.o, v= !u;
 R g= norm2(u), b;
 if(!sgn(g)) return {};
 std::vector<Line<R>> ls;
 for(R a: {c.r - d.r, c.r + d.r}) {
  if(int s= sgn(b= g - a * a); !s) ls.emplace_back(Line<R>(c.o + c.r * a / g * u, v));
  else if(s > 0) {
   Point<R> x= a / g * u, y= std::sqrt(b) / g * v, e= x + y, f= x - y;
   ls.emplace_back(Line<R>(c.o + c.r * e, !e)), ls.emplace_back(Line<R>(c.o + c.r * f, !f));
  }
 }
 return ls;
}
