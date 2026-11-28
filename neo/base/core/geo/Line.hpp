#pragma once
#include <cassert>
#include <cmath>
#include <tuple>
#include <vector>
#include "neo/base/core/geo/Point.hpp"
// 直線 p + t d。where は左側が +1、上 (直線上) が 0、右側が -1。
template <class K> struct Line {
 using P= Point<K>;
 P p, d;
 Line() {}
 Line(const P& p, const P& d): p(p), d(d) { assert(sgn(norm2(d))); }
 // ax + by + c = 0。ax + by + c > 0 が左側。
 Line(K a, K b, K c) {
  int sa= sgn(a), sb= sgn(b);
  assert(sa || sb);
  d= P{b, -a};
  p= sb ? P{0, -c / b} : P{-c / a, 0};
 }
 bool operator==(const Line& l) const { return !sgn(cross(d, l.d)) && !where(l.p); }
 bool operator!=(const Line& l) const { return sgn(cross(d, l.d)) || where(l.p); }
 int where(const P& q) const { return sgn(cross(d, q - p)); }
 P project(const P& q) const { return p + dot(q - p, d) / norm2(d) * d; }
 std::tuple<K, K, K> coef() const { return {-d.y, d.x, cross(p, d)}; }
 friend std::ostream& operator<<(std::ostream& os, const Line& l) { return os << l.p << " + t" << l.d; }
};
template <class K> Line<K> line_through(const Point<K>& p, const Point<K>& q) { return Line<K>(p, q - p); }
template <class K> bool is_parallel(const Line<K>& l, const Line<K>& m) { return !sgn(cross(l.d, m.d)); }
template <class K> bool is_orthogonal(const Line<K>& l, const Line<K>& m) { return !sgn(dot(l.d, m.d)); }
// 1 点で交わるならその点、平行で交わらないなら空、同じ直線なら 2 点。
template <class K> std::vector<Point<K>> cross_points(const Line<K>& l, const Line<K>& m) {
 K a= cross(m.d, l.d), b= cross(l.p - m.p, l.d);
 if(sgn(a)) return {m.p + b / a * m.d};
 if(sgn(b)) return {};
 return {m.p, m.p + m.d};
}
// p と q の垂直二等分線。
template <class K> Line<K> bisector(const Point<K>& p, const Point<K>& q) { return Line<K>((p + q) / 2, !(q - p)); }
template <class K> wide_t<K> dist2(const Line<K>& l, const Point<K>& p) {
 auto a= cross(l.d, p - l.p);
 return a * a / norm2(l.d);
}
template <class K> wide_t<K> dist2(const Point<K>& p, const Line<K>& l) { return dist2(l, p); }
template <class K> wide_t<K> dist2(const Line<K>& l, const Line<K>& m) { return is_parallel(l, m) ? dist2(l, m.p) : 0; }
template <class K> long double dist(const Line<K>& l, const Point<K>& p) { return std::sqrt((long double)dist2(l, p)); }
template <class K> long double dist(const Point<K>& p, const Line<K>& l) { return dist(l, p); }
// 直線 l に対する折り返し。
template <class K> Affine<K> reflect(const Line<K>& l) {
 K a= l.d.x * l.d.x, b= l.d.x * l.d.y * 2, c= l.d.y * l.d.y, d= a + c;
 a/= d, b/= d, c/= d, d= a - c;
 return {d, b, b, -d, Point<K>{c * 2 * l.p.x - b * l.p.y, a * 2 * l.p.y - b * l.p.x}};
}
