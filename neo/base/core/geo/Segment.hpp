#pragma once
#include <algorithm>
#include <cmath>
#include <vector>
#include "neo/base/core/geo/Line.hpp"
// 線分 p q。
template <class K> struct Segment {
 using P= Point<K>;
 P p, q;
 Segment() {}
 Segment(const P& p, const P& q): p(p), q(q) {}
 bool operator==(const Segment& s) const { return (p == s.p && q == s.q) || (p == s.q && q == s.p); }
 bool operator!=(const Segment& s) const { return !(*this == s); }
 bool on(const P& r) const { return ccw(p, q, r) == Ccw::on_segment; }
 P& operator[](int i) { return i ? q : p; }
 const P& operator[](int i) const { return i ? q : p; }
 long double length() const { return std::sqrt((long double)dist2(p, q)); }
 P closest_point(const P& r) const {
  P d= q - p;
  K a= dot(r - p, d), b;
  return sgn(a) > 0 ? sgn(a - (b= norm2(d))) < 0 ? p + a / b * d : q : p;
 }
 friend std::ostream& operator<<(std::ostream& os, const Segment& s) { return os << s.p << "---" << s.q; }
};
// 1 点で交わるならその点、交わらないなら空、同じ直線なら 2 点。
template <class K> std::vector<Point<K>> cross_points(const Segment<K>& s, const Line<K>& l) {
 Point<K> d= s.q - s.p;
 K a= cross(d, l.d), b= cross(l.p - s.p, l.d);
 if(sgn(a)) {
  if(b/= a; sgn(b) < 0 || sgn(b - 1) > 0) return {};
  else return {s.p + b * d};
 }
 if(sgn(b)) return {};
 return {s.p, s.q};
}
template <class K> std::vector<Point<K>> cross_points(const Line<K>& l, const Segment<K>& s) { return cross_points(s, l); }
// 1 点で交わるならその点、交わらないなら空、同じ直線なら重なりの端点。
template <class K> std::vector<Point<K>> cross_points(const Segment<K>& s, const Segment<K>& t) {
 Point<K> d= s.q - s.p, e= t.q - t.p;
 K a= cross(d, e), b= cross(t.p - s.p, e);
 if(sgn(a)) {
  if(b/= a; sgn(b) < 0 || sgn(b - 1) > 0) return {};
  if(b= cross(d, s.p - t.p) / a; sgn(b) < 0 || sgn(b - 1) > 0) return {};
  return {t.p + b * e};
 }
 if(sgn(b)) return {};
 std::vector<Point<K>> ps;
 auto insert_if_possible= [&](const Point<K>& p) {
  for(auto q: ps)
   if(p == q) return;
  ps.emplace_back(p);
 };
 if(sgn(dot(t.p - s.p, t.q - s.p)) <= 0) insert_if_possible(s.p);
 if(sgn(dot(t.p - s.q, t.q - s.q)) <= 0) insert_if_possible(s.q);
 if(sgn(dot(s.p - t.p, s.q - t.p)) <= 0) insert_if_possible(t.p);
 if(sgn(dot(s.p - t.q, s.q - t.q)) <= 0) insert_if_possible(t.q);
 return ps;
}
enum class Intersection { crossing, touching, disjoint, overlap };
template <class K> Intersection intersection(const Segment<K>& s, const Segment<K>& t) {
 auto cp= cross_points(s, t);
 if(cp.size() == 0) return Intersection::disjoint;
 if(cp.size() == 2) return Intersection::overlap;
 // 端点かどうかは距離で見る。浮動小数では、交点が端点と計算誤差だけずれるので、座標の == では足りない。
 auto at_end= [&](const Point<K>& q) { return !sgn(dist2(cp[0], q)); };
 return at_end(s.p) || at_end(s.q) || at_end(t.p) || at_end(t.q) ? Intersection::touching : Intersection::crossing;
}
template <class K> wide_t<K> dist2(const Segment<K>& s, const Point<K>& p) { return dist2(p, s.closest_point(p)); }
template <class K> wide_t<K> dist2(const Point<K>& p, const Segment<K>& s) { return dist2(s, p); }
template <class K> wide_t<K> dist2(const Segment<K>& s, const Line<K>& l) { return cross_points(s, l).size() ? 0 : std::min(dist2(s.p, l), dist2(s.q, l)); }
template <class K> wide_t<K> dist2(const Line<K>& l, const Segment<K>& s) { return dist2(s, l); }
template <class K> wide_t<K> dist2(const Segment<K>& s, const Segment<K>& t) { return cross_points(s, t).size() ? 0 : std::min({dist2(s, t.p), dist2(s, t.q), dist2(t, s.p), dist2(t, s.q)}); }
template <class K> long double dist(const Segment<K>& s, const Point<K>& p) { return std::sqrt((long double)dist2(s, p)); }
template <class K> long double dist(const Point<K>& p, const Segment<K>& s) { return dist(s, p); }
template <class K> long double dist(const Segment<K>& s, const Segment<K>& t) { return std::sqrt((long double)dist2(s, t)); }
template <class K> Segment<K> Affine<K>::operator()(const Segment<K>& s) const { return {(*this)(s.p), (*this)(s.q)}; }
