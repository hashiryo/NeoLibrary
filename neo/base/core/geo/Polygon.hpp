#pragma once
#include <algorithm>
#include <cassert>
#include <limits>
#include <vector>
#include "neo/base/core/geo/Segment.hpp"
// 反時計回りの多角形。area は 2 倍の面積を 2 で割ったもので、整数なら切り捨て、area2 は 2 倍の面積。
template <class K> class Polygon {
 using P= Point<K>;
 wide_t<K> a2= 0;
protected:
 std::vector<P> dat;
 void build() {
  if(dat.empty()) return;
  a2= cross(dat.back(), dat[0]);
  for(int i= this->size(); --i;) a2+= cross(dat[i - 1], dat[i]);
  if(a2 < 0) std::reverse(dat.begin(), dat.end()), a2= -a2;
 }
 template <int opp> bool contain(const Segment<K>& s) const {
  assert(dat.size());
  if(where(s.p) == opp || where(s.q) == opp) return false;
  std::vector<P> ps, qs;
  for(const auto& e: edges())
   if(auto cp= cross_points(s, e); cp.size()) ps.insert(ps.end(), cp.begin(), cp.end());
  if(ps.empty()) return true;
  std::sort(ps.begin(), ps.end()), ps.erase(std::unique(ps.begin(), ps.end()), ps.end());
  for(int i= ps.size(); --i;) qs.emplace_back((ps[i] + ps[i - 1]) / 2);
  for(const auto& q: qs)
   if(where(q) == opp) return false;
  return true;
 }
public:
 Polygon() {}
 Polygon(const std::vector<P>& ps): dat(ps) { build(); }
 int prev(int i) const { return i ? i - 1 : (int)this->size() - 1; }
 int next(int i) const { return i + 1 >= (int)this->size() ? 0 : i + 1; }
 const P& operator[](int i) const { return dat[i]; }
 auto begin() const { return dat.begin(); }
 auto end() const { return dat.end(); }
 std::size_t size() const { return dat.size(); }
 std::vector<Segment<K>> edges() const {
  std::vector<Segment<K>> es;
  for(int i= 0, e= dat.size(); i < e; ++i) es.emplace_back(dat[i], dat[next(i)]);
  return es;
 }
 bool is_convex() const {
  assert(dat.size());
  for(int i= dat.size(); i--;)
   if(P a= dat[i], b= a - dat[prev(i)], c= dat[next(i)] - a; sgn(cross(b, c)) < 0) return false;
  return true;
 }
 wide_t<K> area() const { return a2 / 2; }
 wide_t<K> area2() const { return a2; }
 // +1: 内, 0: 周上, -1: 外
 int where(const P& p) const {
  assert(dat.size());
  bool in= false;
  for(int i= dat.size(); i--;) {
   Point<K> a= dat[i] - p, b= dat[next(i)] - p;
   if(a.y > b.y) std::swap(a, b);
   int s= sgn(cross(a, b));
   if(!s && sgn(dot(a, b)) <= 0) return 0;
   if(s < 0 && sgn(a.y) <= 0 && 0 < sgn(b.y)) in= !in;
  }
  return in ? 1 : -1;
 }
 bool in(const Segment<K>& s) const { return contain<-1>(s); }
 bool out(const Segment<K>& s) const { return contain<+1>(s); }
};
template <class K> wide_t<K> dist2(const Polygon<K>& g, const Point<K>& p) {
 if(g.where(p) != -1) return 0;
 auto ret= std::numeric_limits<wide_t<K>>::max();
 for(const auto& e: g.edges()) ret= std::min(ret, dist2(e, p));
 return ret;
}
template <class K> wide_t<K> dist2(const Point<K>& p, const Polygon<K>& g) { return dist2(g, p); }
template <class K> wide_t<K> dist2(const Polygon<K>& g, const Line<K>& l) {
 auto ret= std::numeric_limits<wide_t<K>>::max();
 for(const auto& e: g.edges()) ret= std::min(ret, dist2(e, l));
 return ret;
}
template <class K> wide_t<K> dist2(const Line<K>& l, const Polygon<K>& g) { return dist2(g, l); }
template <class K> wide_t<K> dist2(const Polygon<K>& g, const Segment<K>& s) {
 if(g.where(s.p) != -1 || g.where(s.q) != -1) return 0;
 auto ret= std::numeric_limits<wide_t<K>>::max();
 for(const auto& e: g.edges()) ret= std::min(ret, dist2(e, s));
 return ret;
}
template <class K> wide_t<K> dist2(const Segment<K>& s, const Polygon<K>& g) { return dist2(g, s); }
template <class K> wide_t<K> dist2(const Polygon<K>& g, const Polygon<K>& h) {
 auto ret= std::numeric_limits<wide_t<K>>::max();
 for(const auto& e: g.edges()) ret= std::min(ret, dist2(h, e));
 return ret;
}
template <class K> long double dist(const Polygon<K>& g, const Point<K>& p) { return std::sqrt((long double)dist2(g, p)); }
template <class K> long double dist(const Point<K>& p, const Polygon<K>& g) { return dist(g, p); }
