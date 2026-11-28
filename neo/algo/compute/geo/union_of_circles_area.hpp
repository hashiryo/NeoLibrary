#pragma once
#include <algorithm>
#include <cmath>
#include <numbers>
#include <utility>
#include <vector>
#include "neo/base/core/geo/Circle.hpp"
// 円の和集合の面積。各円の周のうち、他のどの円の内側にも入らない弧を集め、Green の定理で面積にする。
// O(N^2 log N)。R は浮動小数。同じ円は先に 1 つにする。
namespace union_of_circles_area_internal {
// 円 c の弧 [a, b] (中心からの偏角) が作る面積 (1/2 ∫ (x dy - y dx))。
template <class R> R arc_area(const Circle<R>& c, R a, R b) {
 return (c.r * c.r * (b - a) + c.r * (c.o.x * (std::sin(b) - std::sin(a)) - c.o.y * (std::cos(b) - std::cos(a)))) / 2;
}
}
template <class R> R union_of_circles_area(std::vector<Circle<R>> cs) {
 using namespace union_of_circles_area_internal;
 const R tau= 2 * std::numbers::pi_v<R>;
 std::sort(cs.begin(), cs.end(), [](const Circle<R>& a, const Circle<R>& b) {
  if(a.o.x != b.o.x) return a.o.x < b.o.x;
  if(a.o.y != b.o.y) return a.o.y < b.o.y;
  return a.r < b.r;
 });
 cs.erase(std::unique(cs.begin(), cs.end(), [](const Circle<R>& a, const Circle<R>& b) { return a.o == b.o && a.r == b.r; }), cs.end());
 R ans= 0;
 for(int i= 0; i < (int)cs.size(); ++i) {
  const R ri= cs[i].r;
  std::vector<std::pair<R, R>> iv;
  bool covered= false;
  for(int j= 0; j < (int)cs.size(); ++j) {
   if(i == j) continue;
   const R rj= cs[j].r, d= dist(cs[i].o, cs[j].o);
   // i が j の中にあるなら、i の周は 1 点を除いて j の内側。
   if(sgn(d + ri - rj) <= 0) {
    covered= true;
    break;
   }
   // j が i の中、または外で離れるなら、i の周は j に覆われない。
   if(sgn(d + rj - ri) <= 0 || sgn(d - (ri + rj)) >= 0) continue;
   // i の周のうち j の内側に入る偏角の区間。
   R c= (ri * ri + d * d - rj * rj) / (2 * ri * d);
   c= std::clamp(c, R(-1), R(1));
   R a= std::atan2(cs[j].o.y - cs[i].o.y, cs[j].o.x - cs[i].o.x) - std::acos(c);
   R b= a + 2 * std::acos(c);
   a= std::fmod(a, tau), a+= a < 0 ? tau : 0;
   b= std::fmod(b, tau), b+= b < 0 ? tau : 0;
   if(a <= b) iv.emplace_back(a, b);
   else iv.emplace_back(a, tau), iv.emplace_back(0, b);
  }
  if(covered || iv.empty()) {
   if(!covered) ans+= std::numbers::pi_v<R> * ri * ri;
   continue;
  }
  std::sort(iv.begin(), iv.end());
  R prev= 0;
  for(auto [a, b]: iv) {
   if(b <= prev) continue;
   if(a > prev) ans+= arc_area(cs[i], prev, a);
   prev= std::max(prev, b);
  }
  if(prev < tau) ans+= arc_area(cs[i], prev, tau);
 }
 return ans;
}
