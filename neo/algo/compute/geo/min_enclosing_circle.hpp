#pragma once
#include <utility>
#include <vector>
#include "neo/base/core/geo/Point.hpp"
// 最小包含円を、円を決める 1〜3 点 (支援集合) として返す。点をランダムな順に足し、いまの円の外に出た点があれば、
// その点を円周に載せた形で作り直す。整数の座標なら、直径と外接円の内外を 1 段広い型の行列式で exact に判定する。
// 浮動小数なら sgn の EPS で判定する。乱数は splitmix64 の固定の種を使うので、同じ入力なら同じ結果になる。
// 円の中心と半径は、支援集合から sqrt で作る (中心と半径は浮動小数になる)。
namespace mec_internal {
using u64= unsigned long long;
constexpr u64 xs64(u64& s) { return s^= s << 7, s^= s >> 9; }
template <class K> void shuffle(std::vector<Point<K>>& a, u64 seed= 0x9E3779B97F4A7C15ULL) {
 if(a.size() <= 1) return;
 u64 s= seed;
 for(int i= (int)a.size() - 1; i >= 1; --i) {
  u64 r= xs64(s);
  int j= (int)(r % (u64)(i + 1));
  if(i != j) std::swap(a[i], a[j]);
 }
}
// +1: 内, 0: 上, -1: 外
template <class K> constexpr int diameter_side(const Point<K>& a, const Point<K>& b, const Point<K>& p) {
 using W= wide_t<K>;
 W ux= W(p.x) * 2 - a.x - b.x, uy= W(p.y) * 2 - a.y - b.y;
 W lhs= ux * ux + uy * uy, rhs= dist2(a, b);
 return lhs < rhs ? 1 : lhs == rhs ? 0 : -1;
}
// a, b, c の外接円に対する p の、向きを直した incircle 行列式。
template <class K> constexpr wide_t<K> incircle_det(const Point<K>& a, const Point<K>& b, const Point<K>& c, const Point<K>& p) {
 using W= wide_t<K>;
 W ax= a.x - p.x, ay= a.y - p.y;
 W bx= b.x - p.x, by= b.y - p.y;
 W cx= c.x - p.x, cy= c.y - p.y;
 W A= ax * ax + ay * ay, B= bx * bx + by * by, C= cx * cx + cy * cy;
 return ax * (by * C - cy * B) - ay * (bx * C - cx * B) + A * (bx * cy - by * cx);
}
// +1: 内, 0: 上, -1: 外。a, b, c が一直線なら、最も遠い 2 点を直径とする円で判定する。
template <class K> constexpr int circumcircle_side(const Point<K>& a, const Point<K>& b, const Point<K>& c, const Point<K>& p) {
 int s= sgn(cross(b - a, c - a));
 if(s == 0) {
  auto dab= dist2(a, b), dbc= dist2(b, c), dac= dist2(a, c);
  if(dab >= dbc && dab >= dac) return diameter_side(a, b, p);
  if(dbc >= dab && dbc >= dac) return diameter_side(b, c, p);
  return diameter_side(a, c, p);
 }
 auto det= incircle_det(a, b, c, p);
 if(s < 0) det= -det;
 return det > 0 ? 1 : det == 0 ? 0 : -1;
}
}
// 円を決める 1〜3 点。m は点の数で、m = 1 なら a、m = 2 なら a と b を直径とする円、m = 3 なら a, b, c の外接円。
template <class K> struct MecSupport {
 int m= 0;
 Point<K> a, b, c;
 // +1: 円の内, 0: 円の上, -1: 円の外
 constexpr int side(const Point<K>& p) const {
  if(m == 0) return -1;
  if(m == 1) return p == a ? 0 : -1;
  if(m == 2) return mec_internal::diameter_side(a, b, p);
  return mec_internal::circumcircle_side(a, b, c, p);
 }
};
// a, b, c が一直線なら、最も遠い 2 点を直径とする形に直す。
template <class K> void mec_reduce_if_collinear(MecSupport<K>& s) {
 if(s.m != 3) return;
 if(sgn(cross(s.b - s.a, s.c - s.a)) != 0) return;
 auto dab= dist2(s.a, s.b), dbc= dist2(s.b, s.c), dac= dist2(s.a, s.c);
 if(dab >= dbc && dab >= dac) s.m= 2;
 else if(dbc >= dab && dbc >= dac) s.m= 2, s.a= s.b, s.b= s.c;
 else s.m= 2, s.b= s.c;
}
template <class K> MecSupport<K> min_enclosing_circle(std::vector<Point<K>> ps) {
 mec_internal::shuffle(ps);
 MecSupport<K> s;
 for(int i= 0; i < (int)ps.size(); ++i) {
  if(s.side(ps[i]) >= 0) continue;
  s.m= 1, s.a= ps[i];
  for(int j= 0; j < i; ++j) {
   if(s.side(ps[j]) >= 0) continue;
   s.m= 2, s.a= ps[i], s.b= ps[j];
   for(int k= 0; k < j; ++k) {
    if(s.side(ps[k]) >= 0) continue;
    s.m= 3, s.a= ps[i], s.b= ps[j], s.c= ps[k];
    mec_reduce_if_collinear(s);
   }
  }
 }
 return s;
}
