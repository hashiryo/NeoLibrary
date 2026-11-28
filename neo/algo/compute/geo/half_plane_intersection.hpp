#pragma once
#include <algorithm>
#include <deque>
#include <vector>
#include "neo/base/core/geo/Convex.hpp"
// 半平面 {x : cross(d, x - p) >= 0} (Line の左側) の共通部分を O(N log N) で作る。半平面を向きで並べ、
// 両端から要らなくなった半平面を落とす。K は浮動小数 (交点が有理数になるため)。共通部分の面積が 0 なら
// 空の Convex を返す。共通部分が有界でないときは閉じた多角形にならないので、呼ぶ側で有界にしてから渡す。
namespace half_plane_intersection_internal {
// 向きの半周。atan2 の (-pi, pi] の代わりに、[0, 2pi) の順に並べるための 0/1。
template <class K> int half(const Point<K>& d) { return sgn(d.y) < 0 || (sgn(d.y) == 0 && sgn(d.x) < 0); }
// 同じ向きか (平行で、向きが逆でない)。
template <class K> bool same_direction(const Line<K>& l, const Line<K>& m) { return !sgn(cross(l.d, m.d)) && sgn(dot(l.d, m.d)) > 0; }
// 直線 l の右側にあるか。
template <class K> bool on_right(const Line<K>& l, const Point<K>& p) { return sgn(cross(l.d, p - l.p)) < 0; }
// 同じ向きなら、相手の境界の点が自分の外側にある方がきつい (係数の大きさによらない)。
template <class K> bool tighter(const Line<K>& l, const Line<K>& m) { return on_right(l, m.p); }
// 平行でない 2 直線の交点。
template <class K> Point<K> cross_point(const Line<K>& l, const Line<K>& m) {
 K a= cross(l.d, m.d), b= cross(m.p - l.p, m.d);
 return l.p + b / a * l.d;
}
}
template <class K> Convex<K> half_plane_intersection(std::vector<Line<K>> ls) {
 using P= Point<K>;
 using L= Line<K>;
 std::sort(ls.begin(), ls.end(), [](const L& a, const L& b) {
  if(int ha= half_plane_intersection_internal::half(a.d), hb= half_plane_intersection_internal::half(b.d); ha != hb) return ha < hb;
  return sgn(cross(a.d, b.d)) > 0;
 });
 std::vector<L> us;
 for(const L& l: ls) {
  if(!us.empty() && half_plane_intersection_internal::same_direction(us.back(), l)) {
   if(half_plane_intersection_internal::tighter(l, us.back())) us.back()= l;
  } else us.push_back(l);
 }
 std::deque<L> q;
 for(const L& l: us) {
  while(q.size() >= 2 && half_plane_intersection_internal::on_right(l, half_plane_intersection_internal::cross_point(q[q.size() - 2], q.back()))) q.pop_back();
  while(q.size() >= 2 && half_plane_intersection_internal::on_right(l, half_plane_intersection_internal::cross_point(q[0], q[1]))) q.pop_front();
  q.push_back(l);
 }
 while(q.size() >= 3 && half_plane_intersection_internal::on_right(q[0], half_plane_intersection_internal::cross_point(q[q.size() - 2], q.back()))) q.pop_back();
 while(q.size() >= 3 && half_plane_intersection_internal::on_right(q.back(), half_plane_intersection_internal::cross_point(q[0], q[1]))) q.pop_front();
 if(q.size() < 3) return Convex<K>();
 std::vector<P> ps;
 for(int i= q.size(); i--;) ps.emplace_back(half_plane_intersection_internal::cross_point(q[i], q[(i + 1) % q.size()]));
 return Convex<K>(ps);
}
