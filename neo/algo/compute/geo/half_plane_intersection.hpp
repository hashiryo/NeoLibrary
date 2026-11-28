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
// 向きの順。同じ半周なら cross > 0 が前。
template <class K> bool angle_less(const Point<K>& a, const Point<K>& b) {
 if(int ha= half(a), hb= half(b); ha != hb) return ha < hb;
 return sgn(cross(a, b)) > 0;
}
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
// 向きの順に並んだ半平面。隣り合う同じ向きはきつい方だけ残す。
template <class K> std::vector<Line<K>> merge_same_direction(std::vector<Line<K>> ls) {
 std::vector<Line<K>> us;
 for(const auto& l: ls) {
  if(!us.empty() && same_direction(us.back(), l)) {
   if(tighter(l, us.back())) us.back()= l;
  } else us.push_back(l);
 }
 return us;
}
// 向きの順に並んだ半平面の共通部分。両端から要らない半平面を落とす。
template <class K> Convex<K> solve(std::vector<Line<K>> ls) {
 using P= Point<K>;
 using L= Line<K>;
 ls= merge_same_direction(std::move(ls));
 std::deque<L> q;
 for(const L& l: ls) {
  while(q.size() >= 2 && on_right(l, cross_point(q[q.size() - 2], q.back()))) q.pop_back();
  while(q.size() >= 2 && on_right(l, cross_point(q[0], q[1]))) q.pop_front();
  q.push_back(l);
 }
 while(q.size() >= 3 && on_right(q[0], cross_point(q[q.size() - 2], q.back()))) q.pop_back();
 while(q.size() >= 3 && on_right(q.back(), cross_point(q[0], q[1]))) q.pop_front();
 if(q.size() < 3) return Convex<K>();
 std::vector<P> ps;
 for(int i= q.size(); i--;) ps.emplace_back(cross_point(q[i], q[(i + 1) % q.size()]));
 return Convex<K>(ps);
}
}
template <class K> Convex<K> half_plane_intersection(std::vector<Line<K>> ls) {
 std::sort(ls.begin(), ls.end(), [](const Line<K>& a, const Line<K>& b) { return half_plane_intersection_internal::angle_less(a.d, b.d); });
 return half_plane_intersection_internal::solve(std::move(ls));
}
