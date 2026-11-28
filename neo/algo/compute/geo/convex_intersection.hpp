#pragma once
#include <algorithm>
#include <iterator>
#include <vector>
#include "neo/algo/compute/geo/half_plane_intersection.hpp"
// 2 つの凸多角形 (反時計回り) の共通部分を O(N + M) で作る。凸多角形の辺の向きは角度の順に並んでいるので、
// 2 つの辺の列を角度で併合して、半平面交差の deque にそのまま渡す。K は浮動小数 (交点が有理数になるため)。
// 面積が 0 のとき (空、点、線分) は空の Convex を返す。
namespace convex_intersection_internal {
// 最小の向きの辺から始まるように回した辺の列。凸多角形の辺の向きは、この回し方で角度の順になる。
template <class K> std::vector<Line<K>> sorted_edges(const Convex<K>& c) {
 std::vector<Line<K>> es;
 for(const auto& e: c.edges()) es.emplace_back(e.p, e.q - e.p);
 if(es.empty()) return es;
 int st= 0;
 for(int i= 1; i < (int)es.size(); ++i)
  if(half_plane_intersection_internal::angle_less(es[i].d, es[st].d)) st= i;
 std::rotate(es.begin(), es.begin() + st, es.end());
 return es;
}
}
template <class K> Convex<K> convex_intersection(const Convex<K>& a, const Convex<K>& b) {
 using L= Line<K>;
 if(a.size() < 3 || b.size() < 3) return Convex<K>();
 auto ea= convex_intersection_internal::sorted_edges(a), eb= convex_intersection_internal::sorted_edges(b);
 std::vector<L> ls;
 ls.reserve(ea.size() + eb.size());
 std::merge(ea.begin(), ea.end(), eb.begin(), eb.end(), std::back_inserter(ls), [](const L& x, const L& y) { return half_plane_intersection_internal::angle_less(x.d, y.d); });
 return half_plane_intersection_internal::solve(std::move(ls));
}
