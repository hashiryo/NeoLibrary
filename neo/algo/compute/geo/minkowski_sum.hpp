#pragma once
#include <algorithm>
#include <iterator>
#include <vector>
#include "neo/algo/compute/geo/half_plane_intersection.hpp"
// 2 つの凸多角形 (反時計回り) の Minkowski 和 {x + y : x in a, y in b} を O(N + M) で作る。
// 辺のベクトルを角度の順に併合し、いちばん下 (同じなら左) の頂点の和から順に足す。
namespace minkowski_sum_internal {
// いちばん下、同じ高さなら左の頂点。
template <class K> int lowest(const Convex<K>& c) {
 int k= 0;
 for(int i= 1; i < (int)c.size(); ++i)
  if(c[i].y < c[k].y || (c[i].y == c[k].y && c[i].x < c[k].x)) k= i;
 return k;
}
// 頂点 k から反時計回りにたどった辺のベクトル。向きは角度の順になる。
template <class K> std::vector<Point<K>> edge_vectors(const Convex<K>& c, int k) {
 std::vector<Point<K>> es;
 int n= c.size();
 for(int i= 0; i < n; ++i) es.emplace_back(c[(k + i + 1) % n] - c[(k + i) % n]);
 return es;
}
}
template <class K> Convex<K> minkowski_sum(const Convex<K>& a, const Convex<K>& b) {
 using P= Point<K>;
 if(a.size() == 0 || b.size() == 0) return Convex<K>();
 int ka= minkowski_sum_internal::lowest(a), kb= minkowski_sum_internal::lowest(b);
 auto ea= minkowski_sum_internal::edge_vectors(a, ka), eb= minkowski_sum_internal::edge_vectors(b, kb);
 std::vector<P> es;
 es.reserve(ea.size() + eb.size());
 std::merge(ea.begin(), ea.end(), eb.begin(), eb.end(), std::back_inserter(es), [](const P& x, const P& y) { return half_plane_intersection_internal::angle_less(x, y); });
 std::vector<P> ps;
 P cur= a[ka] + b[kb];
 ps.push_back(cur);
 for(const auto& e: es) cur+= e, ps.push_back(cur);
 return Convex<K>(ps);
}
