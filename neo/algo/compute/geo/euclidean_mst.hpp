#pragma once
#include <algorithm>
#include <utility>
#include <vector>
#include "neo/algo/compute/geo/delaunay.hpp"
#include "neo/base/core/geo/Point.hpp"
#include "neo/base/ds/set/UnionFind.hpp"
// ユークリッド距離の MST。最小全域木の辺はドロネー三角形分割の辺に入るので、ドロネーの辺だけを
// 距離の 2 乗で並べて Kruskal する。距離の 2 乗の順と距離の順は同じなので、2 乗のまま比べてよい。
template <class K> std::vector<std::pair<int, int>> euclidean_mst(const std::vector<Point<K>>& ps) {
 auto es= delaunay_edges(ps);
 std::sort(es.begin(), es.end(), [&](auto e, auto f) { return dist2(ps[e.first], ps[e.second]) < dist2(ps[f.first], ps[f.second]); });
 UnionFind uf(ps.size());
 std::vector<std::pair<int, int>> res;
 for(auto [u, v]: es)
  if(!uf.same(u, v)) uf.unite(u, v), res.emplace_back(u, v);
 return res;
}
