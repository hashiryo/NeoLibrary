#pragma once
#include <algorithm>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>
#include "neo/base/core/geo/Point.hpp"
#include "neo/base/ds/range/Fenwick.hpp"
#include "neo/base/ds/set/UnionFind.hpp"
// マンハッタン距離の MST。点を 45 度回す操作を 4 通り試し、各点について x と y の両方が大きくない点のうち
// x - y が最大のものを Fenwick で探して候補の辺にする。候補の辺を距離で並べ、Kruskal で MST を作る。
// 候補の辺の本数は 4n 本なので、全体は O(n log n)。
namespace manhattan_mst_internal {
// Fenwick に載せる値は x - y。座標の差が K に収まることを前提にする (座標が 1e9 までなら long long で足りる)。
// 128 bit で持つと要素が 24 byte になり、x64 で Fenwick の読み書きが重くなる。
template <class K> struct MaxPoint {
 using T= MaxPoint;
 K v= 0;
 int id= -1;
 static constexpr T ti() { return {std::numeric_limits<K>::lowest(), -1}; }
 static constexpr T op(T a, T b) {
  if(a.id == -1) return b;
  if(b.id == -1) return a;
  return a.v < b.v ? b : a;
 }
};
template <class K> constexpr wide_t<K> abs_diff(K a, K b) { return a < b ? wide_t<K>(b) - a : wide_t<K>(a) - b; }
}
template <class K> struct ManhattanMst {
 wide_t<K> weight= 0;
 std::vector<std::pair<int, int>> edges;
};
template <class K> ManhattanMst<K> manhattan_mst(std::vector<Point<K>> ps) {
 int n= ps.size();
 std::vector<std::pair<int, int>> cand;
 std::vector<int> ids(n);
 std::iota(ids.begin(), ids.end(), 0);
 for(int ph= 0; ph < 4; ++ph) {
  std::sort(ids.begin(), ids.end(), [&](int i, int j) {
   auto a= wide_t<K>(ps[i].x) + ps[i].y, b= wide_t<K>(ps[j].x) + ps[j].y;
   if(a != b) return a > b;
   return ps[i].y > ps[j].y;
  });
  std::vector<K> xs;
  xs.reserve(n);
  for(auto& p: ps) xs.push_back(p.x);
  std::sort(xs.begin(), xs.end());
  xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
  int m= xs.size();
  Fenwick<manhattan_mst_internal::MaxPoint<K>> fw(m);
  for(int id: ids) {
   int xi= std::lower_bound(xs.begin(), xs.end(), ps[id].x) - xs.begin();
   auto ma= fw.prod(xi + 1);
   if(ma.id != -1) cand.emplace_back(id, ma.id);
   fw.add(xi, {ps[id].x - ps[id].y, id});
  }
  for(auto& p: ps) std::swap(p.x, p.y);
  if(ph == 1)
   for(auto& p: ps) p.y= -p.y;
 }
 auto w= [&](int i, int j) { return manhattan_mst_internal::abs_diff(ps[i].x, ps[j].x) + manhattan_mst_internal::abs_diff(ps[i].y, ps[j].y); };
 std::sort(cand.begin(), cand.end(), [&](auto e, auto f) { return w(e.first, e.second) < w(f.first, f.second); });
 UnionFind uf(n);
 ManhattanMst<K> res;
 for(auto [u, v]: cand)
  if(!uf.same(u, v)) uf.unite(u, v), res.weight+= w(u, v), res.edges.emplace_back(u, v);
 return res;
}
