#pragma once
#include <utility>
#include <vector>
#include "neo/data_structure/RangeMin.hpp"
#include "neo/graph/HeavyLightDecomposition.hpp"
// 木の lca を O(1) で答える。HeavyLightDecomposition の位置 (行きがけ順) で、位置 i の頂点の親の位置を A[i] に並べ、区間最小
// (RangeMin) を載せる。u と v の位置を x < y とすると、A[x + 1, y] の最小の位置の頂点が lca になる。u が v の祖先なら u の
// 位置が、そうでなければ lca から v へ向かう子の親として lca の位置が範囲に現れ、範囲の頂点はどれも lca の真の子孫なので、
// 親の位置が lca より前に行くことはない。表は N = 5 * 10^5 で 8 MB ほど。HLD の lca は鎖を上がる O(log N) で、二分木のような
// 形では 3.8 倍ほど遅く、一直線に近い木ではこちらより速い (記録は algo-notes の notes/rmq.md と notes/tree-queries.md)。
class LCA {
 std::vector<int> pos_, vert_;
 RangeMin rm_;
 static std::vector<unsigned> parents(const HeavyLightDecomposition& h) {
  const int n= int(h.size());
  std::vector<unsigned> a(n);
  for(int i= 0; i < n; ++i) {
   const int p= h.parent(h.to_vertex(i));
   a[i]= p < 0 ? unsigned(i) : unsigned(h.to_seq(p));
  }
  return a;
 }
public:
 LCA()= default;
 // h は作ったあとに捨ててよい。
 explicit LCA(const HeavyLightDecomposition& h): pos_(h.size()), vert_(h.size()), rm_(parents(h)) {
  for(int i= 0, n= int(h.size()); i < n; ++i) vert_[i]= h.to_vertex(i), pos_[vert_[i]]= i;
 }
 // 同じ木の 2 頂点の lca。
 int operator()(int u, int v) const {
  int x= pos_[u], y= pos_[v];
  if(x == y) return u;
  if(x > y) std::swap(x, y);
  return vert_[int(rm_.min(x + 1, y + 1))];
 }
};
