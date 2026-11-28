#pragma once
#include <algorithm>
#include <span>
#include <utility>
#include <vector>
#include "neo/base/core/geo/Segment.hpp"
#include "neo/base/core/geo/angle.hpp"
// 線分の集合を交点で分けて平面グラフを作る。頂点、辺、面を持ち、面の面積と外側かを返す。
// 交差は線分の全ペアで調べる O(n^2)。辺は両方向に持ち、e と e ^ 1 が同じ辺の逆方向。
// 頂点は、整数なら座標の exact な一致で、浮動小数なら sgn の EPS で同じ点とみなしてまとめる。同じ交点でも、
// 違う線分の組から計算すると末尾の bit が違うので、exact の一致だけでまとめると頂点が分かれて面が壊れる。
template <class K> class SegmentArrangement {
 std::vector<Point<K>> ps;
 std::vector<int> dst, pr, nx, icf, cmp, g, pos;
 static bool same_point(const Point<K>& a, const Point<K>& b) { return !sgn(a.x - b.x) && !sgn(a.y - b.y); }
public:
 explicit SegmentArrangement(const std::vector<Segment<K>>& ss) {
  for(const auto& s: ss) ps.push_back(s.p), ps.push_back(s.q);
  for(int i= ss.size(); i--;)
   for(int j= i; j--;)
    if(auto cp= cross_points(ss[i], ss[j]); cp.size()) ps.insert(ps.end(), cp.begin(), cp.end());
  std::sort(ps.begin(), ps.end());
  ps.erase(std::unique(ps.begin(), ps.end(), [](const Point<K>& a, const Point<K>& b) { return same_point(a, b); }), ps.end());
  const int n= ps.size();
  std::vector<std::pair<int, int>> es;
  for(const auto& s: ss) {
   std::vector<int> l;
   for(int i= n; i--;)
    if(s.on(ps[i])) l.push_back(i);
   std::sort(l.begin(), l.end(), [&](int i, int j) { return ps[i] < ps[j]; });
   for(int i= l.size(); --i;) {
    auto [u, v]= std::minmax(l[i], l[i - 1]);
    es.emplace_back(u, v);
   }
  }
  std::sort(es.begin(), es.end()), es.erase(std::unique(es.begin(), es.end()), es.end());
  for(auto [u, v]: es) dst.push_back(u), dst.push_back(v);
  const int m= dst.size();
  pos.assign(n + 1, 0), g.resize(m), icf.assign(m, -1), pr.resize(m), nx.resize(m);
  for(int v: dst) ++pos[v];
  for(int v= 0; v < n; ++v) pos[v + 1]+= pos[v];
  for(int e= m; e--;) g[--pos[dst[e ^ 1]]]= e;
  std::vector<Point<K>> qs(m);
  for(int e= m; e--;) qs[e]= ps[dst[e]] - ps[dst[e ^ 1]];
  AngleComp<K> ac;
  auto cm= [&](int i, int j) { return ac(qs[i], qs[j]); };
  std::vector<int> nxi(m), pri(m);
  for(int v= n; v--;) {
   std::sort(g.begin() + pos[v], g.begin() + pos[v + 1], cm);
   for(int i= pos[v] + 1; i < pos[v + 1]; ++i) nxi[g[i - 1]]= g[i], pri[g[i]]= g[i - 1];
   nxi[g[pos[v + 1] - 1]]= g[pos[v]], pri[g[pos[v]]]= g[pos[v + 1] - 1];
  }
  for(int e= m; e--;) nx[e]= pri[e ^ 1], pr[e]= nxi[e] ^ 1;
  for(int e= m; e--;) {
   if(icf[e] != -1) continue;
   int f= cmp.size(), x= e;
   cmp.push_back(e);
   do {
    icf[x]= f, x= nx[x];
   } while(x != e);
  }
 }
 std::size_t vertex_size() const { return ps.size(); }
 std::size_t face_size() const { return cmp.size(); }
 std::span<const int> out_edges(int v) const { return {g.data() + pos[v], std::size_t(pos[v + 1] - pos[v])}; }
 Point<K> point(int v) const { return ps[v]; }
 int vertex(const Point<K>& p) const {
  int i= std::lower_bound(ps.begin(), ps.end(), p) - ps.begin();
  if(i < (int)ps.size() && same_point(ps[i], p)) return i;
  if(i > 0 && same_point(ps[i - 1], p)) return i - 1;
  if(i + 1 < (int)ps.size() && same_point(ps[i + 1], p)) return i + 1;
  return i;
 }
 int to_v(int e) const { return dst[e]; }
 int from_v(int e) const { return dst[e ^ 1]; }
 int twin(int e) const { return e ^ 1; }
 int next(int e) const { return nx[e]; }
 int prev(int e) const { return pr[e]; }
 int incident_face(int e) const { return icf[e]; }
 int component_e(int f) const { return cmp[f]; }
 wide_t<K> area2(int f) const {
  wide_t<K> a= 0;
  for(int e= cmp[f], ne= nx[e];; ne= nx[e= ne]) {
   a+= cross(ps[dst[e]], ps[dst[ne]]);
   if(ne == cmp[f]) break;
  }
  return a;
 }
 wide_t<K> area(int f) const { return area2(f) / 2; }
 bool is_outside(int f) const { return area2(f) < 0; }
};
