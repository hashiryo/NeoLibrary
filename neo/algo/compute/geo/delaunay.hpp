#pragma once
#include <algorithm>
#include <numeric>
#include <tuple>
#include <utility>
#include <vector>
#include "neo/base/core/geo/Point.hpp"
// ドロネー三角形分割 (Guibas–Stolfi の分割統治)。点を x の順に並べ、左右の三角形分割を共通接線で繋ぐ。
// 辺は quad-edge 構造で両方向に持ち、辺を消すときは両方向を消して番号を使い回す。
// 座標が整数なら、外接円の内側かの判定は wide_t<K> の行列式で exact。同じ座標の点は 1 つにまとめ、
// 元の添字との辺 (長さ 0) も返す。返すのは三角形分割の辺 (点の添字の組)。
namespace delaunay_internal {
template <class K> struct Builder {
 using P= Point<K>;
 using W= wide_t<K>;
 struct Edge {
  int to, ccw, cw, rev;
  bool enabled= false;
 };
 std::vector<int> open_address;
 std::vector<P> pos;
 std::vector<Edge> edges;
 std::vector<int> mappings;
 int outer_one_edge= -1;
 // d から見て a, b, c の外接円の内側に d があるか。a, b, c は反時計回り。
 static bool inside_circumcircle(P a, P b, P c, P d) {
  a= a - d, b= b - d, c= c - d;
  W val= cross(b, c) * W(norm2(a)) + cross(c, a) * W(norm2(b)) + cross(a, b) * W(norm2(c));
  return val > 0;
 }
 int new_edge_index() {
  if(open_address.empty()) {
   edges.push_back({});
   return (int)edges.size() - 1;
  }
  int res= open_address.back();
  open_address.pop_back();
  return res;
 }
 std::pair<int, int> new_edge(int u, int v) {
  int euv= new_edge_index(), evu= new_edge_index();
  edges[euv].ccw= edges[euv].cw= euv;
  edges[evu].ccw= edges[evu].cw= evu;
  edges[euv].to= v, edges[evu].to= u;
  edges[euv].rev= evu, edges[evu].rev= euv;
  edges[euv].enabled= edges[evu].enabled= true;
  return {euv, evu};
 }
 void erase_single_edge(int e) {
  int eccw= edges[e].ccw, ecw= edges[e].cw;
  edges[eccw].cw= ecw;
  edges[ecw].ccw= eccw;
  edges[e].enabled= false;
 }
 void erase_edge_pair(int e) {
  int ex= edges[e].rev;
  erase_single_edge(e), erase_single_edge(ex);
  open_address.push_back(e), open_address.push_back(ex);
 }
 void insert_ccw_after(int e, int x) {
  int xccw= edges[x].ccw;
  edges[e].ccw= xccw;
  edges[xccw].cw= e;
  edges[e].cw= x;
  edges[x].ccw= e;
 }
 void insert_cw_after(int e, int x) {
  int xcw= edges[x].cw;
  edges[e].cw= xcw;
  edges[xcw].ccw= e;
  edges[e].ccw= x;
  edges[x].cw= e;
 }
 int orientation(int a, int b, int c) const {
  auto ab= pos[b] - pos[a], ac= pos[c] - pos[a];
  return sgn(cross(ab, ac));
 }
 std::pair<int, int> go_next(int, int ea) const {
  int ap= edges[ea].to;
  int eap= edges[edges[ea].rev].ccw;
  return {ap, eap};
 }
 std::pair<int, int> go_prev(int, int ea) const {
  int ap= edges[edges[ea].cw].to;
  int eap= edges[edges[ea].cw].rev;
  return {ap, eap};
 }
 std::tuple<int, int, int, int> go_bottom(int a, int ea, int b, int eb) const {
  while(true) {
   auto [ap, eap]= go_prev(a, ea);
   if(orientation(b, a, ap) > 0) {
    a= ap, ea= eap;
    continue;
   }
   auto [bp, ebp]= go_next(b, eb);
   if(orientation(a, b, bp) < 0) {
    b= bp, eb= ebp;
    continue;
   }
   break;
  }
  return {a, ea, b, eb};
 }
 std::pair<int, int> get_maximum(int a, int ea, bool to_min) const {
  std::pair<int, int> ans= {a, ea};
  int p= a, ep= ea;
  do {
   std::tie(p, ep)= go_next(p, ep);
   if(to_min) ans= std::min(ans, std::make_pair(p, ep));
   else ans= std::max(ans, std::make_pair(p, ep));
  } while(ep != ea);
  return ans;
 }
 std::pair<int, int> merge(int a, int ea, int b, int eb) {
  std::tie(a, ea)= get_maximum(a, ea, false);
  std::tie(b, eb)= get_maximum(b, eb, true);
  auto [al, eal, bl, ebl]= go_bottom(a, ea, b, eb);
  auto [bu, ebu, au, eau]= go_bottom(b, eb, a, ea);
  ebl= edges[ebl].cw;
  ebu= edges[ebu].cw;
  auto [abl, bal]= new_edge(al, bl);
  insert_cw_after(abl, eal);
  insert_ccw_after(bal, ebl);
  if(al == au) eau= abl;
  if(bl == bu) ebu= bal;
  int ap= al, eap= eal;
  int bp= bl, ebp= ebl;
  while(ap != au || bp != bu) {
   int a2= edges[eap].to;
   int b2= edges[ebp].to;
   int nxeap= edges[eap].ccw;
   int nxebp= edges[ebp].cw;
   if(eap != eau && nxeap != abl) {
    int a1= edges[nxeap].to;
    if(inside_circumcircle(pos[ap], pos[bp], pos[a2], pos[a1])) {
     erase_edge_pair(eap);
     eap= nxeap;
     continue;
    }
   }
   if(ebp != ebu && nxebp != bal) {
    int b1= edges[nxebp].to;
    if(inside_circumcircle(pos[b2], pos[ap], pos[bp], pos[b1])) {
     erase_edge_pair(ebp);
     ebp= nxebp;
     continue;
    }
   }
   bool choose_a= ebp == ebu;
   if(eap != eau && ebp != ebu) {
    if(orientation(ap, bp, b2) < 0) choose_a= true;
    else if(orientation(a2, ap, bp) < 0) choose_a= false;
    else choose_a= inside_circumcircle(pos[ap], pos[bp], pos[b2], pos[a2]);
   }
   if(choose_a) {
    nxeap= edges[edges[eap].rev].ccw;
    auto [hab, hba]= new_edge(a2, bp);
    insert_cw_after(hab, nxeap);
    insert_ccw_after(hba, ebp);
    eap= nxeap, ap= a2;
   } else {
    nxebp= edges[edges[ebp].rev].cw;
    auto [hba, hab]= new_edge(b2, ap);
    insert_ccw_after(hba, nxebp);
    insert_cw_after(hab, eap);
    ebp= nxebp, bp= b2;
   }
  }
  return {al, abl};
 }
 std::pair<int, int> solve_range(int l, int r) {
  if(r - l == 2) {
   int u= l, v= l + 1;
   auto [uv, vu]= new_edge(u, v);
   return {u, uv};
  }
  if(r - l == 3) {
   int u= l, v= l + 1, w= l + 2;
   auto [uv, vu]= new_edge(u, v);
   auto [vw, wv]= new_edge(v, w);
   int ccw= orientation(u, v, w);
   if(ccw == 0) insert_ccw_after(vu, vw);
   if(ccw > 0) {
    auto [uw, wu]= new_edge(u, w);
    insert_cw_after(uv, uw);
    insert_cw_after(vw, vu);
    insert_cw_after(wu, wv);
    return {u, uv};
   }
   if(ccw < 0) {
    auto [uw, wu]= new_edge(u, w);
    insert_ccw_after(uv, uw);
    insert_ccw_after(vw, vu);
    insert_ccw_after(wu, wv);
    return {v, vu};
   }
   return {u, uv};
  }
  int m= (l + r) / 2;
  auto [a, ea]= solve_range(l, m);
  auto [b, eb]= solve_range(m, r);
  return merge(a, ea, b, eb);
 }
 void solve() {
  int sz= pos.size();
  if(sz <= 1) return;
  std::vector<int> pi(sz);
  std::iota(pi.begin(), pi.end(), 0);
  std::stable_sort(pi.begin(), pi.end(), [&](int l, int r) { return pos[l].x != pos[r].x ? pos[l].x < pos[r].x : pos[l].y < pos[r].y; });
  auto posbuf= pos;
  int posptr= 0;
  mappings.assign(sz, 0);
  for(int i= 0; i < sz; ++i) {
   int v= pi[i];
   if(i == 0 || !(posbuf[pi[posptr - 1]] == posbuf[v])) {
    pi[posptr]= v;
    pos[posptr++]= posbuf[v];
    mappings[v]= v;
   } else mappings[v]= pi[posptr - 1];
  }
  if(posptr >= 2) outer_one_edge= solve_range(0, posptr).second;
  std::swap(pos, posbuf);
  for(auto& e: edges) e.to= pi[e.to];
 }
 explicit Builder(std::vector<P> x_points): pos(std::move(x_points)) { solve(); }
 std::vector<std::pair<int, int>> get_edges() const {
  std::vector<std::pair<int, int>> res;
  for(int e= 0; e < (int)edges.size(); ++e)
   if(edges[e].enabled) {
    int re= edges[e].rev;
    if(e < re) continue;
    res.push_back({edges[e].to, edges[re].to});
   }
  for(int v= 0; v < (int)mappings.size(); ++v)
   if(mappings[v] != v) res.push_back({v, mappings[v]});
  return res;
 }
};
}
template <class K> std::vector<std::pair<int, int>> delaunay_edges(const std::vector<Point<K>>& ps) {
 delaunay_internal::Builder<K> b(ps);
 return b.get_edges();
}
