#pragma once
#include <algorithm>
#include <cassert>
#include <numeric>
#include <vector>
#include "neo/base/core/geo/Point.hpp"
// 凸包を剥がして、各点が何層目で消えるかを返す (凸の層、onion decomposition)。
// 点を y の順に並べ、区間ごとの左側の凸包の鎖を持つセグメント木を作る。点を消すと、根までの節の鎖を
// 共通接線で繋ぎ直す (Overmars–van Leeuwen の削除つき動的凸包)。右側の鎖は、点を反転して同じ木で持つ。
// 1 層ごとに両方の木の根の鎖を集め、その点を消す。O(n log^2 n)。
// 点はすべて異なることを前提にする。
namespace convex_layers_internal {
template <class K> struct LeftHull {
 using P= Point<K>;
 using W= wide_t<K>;
 const std::vector<P>& ps;
 struct Node {
  int bl, br, L, R, lchd, rchd;
 };
 std::vector<Node> nodes;
 int root;
 bool isleaf(int w) const { return nodes[w].lchd == -1 && nodes[w].rchd == -1; }
 W cross3(int a, int b, int c) const { return cross(ps[b] - ps[a], ps[c] - ps[a]); }
 void pull(int w) {
  int l= nodes[w].lchd, r= nodes[w].rchd;
  W split_y= ps[nodes[r].L].y;
  while(!isleaf(l) || !isleaf(r)) {
   int a= nodes[l].bl, b= nodes[l].br, c= nodes[r].bl, d= nodes[r].br;
   if(a != b && cross3(a, b, c) > 0) l= nodes[l].lchd;
   else if(c != d && cross3(b, c, d) > 0) r= nodes[r].rchd;
   else if(a == b) r= nodes[r].lchd;
   else if(c == d) l= nodes[l].rchd;
   else {
    W s1= cross3(a, b, c), s2= cross3(b, a, d);
    assert(s1 + s2 >= 0);
    if(s1 + s2 == 0 || s1 * ps[d].y + s2 * ps[c].y < split_y * (s1 + s2)) l= nodes[l].rchd;
    else r= nodes[r].lchd;
   }
  }
  nodes[w].bl= nodes[l].L;
  nodes[w].br= nodes[r].L;
 }
 void build(int w, int L, int R) {
  nodes[w].L= L, nodes[w].R= R;
  if(R - L == 1) {
   nodes[w].lchd= nodes[w].rchd= -1;
   nodes[w].bl= nodes[w].br= L;
  } else {
   int M= (L + R) / 2;
   nodes[w].lchd= w + 1;
   nodes[w].rchd= w + 2 * (M - L);
   build(nodes[w].lchd, L, M);
   build(nodes[w].rchd, M, R);
   pull(w);
  }
 }
 int erase(int w, int L, int R) {
  if(R <= nodes[w].L || L >= nodes[w].R) return w;
  if(L <= nodes[w].L && R >= nodes[w].R) return -1;
  nodes[w].lchd= erase(nodes[w].lchd, L, R);
  nodes[w].rchd= erase(nodes[w].rchd, L, R);
  if(nodes[w].lchd == -1) return nodes[w].rchd;
  if(nodes[w].rchd == -1) return nodes[w].lchd;
  pull(w);
  return w;
 }
 void get_hull(int w, int l, int r, std::vector<int>& res) const {
  if(isleaf(w)) res.push_back(nodes[w].L);
  else if(r <= nodes[w].bl) get_hull(nodes[w].lchd, l, r, res);
  else if(l >= nodes[w].br) get_hull(nodes[w].rchd, l, r, res);
  else {
   get_hull(nodes[w].lchd, l, nodes[w].bl, res);
   get_hull(nodes[w].rchd, nodes[w].br, r, res);
  }
 }
 explicit LeftHull(const std::vector<P>& ps): ps(ps), nodes(ps.size() * 2), root(0) { build(0, 0, ps.size()); }
 std::vector<int> get_hull() const {
  if(root == -1) return {};
  std::vector<int> res;
  get_hull(root, 0, (int)ps.size() - 1, res);
  return res;
 }
 void erase(int L) { root= erase(root, L, L + 1); }
};
}
// ps はすべて異なる点。返り値の i 番目は、ps[i] が消える層 (1 始まり)。
template <class K> std::vector<int> convex_layers(const std::vector<Point<K>>& ps) {
 int n= ps.size();
 std::vector<int> ord(n);
 std::iota(ord.begin(), ord.end(), 0);
 std::sort(ord.begin(), ord.end(), [&](int i, int j) { return ps[i].y != ps[j].y ? ps[i].y < ps[j].y : ps[i].x < ps[j].x; });
 std::vector<Point<K>> sp(n);
 for(int i= 0; i < n; ++i) sp[i]= ps[ord[i]];
 convex_layers_internal::LeftHull<K> left(sp);
 std::vector<Point<K>> rp(n);
 for(int i= 0; i < n; ++i) rp[i]= -sp[n - 1 - i];
 convex_layers_internal::LeftHull<K> right(rp);
 std::vector<int> layer(n, 0);
 int cnt= 0;
 for(int l= 1; cnt < n; ++l) {
  std::vector<int> hull;
  for(int i: left.get_hull()) hull.push_back(i);
  for(int i: right.get_hull()) hull.push_back(n - 1 - i);
  std::sort(hull.begin(), hull.end());
  hull.erase(std::unique(hull.begin(), hull.end()), hull.end());
  for(int i: hull) {
   if(layer[i]) continue;
   layer[i]= l, ++cnt;
   left.erase(i), right.erase(n - 1 - i);
  }
 }
 std::vector<int> ans(n);
 for(int i= 0; i < n; ++i) ans[ord[i]]= layer[i];
 return ans;
}
