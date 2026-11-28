#pragma once
#include <algorithm>
#include <iterator>
#include <set>
#include <vector>
#include "neo/base/core/geo/Point.hpp"
// 点を足すだけの動的凸包。下側の鎖を std::set で持ち、点を足したときに要らなくなった点を左右から消す。
// 上側の鎖は、点を反転して同じ構造で持つ。where は +1: 内, 0: 周上, -1: 外。
template <class K> class IncrementalConvexHull {
 using P= Point<K>;
 struct Lower {
  std::set<P> S;
  template <class A, class R> void insert(const P& p, const A& ad, const R& rm) {
   if(where(p) >= 0) return;
   S.insert(p);
   std::vector<P> l, r;
   auto st= S.find(p);
   for(auto it= st; it != S.begin();) {
    if(--it; l.empty()) {
     l.emplace_back(*it);
     continue;
    }
    if(sgn(cross(*it - p, l.back() - p)) > 0) break;
    rm(*it, l.back()), l.emplace_back(*it);
   }
   for(auto it= st; ++it != S.end();) {
    if(r.empty()) {
     r.emplace_back(*it);
     continue;
    }
    if(sgn(cross(r.back() - p, *it - p)) > 0) break;
    rm(r.back(), *it), r.emplace_back(*it);
   }
   if(l.size() > 1) S.erase(std::next(S.find(l.back())), S.find(p));
   if(l.size()) ad(l.back(), p);
   if(r.size() > 1) S.erase(std::next(S.find(p)), S.find(r.back()));
   if(r.size()) ad(p, r.back());
   if(l.size() && r.size()) rm(l[0], r[0]);
  }
  int where(const P& p) const {
   auto r= S.lower_bound(p);
   if(r == S.begin()) return S.size() && *r == p ? 0 : -1;
   if(r == S.end()) return -1;
   return sgn(cross(*std::prev(r) - p, *r - p));
  }
 } L, U;
 std::size_t m;
 wide_t<K> a;
public:
 IncrementalConvexHull(): m(0), a(0) {}
 std::size_t edge_size() const { return m; }
 wide_t<K> area() const { return a / 2; }
 // 2 倍の面積。
 wide_t<K> area2() const { return a; }
 // +1: 内, 0: 周上, -1: 外
 int where(const P& p) const {
  int l= L.where(p), u= U.where(-p);
  return !l || !u ? 0 : std::min(l, u);
 }
 void insert(const P& p) {
  auto ad= [&](const P& q, const P& r) { a+= cross(q, r), ++m; };
  auto rm= [&](const P& q, const P& r) { a-= cross(q, r), --m; };
  L.insert(p, ad, rm), U.insert(-p, ad, rm);
 }
};
