#pragma once
#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>
#include "neo/base/core/geo/Point.hpp"
// 最近点対。y で並べた列を分割統治し、統合のときに x の差が best 未満の点だけを見る。
template <class K> std::pair<Point<K>, Point<K>> closest_pair(std::vector<Point<K>> ps) {
 int n= ps.size();
 assert(n >= 2);
 std::sort(ps.begin(), ps.end(), [](const Point<K>& p, const Point<K>& q) { return p.y < q.y; });
 std::vector<Point<K>> memo(n);
 std::pair<Point<K>, Point<K>> ret= {ps[0], ps[1]};
 auto best= dist2(ps[0], ps[1]);
 auto rec= [&](auto& rec, int l, int r) -> void {
  if(r - l == 1) return;
  int m= (l + r) / 2;
  K y= ps[m].y, d;
  wide_t<K> tmp;
  rec(rec, l, m), rec(rec, m, r), std::inplace_merge(ps.begin() + l, ps.begin() + m, ps.begin() + r);
  for(int i= l, cnt= 0; i < r; ++i) {
   if(d= ps[i].y - y; wide_t<K>(d) * d >= best) continue;
   for(int j= cnt; j--;) {
    if(d= ps[i].x - memo[j].x, tmp= wide_t<K>(d) * d; tmp >= best) break;
    if(d= ps[i].y - memo[j].y, tmp+= wide_t<K>(d) * d; best > tmp) best= tmp, ret= {ps[i], memo[j]};
   }
   memo[cnt++]= ps[i];
  }
 };
 return rec(rec, 0, n), ret;
}
