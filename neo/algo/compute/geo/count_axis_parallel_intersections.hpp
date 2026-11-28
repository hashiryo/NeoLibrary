#pragma once
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <vector>
#include "neo/base/core/geo/Segment.hpp"
#include "neo/base/ds/range/Fenwick.hpp"
// 軸平行な線分の集合で、縦の線分と横の線分が交わる組の数を数える。x の小さい順に、横の線分の y を
// 左端で Fenwick に足し、右端で引き、縦の線分の x で y の範囲の和を取る。同じ x では足すのが先で、
// 縦の問い合わせのあとに引くので、端点で接する組も 1 回だけ数える。
// 同じ直線上の線分が端点を共有すると、1 つの交点を 2 組が数えるので、その形の入力には使えない。
namespace count_axis_parallel_intersections_internal {
struct Sum {
 using T= long long;
 static constexpr T ti() { return 0; }
 static constexpr T op(T a, T b) { return a + b; }
};
}
template <class K> long long count_axis_parallel_intersections(const std::vector<Segment<K>>& ss) {
 struct Ev {
  K x;
  int type;  // 0: 横の左端, 1: 縦, 2: 横の右端
  K a, b;    // 縦なら y の範囲 [a, b]、横なら a = b = y
 };
 std::vector<Ev> ev;
 std::vector<K> ys;
 ev.reserve(ss.size() * 2), ys.reserve(ss.size());
 for(const auto& s: ss) {
  assert(s.p.x == s.q.x || s.p.y == s.q.y);
  if(s.p.y == s.q.y) {
   K l= std::min(s.p.x, s.q.x), r= std::max(s.p.x, s.q.x);
   ev.push_back({l, 0, s.p.y, s.p.y}), ev.push_back({r, 2, s.p.y, s.p.y});
   ys.push_back(s.p.y);
  } else {
   K l= std::min(s.p.y, s.q.y), r= std::max(s.p.y, s.q.y);
   ev.push_back({s.p.x, 1, l, r});
  }
 }
 std::sort(ev.begin(), ev.end(), [](const Ev& a, const Ev& b) { return a.x != b.x ? a.x < b.x : a.type < b.type; });
 std::sort(ys.begin(), ys.end()), ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
 Fenwick<count_axis_parallel_intersections_internal::Sum> fw(ys.size());
 long long ans= 0;
 for(const auto& e: ev) {
  if(e.type == 0) fw.add(int(std::lower_bound(ys.begin(), ys.end(), e.a) - ys.begin()), 1);
  else if(e.type == 2) fw.add(int(std::lower_bound(ys.begin(), ys.end(), e.a) - ys.begin()), -1);
  else {
   int l= int(std::lower_bound(ys.begin(), ys.end(), e.a) - ys.begin());
   int r= int(std::upper_bound(ys.begin(), ys.end(), e.b) - ys.begin());
   ans+= fw.prod(r) - fw.prod(l);
  }
 }
 return ans;
}
