#pragma once
#include <array>
#include <cstdint>
#include <utility>
#include <vector>
#include "neo/base/core/geo/Point.hpp"
// 三角形の内部 (境界を除く) にある B の点の数を、A の 3 点の添字のクエリごとに返す。
// 各 (A_x, A_y) について、A_x -> A_y の左側にある B の点をビットセットで持ち、クエリでは三角形の 3 辺の
// ビットセットの積を数える。N, M が数百のときの形で、O(N^2 M + Q M / 64) 時間、O(N^2 M / 64) メモリ。
// N, M がもっと大きいときは、各 A_x で B を角度で並べて Fenwick で数える O(NM log(NM) + Q) の形にする。
template <class K> std::vector<int> count_points_in_triangle(const std::vector<Point<K>>& a, const std::vector<Point<K>>& b, const std::vector<std::array<int, 3>>& qs) {
 int n= a.size(), m= b.size(), w= (m + 63) >> 6;
 std::vector<std::uint64_t> left((std::size_t)n * n * w, 0);
 for(int x= 0; x < n; ++x)
  for(int y= 0; y < n; ++y) {
   if(x == y) continue;
   auto d= a[y] - a[x];
   std::uint64_t* row= &left[((std::size_t)x * n + y) * w];
   for(int j= 0; j < m; ++j)
    if(sgn(cross(d, b[j] - a[x])) > 0) row[j >> 6]|= std::uint64_t(1) << (j & 63);
  }
 std::vector<int> ans;
 ans.reserve(qs.size());
 for(auto& q: qs) {
  int x= q[0], y= q[1], z= q[2], cnt= 0;
  if(int s= sgn(cross(a[y] - a[x], a[z] - a[x])); s) {
   if(s < 0) std::swap(y, z);
   const std::uint64_t* p= &left[((std::size_t)x * n + y) * w];
   const std::uint64_t* r= &left[((std::size_t)y * n + z) * w];
   const std::uint64_t* t= &left[((std::size_t)z * n + x) * w];
   for(int i= 0; i < w; ++i) cnt+= __builtin_popcountll(p[i] & r[i] & t[i]);
  }
  ans.push_back(cnt);
 }
 return ans;
}
