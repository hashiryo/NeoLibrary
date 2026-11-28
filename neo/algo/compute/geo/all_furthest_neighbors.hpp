#pragma once
#include <vector>
#include "neo/base/core/geo/Point.hpp"
// 狭義凸多角形の頂点 p_0, ..., p_{n-1} (反時計回り) を受け、各 i について最も遠い頂点 f(i) を返す。
// 距離の行列 D(i, j) = |p_i - p_{j mod n}|^2 (i < j < i + n) は単調なので、行ごとの最大を分割統治で求める。
// 単純な回転カリパスは、頂点からの距離が境界に沿って単峰でない多角形で誤る (Library Checker の non_unimodal)。
namespace all_furthest_neighbors_internal {
// 各行の最大の位置が単調な行列で、行ごとの最大 (いちばん左) を求める。
template <class F> std::vector<int> monotone_maxima(int h, int w, const F& value) {
 std::vector<int> max_col(h);
 auto dfs= [&](auto&& self, int row_l, int row_r, int col_l, int col_r) -> void {
  if(row_l == row_r) return;
  int row= (row_l + row_r) / 2, best= col_l;
  for(int col= col_l + 1; col < col_r; ++col)
   if(value(row, best) < value(row, col)) best= col;
  max_col[row]= best;
  self(self, row_l, row, col_l, best + 1);
  self(self, row + 1, row_r, best, col_r);
 };
 return dfs(dfs, 0, h, 0, w), max_col;
}
}
// ps は狭義凸多角形の頂点を反時計回りに並べたもの。n >= 3 を前提にする。
template <class K> std::vector<int> all_furthest_neighbors(const std::vector<Point<K>>& ps) {
 int n= ps.size();
 auto value= [&](int i, int j) -> wide_t<K> {
  if(j <= i) return j - i;
  if(j >= i + n) return -1;
  return dist2(ps[i], ps[j % n]);
 };
 std::vector<int> ans= all_furthest_neighbors_internal::monotone_maxima(n, 2 * n - 1, value);
 for(auto& x: ans) x%= n;
 return ans;
}
