#pragma once
#include <utility>
#include <vector>
#include "neo/data_structure/WaveletMatrixSum.hpp"
#include "neo/sort/argsort.hpp"
// 静的な 2 次元の重み付きの点の集合について、長方形 [xl, xr) × [yl, yr) に入る点の重みの和を引く。点を x の順に並べ
// (neo/sort/argsort.hpp)、並べた順の y の列と重みの列から WaveletMatrixSum を作る。x の範囲から位置の範囲への変換は、y の値から
// 順位への変換と同じ形の表 (値の範囲を最大 65536 個の区間に分けた表と cmov の二分探索) で引く。procon-judge の
// yosupo-rectangle-sum で書き比べた wm_sep_brtbl の形。
template <class T= int, class W= long long> class RectangleSum {
 wavelet_matrix_internal::Lookup<T> xs;  // 点の x の昇順 (重複を含む) と、x 未満の点の数を引く表
 WaveletMatrixSum<T, W> wm;              // x の順に並べた点の y と重み
public:
 RectangleSum(): RectangleSum(std::vector<T>(), std::vector<T>(), std::vector<W>()) {}
 // 重み w[i] の点 (x[i], y[i]) の集合から作る。
 RectangleSum(const std::vector<T>& x, const std::vector<T>& y, const std::vector<W>& w) {
  const int n= x.size();
  const std::vector<int> o= argsort(x);
  std::vector<T> ys(n);
  std::vector<W> ws(n);
  xs.v.resize(n);
  for(int j= 0; j < n; ++j) xs.v[j]= x[o[j]], ys[j]= y[o[j]], ws[j]= w[o[j]];
  xs.build();
  wm= WaveletMatrixSum<T, W>(std::move(ys), std::move(ws));
 }
 int size() const { return wm.size(); }
 // x が [xl, xr)、y が [yl, yr) に入る点の重みの和。xl >= xr や yl >= yr なら 0。
 [[gnu::always_inline]] W sum(T xl, T xr, T yl, T yr) const {
  if(!(xl < xr)) return W();
  return wm.count_sum(xs.template rank<false>(xl), xs.template rank<false>(xr), yl, yr).second;
 }
};
