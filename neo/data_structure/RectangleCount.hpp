#pragma once
#include <utility>
#include <vector>
#include "neo/data_structure/WaveletMatrix.hpp"
#include "neo/sort/argsort.hpp"
// 静的な 2 次元の点の集合について、長方形 [xl, xr) × [yl, yr) に入る点の数を引く。点を x の順に並べ (neo/sort/argsort.hpp)、
// 並べた順の y の列から WaveletMatrix を作る。x の範囲から位置の範囲への変換は、y の値から順位への変換と同じ形の表 (値の範囲を
// 最大 65536 個の区間に分けた表と cmov の二分探索) で引く。procon-judge の aoj-2426 で書き比べた wm4 の形。
template <class T= int> class RectangleCount {
 wavelet_matrix_internal::Lookup<T> xs;  // 点の x の昇順 (重複を含む) と、x 未満の点の数を引く表
 WaveletMatrix<T> wm;                    // x の順に並べた点の y
public:
 RectangleCount(): RectangleCount(std::vector<T>(), std::vector<T>()) {}
 // 点 (x[i], y[i]) の集合から作る。
 RectangleCount(const std::vector<T>& x, const std::vector<T>& y) {
  const int n= x.size();
  const std::vector<int> o= argsort(x);
  std::vector<T> ys(n);
  xs.v.resize(n);
  for(int j= 0; j < n; ++j) xs.v[j]= x[o[j]], ys[j]= y[o[j]];
  xs.build();
  wm= WaveletMatrix<T>(std::move(ys));
 }
 int size() const { return wm.size(); }
 // x が [xl, xr)、y が [yl, yr) に入る点の数。xl >= xr や yl >= yr なら 0。
 [[gnu::always_inline]] int count(T xl, T xr, T yl, T yr) const {
  if(!(xl < xr)) return 0;
  return wm.count(xs.template rank<false>(xl), xs.template rank<false>(xr), yl, yr);
 }
};
