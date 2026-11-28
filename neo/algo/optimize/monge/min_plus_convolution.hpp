#pragma once
#include <algorithm>
#include <limits>
#include <vector>
#include "neo/algo/optimize/monge/monotone_minima.hpp"
// 凸な列どうしの min-plus convolution。c[k] = min_{i+j=k} a[i] + b[j]。
// 差分の列 a[i+1]-a[i] と b[j+1]-b[j] はどちらも広義単調増加なので、2 つを併合した列の前からの累積和が
// c の隣り合う差になる。O(n + m)。
template <class T> std::vector<T> min_plus_convolution_convex_convex(const std::vector<T>& a, const std::vector<T>& b) {
 const int n= (int)a.size(), m= (int)b.size();
 std::vector<T> c;
 if(n == 0 || m == 0) return c;
 c.resize(n + m - 1), c[0]= a[0] + b[0];
 int i= 0, j= 0, k= 1;
 while(i + 1 < n && j + 1 < m) {
  if(a[i + 1] - a[i] <= b[j + 1] - b[j]) c[k]= c[k - 1] + (a[i + 1] - a[i]), ++i;
  else c[k]= c[k - 1] + (b[j + 1] - b[j]), ++j;
  ++k;
 }
 while(i + 1 < n) c[k]= c[k - 1] + (a[i + 1] - a[i]), ++i, ++k;
 while(j + 1 < m) c[k]= c[k - 1] + (b[j + 1] - b[j]), ++j, ++k;
 return c;
}
// a が凸で b が任意の min-plus convolution。c[k] = min_{i+j=k} a[i] + b[j]。
// 行 k 列 j の値を a[k-j] + b[j] にすると (j > k や k-j >= n は最大値)、行の最小の列が行について
// 単調になるので、分割統治 (monotone_minima) で作れる。SMAWK でも作れるが、計測では 1 候補 1 評価の
// 分割統治のほうが x64 の gcc と clang の両方で 1.6 倍速かった (algo-notes の notes/dp-optimization.md)。
// ここでは monotone_minima に値を渡さず、行 k の有効な列 [max(0, k-(n-1)), min(m-1, k)] を
// 分割統治の候補の範囲に切り落とす。有効でない (k, j) を評価しないぶん、値の関数を渡す形より
// x64-gcc で 3 割ほど速い。a[i] + b[j] が T に収まること。
template <class T> std::vector<T> min_plus_convolution_convex(const std::vector<T>& a, const std::vector<T>& b) {
 const int n= (int)a.size(), m= (int)b.size();
 std::vector<T> c;
 if(n == 0 || m == 0) return c;
 const int h= n + m - 1;
 c.resize(h);
 auto rec= [&](auto& rec, int h1, int h2, int w1, int w2) -> void {
  if(h1 == h2) return;
  const int hm= (h1 + h2) / 2;
  const int wl= std::max(w1, std::max(0, hm - (n - 1))), wr= std::min(w2, std::min(m - 1, hm) + 1);
  int bj= wl;
  T bv= a[hm - wl] + b[wl];
  for(int j= wl + 1; j < wr; ++j) {
   const T v= a[hm - j] + b[j];
   if(v < bv) bv= v, bj= j;
  }
  c[hm]= bv;
  rec(rec, h1, hm, w1, bj + 1);
  rec(rec, hm + 1, h2, bj, w2);
 };
 rec(rec, 0, h, 0, m);
 return c;
}
// a が凹で b が任意の min-plus convolution。c[k] = min_{i+j=k} a[i] + b[j]。
// 行 k 列 j の有効な列は [max(0, k-(n-1)), min(m-1, k)] の階段になる。この行列は反 Monge なので、行の
// 最小の列は階段の上では行について単調にならない (a = [-19,-16,-13,-11,-10,-14]、b = [-3,-7,-6,-3] で、
// 行 4 の最小は列 2 か 3、行 5 は 0 か 1 か 2、行 6 は 1)。そこで、階段を長いほうの辺で半分に切ることを
// 繰り返し、切り口の階段の端で子の範囲を狭める。範囲が全部有効になった矩形では、列を逆順にすると
// 反 Monge が Monge になるので、monotone_minima で各行の最小を引ける。Library Checker の
// min_plus_convolution_concave_arbitrary の想定解 (ac_func.cpp) と同じ形。a[i] + b[j] が T に収まること。
template <class T> std::vector<T> min_plus_convolution_concave(const std::vector<T>& a, const std::vector<T>& b) {
 const int n= (int)a.size(), m= (int)b.size();
 std::vector<T> c;
 if(n == 0 || m == 0) return c;
 const int h= n + m - 1;
 // 行 x の有効な列 [y_min[x], y_max[x]] と、列 y の有効な行 [x_min[y], x_max[y]]。
 std::vector<int> y_min(h), y_max(h), x_min(m), x_max(m);
 for(int x= 0; x < h; ++x) y_min[x]= std::max(0, x - (n - 1)), y_max[x]= std::min(m - 1, x);
 for(int y= 0; y < m; ++y) x_min[y]= y, x_max[y]= y + n - 1;
 c.assign(h, std::numeric_limits<T>::max());
 auto rf= [&](auto& rf, int x1, int x2, int y1, int y2) -> void {
  // 矩形が全部有効なら、列を逆順にして Monge にし、monotone_minima で引く。
  if(y_max[x1] >= y2 && y1 >= y_min[x2]) {
   auto value= [&](int i, int jj) -> T {
    const int y= y2 - jj;
    return a[(x1 + i) - y] + b[y];
   };
   auto best= monotone_minima<T>(x2 - x1 + 1, y2 - y1 + 1, value);
   for(int i= 0, ie= x2 - x1; i <= ie; ++i)
    if(const T v= value(i, best[i]); v < c[x1 + i]) c[x1 + i]= v;
   return;
  }
  // 小さい矩形は素朴に見る。
  if((long long)(x2 - x1) * (y2 - y1) < 1000) {
   for(int x= x1; x <= x2; ++x)
    for(int y= std::max(y_min[x], y1), ye= std::min(y_max[x], y2); y <= ye; ++y)
     if(const T v= a[x - y] + b[y]; v < c[x]) c[x]= v;
   return;
  }
  // 長いほうの辺で半分に切る。切り口の階段の端で、子の範囲を狭める。
  if(x2 - x1 > y2 - y1) {
   const int xm= (x1 + x2) / 2;
   const int ny2= std::min(y_max[xm], y2);
   if(y1 <= ny2) rf(rf, x1, xm, y1, ny2);
   const int ny1= std::max(y_min[xm], y1);
   if(ny1 <= y2) rf(rf, xm + 1, x2, ny1, y2);
  } else {
   const int ym= (y1 + y2) / 2;
   const int nx2= std::min(x_max[ym], x2);
   if(x1 <= nx2) rf(rf, x1, nx2, y1, ym);
   const int nx1= std::max(x_min[ym], x1);
   if(nx1 <= x2) rf(rf, nx1, x2, ym + 1, y2);
  }
 };
 rf(rf, 0, h - 1, 0, m - 1);
 return c;
}
