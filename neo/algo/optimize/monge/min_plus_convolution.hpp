#pragma once
#include <limits>
#include <vector>
#include "neo/algo/optimize/monge/smawk.hpp"
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
// 単調になるので、SMAWK で O(n + m) 回の評価で作れる。分割統治 (monotone_minima) でも作れ、そちらが
// 速いコンパイラもある。a[i] + b[j] が T に収まること。
template <class T> std::vector<T> min_plus_convolution_convex(const std::vector<T>& a, const std::vector<T>& b) {
 const int n= (int)a.size(), m= (int)b.size();
 std::vector<T> c;
 if(n == 0 || m == 0) return c;
 const T inf= std::numeric_limits<T>::max() / 4;
 auto value= [&](int k, int j) -> T {
  const unsigned i= (unsigned)(k - j);
  return i < (unsigned)n ? a[i] + b[j] : inf;
 };
 const int h= n + m - 1;
 auto r= smawk<T>(h, m, value);
 c.resize(h);
 for(int k= 0; k < h; ++k) c[k]= a[k - r[k]] + b[r[k]];
 return c;
}
