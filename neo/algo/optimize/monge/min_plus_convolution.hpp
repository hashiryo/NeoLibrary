#pragma once
#include <algorithm>
#include <vector>
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
