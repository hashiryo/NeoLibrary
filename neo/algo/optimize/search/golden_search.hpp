#pragma once
#include <cassert>
#include <utility>
// [l, r] の実数の上で、単峰な f の最大 (MAX=true) か最小を黄金分割探索で求める。区間を iter 回
// 縮めるので、f の評価は 2 + iter 回。返すのは (x, f(x))。
// 今の Library の golden_search から、MinMaxEnum を bool MAX に変えたもの。
template <bool MAX= false, class F> auto golden_search(const F& f, long double l, long double r, int iter= 100) {
 using T= std::decay_t<decltype(f(l))>;
 static constexpr long double c= 0.38196601125;
 assert(l <= r);
 long double x= l + (r - l) * c, y= r - (r - l) * c;
 T fx= f(x), fy= f(y);
 for(bool g; iter--;) {
  if constexpr(MAX) g= fx > fy;
  else g= fx < fy;
  if(g) r= y, y= x, fy= fx, fx= f(x= l + (r - l) * c);
  else l= x, x= y, fx= fy, fy= f(y= r - (r - l) * c);
 }
 return std::pair<long double, T>{x, fx};
}
