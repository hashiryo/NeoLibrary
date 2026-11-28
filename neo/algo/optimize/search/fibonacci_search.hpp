#pragma once
#include <cassert>
#include <utility>
// [l, r] の整数の上で、単峰な f の最大 (MAX=true) か最小を Fibonacci 探索で求める。f の評価は
// O(log(r - l)) 回。返すのは (x, f(x))。単峰でない f では答えは保証されない。
// 今の Library の fibonacci_search から、MinMaxEnum を bool MAX に変えたもの。
// 記録は algo-notes の notes/dp-optimization.md にある。
template <bool MAX= false, class F> auto fibonacci_search(const F& f, long long l, long long r) {
 using T= std::decay_t<decltype(f(l))>;
 assert(l <= r);
 long long s= 1, t= 2, a= l - 1, x, b, y;
 for(long long e= r - l + 2; t < e;) std::swap(s+= t, t);
 b= a + t, x= b - s;
 T fx= f(x), fy;
 for(bool g; a + b != 2 * x;) {
  if(y= a + b - x; r < y) b= a, a= y;
  else {
   if constexpr(MAX) g= fx > (fy= f(y));
   else g= fx < (fy= f(y));
   if(g) b= a, a= y;
   else a= x, x= y, fx= fy;
  }
 }
 return std::pair<long long, T>{x, fx};
}
