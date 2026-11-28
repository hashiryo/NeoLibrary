#pragma once
#include <limits>
#include <vector>
// dp[i] = min_{j < i} (dp[j] + w(i, j)) を、w が Monge (行の最小の列が行について単調) のときに解く。
// dp[0] = 0 で、返すのは dp[0..n]。w の評価は O(n log n) 回。
// 行の dp を確定させる順を分割統治で決める (simplified LARSCH)。行 m の最良の列は、行 l と行 r の
// 最良の列の間に絞れる。先に m を確定させ、[l, m] を解いてから、r を [l+1, m] で引き直し、[m, r] を解く。
// 添字は rows = dp の添字 (0..n)、cols = 直前に選んだ列。w(i, j) は j < i でだけ呼ばれる。
// procon-judge の yukicoder 409、703〜705、AOJ 3086、2603 で使う今の Library の
// simplified_larsch_dp.hpp を写したもの。記録は algo-notes の notes/dp-optimization.md にある。
template <class T, class F> std::vector<T> larsch_dp(int n, const F& w) {
 std::vector<T> dp(n + 1, std::numeric_limits<T>::max());
 std::vector<int> best(n + 1);
 dp[0]= 0;
 auto relax= [&](int i, int j) {
  const T v= dp[j] + w(i, j);
  if(v < dp[i]) dp[i]= v, best[i]= j;
 };
 relax(n, 0);
 auto rec= [&](auto& rec, int l, int r) -> void {
  if(r - l <= 1) return;
  const int m= (l + r) / 2;
  for(int j= best[l]; j <= best[r]; ++j) relax(m, j);
  rec(rec, l, m);
  for(int j= l + 1; j <= m; ++j) relax(r, j);
  rec(rec, m, r);
 };
 rec(rec, 0, n);
 return dp;
}
