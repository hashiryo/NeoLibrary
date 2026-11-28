#pragma once
#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>
// 区間 DP dp[l][r] = min_{l <= k < r} (dp[l][k] + dp[k+1][r]) + w(l, r)、dp[i][i] = 0 を、
// Knuth の最適化で O(n^2) で解く。最適な k を opt[l][r] とすると、費用 w が
//   quadrangle inequality: w(a, c) + w(b, d) <= w(a, d) + w(b, c)  (a <= b <= c <= d)
//   単調性: a <= b <= c <= d なら w(b, c) <= w(a, d)
// を満たすとき opt[l][r-1] <= opt[l][r] <= opt[l+1][r] が成り立つので、各区间の k を
// その範囲に絞れる。返すのは dp[0][n-1]。w(l, r) の呼び出しは O(n^2) 回。
// dp は三角配列 (n(n+1)/2 個) に置き、opt は直前の長さの層だけ持つ。
// 記録は algo-notes の notes/dp-optimization.md にある。
template <class T, class F> T knuth_dp(int n, const F& w) {
 if(n <= 1) return 0;
 // off[l] は dp の行 l の先頭。dp[l][r] は off[l] + (r - l)。
 std::vector<int> off(n + 1);
 for(int l= 0; l < n; ++l) off[l + 1]= off[l] + (n - l);
 std::vector<T> dp((std::size_t)off[n]);
 // opt[l] は直前の長さの層の opt[l][l + len - 2] (len = 今の長さ)。
 std::vector<int> opt(n), nopt(n);
 for(int l= 0; l < n; ++l) opt[l]= l;
 for(int len= 2; len <= n; ++len) {
  for(int l= 0, r= len - 1; r < n; ++l, ++r) {
   const int lo= opt[l], hi= std::min(opt[l + 1], r - 1);
   const int base_l= off[l], base_lo= off[lo + 1];
   int bk= lo;
   T bv= dp[base_l + (lo - l)] + dp[base_lo + (r - lo - 1)];
   for(int k= lo + 1; k <= hi; ++k) {
    const T v= dp[base_l + (k - l)] + dp[off[k + 1] + (r - k - 1)];
    if(v < bv) bv= v, bk= k;
   }
   dp[base_l + (r - l)]= bv + w(l, r);
   nopt[l]= bk;
  }
  std::swap(opt, nopt);
 }
 return dp[n - 1];
}
