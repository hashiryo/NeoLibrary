---
title: LARSCH (オンラインの行最小)
documentation_of: ../../neo/algo/optimize/monge/larsch.hpp
---

`dp[i] = min_{j < i} (dp[j] + w(i, j))` を、w が Monge (行の最小の列が行について単調) のときに解きます。`dp[0] = 0` で、返すのは `dp[0..n]` です。w の評価は O(n log n) 回です。

```cpp
template <class T, class F> std::vector<T> larsch_dp(int n, const F& w);
```

`w(i, j)` は `j < i` でだけ呼ばれます。dp の行を確定させる順を分割統治で決めます。行 m の最良の列は、行 l と行 r の最良の列の間に絞れるので、先に m を確定させ、[l, m] を解いてから、r を [l+1, m] で引き直し、[m, r] を解きます。

## 使い方

```cpp
// dp[i] = min_{j < i} dp[j] + (S_i - S_j)^2 (S は広義単調増加)
auto w = [&](int i, int j) {
  const long long d = s[i] - s[j];
  return d * d;
};
std::vector<long long> dp = larsch_dp<long long>(n, w);
```

行を 1 つずつ確定させるので、w が dp[j] を使う形でも解けます。行 i の dp を全部の j について引く [分割統治](monotone_minima.md) は、行の順に dp が決まっている必要があります。

## 中身

今の Library の simplified_larsch_dp を写したものです。procon-judge の yukicoder 409、703〜705、AOJ 3086、2603 で書き比べた記録は、algo-notes の notes/dp-optimization.md にあります。
