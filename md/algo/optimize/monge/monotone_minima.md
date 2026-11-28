---
title: 単調な行最小 (分割統治)
documentation_of: ../../neo/algo/optimize/monge/monotone_minima.hpp
---

行の最小の列が行について単調な行列の、各行の最小の列を分割統治で求めます。`value` の評価は O(h log w) 回です。

```cpp
template <class T, class F> std::vector<int> monotone_minima(int h, int w, const F& value);
template <class T, class F, class Lo, class Hi> std::vector<int> monotone_minima_bounded(int h, const F& value, const Lo& lo, const Hi& hi);
template <class T, class F, class Lo, class Hi> std::vector<T> monotone_minima_bounded_values(int h, const F& value, const Lo& lo, const Hi& hi);
```

`value(r, c)` は行 r 列 c の値を `T` で返します。有効でない (r, c) は `T` の最大値を返してください。行の最小の列が行について単調になるのは、行列が totally monotone のとき (Monge や逆 Monge を含む) です。

## 使い方

```cpp
// c[k] = min_{i+j=k} a[i] + b[j] (a が凸、b は任意)
const long long INF = std::numeric_limits<long long>::max() / 4;
auto r = monotone_minima<long long>(n + m - 1, m, [&](int k, int j) {
  const unsigned i = (unsigned)(k - j);
  return i < (unsigned)n ? a[i] + b[j] : INF;
});
for (int k = 0; k < n + m - 1; ++k) c[k] = a[k - r[k]] + b[r[k]];
```

## 行の有効な範囲が階段になるとき

`monotone_minima_bounded` は、行 r の有効な列が lo(r) 以上 hi(r) 未満のときに使います。lo と hi は広義単調増加で、どの行も lo(r) < hi(r) です。有効でないセルを評価しないので、値の関数の中で範囲を判定して `T` の最大値を返す形より速いです。

`monotone_minima_bounded_values` は同じ計算で、各行の最小の値だけを返します。呼ぶ側が添字から値を引き直さないので、dp の各層のように最小の値がそのまま次の層に要るときはこちらが速いです (procon-judge の cses-2086 で x64-gcc 56.4 対 59.2 ms)。

```cpp
// dp_t[i] = min_{j < i} dp_{t-1}[j] + (S_i - S_j)^2 の 1 層。行 r は i = t + r。
auto value = [&](int r, int j) { const int i = t + r; const long long d = s[i] - s[j]; return dp[j] + d * d; };
auto lo = [&](int) { return t - 1; };
auto hi = [&](int r) { return t + r; };
std::vector<long long> next = monotone_minima_bounded_values<long long>(n + 1 - t, value, lo, hi);
```

## 中身

行の真ん中で列を走査して最良の列を選び、上の行はその列まで、下の行はその列から、と範囲を狭めて再帰します。走査では最良の値を覚え、候補ごとに `value` を 1 回だけ評価します。比較の述語 (`select(r, c1, c2)`) を受ける形だと 1 候補あたり 2 回評価になり、計測で 1.8 倍遅かったため、値の関数を受ける形にしています。

procon-judge の yosupo-min-plus-convolution-convex-arbitrary で書き比べた記録は、algo-notes の notes/dp-optimization.md にあります。
