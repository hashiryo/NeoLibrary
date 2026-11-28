---
title: 黄金分割探索 (実数)
documentation_of: ../../neo/algo/optimize/search/golden_search.hpp
---

区間 [l, r] の実数の上で、単峰な f の最大か最小を求めます。区間を `iter` 回縮めるので、f の評価は 2 + iter 回です。

```cpp
template <bool MAX = false, class F> auto golden_search(const F& f, long double l, long double r, int iter = 100);
```

`MAX = true` で最大、既定は最小です。返すのは (x, f(x)) の組です。

## 使い方

```cpp
auto [x, fx] = golden_search<false>(f, -1000.0L, 1000.0L);
```

整数の上で探すときは [fibonacci_search](fibonacci_search.md) を使います。
