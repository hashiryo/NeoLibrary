---
title: Fibonacci 探索 (整数)
documentation_of: ../../neo/algo/optimize/search/fibonacci_search.hpp
---

区間 [l, r] の整数の上で、単峰な f の最大か最小を求めます。f の評価は O(log(r - l)) 回です。

```cpp
template <bool MAX = false, class F> auto fibonacci_search(const F& f, long long l, long long r);
```

`MAX = true` で最大、既定は最小です。返すのは (x, f(x)) の組です。単峰でない f では答えは保証されません。

## 使い方

```cpp
// 凹な f の最大を [0, 2e9] で探す
auto [x, fx] = fibonacci_search<true>(f, 0, 2000000000LL);
```

実数の上で探すときは [golden_search](golden_search.md) を使います。
