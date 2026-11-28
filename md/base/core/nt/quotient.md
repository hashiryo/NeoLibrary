---
title: 商の列挙と floor(N/i) の配列
documentation_of: ../../neo/base/core/nt/quotient.hpp
---

$\lfloor N/i \rfloor$ が同じになる $i$ の区間を列挙する関数と、$F(n) = \sum_{i=1}^n f(i)$ を $n \in S = \{\lfloor N/x \rfloor \mid 1 \le x \le N\}$ の所だけ持つ配列です。

```cpp
std::vector<std::tuple<unsigned long long, unsigned long long, unsigned long long>>
enumerate_quotients(unsigned long long N);
// (q, l, r): i in (l, r] のとき floor(N/i) = q

template <class T> struct CumSumQuotient {
 explicit CumSumQuotient(unsigned long long N);
 T& operator[](unsigned long long n);  // n in S
 T operator()(unsigned long long n) const;
 void add(unsigned long long i, T v);  // f(i) <- f(i) + v
 T sum(unsigned long long n) const;    // F(n)
 T sum() const;                        // F(N)
 // +=、-=、単項 -、スカラー倍、2 つの + と -
};
```

`enumerate_quotients` は区間を $q$ の昇順で返します。$N = 10$ なら `(1, 5, 10)`、`(2, 3, 5)`、`(3, 2, 3)`、`(5, 1, 2)`、`(10, 0, 1)` の 5 つです。区間の数は $\lfloor 2\sqrt{N} \rfloor$ ほどで、計算量は $O(\sqrt{N})$ です。

`CumSumQuotient` は $n \le K = \lfloor\sqrt N\rfloor$ を `X[n]`、$n > K$ を `X[K + N/n]` に置きます。`add(i, v)` は $F(n)$ を $n \ge i$ の所へ足します。今の Library の `add` は $i > K + N/i$ のときに何もしないため、$i$ が $N/j$ の形の大きい値だと $F(N)$ が変わらない不具合がありました。ここでは `X[K + j]` ($j \le N/i$) と `X[i..K]` の両方を足す形にして、どの $i$ でも正しく動きます。`sum()` も $N = 1$ で `X[K + 1]` が $F(N)$ と違う位置になるので、$N$ の位置を引きます。

$N = 10^{12}$ で 1999999 個の区間が出ます。`add` は $O(\sqrt N)$、$n$ を引くのは $O(1)$、2 つの配列の加減算は $O(\sqrt N)$ です。

procon-judge の `yosupo-enumerate-quotients` ($N = 10^{12}$) の手元の Apple M2 Max で、新しい実装は最大 300 ms、95 MB、今の Library は 368 ms、98 MB でした。新しい実装の `CumSumQuotient` は `self-nt-cum-sum-quotient` で、$N \le 200$ の素朴な累積和と、2 つの配列の加減算を突き合わせています。
