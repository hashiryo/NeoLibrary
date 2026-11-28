---
title: 階乗と二項係数 (素数 mod)
documentation_of: ../../neo/base/core/nt/factorial.hpp
---

素数 $p$ を法とする階乗と、その逆元の表を持ちます。

```cpp
template <class Z> struct Factorial {
 explicit Factorial(int n = 0);
 void ensure(int n);
 Z fact(int n);
 Z finv(int n);
 Z inv(int n);  // 1 <= n
 Z nPr(int n, int r);
 Z nCr(int n, int r);
 Z nHr(int n, int r);
};
```

`Z` は `ZMod` です。法が素数で、$0 \le n < p$ の範囲で使います。$n \ge p$ では階乗が 0 になり、`nCr` は Lucas の定理が要ります。`r < 0` と `n < r` では 0 を返します。

表は `fac` と `ifac` の 2 本です。`ifac[n] = fac[n]^{-1}` を 1 回の逆元で作り、下向きの掛け算で `ifac[0]` まで埋めます。`ensure` は足りない分だけを延ばし、前の表はそのまま使います。今の Library の `FactorialPrecalculation` は `iv`、`fct`、`fiv` の 3 本を持ち、`fiv` を 1 つずつ `inv(k)` を掛けて作ります。2 本にするとメモリが 3 分の 2 になり、`fiv` の作り方も掛け算だけになります。

`ensure(n)` は $O(n)$、各問いは $O(1)$ です。

procon-judge の `yosupo-binomial-coefficient-prime-mod` ($T = 10^6$、$n \le 10^7$) の手元の Apple M2 Max の最大ケースで、新しい実装は 3.40 s、87 MB、今の Library は 4.33 s、281 MB でした (raw のプロセス時間とピーク RSS)。CI の値はまだ取っていません。
