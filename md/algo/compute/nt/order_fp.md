---
title: 乗法の位数 (素数 mod)
documentation_of: ../../neo/algo/compute/nt/order_fp.hpp
---

素数 $p$ を法とする元の乗法の位数を返します。

```cpp
class OrderFp {
 explicit OrderFp(unsigned long long p);  // p は素数、2 <= p < 2^64
 unsigned long long operator()(unsigned long long x) const;  // x^k ≡ 1 となる最小の k。x ≡ 0 なら 0
 unsigned long long primitive_root() const;
};
```

$p - 1$ を `factorize_pairs` で分解し、位数の候補 $p - 1$ から素因数 $q$ を 1 つずつ落とせるか試します。$x^{(p-1)/q} \equiv 1$ なら位数を $q$ で割れるので、割れなくなるまで繰り返します。累乗は Montgomery で、$p < 2^{64}$ でも桁あふれしないように和を 128 bit で取ります。`primitive_root()` は位数 $p - 1$ の元を小さい方から探します。

前計算は $p - 1$ の因数分解で、1 回の問いは素因数の指数の和だけ累乗します ($O(\log^2 p)$ ほど)。今の Library の `OrderFp` は constexpr で、位数を素因数の積の二分木でまとめて求めます。ここでは `factorize` を使う実行時の形にしました。

Library Checker に位数そのものの問題が無いため、順位表では測っていません。procon-judge の `self-nt-order-fp` で、$p < 500$ の全元を素朴な掛け算と比べ、64 bit の素数では位数が $p - 1$ を割ること、$x^{位数} \equiv 1$、$g^k$ の位数が $(p-1)/\gcd(p-1,k)$ になることを確かめています (手元の M2 Max で 76 ms)。
