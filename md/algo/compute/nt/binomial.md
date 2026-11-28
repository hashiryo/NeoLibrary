---
title: 二項係数 (合成数の法)
documentation_of: ../../neo/algo/compute/nt/binomial.hpp
---

法 $m \le 10^6$ (合成数でもよい) の二項係数を返します。

```cpp
class BinomialCoefficient {
 explicit BinomialCoefficient(int m, std::size_t pre_size = 1 << 14);
 int nCr(long long n, long long r) const;  // 0 <= n <= 1e18、0 <= r
 int nHr(long long n, long long r) const;
};
```

$m = \prod p^e$ と分け、素数べき $p^e$ ごとに Granville の方法で $C(n, r) \bmod p^e$ を求めて CRT で束ねます。$p$ を取り除いた階乗の表 `facts` と、$p$ 進の桁和の表 `ds` を持ち、$n$ を $pp = p^{e+k}$ で割りながら桁を見ます。$pp$ は `pre_size` 以下の最大の $p$ の冪 ($p^e$ の倍数) です。

法が素数 1 つのときは `pre_size` を $2^{20}$ に、それ以外は $2^{14}$ にします。前計算は $m$ に対して $O(m)$ ほど、1 問は $O(\log_m n)$ です。$m = 1$ はすべて 0、$n < r$ も 0 です。

今の Library の `BinomialCoefficient` と同じ算法で、因数分解に `factorize_pairs`、逆元に `inv_mod` を使います。`yosupo-binomial-coefficient` ($T = 2 \times 10^5$、$n \le 10^{18}$、$m \le 10^6$) の手元の Apple M2 Max で、新しい実装は最大 1363 ms、9.3 MB、今の Library は 1370 ms、9.3 MB でした。CI の値はまだ取っていません。
