---
title: 多項式の因数分解
documentation_of: ../../neo/polynomial/poly_factorize.hpp
---

素数 $p$ を法とする多項式を、monic な既約多項式の積に分解します。

```cpp
template <class Z> std::vector<std::pair<std::vector<Z>, int>> poly_factorize(std::vector<Z> f);
```

返すのは (既約因子、重複度) の組で、順は決めていません。$f$ の最高次の係数は落とし、monic にしてから分解します。$f$ が 0 か定数なら空を返します。`Z` は `ZMod` で、法は素数であることが要ります。FFT を使わない素朴な形で、次数 100 程度までを想定しています。

## 中身

3 つの段で分けます。

1. 平方因子なしの分解は、$f'$ と gcd を取る Yun の方法です。$f' = 0$ の部分は $f = h(x^p)$ なので、$h$ に直して重複度を $p$ 倍します。
2. 次数ごとの分解は、法 $s$ での Frobenius $h \mapsto h^p$ を、行 $x^{ip} \bmod s$ を並べた行列で $O(n^2)$ で当て、$\gcd(s, x^{p^d} - x)$ を取ります。
3. 同じ次数 $d$ の因子の積は、乱択の $a$ で分けます。$p$ が奇数なら $\gcd(g, a^{(p^d - 1)/2} - 1)$、$p = 2$ なら $\gcd(g, a + a^2 + \cdots + a^{2^{d-1}})$ を取ります。

procon-judge の yosupo-factorization-of-polynomials ($N \le 100$、$p \le 998244353$) で書き比べた naive を写しました。Xeon 6975P-C の x64-gcc (run 44231884007) で、naive の最大は 23.77 ms でした。
