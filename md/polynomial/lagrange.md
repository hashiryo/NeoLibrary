---
title: 多項式の補間
documentation_of: ../../neo/polynomial/lagrange.hpp
---

連続する点の値からの 1 点の値、任意の点からの係数、多項式と指数の積の和を求めます。

```cpp
template <class K> K interpolate_iota(const std::vector<K>& y, unsigned long long x);
template <class K> std::vector<K> interpolate(const std::vector<K>& x, const std::vector<K>& y);
template <class K> K sum_exp_poly(const std::vector<K>& y, K r, unsigned long long n);
template <class K> K sum_exp_poly_limit(const std::vector<K>& y, K r);
```

- `interpolate_iota(y, x)` は、$y_i = f(i)$ ($0 \le i < n$) を満たす $n - 1$ 次以下の $f$ の $f(x)$ を $O(n)$ で返します。今の Library の `sample_points_shift` にあたります。
- `interpolate(x, y)` は、点 $(x_i, y_i)$ を通る $n - 1$ 次以下の多項式の係数を $O(n^2)$ で返します。$x_i$ は相異なることが要ります。
- `sum_exp_poly(y, r, n)` は、$y_i = f(i)$ ($0 \le i \le d$) を満たす $d$ 次以下の $f$ について $\sum_{i=0}^{n-1} r^i f(i)$ を $O(d)$ で返します。$0^0 = 1$ とします。
- `sum_exp_poly_limit(y, r)` は $\sum_{i=0}^{\infty} r^i f(i)$ ($r \ne 1$) です。和は分母が $(1 - r)^{d+1}$ の有理式として定めます。

`K` は体で、`ZMod` なら法は素数であることが要ります。べき乗和 $\sum_{i<n} i^d$ は、`pow_table<K>(d, d)` (`neo/number_theory/pow_table.hpp`) を $r = 1$ で `sum_exp_poly` に渡せば求まります。

## 中身

`interpolate_iota` は、$P = \prod_j (x - j)$ として次の式で求めます。

$$f(x) = \frac{P}{(n-1)!} \sum_k (-1)^{n-1-k} \binom{n-1}{k} \frac{y_k}{x - k}$$

`ZMod` では、$1/k$ と $1/(x - k)$ を 8 本の鎖でまとめて求めます。$i \equiv b \pmod 8$ の数の累積積をレーン $b$ が受け持ち、AVX2 の Montgomery の積で前向きに置いてから、レーンごとの総積の逆元から後ろ向きに戻します。$1/(x - k)$ は後ろ向きの段でそのまま和にためて、配列に置きません。$\binom{n-1}{k}$ は比 $(n - k)/k$ の累積積で、8 個の区間の鎖に分けてから前の区間の積を掛けます。2 MB 以上の作業の配列は、Linux では 2 MB 境界の mmap に置いて huge page を頼みます。

`sum_exp_poly` は、$S(x) = \sum_{i<x} r^i f(i)$ が $c + r^x g(x)$ ($g$ は $d$ 次以下) と書けることを使います (min_25 の方法)。$g$ の $d + 1$ 階差分が 0 なので、$c$ は次の式で求まります。

$$c = \frac{r^{d+1}}{(1 - r)^{d+1}} \sum_k \binom{d+1}{k} (-1)^{d+1-k} r^{-k} S(k)$$

そのあと $g(k) = (S(k) - c) r^{-k}$ ($0 \le k \le d$) から $g(n)$ を補間します。

procon-judge の yosupo-sum-of-exponential-times-polynomial と -limit ($d \le 10^7$) で書き比べた se1 を写し、作業の配列を減らしました。Xeon 6975P-C の x64-gcc (run 44231884007) で、se1 の最大は 409.99 ms と 262.99 ms でした。今の Library は 712.97 ms と 545.93 ms です。
