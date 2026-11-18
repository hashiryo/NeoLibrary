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
- `sum_exp_poly(y, r, n)` は $\sum_{i=0}^{n-1} r^i f(i)$ を $O(d)$ で返します。$f$ は $y_i = f(i)$ ($0 \le i \le d$) を満たす $d$ 次以下の多項式で、$0^0 = 1$ とします。
- `sum_exp_poly_limit(y, r)` は $\sum_{i=0}^{\infty} r^i f(i)$ ($r \ne 1$) です。和は分母が $(1 - r)^{d+1}$ の有理式として定めます。

`K` は体で、`ZMod` なら法は素数であることが要ります。べき乗和 $\sum_{i<n} i^d$ は、`pow_table<K>(d, d)` (`neo/number_theory/pow_table.hpp`) を $r = 1$ で `sum_exp_poly` に渡せば求まります。

## 中身

`interpolate_iota` は、ラグランジュの基底 $L_k(x) = \prod_{j \ne k} (x - j)/(k - j)$ で $f(x) = \sum_k y_k L_k(x)$ を求めます。隣どうしの比は次の形です。

$$\frac{L_{k+1}(x)}{L_k(x)} = -\frac{(n-1-k)(x-k)}{(k+1)(x-k-1)}$$

`ZMod` (法は奇数) では、分母 $(k+1)(x-k-1)$ の逆元を 8 本の鎖でまとめて求めます。$i \equiv b \pmod 8$ の値の累積積をレーン $b$ が受け持ち、AVX2 の Montgomery の積で前向きに置いてから、レーンごとの総積の逆元から後ろ向きに戻します。分母の値はレジスタの中で作り、配列に置きません。比の累積積は、列を 4 つの区間に分けて 4 本の鎖で作ります。ベクトルの中の累積積は 3 段の並べ替えと積で作り、前のベクトルからの繰り上げを掛けます。累積積は配列に置かず、その場で $y_k$ と掛けて和にためます。区間の和には最後に前の区間の総積を掛けます。$L_0(x)$ は、分母の総積 $Q = (n-1)! \prod_{j \ge 1} (x - j)$ から求めます。作業の配列は、Linux では 2 MB 境界の mmap に置いて huge page を頼みます。

`interpolate` は $f = \sum_i w_i P(X)/(X - x_i)$ ($P = \prod_j (X - x_j)$、$w_i = y_i / P'(x_i)$) です。`ZMod` (法は 2 より大きい素数) では、$P/(X - x_i)$ の $k$ 次の係数が $\sum_{t>k} P_t x_i^{t-k-1}$ であることを使います。$s_m = \sum_i w_i x_i^m$ とすると、$f_k = \sum_m P_{k+1+m} s_m$ です。$P$ は Shoup の方法で 8 個ずつ作り、$P'(x_i)$ と $s_m$ は点を 8 個ずつレーンに載せた Montgomery の積 (4 組を並べて待ちを隠す) で求めます。$f$ は、逆に並べた $P$ と $s$ の畳み込みの下半分を `convolve_naive` の核で求めます。ほかの型は、点ごとに組立除法で $P/(X - x_i)$ を作って足します。

`sum_exp_poly` は、$S(x) = \sum_{i<x} r^i f(i)$ が $c + r^x g(x)$ ($g$ は $d$ 次以下) と書けることを使います (min_25 の方法)。$g$ の $d + 1$ 階差分が 0 なので、$c$ は次の式で求まります。

$$c = \frac{r^{d+1}}{(1 - r)^{d+1}} \sum_k \binom{d+1}{k} (-1)^{d+1-k} r^{-k} S(k)$$

そのあと $g(k) = (S(k) - c) r^{-k}$ ($0 \le k \le d$) から $g(n)$ を補間します。`ZMod` では、$S$ を 8 レーンの $r^i$ の鎖とベクトルの中の累積和で作ります。$c$ の重みは比 $-(d+2-j)/(jr)$ の累積積で、$g(n)$ の重みは $L_k(x) r^{-k}$ の比の累積積です。どちらも `interpolate_iota` と同じ 2 つの段で求めるので、$g$ は配列に作りません。

procon-judge の yosupo-sum-of-exponential-times-polynomial と -limit ($d \le 10^7$) で書き比べた se2 を、法を定数にしない形で写しました。EPYC 7763 の x64-gcc (run 44256846647) で、この形 (`pow_table` と `sum_exp_poly`) の最大は 224.01 ms と 148.50 ms でした。写し元の se2 は 187.48 ms と 116.16 ms、今の Library は 864.07 ms と 669.25 ms です。se2 は $i^d$ の表を作業の配列に直接作るのに対し、この形は `pow_table` の答えを写すので、$10^7$ 要素の配列 1 本ぶん遅くなります。係数の補間は procon-judge の self-poly-interpolation (5000 点以下) で書き比べた interp_ps の写しで、Xeon 6975P-C で 16.75 ms (組立除法の参照実装は 197.95 ms) でした。
