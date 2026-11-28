---
title: 多項式係数の線形漸化式の第 N 項
documentation_of: ../../neo/algo/compute/poly/p_recursive.hpp
---

多項式係数の線形漸化式 $\sum_{i=0}^{d} p_i(n) a_{n+i} = 0$ の第 $N$ 項 $a_N$ を求めます。

```cpp
template <class Conv, class K> K nth(const std::vector<std::vector<K>>& p, const std::vector<K>& a0, unsigned long long N);
```

`p[i][j]` は $p_i$ の $n^j$ の係数で、`p` は $p_0, \ldots, p_d$ の $d+1$ 本です。`a0` には初項 $a_0, \ldots, a_{d-1}$ を渡します。$N < d$ なら $a_N$ をそのまま返します。$p_d(n) \ne 0$ ($0 \le n \le N - d$) を仮定します。

`K` は体です。法 (標数) は $2 e k + 1$ より大きいことが要ります ($e$ は係数の次数の上限、$k$ は $\sqrt{N}$ くらいの 2 冪)。$1/j!$ を使うためです。

$d$ と $e$ を固定すると、手間は $N$ の平方根の準線形で、メモリは $O(d^2 \sqrt{eN})$ です。$N$ が法の 4 倍より長く、`K` が `ZMod` のときは、遷移行列が法を周期に持つので $O(\sqrt{\text{法}})$ になります。

$e = 0$ (定数係数) のときは、行列の冪で $O(d^3 \log N)$ です。定数係数なら `linear_recurrence` の方が速いです。

## 中身

$v_n = (a_n, \ldots, a_{n+d-1})$ とすると、$p_d(n) v_{n+1} = M(n) v_n$ で、$M(n)$ の最後の行が $-p_i(n)$、行 $r < d-1$ が $e_{r+1}$ です。割り算を避けるため $y_n = \left(\prod_{i<n} p_d(i)\right) v_n$ と置くと、$y_{n+1} = B(n) y_n$ で、$B(n)$ のずらしの行にも $p_d(n)$ が掛かります。よって

$$a_N = \frac{e_{d-1}^\top B(N-d) \cdots B(0)\, v_0}{\prod_{n=0}^{N-d} p_d(n)}$$

で、割り算は最後の 1 回だけです。

積は Bostan–Gaudry–Schost の Theorem 8 の形で取ります。ブロック多項式 $M_t(X) = B(X) B(X+1) \cdots B(X+t-1)$ を係数では作らず、等差数列 $X = jk$ 上の値 $M_t(jk)$ を持ち、$M_{2t}(jk) = M_t(jk + t) M_t(jk)$ で倍々にします。$X$ で $t$ ずらすのは添字で $t/k$ ずらすことになるので、標本点のずらし (`sample_points_shift`) 1 回で作れます。倍々の各段では、$M_t(jk)$ が $j$ の $e t$ 次以下の多項式なので、次数を超える先の値もずらしで延ばします。最後に $M_k(0), M_k(k), \ldots, M_k((q-1)k)$ を順に当てます。

法を周期にする道は、$B(n + M) = B(n)$ から $\prod_{i<qM+r} B(i) = \left(\prod_{i<r} B(i)\right)\left(\prod_{i<M} B(i)\right)^q$ を使います。1 周期分の積は同じ BSGS で行列として作り、行列の冪で $q$ 乗します。

procon-judge の self-poly-p-recursive (d 1〜4、e 0〜3、N は 1e10 まで) で測りました。AMD EPYC 9V74 の x64-gcc (run 46106034973) で、合計 2465 ms (最大 969 ms、29 MB) です。多点評価で係数を作る版は 4949 ms、AVX2 の核では 1270 ms でした。
