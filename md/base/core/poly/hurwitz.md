---
title: Hurwitz 級数環 (二項畳み込み)
documentation_of: ../../neo/base/core/poly/hurwitz.hpp
---

可換環 $R$ の列 $(a_n)$ に、項ごとの和と二項畳み込み
$$(a b)_n = \sum_{k=0}^{n} \binom{n}{k} a_k b_{n-k}$$
を入れた環を、Hurwitz 級数環 $H(R)$ と呼びます (Keigher 1997)。数え上げでは EGF の環と呼びます。微分は添字のずらし $(\partial a)_n = a_{n+1}$ で、`inv`、`log`、`exp` に整数の割り算が要らないので、合成数の法でも $2^{64}$ でも定義できます。

```cpp
template <class Conv, class Crt, class K> std::vector<K> hurwitz::convolve(const std::vector<K>& a, const std::vector<K>& b);
template <class Conv, class Crt, class K> std::vector<K> hurwitz::inv(const std::vector<K>& f, int n);
template <class Conv, class Crt, class K> std::vector<K> hurwitz::log(const std::vector<K>& f, int n);
template <class Conv, class Crt, class K> std::vector<K> hurwitz::exp(const std::vector<K>& f, int n);
```

- `hurwitz::convolve<Conv, Crt>(a, b)`: 二項畳み込み $c_n = \sum_k \binom{n}{k} a_k b_{n-k}$ の全体。どちらかが空なら空。
- `hurwitz::inv<Conv, Crt>(f, n)`: 二項畳み込みの環で $f g \equiv 1 \pmod{x^n}$ を満たす $g$ の先頭 $n$ 項。`f[0]` は逆元を持つことが要ります。
- `hurwitz::log<Conv, Crt>(f, n)`: $\log f$ の先頭 $n$ 項。`f[0] = 1` が要ります。
- `hurwitz::exp<Conv, Crt>(f, n)`: $\exp f$ の先頭 $n$ 項。`f[0] = 0` が要ります。

`Conv` は普通の畳み込みの核、`Crt` は多素数の核です。どちらも省くと、include した核の既定になります (`convolve.hpp` の `kernel_t` と、`crt.hpp` の `default_crt_engine`)。

## 道

法 $M$ の素因数がすべて $n + m - 1$ 以上なら、$1/k!$ が法で割り切れないので、$a_k \mapsto a_k / k!$ で普通の畳み込みに直し、$n!$ を掛けて戻します (1 回の畳み込みで済みます)。

そうでなければ、EntropyIncreaser「任意模数二项卷积」の方法を使います。$M$ の素因数のうち $n + m - 1$ 未満のものを $S$ として、$k! = u_k \prod_{p \in S} p^{v_p(k!)}$ と分け ($u_k$ は $M$ と互いに素)、
$$c_n = u_n \sum_k \Bigl(\prod_{p \in S} p^{d_{p,k}}\Bigr) (a_k u_k^{-1}) (b_{n-k} u_{n-k}^{-1}), \qquad d_{p,k} = v_p(n!) - v_p(k!) - v_p((n-k)!)$$
の和 $S_n$ を、NTT の素数 $q$ の下で $A_k = (a_k u_k^{-1}) \prod_p p^{-v_p(k!)}$ と $B_j$ を畳み込んで $\prod_p p^{v_p(n!)}$ を掛けて求め、上界 $N^{1+|S|} M^2$ を覆う本数の素数で Garner して戻します。上界が素数の積で覆えない法 (小さい素因数が多い大きい法) は、$O(n m)$ の素朴な形に落ちます。

## 中身

- `inv` は Newton 法 ($g \leftarrow g (2 - f g)$) を倍々に回します。二項畳み込みでは $x^m$ を掛けても添字がずれるだけにならないので、1 段は $(f g)$ の $m$ 次以上の尾を次数 $m$ に置き直し、$g$ の $m$ 次から $2m - 1$ 次までを $-(g t)$ の同じ次数で直します。
- `log` は $\partial(\log f) = (\partial f) / f$ を積分します。Hurwitz 級数環の積分は添字をずらすだけで、割り算は要りません。
- `exp` は Newton 法 $g \leftarrow g (1 + f - \log g)$ を倍々に回します。

ふつうの冪級数 (`fps.hpp`) と共通のフレームワークは作りません。共通にできるのは Newton 法の骨組みだけで、商と余りや多点評価のような多項式の道具は Hurwitz 級数環で使う場面がありません (Hurwitz 級数環は $R[[x]]$ と同型でないことがあり、$\mathbb{F}_p$ では多変数の切り詰めの環になります)。
