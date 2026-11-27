---
title: 形式的冪級数
documentation_of: ../../neo/base/core/poly/fps.hpp
---

形式的冪級数の先頭 $n$ 項を求める道具です。

```cpp
template <class Conv, class K> std::vector<K> fps::inv(const std::vector<K>& f, int n);
template <class Conv, class K> std::vector<K> fps::log(const std::vector<K>& f, int n);
template <class Conv, class K> std::vector<K> fps::exp(const std::vector<K>& f, int n);
template <class Conv, class K> std::vector<K> fps::pow(const std::vector<K>& f, unsigned long long k, int n);
template <class Conv, class K> std::vector<K> fps::sqrt(const std::vector<K>& f, int n);
```

- `fps::inv<Conv>(f, n)`: $f g \equiv 1 \pmod{x^n}$ を満たす $g$ の先頭 $n$ 項。`f[0]` は逆元を持つことが要ります。
- `fps::log<Conv>(f, n)`: $\log f$ の先頭 $n$ 項。`f[0] = 1` が要ります。
- `fps::exp<Conv>(f, n)`: $\exp f$ の先頭 $n$ 項。`f[0] = 0` が要ります。
- `fps::pow<Conv>(f, k, n)`: $f^k$ の先頭 $n$ 項 ($k \ge 0$)。`f[0] = 0` のときは、0 でない最初の次数 $v$ として $x^{vk}$ を外に出します。$0^0$ は 1 です。
- `fps::sqrt<Conv>(f, n)`: $\sqrt f$ の先頭 $n$ 項。平方根が無ければ空を返します。`f[0]` は 1 でなくてもかまいません (0 のときは、0 でない最初の次数 $v$ が偶数なら $x^{v/2}$ を外に出し、定数 $c$ の平方根は `ZMod` なら `sqrt_mod` で求めます)。

`f` は $n$ 項より短くても長くてもかまいません (足りない係数は 0 と見ます)。核はテンプレート引数 `Conv` で受け、畳み込みの共通の口 (`neo/base/core/convolve/convolve.hpp`) の `convolve_range` に渡します。省くと、include した核のうち `K` を受け持つ既定の核になります。

## 中身

- `inv` は Newton 法で、長さ $m$ から $2m$ ずつ倍にします。$f g = 1 + x^m e$ のとき $g (2 - f g) = g - x^m g e$ なので、1 段は、$e$ を $f g$ の $m$ 次から $2m - 1$ 次までで求め、$g$ の $m$ 次から $2m - 1$ 次までを $-g e$ の先頭 $m$ 項で直します。
- `log` は $\log f = \int f' / f$ です。`inv` で $1/f$ を作り、$f'$ との積を 1 回の範囲の積で求めて積分します。
- `exp` は $g \leftarrow g (1 - \log g + f)$ を倍々に回します。段の入口で $g = \exp f$ の $x^m$ までと $h = 1/g$ の $x^m$ までを持ち、1 段で $h$ を Newton 1 段で $2m$ 桁に直し、$(\log g)' = g' h$ から $w = f - \log g$ の $x^m$ 以上を作って $g$ を直し、$h$ も $g$ が変わった分だけ直します。変換の長さは 4 回とも $2m$ で済みます。
- `pow` は、0 でない最初の次数を外に出し、先頭を 1 にした級数を `exp` と `log` で $k$ 乗し、定数の $k$ 乗を掛けて戻します。
- `sqrt` は、0 でない最初の次数 $v$ を外に出し、残りを Newton 法 ($h \leftarrow (h + (f / x^v) / h) / 2$) で作ります。$v$ が奇数か、定数の平方根が無ければ空を返します。

範囲の積は、核が巡回の長さを必要な分だけに詰めるので、`inv` の 1 段の変換の長さは $2m$ で済みます。核を差し替えると、998244353 でも $10^9 + 7$ でも動きます (変換を使わない素朴な核でも同じ関数が動きます)。
