---
title: 形式的冪級数
documentation_of: ../../neo/base/core/poly/fps.hpp
---

形式的冪級数の先頭 $n$ 項を求める道具です。今は逆数だけがあります。

```cpp
template <class Conv, class K> std::vector<K> fps::inv(const std::vector<K>& f, int n);
```

`fps::inv<Conv>(f, n)` は $f g \equiv 1 \pmod{x^n}$ を満たす $g$ の先頭 $n$ 項を返します。`f[0]` は逆元を持つことが要ります。`f` は $n$ 項より短くても長くてもかまいません (足りない係数は 0 と見ます)。

核はテンプレート引数 `Conv` で受け、畳み込みの共通の口 (`neo/base/core/convolve/convolve.hpp`) の `convolve_range` に渡します。省くと、include した核のうち `K` を受け持つ既定の核になります。

## 中身

Newton 法で、長さ $m$ から $2m$ ずつ倍にします。$f g = 1 + x^m e$ のとき $g (2 - f g) = g - x^m g e$ なので、1 段は、$e$ を $f g$ の $m$ 次から $2m - 1$ 次までで求め、$g$ の $m$ 次から $2m - 1$ 次までを $-g e$ の先頭 $m$ 項で直します。$e$ と $g e$ はどちらも範囲の積 (`convolve_range`) で、核が巡回の長さを必要な分だけに詰めるので、1 段の変換の長さは $2m$ で済みます。

核を差し替えると、998244353 でも $10^9 + 7$ でも動きます (変換を使わない素朴な核でも同じ関数が動きます)。
