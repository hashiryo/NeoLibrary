---
title: relaxed (オンライン) の畳み込み
documentation_of: ../../neo/base/core/poly/relaxed.hpp
---

係数を 1 つずつ受け取り、積の係数を次数の順に返す畳み込みです。計算量は $O(n \log^2 n)$ ($n$ は `next` を呼んだ回数) です。

```cpp
template <class K, class Conv> struct relaxed::SemiRelaxedConvolution;
template <class K, class Conv> struct relaxed::Convolution;
template <class Conv, class K> std::vector<K> relaxed::exp(const std::vector<K>& f, int n);
template <class Conv, class K> std::vector<K> relaxed::inv(const std::vector<K>& f, int n);
```

- `relaxed::SemiRelaxedConvolution<K, Conv> s(g)`: $g$ を構築で全部受け取り (オフライン)、`s.next(f_i)` が $(f g)_i$ を返します。
- `relaxed::Convolution<K, Conv> s`: $f$ と $g$ の両方を `s.next(f_i, g_i)` で受け取り、$(f g)_i$ を返します。
- `relaxed::exp<Conv>(f, n)` と `relaxed::inv<Conv>(f, n)`: `fps::exp` と `fps::inv` と同じ値の先頭 $n$ 項。どちらも上の畳み込みの上で係数を 1 つずつ決めます。

どちらもまだ渡していない次数の係数は 0 として扱います。`next` は何回でも呼べて、$i$ 回目 (0 始まり) の呼び出しが $i$ 次の係数を返します。長さ 0 は `next` を 1 回も呼ばないことで、`SemiRelaxedConvolution` の $g$ は空でもかまいません。$f$ と $g$ の長さが違う場合も、短い側を 0 で埋めた積の係数になります。あとから渡す係数が前の次数に影響することはありません。

係数の型 `K` が先頭のテンプレート引数です。クラステンプレートは既定の引数のあとに既定の無い引数を置けないためです。`SemiRelaxedConvolution<Z>` が既定の核で、`SemiRelaxedConvolution<Z, NttR4>` のように核を選べます。

`t` 回目の呼び出しで $2^k \mid t$ となる大きさ $2^k$ の $f$ のブロック ($t - 2^k$ から $t$ の手前まで) と、$g$ のブロック ($2^k$ から $2^{k+1}$ の手前まで) の積を 1 回の畳み込みで足します。両方がオンラインのときは転置の対も $t \ge 2^{k+1}$ で取ります。$t = 2^{k+1}$ のときは同じ正方形になるので 1 回だけです。0 次の行と列 ($f_0 g_j$ と $f_i g_0$) は `next` の中で直接足します。分割は hly1204 さんの記事 (UOJ のブログ 7319) の平方分割を、再帰でなく `next` の反復に落としたものです。

内部では核の `convolve_range` だけを呼ぶので、`dft` や `idft` を持たない核 (`Naive` など) でも動きます。速さは核の呼び出しごとの確保に縛られます。

`exp` と `inv` は、オフライン側を $(k f_k)$ と $(f_{k+1})$ にできるので、半オンラインの 1 本で書けます。両方オンラインの版より 1.6〜1.7 倍速いです。
