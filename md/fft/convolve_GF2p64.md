---
title: 畳み込み ($\mathbb{F}_{2^{64}}$)
documentation_of: ../../neo/fft/convolve_GF2p64.hpp
---

$\mathbb{F}_{2^{64}}$ 上の多項式の積を求めます。

```cpp
std::vector<GF2p64> convolve(const std::vector<GF2p64>& a, const std::vector<GF2p64>& b);
```

長さ $n$ の `a` と長さ $m$ の `b` から、$c_k = \sum_{i + j = k} a_i b_j$ ($0 \le k \le n + m - 2$) を並べた長さ $n + m - 1$ の vector を返します。どちらかが空なら空の vector を返します。

Cantor 基底の加法的 FFT で計算します。$N = n + m$ として、積は $O(N \log N)$ 回、基底変換の XOR は $O(N \log N \log \log N)$ 回です。$n + m - 1 \le 64$ のときは素朴に掛けます。x86-64 では VPCLMULQDQ が使えるかを実行時に調べ、使えれば 256 bit の clmul で計算します。
