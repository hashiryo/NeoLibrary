---
title: 多変数の畳み込み (truncated)
documentation_of: ../../neo/base/core/convolve/multivariate.hpp
---

`dim` で各次元の大きさを受け取り、$n = \prod_i \mathrm{dim}_i$ 個の配列 2 本の、各次元で $\mathrm{dim}_i$ 以上の次数を落とした積を求めます。

```cpp
template <class Ntt, class K> std::vector<K> multivariate_convolve(const std::vector<int>& dim, const std::vector<K>& f, const std::vector<K>& g);
```

ランクの変数 $\chi_i = \sum_j \lfloor i / (\mathrm{dim}_0 \cdots \mathrm{dim}_{j-1}) \rfloor \bmod k$ を使います。$\chi_i + \chi_j \equiv \chi_{i+j} \pmod k$ で、桁上がりは鎖が $k$ だけ進むことに吸収されるので、クラスごとに長さ $2n$ 以上の変換をして、$(i, j)$ のクラス $i + j \bmod k$ に足せば truncated の積が戻ります。変換の回数は $2k$ (順) $+ k$ (逆) で、点ごとの積が $k^2$ 本です。

変換は単独の NTT の核の `dft` と `idft` を使うので、核を差し替えられます (`Ntt` でも実験用の `NttAvx2` でも同じ値になることを確かめてあります)。係数の型は法がコンパイル時に決まった `ZMod` だけです (`ZMod<0>` は使えません)。
