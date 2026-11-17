---
title: 半環の部分集合畳み込みの AVX2 の核 (内部)
documentation_of: ../../neo/set_power_series/semiring_simd.hpp
---

引き算の無い半環の `sps::convolve` と `sps::exp` が呼ぶ内部の核です。直接は使いません。型は `MinPlus`、`MaxPlus`、`MaxMin` ([semiring.hpp](../algebra/semiring.md)) で、値が 32 bit か 64 bit の符号付き整数のものです。ほかの半環は、汎用の道 ([generic.hpp](generic.md)) の $3^n$ で行います。

## 部分集合畳み込み

$h[S] = \bigoplus_{T \subseteq S} f[T] \otimes g[S \setminus T]$ を求めます。添字の下の $B$ bit (値が 4 byte なら 3、8 byte なら 2) は、1 本の `__m256i` のレーンに置きます。

- 上の bit は 2 つで回す。1 つは $T$ の上の部分 $t$ で、もう 1 つは $t$ の補集合の部分集合 $v$ (S の上の部分は $t \cup v$)。組は $3^{n-B}$ 個で、1 組で $3^B$ 項になる。
- $t$ ごとに、$f$ の塊をレーンの xor で並べ替えた $2^B$ 本を作っておく。$u \subseteq s$ のレーン $s$ には $f[t, s \oplus u]$ を置き、ほかのレーンは掛けると和の単位元となる値で埋める。
- $v$ ごとに $g$ の塊の $2^B$ 個を 1 つずつ配って掛け、木の形で畳んでから $h$ の塊へ畳む。

`MinPlus` と `MaxPlus` の積は INF で止める足し算です。核では止めずに足し、最後に INF の側へ外れた値を INF に直します (行列の積の速い道と同じ)。値がすべて 0 以上 INF 以下 (`MaxPlus` は -INF 以上 0 以下) か、すべて ±INF か絶対値が INF / 4 以下なら、厳密な積と同じ値になります。どちらでもなければ汎用の道で行います。`MaxMin` は min と max だけなので、いつも核で行います。

## exp

最大の要素を含む塊を選ぶ形にします。$h[2^i + S] = \bigoplus_{T \subseteq S} f[2^i + T] \otimes h[S \setminus T]$ を求める部分集合畳み込みを、$i = 0, \ldots, n-1$ と続けます。どの段も上の核で、段ごとに値の範囲を確かめます。範囲の外の段が出たら、全体を汎用の道 (最小の要素を含む塊を選ぶ形) でやり直します。値が INF に近いときの止め方を、今までの形と同じにするためです。
