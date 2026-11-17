---
title: 集合冪級数の速い核 (内部)
documentation_of: ../../neo/set_power_series/ranked_xor.hpp
---

法が奇数の `ZMod` (法 $p < 2^{30}$) のときに、`sps` の各関数が呼ぶ内部の核です。直接は使いません。集合の大きさ $n$ が 4 (exp などは 5) 以上 22 以下で使い、それ以外の型と大きさは汎用の道 (generic.hpp) に回します。

## rank 付きの XOR 変換

$f[S]$ を $f[S]\, t^{|S|}$ に置いて XOR 畳み込みをすると、$S \oplus T = U$ の組は $t^{|S|+|T|} = t^{|U| + 2|S \cap T|}$ の項に入ります。そのため $t^{|U|}$ の係数は、互いに素な組だけの和、すなわち部分集合畳み込みになります。

部分集合の大きさより rank の大きい項 (重なる組から出る項) を余りの項と呼びます。余りの項は、足し算と掛け算のどちらを通しても余りの項のままです。(rank $r_1$, 集合 $U_1$) と (rank $r_2$, 集合 $U_2$) の積は (rank $r_1 + r_2$, 集合 $U_1 \oplus U_2$) です。$r_1 > |U_1|$ なら $r_1 + r_2 > |U_1| + |U_2| \ge |U_1 \oplus U_2|$ となります。そのため、大きさと等しい rank の係数は途中で取り出さなくても正しいままで、逆変換と取り出しは最後に 1 回だけで済みます。余りの項はいつ捨ててもよいので、$i$ bit の集合冪級数は rank $i$ までで切れます。

これを使って、上の要素から 1 つずつ足す再帰を、変換の領域の各点の計算として書きます。$f_i[S'] = f[S' \cup \{i\}]$ とします。

- exp: $G^{(i)} = \exp(f|_{[0, i)})$ は、変換の領域で $\hat{G}^{(i+1)}[y\pm] = \hat{G}^{(i)}[y] (1 \pm t \hat{f}_i[y])$ になる。
- log: 変換の領域で $\hat{f}_i = \hat{h}_i / \hat{h}_{<i}$ を各点の冪級数の割り算 (定数項は 1) で求め、$f_i$ ごとに $i$ bit の逆変換をする。
- 合成: 係数 $e_k$ の EGF の $j$ 階の導関数に $g$ を入れたものを $G_j$ とする。変換の領域で $\hat{G}_j^{(i+1)}[y\pm] = \hat{G}_j^{(i)}[y] \pm t \hat{G}_{j+1}^{(i)}[y] \hat{g}_i[y]$ になる。
- power projection は合成の転置で、$w$ を順変換してから上の bit から双対を縮める。掛け算の転置は相関になる。
- 割り算 $f/g$ は、両方を変換して各点で割り、1 回だけ逆変換する。

どれも割り算 (整数での割り算) を使いません。逆変換で $2^{1-n}$ を掛けるので、法は奇数に限ります。

## 持ち方

変換の添字 $y$ を補集合に移すと $t$ が $-t$ に変わるので、$y$ の bit 0 が 0 の半分だけを求めます。bit 0 の変換は、$f[2S']$ と $f[2S'+1]$ を同じ位置の rank 別の値に置くだけです。値は u32 で $[0, 2p]$ に持ち、位置 8 個を 1 本の `__m256i` に並べます (塊は位置 8 個と rank の数の積)。変換の 1 段は、足し算と引き算と、$2p$ を引いて min を取る還元です。

下の bit から変換すると、bit $b$ まで済んだ位置の値は、rank が [残りの上の bit の数, それ + $b$ + 1] の範囲にしかありません。その範囲 (窓) だけを変換します。逆変換は上の bit からにして、最後の数段は答えに要る rank だけを持ちます。

各点の積は、値を $[0, p]$ に直してから 12 項までずつ u64 にためて Montgomery で還元します。塊を $2^9$ 個ずつの区画に分けて下の段を区画の中で済ませ、区画をまたぐ上の段は列ごとに作業の配列に写して行います。大きい作業の配列には huge page を頼みます。

## 速さ

procon-judge の Library Checker の 5 問で書き比べた提出 (xw1、xe1、xl1、xc1、xp1) を、法を実行時の奇数に広げて写しました。CI の x64-gcc の値は次のとおりです (algo の最大)。

| 問題 | 書き比べの提出 | 首位の写し | 今の Library | CPU |
| --- | ---: | ---: | ---: | --- |
| subset_convolution | 41.62 ms | 138.44 ms | 452.86 ms | EPYC 9V74 |
| exp_of_set_power_series | 39.25 ms | 135.54 ms | 386.65 ms | EPYC 9V74 |
| log_of_set_power_series | 45.48 ms | 171.73 ms | 569.65 ms | EPYC 7763 |
| polynomial_composite_set_power_series | 68.91 ms | 291.20 ms | 670.39 ms | EPYC 7763 |
| power_projection_of_set_power_series | 74.24 ms | 295.35 ms | なし | EPYC 7763 |

この核 (NeoLibrary 5c89acb で、積和のループと法の定数の持ち方を直したあと) は、写し元と並びます。x64-gcc の EPYC 9V45 の値は次のとおりです。subset_convolution は run 44047217621、ほかの 4 問は Solver の形にしたハーネスで測った回 (judge 5c073dd4) です。

| 問題 | この核 | 写し元 | 首位の写し |
| --- | ---: | ---: | ---: |
| subset_convolution | 37.74 ms | 45.48 ms | 109.52 ms |
| exp_of_set_power_series | 34.79 ms | 33.89 ms | 93.53 ms |
| log_of_set_power_series | 34.37 ms | 33.44 ms | 113.42 ms |
| polynomial_composite_set_power_series | 51.53 ms | 52.39 ms | 223.11 ms |
| power_projection_of_set_power_series | 52.65 ms | 58.04 ms | 186.59 ms |

記録は algo-notes の notes/set-power-series.md にあります。
