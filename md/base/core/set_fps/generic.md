---
title: 集合冪級数の汎用の道 (内部)
documentation_of: ../../neo/base/core/set_fps/generic.hpp
---

`sps` の各関数が、法が奇数の `ZMod` 以外の型と、速い核で受けない大きさのときに呼ぶ内部の道です。直接は使いません。

値の型 `T` の演算子だけを使い、割り算を使いません (ランク付きのゼータ変換)。そのため法が偶数の `ZMod`、整数、mod $2^{64}$、実数でも動きます。今の Library (hashiryo/Library) の mylib/algebra/set_power_series.hpp の写しです。名前空間と、後ろに足した 2 つの関数のほかは変えていません。足したのは、引き算の無い半環の部分集合畳み込み (`semiring_convolve<M>`、$3^n$) と、半環の exp (`semiring_exp<M>`、最小の要素を含む塊を選ぶ $3^n$) です。どちらも構造 `M` の `add`、`mul`、`one` で書いてあり、値の型は `M::T` です。
