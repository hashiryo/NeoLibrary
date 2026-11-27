---
title: 多点評価と補間
documentation_of: ../../neo/base/core/poly/multipoint.hpp
---

点での値の列と、多項式の係数を相互に変換します。

```cpp
template <class Conv, class K>
std::vector<K> multipoint_eval(const std::vector<K>& f, const std::vector<K>& x);
template <class Conv, class K>
std::vector<K> interpolate(const std::vector<K>& x, const std::vector<K>& y);
```

- `multipoint_eval<Conv>(f, x)`: $f$ の各点 $x_i$ での値。手間は $O(M(n) \log n)$。
- `interpolate<Conv>(x, y)`: 点 $(x_i, y_i)$ を通る $n - 1$ 次以下の多項式の係数。手間は $O(M(n) \log n)$。

`K` は体で、`ZMod` なら法は素数です。点は相異なることが要ります。小さい $n$ では、同じ約束の `lagrange::interpolate` (`neo/algo/compute/poly/lagrange.hpp`、$O(n^2)$) の方が速いことがあります (手元の M2 で交差は $n \approx 1000$)。

核はテンプレート引数 `Conv` で受け、畳み込みの共通の口の `convolve_range` に渡します。省くと、include した核のうち `K` を受け持つ既定の核になります。

## 中身

どちらも部分積木 (葉が $X - x_i$、上が子の積) の上で、転置の形で計算します。木の多項式は反転して持ち ($P^* = X^{\deg P} P(1/X)$、定数項が 1)、内側の節は子の最高次の 1 を外して掛けます。親の反転は $1 + (u + w + uw)$ なので、積 $uw$ は $2k - 1$ 項で済み、巡回の長さが子の 2 倍になりません。

多点評価は、根で $S = f \bmod P$ を作り、$v = \operatorname{rev}(S) / P^*$ を 1 回だけ求めてから、各節のブロック (節の次数と同じ長さ) を

- 左の子: $v \cdot P_R^*$ の $\deg P_R$ 次から $\deg P$ 次の窓
- 右の子: $v \cdot P_L^*$ の $\deg P_L$ 次から $\deg P$ 次の窓

で子へ降ろします。葉のブロックが $f(x_i)$ です。割り算は根の 1 回だけで、あとは範囲の積 (`convolve_range`) だけです。子のブロックは親のブロックと同じ場所に入るので、作業の配列は 1 本で済みます。

補間は、$P'$ を同じ降下で評価して $w_i = y_i / P'(x_i)$ を作り、下から

$$G = G_L P_R^* + G_R P_L^*$$

でまとめます。反転したまままとまるので、最後に 1 回反転して係数に戻します。

点の無い葉は多項式 1 にして積で消します ($X$ にすると $x_i = 0$ の点と根が重なり、補間が割れます)。
