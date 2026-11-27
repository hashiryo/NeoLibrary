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

多点評価は、部分積木 (葉が $X - x_i$、上が子の積) を作り、根で $f$ を根の多項式で割った余りから、子へ降りながら割っていきます。葉の余り (0 次) が $f(x_i)$ です。余りは `poly_divmod_fast` で作ります。

補間は、根の多項式 $P$ の導関数 $P'$ を多点評価して $P'(x_i)$ を作り、$w_i = y_i / P'(x_i)$ から $f = \sum_i w_i P/(X - x_i)$ を、下から $f_{\mathrm{left}} P_{\mathrm{right}} + f_{\mathrm{right}} P_{\mathrm{left}}$ でまとめます。
