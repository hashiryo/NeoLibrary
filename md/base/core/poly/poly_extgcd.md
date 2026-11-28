---
title: 多項式の拡張 gcd (half-gcd)
documentation_of: ../../neo/base/core/poly/poly_extgcd.hpp
---

$f x + g y = \gcd(f, g)$ を満たす $(d, x, y)$ を $O(M(n) \log n)$ で求めます。

```cpp
template <class Conv, class K> std::vector<K> poly_extgcd(const std::vector<K>& f, const std::vector<K>& g, std::vector<K>& x, std::vector<K>& y);
```

`f` と `g` は末尾の 0 を取り除いた係数の vector です (0 は空)。`x` と `y` に係数が入り、戻り値は $d = \gcd(f, g)$ です。$d$ は monic にしません ($f$ と $g$ の線形結合のまま)。`f` と `g` は空でないことが要ります。

## 中身

half-gcd で、2 つの再帰から成ります。

- `hgcd(p0, p1)` は、$\deg p_0 > \deg p_1$ のときに、$M (p_0, p_1)^T$ の 2 行目が $\deg p_0 / 2$ 未満になる 2×2 の多項式行列 $M$ を返します。高次の半分を切り出して再帰し、その行列で 1 段進めてから、残りをまた半分ずつ進めます。途中で余りが目標の次数より短くなったら、その時点で打ち切ります。
- `cogcd(p0, p1)` は、$p_0 x + p_1 y = \gcd(p_0, p_1)$ の $(x, y)$ を返します。小さい方は商と余りを直接取り、大きい方は `hgcd` の行列で進めます。

商と余りは `poly_divmod_fast` (Newton 法の逆数)、多項式の積は核の `convolve_range` を使うので、$M(n)$ は選んだ核の速さになります。今の Library の `extgcd` (`mylib/fft/extgcd.hpp`) と同じ組み立てです。
