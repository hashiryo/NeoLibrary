---
title: 余因子行列
documentation_of: ../../neo/linalg/adjugate.hpp
---

正方行列の余因子行列 $\mathrm{adj}(A)$ を求めます。$\mathrm{adj}(A)_{ij}$ は、$A$ から行 $j$ と列 $i$ を除いた行列の行列式に $(-1)^{i+j}$ を掛けたものです。$A\,\mathrm{adj}(A) = \mathrm{adj}(A)\,A = \det(A) I$ を満たします。

```cpp
template <class T> Matrix<T> adjugate(const Matrix<T>& A);
```

`T` は体で、`ZMod` なら法は素数でなければなりません。浮動小数点数では、絶対値が $10^{-8}$ 以下の要を 0 とみなします。$n \times n$ の行列に $O(n^3)$ の時間がかかります。乱数は使いません。

## 中身

階数で分けます。階数が $n$ なら $\mathrm{adj}(A) = \det(A) A^{-1}$ で、$n - 2$ 以下なら 0 です。$n - 1$ なら、$Ax = 0$ と $y^T A = 0$ の 0 でない解で $\mathrm{adj}(A) = c\,x y^T$ と書けます。

$[A \mid I]$ を前進消去すると (要のない列は 1 本まで飛ばし、2 本目で止めて 0 を返します)、要にならなかった行 $r_0$ の右半分が $y$ ($y_{r_0} = 1$) になります。要のない列 $f$ からは、後退代入で $x$ ($x_f = 1$) が求まります。$c = \mathrm{adj}(A)_{f r_0}$ は、$A$ から行 $r_0$ と列 $f$ を除いた行列の行列式に $(-1)^{f + r_0}$ を掛けたものです。その行列式は、要の値の積に、要の行の並びの置換の符号を掛ければ求まります。前進消去と後退代入は `inverse` の核 (右半分の付け替えつき) をそのまま使います。

procon-judge の yosupo-adjugate-matrix の首位だった zmod_adj_blk8t は、乱択の行と列で 1 回り大きくした行列の逆行列から Jacobi の公式で求める形でした。こちらは付け替えつきの逆行列の核を使えます。手元の M2 では、zmod_adj_blk8t の 0.78 倍 (最大のケース) と 0.69 倍 (全ケースの和) の時間でした。法が小さい素数でも、乱択がはずれて答えを誤ることがありません。記録は algo-notes の notes/linear-algebra.md にあります。
