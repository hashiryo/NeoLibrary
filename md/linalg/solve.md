---
title: 連立一次方程式と核
documentation_of: ../../neo/linalg/solve.hpp
---

連立一次方程式 $Ax = b$ の解と、$Ax = 0$ の解の空間 (核) の基底を求めます。

```cpp
template <class T> std::optional<std::vector<T>> solve(const Matrix<T>& A, const std::vector<T>& b);
template <class T> std::optional<std::pair<std::vector<T>, Matrix<T>>> solve_with_kernel(const Matrix<T>& A, const std::vector<T>& b);
template <class T> Matrix<T> kernel(const Matrix<T>& A);
```

`solve` は解を 1 つ返し、解がなければ `std::nullopt` を返します。`solve_with_kernel` は、その解と核の基底の組を返します。`kernel` は核の基底だけを返します。基底は $(w - \mathrm{rank}\,A) \times w$ の行列で、各行が基底の 1 本です。

左から順に掃き出したとき、要にならなかった列を自由な変数と呼びます。返す解は、自由な変数をすべて 0 にしたものです。基底の $k$ 本目は、$k$ 番目の自由な変数を 1、ほかの自由な変数を 0 にしたものです。

`T` は体で、`ZMod` なら法は素数でなければなりません。浮動小数点数では、列の中で絶対値が最大の値を要に選び、絶対値が $10^{-8}$ 以下なら 0 とみなします (解があるかの判定にも使います)。$h \times w$ の行列に $O(hw \min(h, w))$ の時間がかかり、基底を返すときは後退代入に $O(\mathrm{rank}^2 (w - \mathrm{rank}))$ が加わります。`solve` は基底を作らないので、その分を払いません。

## 中身

`ZMod` は、$[A \mid b]$ を前進消去し、要の列が作る上三角の行列 $T$ で、自由な列と $b$ の列をまとめて後退代入します ($TX = Y$)。前進消去は `det` と同じ遅らせた更新で、行を u64 で持ち、要を 8 本ずつ貯めます。後退代入は下から 8 行ずつの組にして、組より下の行の分を 8 行まとめて足します。$X$ の 4 列ずつ 8 行分の和をレジスタに持ち、下の行の 4 列を 1 回読んで 8 行に足す形です。$\max(h, w) \le 16$ の `ZMod` とほかの体は、素朴な掃き出しで簡約した階段の形にしてから読みます。

procon-judge の yosupo-system-of-linear-equations で書き比べた zmod_sle_blk8t を写しました。写し元は x64-gcc の Xeon 8370C で 7.02 ms と首位で、今の Library の 9 倍の速さです。lib-neo の提出は `solve_with_kernel` を呼び、結果をハーネスの形に直す分も含めて、手元の M2 で写し元と 1% から 4% の差でした。記録は algo-notes の notes/linear-algebra.md にあります。
