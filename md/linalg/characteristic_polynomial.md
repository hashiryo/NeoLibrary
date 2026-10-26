---
title: 特性多項式
documentation_of: ../../neo/linalg/characteristic_polynomial.hpp
---

正方行列の特性多項式 $\det(xI - A)$ を求めます。

```cpp
template <class T> std::vector<T> characteristic_polynomial(Matrix<T> A);
```

返す値は係数 $p_0, p_1, \ldots, p_n$ の列 ($p_n = 1$) です。`T` は体で、`ZMod` なら法は素数でなければなりません。浮動小数点数では、列の中で絶対値が最大の値を要に選びます。$n \times n$ の行列に $O(n^3)$ の時間がかかります。

## 中身

相似変換で上 Hessenberg の形 $H$ に直してから、左上 $k \times k$ の特性多項式 $p_k$ を次の漸化式で順に求めます。

$$p_{k+1} = (x - H_{kk})\,p_k - \sum_{i<k} H_{ik} \Bigl(\prod_{j=i+1}^{k} H_{j,j-1}\Bigr) p_i$$

列 $k$ の手順では、行 $k + 1$ 以降から列 $k$ の値が 0 でない行を選んで行 $k + 1$ に置きます (行と列を同時に入れ替えます)。$f_i = H_{ik} / H_{k+1,k}$ ($i \ge k + 2$) とします。行 $i$ から行 $k + 1$ の $f_i$ 倍を引き (左から $L^{-1}$)、列 $i$ の $f_i$ 倍を列 $k + 1$ に足します (右から $L$)。

`ZMod` は行列を $2M$ 未満の u32 で持ちます。列の操作では、全部の行で「行の列 $k + 2$ 以降と $f$ の内積」を列 $k + 1$ へ足します。そこで行 $k + 2$ 以降は、行の更新と内積を 1 回の読み書きで済ませます。行の更新は、掛ける数を行ごとに決まった値として Shoup の方法で $2M$ 未満に直し、内積は u64 にためて 4 回ごとに畳みます。漸化式の和は、多項式を u64 で持って 8 本ずつまとめて足します。ほかの体は素朴な形にします。

procon-judge の yosupo-characteristic-polynomial で書き比べた zmod_hess32 を写し、法を定数から変数に変えました。写し元は x64-gcc の Xeon 8573C で 30.72 ms と首位で、前の首位 (116.25 ms) の 3.8 倍の速さです。lib-neo の提出は、手元の M2 で写し元より 2% ほど速くなりました。記録は algo-notes の notes/linear-algebra.md にあります。
