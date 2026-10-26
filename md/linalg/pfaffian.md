---
title: パフィアン
documentation_of: ../../neo/linalg/pfaffian.hpp
---

交代行列 ($A^T = -A$、対角が 0) のパフィアン $\mathrm{Pf}(A)$ を求めます。$\mathrm{Pf}(A)^2 = \det(A)$ です。

```cpp
template <class T> T pfaffian(Matrix<T> A);
```

辺が奇数なら 0 を返します。`T` は体で、`ZMod` なら法は素数でなければなりません。浮動小数点数では、列の中で絶対値が最大の値を要に選びます。$n \times n$ の行列に $O(n^3)$ の時間がかかります。$A$ が交代行列かどうかは確かめません。

## 中身

列 $c = 0, 1, \ldots, n - 2$ の順に、要を行 $c + 1$ 以降から選んで行 $c + 1$ に置き、行 $c + 2$ 以降から行 $c + 1$ の倍数を引いて列 $c$ の値を消します。行を入れ替えるときは同じ番号の列も入れ替え、符号を変えます。列の操作はしません。

消し終えた行列 $H$ は上 Hessenberg の形です。引いた操作を対角が 1 の下三角の $L$ とすると $H = LA$ なので、$LAL^T = HL^T$ は交代行列で、かつ上 Hessenberg の形、つまり 3 重対角になります。その下の対角は $H$ と同じです。$\mathrm{Pf}(LAL^T) = \det(L)\,\mathrm{Pf}(A) = \mathrm{Pf}(A)$ なので、$\mathrm{Pf}(A)$ は $c$ が偶数の列の $-H_{c+1,c}$ の積になります。$c$ が奇数の列はこの積に入らないので、要がなければ何もせずに次の列へ進みます。

`ZMod` は、行を u64 で持ち、要を 8 本ずつ貯めて残りの各行を 1 回だけ読み書きします (遅らせた更新)。列を入れ替えるときは、残りの行と、貯めている要の行の両方で入れ替えます。$n \le 16$ の `ZMod` とほかの体は素朴な形にします。

procon-judge の yosupo-pfaffian-of-matrix で書き比べた zmod_pf_blk8 を写しました。写し元は x64-gcc の Xeon 8370C で 8.70 ms と首位でした。lib-neo の提出は、手元の M2 で写し元と 2% 以内の差でした。写し元と naive と yosupo_fastest (cp-algo の提出) は、奇数の列に要がないときに 0 を返します。そのため、$2 \times 2$ の塊を対角に並べた行列 (パフィアンは塊の値の積) で 0 を返しますが、こちらは正しい値を返します。記録は algo-notes の notes/linear-algebra.md にあります。
