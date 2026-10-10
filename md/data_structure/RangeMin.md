---
title: 区間最小 (静的な列)
documentation_of: ../../neo/data_structure/RangeMin.hpp
---

32 bit の符号なし整数の列について、区間の最小を答えます。列は作ったあと変えられません。

```cpp
explicit RangeMin(const std::vector<unsigned>& a);
unsigned min(int l, int r) const;  // min(a[l], ..., a[r - 1])
int size() const;
```

`min(l, r)` は $0 \le l < r \le N$ で使います。構築は $O(N)$、1 回の `min` は $O(1)$ です。

列を 16 個ずつの塊に分け、塊の中で左から累積した最小と、右から累積した最小を持ちます。塊ごとの最小の列には sparse table を載せます。区間の両端が別の塊なら、両端の累積と、間の塊の sparse table の 2 か所を引いた min が答えです。同じ塊の中の区間は、塊の 16 個 (キャッシュの 1 行) を AVX2 で読んで、区間の外を最大値で埋めた min を返します。表は $N = 5 \times 10^5$ で 8 MB ほどで、1 MB 以上なら 2 MB 境界に置いて huge page を頼みます (Linux のみ)。

procon-judge の yosupo-staticrmq ($N = Q = 5 \times 10^5$) で測ると、構築とクエリを合わせた時間は、EPYC 7763 の x64-gcc で最大 4.2 ms ほどです。今の Library の SparseTable の 5.8 分の 1 にあたります。
