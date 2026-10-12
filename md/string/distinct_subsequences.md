---
title: 相異なる部分列の数
documentation_of: ../../neo/string/distinct_subsequences.hpp
---

列の空でない相異なる部分列の数を、法 `mod` で割った余りとして求めます。

```cpp
template <class T> unsigned distinct_subsequences(const std::vector<T>& a, unsigned mod = 998244353);
unsigned distinct_subsequences(const std::string& s, unsigned mod = 998244353);
```

`T` は整数の型で、`mod` は 1 以上 $2^{31}$ 未満です。文字列の版は、文字を unsigned char として比べます。空の列では 0 を返します。値の幅 (最大 - 最小 + 1) が $2^{32}$ 以下なら $O(n)$ で、それより広い 64 bit の値では $O(n \log n)$ です。

空も含めた `a[0, i)` の相異なる部分列の数を $d_i$ とすると、$d_{i+1} = 2 d_i - d_j$ です。$j$ は $i$ より前で `a[i]` と同じ値が最後に出た位置で、無ければ引きません。答えは $d_n - 1$ です。値の幅が、要素の数の 2 倍と $2^{16}$ の大きい方以下なら、値から直接引く表に最後の $d$ を持って 1 回の走査で済ませます。幅が広ければ、(値 - 最小, 添字) の組を 11 bit ずつの基数ソートで並べて前の位置を求めます。全要素で同じになる桁の段は、値の OR と AND で分かるので飛ばします。幅が $2^{32}$ を超えるときは、`std::sort` で値に番号を振り直してから同じ手順で数えます。

procon-judge の yosupo-number-of-subsequences ($n \le 5 \times 10^5$、値は $10^9$ 以下) で書き比べました。Xeon 8573C の x64-gcc で、最も遅いケースの組の中央値が 10.27 ms でした (基数ソートだけの形は 13.25 ms)。組の列を huge page を頼んだ 1 つの領域に置き、前の位置と $d$ を並べ終えた側でない組の列の領域に置いた分が効いています。記録は algo-notes の notes/number-of-subsequences.md にあります。
