---
title: Z 配列
documentation_of: ../../neo/string/z_algorithm.hpp
---

列の Z 配列を求めます。

```cpp
std::vector<int> z_algorithm(const std::string& s);
template <class T> std::vector<int> z_algorithm(const std::vector<T>& s);
```

返す配列の $i$ 番目は、`s` と `s[i, n)` の最長共通接頭辞の長さで、0 番目は $n$ です。どちらの版も $O(n)$ です。

- `z_algorithm(s)` は文字列を受け、文字を unsigned char として比べます。
- `std::vector<T>` の版は値を `==` で比べます。

右端が最も右にある一致の区間 `[l, r)` を持ち、位置 $i$ がその中で $z_{i - l} < r - i$ なら、比べずに $z_i = z_{i - l}$ とします。そうでなければ $r - i$ から文字を比べて伸ばします。

文字列の版は、$i > 0$ で $z_i$ が 1 以上になるのは $s_i = s_0$ のときだけなのを使います。0 で埋めた出力に、$s_0$ と同じ文字の位置だけを書きます。その位置は、32 byte ずつ AVX2 で比べた印から拾います。一様な文字列では、調べる位置が文字の種類分の 1 になります。伸ばすときは、次の 1 文字を比べてから、8 byte ずつ排他的論理和で一致を数えます。

procon-judge の yosupo-zalgorithm ($n \le 5 \times 10^5$) で書き比べました。EPYC 7763 の x64-gcc では、ケースの組ごとの中央値で最も遅い組が 1.38 ms でした。同じ回に、今の Library の z_algorithm は 4.05 ms でした。このうち 0.88 ms は、出力の 2 MB を確保して 0 で埋める手間です。OS は、確保した領域の 4 KB ごとのページへ初めて書き込まれたところで、メモリを割り当てて 0 で埋めるからです。記録は algo-notes の notes/z-array-palindromes.md にあります。
