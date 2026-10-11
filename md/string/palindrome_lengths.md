---
title: 回文の長さ
documentation_of: ../../neo/string/palindrome_lengths.hpp
---

列の各中心について、最長の回文の長さを Manacher の方法で求めます。

```cpp
std::vector<int> palindrome_lengths(const std::string& s);
template <class T> std::vector<int> palindrome_lengths(const std::vector<T>& s);
```

返す配列の長さは $2n - 1$ で、$c$ 番目は位置 $c / 2$ を中心とする最長の回文の長さです。$c$ が偶数なら文字 $s_{c / 2}$ を中心とする奇数の長さ、奇数なら $s_{(c - 1) / 2}$ と $s_{(c + 1) / 2}$ の間を中心とする偶数の長さです。Library Checker の Enumerate Palindromes の出力と同じ並びです。どちらの版も $O(n)$ です。

- `palindrome_lengths(s)` は文字列を受け、文字を unsigned char として比べます。
- `std::vector<T>` の版は値を `==` で比べます。

中心 $c$ の回文を `[lo, hi)` とすると $lo + hi = c + 1$ なので、回文は右端だけで決まります。$2n - 1$ 個の中心を 1 本のループで回し、右端が最も右にある回文について対称な中心の答えから、その回文の中で決まる分は比べずに済ませます。

文字列の版は、両隣が一致しない中心の答えが自明なのを使います。文字の中心なら長さ 1、文字の間の中心なら 0 です。出力の偶数番目を 1、奇数番目を 0 で埋めてから、両隣が一致する中心だけを中心の順に調べます。両隣の一致は、32 文字ずつ AVX2 で比べた 2 つの印から拾います。伸ばすときは、次の 1 文字を比べてから、右側の 8 byte と、左側の 8 byte を byte の順で反転したものの排他的論理和で一致を数えます。

procon-judge の yosupo-enumerate-palindromes ($n \le 5 \times 10^5$) で書き比べました。EPYC 7763 の x64-gcc では、ケースの組ごとの中央値で最も遅い組 (すべて同じ文字の列) が 3.67 ms、一様な乱数の列が 1.95 ms でした。同じ回に、今の Library のローリングハッシュと二分探索で求める解は 177.60 ms でした。記録は algo-notes の notes/z-array-palindromes.md にあります。
