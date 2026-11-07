---
title: 接尾辞の LCP と部分文字列の比べ
documentation_of: ../../neo/string/SuffixLCP.hpp
---

文字列の 2 つの接尾辞の最長共通接頭辞 (LCP) の長さと、2 つの部分文字列の辞書順の大小を答えます。

```cpp
explicit SuffixLCP(const std::string& s);
int size() const;
int lcp(int i, int j) const;                      // s[i, n) と s[j, n) の LCP の長さ
int compare(int l1, int r1, int l2, int r2) const;  // s[l1, r1) と s[l2, r2) の大小
```

`lcp(i, j)` は $0 \le i, j \le n$ で使い、$i = n$ か $j = n$ なら 0 を返します。`compare` は、前の部分文字列が小さければ負、等しければ 0、大きければ正を返します。構築は $O(n)$、どちらのクエリも $O(1)$ です。文字は unsigned char として比べます。

[接尾辞配列](suffix_array.md)と LCP 配列を作り、LCP 配列に区間最小 (RangeMin) を載せます。`lcp(i, j)` は、2 つの接尾辞の接尾辞配列での順位の間の、LCP 配列の最小です。1 つのクエリで、順位の表を 2 回と区間最小の表を引きます。`compare` は LCP の次の 1 文字を比べます。

procon-judge の self-string-suffix-lcp (長さ $10^6$ の文字列に $10^6$ 個のクエリ) で書き比べました。EPYC 9V45 の x64-gcc では、構築とクエリを合わせた時間がケースごとの最大で 58.7 ms でした (abacabad... の形の列)。同じ回に、今の Library の SuffixArray と LCPArray は 95.5 ms、ローリングハッシュで長さを二分探索する形は 492.5 ms でした。先頭の 32 文字を直接比べ、それで決まらないクエリが初めて来たときに構造を作る形は、一様な文字列では構造を作らずに済んで 3 ms ほどになります。ただし LCP の長いクエリばかりの列では、直接比べた分が上乗せになって 73.4 ms かかりました。記録は algo-notes の notes/string-basics.md にあります。
