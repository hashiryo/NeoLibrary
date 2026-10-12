---
title: 両端の回文木
documentation_of: ../../neo/string/PalindromeDeque.hpp
---

両端に文字を足し引きできる列について、異なる回文の数と、最長の回文接頭辞と接尾辞の長さを保ちます。

```cpp
template <int K = 26> class PalindromeDeque {
 public:
  explicit PalindromeDeque(int n);
  void push_back(int c);
  void push_front(int c);
  void pop_back();
  void pop_front();
  int size() const;
  int distinct() const;
  int longest_prefix() const;
  int longest_suffix() const;
};
```

文字は 0 以上 K 未満の整数で、K は 256 以下です。`n` は `push_back` と `push_front` を合わせて呼ぶ回数の上限です。`pop_back` と `pop_front` は空でない列で呼びます。`distinct()` は列に現れる空でない異なる回文の数で、`longest_prefix()` と `longest_suffix()` は最長の回文接頭辞と回文接尾辞の長さです。空の列ではどれも 0 です。各操作は $O(K)$ で、繰り返しはありません。初めに取る領域は $n$ に比例し、頂点あたり $8K + 24$ バイトと位置あたり 9 バイトです。

Wang、Yang、Zhu の double-ended eertree (arXiv 2210.02292) のうち、表面を記録する方法で作りました。回文の頂点ごとに、文字ごとの遷移と direct link (自分の真の回文接尾辞のうち、自分の中で直前の文字が c であるものの最長) を持ちます。列の各位置には、そこから始まる表面とそこで終わる表面を持ちます。表面は回文の出現のうち、ほかの回文の出現の接頭辞と接尾辞のどちらにもならないものです。最長の回文接頭辞は先頭の位置から始まる表面で、最長の回文接尾辞は末尾の位置で終わる表面です。足すときは、最長の回文接尾辞 (接頭辞) から direct link を 1 回引いて伸ばせる回文を見つけ、無ければ頂点を作ります。除くときは、その回文が列の中で一意なら頂点を消します。

procon-judge の yosupo-palindromes-in-deque ($q \le 5 \times 10^5$) で書き比べました。表を幅 K のまま 1 つの領域に置いて huge page を頼む形が、EPYC 9V74 の x64-gcc でケースの組ごとの中央値の最大 8.82 ms でした。同じ形で配列を 4 KB のページに取ると 21.11 ms で、半分以上が確保したばかりの領域に初めて書く代価でした。記録は algo-notes の notes/palindromes-in-deque.md にあります。
