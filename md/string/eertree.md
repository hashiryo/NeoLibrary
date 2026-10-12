---
title: 回文木
documentation_of: ../../neo/string/eertree.hpp
---

列の回文木 (eertree) を求めます。

```cpp
struct Eertree {
  std::vector<int> len, par, link, last;
};
Eertree eertree(const std::string& s);
template <class T> Eertree eertree(const std::vector<T>& s);
```

頂点 0 は長さ $-1$ の根、頂点 1 は長さ 0 の根で、頂点 2 以降が列に現れる異なる回文です。回文の頂点は、初めて現れて終わる位置の順に並びます。`len[v]` は長さで、`par[v]` は両端の文字を除いた回文の頂点です (長さ 1 の回文なら 0)。`link[v]` は suffix link で、自分より短い最長の回文接尾辞の頂点を指します (長さ 1 なら 1)。根の `par` と `link` は 0 です。`last[i]` は `s[0, i]` の最長の回文接尾辞の頂点です。文字列の版は文字を unsigned char として比べ、$O(n)$ です。`std::vector<T>` の版は値を `==` と `<` で比べ、値に番号を振り直す整列の分だけ $O(n \log n)$ です。

列を左から読み、位置 $i$ の値 $x$ について、直前の最長の回文接尾辞から suffix link をたどって、$s_{i - \ell - 1} = x$ となる回文 $u$ を見つけます ($\ell$ は $u$ の長さ)。$u$ の $x$ の遷移が無ければ頂点を作ります。新しい頂点の suffix link は、$u$ の suffix link から同じようにたどった回文の $x$ の遷移です。たどる回数は全体で $2n$ 以下です。遷移は、値の種類が 32 以下なら頂点ごとの幅の表で持ちます。それより多ければ、1 つ目の子を頂点に直接持ち、2 つ目以降の子をハッシュ表に置きます。

procon-judge の yosupo-eertree ($n \le 10^6$) で書き比べました。EPYC 9V74 の x64-gcc で、ケースの組ごとの中央値の最大が 11.88 ms でした。頂点の列と遷移の表を 1 つの領域に置いて huge page を頼み、表の添字を種類の数との掛け算で求める形です。表の幅を 2 の冪に切り上げて添字をシフトで求めると、x64 では 2 倍ほど遅くなりました。記録は algo-notes の notes/eertree.md にあります。
