---
title: 接頭辞関数 (KMP)
documentation_of: ../../neo/string/prefix_function.hpp
---

列の接頭辞関数 (KMP の失敗関数) を求めます。

```cpp
std::vector<int> prefix_function(const std::string& s);
template <class T> std::vector<int> prefix_function(const std::vector<T>& s);
```

返す配列の $i$ 番目は、`s[0, i]` の、自分より短い最長の border (接頭辞と接尾辞の両方に一致する部分) の長さです。0 番目は 0 です。どちらの版も $O(n)$ です。最小の周期は $n - \pi_{n - 1}$ で求まります。

- `prefix_function(s)` は文字列を受け、文字を unsigned char として比べます。
- `std::vector<T>` の版は値を `==` で比べます。

KMP の手順で求めます。1 つ前の接頭辞の border を伸ばせるかを見て、伸ばせなければ border の border へたどります。一致の判定は分岐で書きます。`k += (s[i] == s[k])` と書くと、GCC は比べた結果を足し算にするので、次に読む `s[k]` の位置が前の読み込みの結果を待つ鎖になります。分岐にすれば、外れない列では先へ進めます。

文字列の版は、使われている文字が 4 種類以下なら、border の長さが 64 未満の間は、KMP のオートマトンの遷移の表 (64 行 4 列) を 1 回引いて次の長さを求めます。2 文字の一様な列では border をたどる分岐が半分ほど外れるので、表の読み込みの方が速いためです。64 以上の長さと、5 種類以上の文字の列は、分岐の KMP で求めます。

procon-judge の self-string-prefix-function (長さ $10^6$) で書き比べました。EPYC 7763 の x64-gcc では、ケースごとの最大が 4.45 ms でした (2 文字の一様な列)。同じ回に、一致の判定を足し算で書いた KMP は 8.20 ms、今の Library の z_algorithm と z_to_kmp で Z 配列から直す形は 16.00 ms でした。記録は algo-notes の notes/string-basics.md にあります。
