---
title: Lyndon 分解
documentation_of: ../../neo/string/lyndon_factorization.hpp
---

列を Lyndon 分解します。Lyndon 語は、自分のどの真の接尾辞よりも辞書順で小さい列です。どの列も、辞書順で減っていかない Lyndon 語の並びへ、ただ 1 通りに分けられます。

```cpp
std::vector<int> lyndon_factorization(const std::string& s);
template <class T> std::vector<int> lyndon_factorization(const std::vector<T>& s);
```

返すのは切れ目の列 $0 = a_0 < a_1 < \cdots < a_k = n$ で、`s[a_{t-1}, a_t)` が $t$ 番目の Lyndon 語です。Library Checker の Lyndon Factorization の出力と同じ形です。どちらの版も $O(n)$ です。

- `lyndon_factorization(s)` は文字列を受け、文字を unsigned char として比べます。
- `std::vector<T>` の版は値を `<` で比べます。

Duval の方法 (1983) で求めます。$i$ から始まる部分を読み、$j$ を進めながら周期の頭 $k$ と比べます。$s_k = s_j$ なら $k$ も進め、$s_k < s_j$ なら $k$ を $i$ に戻します。$s_k > s_j$ で止まったら、周期 $j - k$ の Lyndon 語を $i$ から切り出します。

文字列の版は 3 つの手を足しています。$k$ が $i$ に戻っている間は、$s_j > s_i$ の位置では $j$ が進むだけなので、$s_j \le s_i$ になる最初の $j$ を 32 byte ずつ AVX2 で探して飛びます。$s_k = s_j$ の続く長さは、1 文字目を比べてから 8 byte ずつ排他的論理和で数えます。切れ目の列は $n + 1$ 個分を先に取り、同じ周期の因子が続く分はまとめて書きます。

procon-judge の yosupo-lyndon-factorization ($n \le 5 \times 10^5$) で書き比べました。EPYC 7763 の x64-gcc では、ケースの組ごとの中央値で最も遅い組 (すべて同じ文字の列) が 0.90 ms、一様な乱数の列が 0.31 ms でした。同じ回に、Duval の方法をそのまま書いたものは 3.08 ms でした。すべて同じ文字の列は因子が $n$ 個あるので、その多くは出力の 2 MB に初めて触る手間です。記録は algo-notes の notes/lyndon.md にあります。
