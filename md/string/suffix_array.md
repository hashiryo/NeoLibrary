---
title: 接尾辞配列
documentation_of: ../../neo/string/suffix_array.hpp
---

列の接尾辞配列を SA-IS で求めます。

```cpp
std::vector<int> suffix_array(const std::string& s);
std::vector<int> suffix_array(const std::vector<int>& s, int K);
template <class T> std::vector<int> suffix_array(const std::vector<T>& s);
```

返す配列の $i$ 番目は、辞書順で $i$ 番目の接尾辞の開始位置です。長さを $n$ とします。

- `suffix_array(s)` は文字列を受け、文字を unsigned char として比べます。$O(n)$ です。
- `suffix_array(s, K)` は、値がどれも $0$ 以上 $K$ 未満の列を受けます。$O(n + K)$ です。
- `suffix_array(s)` に `std::vector<T>` を渡すと、値を座標圧縮してから解きます。`T` は `<` で比べられれば何でもよく、$O(n \log n)$ です。

1 回目の induced sorting は、libsais (Ilya Grebnov) の組み立て方で回します。位置を「自分が L 型か」と「1 つ前が L 型か」の 4 種類に分けて文字ごとに数え、次の走査で読む位置だけを区画に詰めて置くので、走査の中に分岐がありません。同じ走査の中で、LMS 部分文字列の順位も付けます。文字の種類が長さに比べて多い段では、この表が大きくなりすぎるので、表を 1 本だけ持つ素直な induced sorting で並べます。最後の induced sorting では、置く値の符号に「1 つ前を今の走査で induce するか」を持たせる sais-lite (Yuta Mori) の手を使います。L 型と S 型の種類を配列に持たず、LMS の位置は 1 bit の印で持ちます。

procon-judge の yosupo-suffixarray ($n \le 5 \times 10^5$) で書き比べました。EPYC 7763 の x64-gcc では、最も遅いケースで 14.74 ms でした。同じ回に、Library の SuffixArray は 26.18 ms、ACL は 41.72 ms、libsais 2.10.4 は 15.30 ms でした。libsais から先読みを抜くと 13.23 ms になり、まだ 1 割ほど及びません。
