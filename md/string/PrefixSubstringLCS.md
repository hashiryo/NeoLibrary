---
title: 接頭辞と部分文字列の LCS
documentation_of: ../../neo/string/PrefixSubstringLCS.hpp
---

列 `s` の接頭辞と列 `t` の部分文字列の最長共通部分列の長さを、構築したあと 1 回あたり $O(1)$ で答えます。

```cpp
class PrefixSubstringLCS {
 public:
  PrefixSubstringLCS(const std::string& s, const std::string& t);
  template <class T> PrefixSubstringLCS(const std::vector<T>& s, const std::vector<T>& t);
  int query(int a, int b, int c) const;
};
```

`query(a, b, c)` は `s[0, a)` と `t[b, c)` の最長共通部分列の長さで、$0 \le a \le |s|$、$0 \le b \le c \le |t|$ で呼びます。要素は `==` で比べます。$|t| < 32767$ です。$m = |s|$、$n = |t|$ として、構築は $O(mn + m (n / 64)^2)$ です。表は $(m + 1)((n / 64 + 1)^2 \cdot 16 + (n / 64 + 2) \cdot 128)$ バイトです。$m = n = 1000$ で 6 MB ほどになります。

Tiskin の seaweed の組み替えで作りました。`s` の要素を 1 つずつ足しながら、`t` の各位置の下端に出る糸の出どころ $h_j$ を求めます。上端の位置 $x$ で入った糸は $x$、左端で入った糸は $-1$ とします。`s[0, a)` を足したあとの $h$ で、答えは `[b, c)` のうち $h_j < b$ となる $j$ の数です。これは、$j < c$ で $h_j < b$ となる数から $b$ を引いた値に等しくなります。$j < b$ なら $h_j \le j < b$ だからです。

各 $a$ の行で、$b$ が 64 の倍数のところの集合 $\{j : h_j < b\}$ を 64 bit の語の列と累積の数で控え、糸が下端に出る位置も持ちます。クエリは、$b$ を 64 の倍数に切り下げた控えから引き、残りの 63 個以下の糸を AVX2 の比較 4 回で数えます。

procon-judge の yosupo-prefix-substring-lcs ($|s|, |t| \le 1000$、クエリ $5 \times 10^5$ 個) で書き比べました。EPYC 9V74 の x64-gcc で、最も遅いケースの組の中央値が 7.75 ms でした。Library Checker の参照解の形 (クエリを $a$ ごとに分け、Fenwick 木で数える) は 26.27 ms でした。arm では AVX2 を SIMDe で真似るので遅く、クエリに SIMD を使わない形の 3 倍ほどかかります。記録は algo-notes の notes/prefix-substring-lcs.md にあります。
