---
title: トライ
documentation_of: ../../neo/string/Trie.hpp
---

文字列の多重集合を持ち、ある文字列で始まる要素の数を答えます。

```cpp
template <int K = 26, char Base = 'a'> class Trie;
Trie();
void insert(const std::string& s);
void erase(const std::string& s);              // s は今の多重集合にあること
int size() const;                              // 要素の数
int count(const std::string& s) const;         // s と等しい要素の数
int count_prefix(const std::string& p) const;  // p で始まる要素の数
```

文字は `Base` から `Base + K - 1` までで、既定は英小文字です ($K \le 255$)。どの操作も、渡す文字列の長さに比例する時間です。空の文字列も要素にできます。

頂点ごとに、その頂点を通る要素の数と、その頂点で終わる要素の数を持ちます。頂点の大きさは子の数で変えます (adaptive radix tree の Node4 と Node256 の考え方)。子が 4 つまでの頂点は 32 byte で、子の文字 4 つを 1 語に詰め、語の中の byte の一致で子を探します。5 つ目の子が来たら $K$ 個の子の欄を持つ頂点に替え、親の欄を書き換えます。根は初めから $K$ 個の欄を持ちます。深い頂点はほとんど子が 1 つなので、$K$ 個の欄を全部の頂点に持つ形より表がずっと小さくなります。頂点の列は 2 MB 以上なら、Linux では huge page を頼みます。

procon-judge の self-string-prefix-count (追加と削除と接頭辞の個数の問い合わせ、文字の合計が 300 万前後) で書き比べました。EPYC 7763 の x64-gcc では、ケースごとの最大が 32.4 ms でした。同じ回に、26 個の子の欄を全部の頂点に持つ形は 230.2 ms、その表を huge page の大きい仮想の領域に置いた形は 59.4 ms でした。接頭辞のハッシュから個数への表で答える形は 76.6 ms です。頂点の列を写し直さずに大きい仮想の領域で伸ばす形 (34.6 ms) は、huge page の vector を倍々に伸ばすこの形とほぼ同じでした。記録は algo-notes の notes/string-basics.md にあります。
