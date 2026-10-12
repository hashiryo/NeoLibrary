---
title: Aho–Corasick
documentation_of: ../../neo/string/aho_corasick.hpp
---

文字列の列のトライと、各頂点の suffix link を求めます。

```cpp
struct AhoCorasick {
  std::vector<int> par, link, word;
};
AhoCorasick aho_corasick(const std::vector<std::string>& ss);
```

頂点 0 が根で、頂点の番号は文字列を 0 番から順にトライへ入れたときにできる順です。`par[v]` は親、`link[v]` は suffix link (頂点 v の文字列の真の接尾辞のうち、トライにある最長のもの) で、`word[i]` は i 番目の文字列の頂点です。`par[0]` と `link[0]` は 0 です。文字は unsigned char として扱います。文字の数の和を $L$、使われている文字の種類を $k$ として $O(Lk)$ です。

作りは入力の形で 2 つから選びます。1 つ目は、遷移の表を幅 $k$ で全部持つ形です。幅優先の順に頂点を取り出し、子の欄を待ち行列に積んでから、空いた欄を suffix link の行で埋めます。子の有無で分岐しないよう、子の無い欄でも待ち行列の末尾に書いて、末尾は子のときだけ進めます。2 つ目は、1 つ目の子を頂点に直接持ち、枝分かれする頂点にだけ幅 $k$ の行を割り当てる形です。suffix link は親の suffix link からたどって求め、途中で行のある頂点に着いたらその欄で止まります。

1 つ目は、suffix link の先が散らばる入力 (でたらめな文字列を多く並べたもの) で速いです。2 つ目は、ほぼ 1 本道のトライ (周期的な文字列やフィボナッチ文字列) と、文字の種類が多い入力で速く、表も小さくなります。文字の種類が 9 以上か、文字列が 2 本以下か、先頭の長さ 12 の部分文字列の種類が少なければ 2 つ目を使います。

procon-judge の yosupo-aho-corasick (文字の数の和が $10^6$ 以下) で書き比べました。EPYC 9V74 の x64-gcc で、ケースの組ごとの中央値の最大が 10.34 ms でした。今の Library の AhoCorasick から同じ答えを組み直す形は 145 ms で、周期的な入力では 1 GB を超えました。記録は algo-notes の notes/aho-corasick.md にあります。
