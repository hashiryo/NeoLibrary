---
title: パターンの出現回数 (FM-index)
documentation_of: ../../neo/string/FMIndex.hpp
---

文字列の中のパターンの出現回数を答えます。

```cpp
explicit FMIndex(const std::string& s);
int size() const;
int count(const std::string& p) const;  // s[i, i + |p|) = p となる i の数
```

構築は接尾辞配列を作るので $O(n)$ です。`count(p)` は、パターンが 32 文字以下なら $O(|p|)$、33 文字以上なら $O(|p| \log n)$ です。空のパターンには $n + 1$ を返します。文字は unsigned char として扱います。

接尾辞配列から、末尾に最小の番兵を足した文字列の BWT を作り、FM-index (Ferragina と Manzini、2000) で答えます。パターンを後ろから 1 文字ずつ読み、パターンの接尾辞で始まる接尾辞の範囲を、文字ごとの rank で絞ります。rank の表は、文字ごとに BWT の位置の bit の列と、64 個ごとのその文字の累積の数を 16 byte の組で持ち、1 回の読み込みと popcount で引きます。

文字は使われているものに番号を振り直し、3 文字の組ごとに、それより小さい接尾辞の数も持ちます。3 文字以下のパターンはこの表だけで答え、それより長いパターンは、最後の 3 文字の範囲から絞り始めます。表は、文字が $k$ 種類のとき $(k + 1)^3$ 個で、$(k + 1)^3 \le 4n$ のときだけ作ります。33 文字以上のパターンは、最初の 3 文字の範囲の中で、接尾辞配列を二分探索します (8 byte ずつ比べます)。長いパターンでは、二分探索の $\log n$ 回の比べは前から続けて読むだけなので、FM-index の 1 文字ごとの飛び飛びの読み込みより速くなります。

procon-judge の self-string-count-occurrences (長さ $5 \times 10^5$ のテキストに多数のパターン) で書き比べました。EPYC 7763 の x64-gcc では、構築とクエリを合わせた時間がケースごとの最大で 33.3 ms でした。同じ回に、今の Library の SuffixArray::pattern_matching は 150.6 ms、26 個の遷移の欄を持つ接尾辞オートマトンは 382.6 ms でした。接尾辞オートマトンが遅いのは、パターンの 1 文字ごとに 104 MB の遷移の表を飛び飛びに読むためで、表を huge page に置いても 251.3 ms でした。3 文字の表を持たずに 1 文字目から FM-index で絞る形とは、最大の時間がほぼ同じです。記録は algo-notes の notes/string-basics.md にあります。
