---
title: パターンの出現位置
documentation_of: ../../neo/string/find_all.hpp
---

テキストの中のパターンの出現位置をすべて求めます。

```cpp
std::vector<int> find_all(const std::string& text, const std::string& pattern);
```

`text[i, i + |pattern|)` が `pattern` と等しい $i$ を昇順に返します。パターンが空なら $0$ から $|\mathit{text}|$ までのすべてを返します。全体で $O(|\mathit{text}| + |\mathit{pattern}|)$ です。文字は byte として比べます。

パターンが 8 文字以下なら、テキストの 32 か所の位置について、パターンの全部の文字を AVX2 で比べて and を取り、そのまま答えにします。9 文字以上なら、先頭の 4 文字と末尾の 1 文字で候補を絞り、候補の位置で残りを 32 byte ずつ比べます。W. Muła の SIMD-friendly algorithms for substring searching の generic な形を広げたものです。一様なテキストでは、候補はほとんど残りません。

同じ文字だけのテキストと長いパターンのように候補が密な入力では、候補ごとに長く比べるので、このままでは $O(nm)$ かかります。そこで、比べた 32 byte の塊の数が $n / 2 + m$ を超えた所から KMP に切り替えます。KMP の一致の判定は分岐で書きます。`k += (t[i] == p[k])` と書くと、次に読む `p[k]` の位置が前の読み込みを待つ鎖になり、出現が密な入力で遅くなるためです。

procon-judge の self-string-pattern-search (テキストの長さ $10^6$) で書き比べました。EPYC 7763 の x64-gcc では、ケースごとの最大が 4.38 ms、合計が 7.70 ms でした。最大は同じ文字だけのテキストに長さ 1000 のパターンを当てるケースで、出現が $10^6$ 個近くあります。同じ回に、KMP だけで探す形は最大 8.44 ms、今の Library の z_algorithm で P と T を繋いだ列の Z 配列を作る形は 14.28 ms でした。C ライブラリの memmem を見つけた位置の次から呼び直す形は、呼ぶたびにパターンの前処理をし直すので、出現が密なケースで 3 秒かかりました。記録は algo-notes の notes/string-basics.md にあります。
