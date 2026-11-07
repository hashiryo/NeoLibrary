---
title: 接尾辞オートマトン
documentation_of: ../../neo/string/SuffixAutomaton.hpp
---

文字を 1 つずつ末尾に足しながら、文字列の接尾辞オートマトンを作ります。文字列ごとに根から足し直すと、複数の文字列の広義の接尾辞オートマトンになります。

```cpp
template <int K = 26, char Base = 'a'> class SuffixAutomaton;
SuffixAutomaton();
explicit SuffixAutomaton(const std::string& s);
explicit SuffixAutomaton(const std::vector<std::string>& ss);  // 広義
void push(char c);                     // c を今の文字列の末尾に足す
void add(const std::string& s);        // s を根から足す
int size() const;                      // 状態の数 (状態 0 は空の文字列)
int last() const;                      // 今の文字列 (最後に足した文字列) 全体を受理する状態
int len(int v) const;                  // 状態 v の文字列の最長の長さ
int link(int v) const;                 // 状態 v の suffix link (根は -1)
int next(int v, char c) const;         // 状態 v から c で移る先 (無ければ -1)
long long distinct_substrings() const; // 足した文字列の、異なる空でない部分文字列の数
```

文字は `Base` から `Base + K - 1` までで、既定は英小文字です ($K \le 255$)。状態は文字の合計の 2 倍以下で、`push` はならして $O(1)$ です。`distinct_substrings` は文字を足すたびに足していくので $O(1)$ です。

`add(s)` は今の文字列を空に戻してから `s` を足します。足す先の状態が既にその文字の遷移を持つとき (前に足した文字列と同じ接頭辞を足すとき) は、新しい状態を作らずに、長さが合えばその先へ進み、合わなければ状態を複製してそこへ進みます。1 本の文字列では、末尾の状態は遷移を持たないので、この判定は成り立ちません。

遷移の持ち方は遷移の数で変えます (adaptive radix tree の Node4 と Node256 の考え方)。状態は 32 byte で、suffix link と最長の長さのほかに、4 本までの遷移を、文字 4 つを詰めた 1 語と行き先 4 つで持ち、語の中の byte の一致で文字を探します。5 本目の遷移が来たら、$K$ 個の欄の表を別の列に取ってそこへ移します。状態の遷移はほとんど 1〜2 本なので、$K$ 個の欄を全部の状態に持つ形 (26 種類で 112 byte) より表がずっと小さくなります。列は 2 MB 以上なら、Linux では huge page を頼みます。

procon-judge の self-string-online-distinct-substrings で書き比べました。長さ $5 \times 10^5$ 以下の文字列へ 1 文字ずつ足し、足すたびに異なる部分文字列の数を答える問題です。EPYC 7763 の x64-gcc では、同じ持ち方で、列を 4 KB のページの vector に置くと最大 63.1 ms、huge page の大きい仮想の領域で写し直さずに伸ばすと 55.7 ms でした。同じ回に、26 個の遷移の欄を全部の状態に持つ形は 191.4 ms、遷移を辺の連結リストで持つ形は 243.3 ms でした。パターンの出現回数を数えるだけなら、[FMIndex](FMIndex.md) の方が速くなります。接尾辞オートマトンは、パターンの 1 文字ごとの遷移の読み込みが飛び飛びになるためです。

広義の形は、procon-judge の self-string-multi-distinct-substrings (複数の文字列の異なる部分文字列の数、文字の合計 $10^6$ 前後) で書き比べました。EPYC 9V45 の x64-gcc で最大 78.8 ms でした。同じ回に、文字列を区切り文字でつないだ文字列の接尾辞配列と LCP 配列で数える形は 41.0 ms で、全体を先に見られるなら接尾辞配列の方が速くなります。記録は algo-notes の notes/string-basics.md にあります。
