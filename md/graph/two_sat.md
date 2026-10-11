---
title: 2-SAT
documentation_of: ../../neo/graph/two_sat.hpp
---

2 個のリテラルの論理和だけからなる論理式を満たす値を求めます。

```cpp
std::vector<int> two_sat(int n, const std::vector<std::array<int, 2>>& clauses);
```

変数は 0 から n - 1 の番号で、変数 i のリテラルは、肯定が `i`、否定が `~i` (= -1 - i) です。`clauses` の各要素 {a, b} が節 a ∨ b です。すべての節を満たせるなら、各変数の値 (0 か 1) の列を返し、満たせなければ空の列を返します。変数と節の数の和に比例する時間で求めます。

## 中身

含意グラフ (節 a ∨ b から辺 ¬a → b と ¬b → a) の強連結成分で解きます。x と ¬x が同じ成分なら満たせず、そうでなければ、トポロジカル順で x の成分が ¬x の成分より後ろにあるとき x を真にします。

含意グラフの辺の列は作りません。リテラル x の出現の並び occ[x] (x を含む節のもう一方のリテラル) の CSR を 1 つだけ組みます。含意グラフでは、x から出る辺の行き先が occ[¬x]、x の入次数が |occ[x]| になるので、節を 1 回なめて出現を数えれば、CSR の区切りと、剥がしに使う入次数の両方が手に入ります。リテラルごとの配列は鍵の 4 byte だけです。強連結成分の求め方は `strongly_connected_components` と同じく、入次数 0 のリテラルを待ち行列で剥がしてから、残りを path-based (Gabow) の反復の DFS で分けます。

作業領域を小さく持つのは、$N = M = 5 \times 10^5$ の乱択の含意グラフ ($10^6$ 頂点) で、強連結成分分解の持ち方のままだと作業領域が L3 から溢れたためです。procon-judge の yosupo-two-sat で書き比べ、x64 の 5 つの CPU と arm のすべてで首位でした。x64-gcc の最大ケースは EPYC 9V74 で 15.4 ms、Xeon Platinum 8370C で 37.6 ms で、今の Library の TwoSatisfiability の 3.1 倍から 6.0 倍の速さです。記録は algo-notes の notes/graph-basics.md にあります。
