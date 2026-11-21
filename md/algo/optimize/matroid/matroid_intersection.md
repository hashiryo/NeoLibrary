---
title: マトロイドの交差
documentation_of: ../../neo/algo/optimize/matroid/matroid_intersection.hpp
---

台集合が同じ 2 つのマトロイド $M_1$ と $M_2$ (要素は $0$ から $n-1$ の番号) について、両方で独立な集合 (共通独立集合) を求めます。

```cpp
template <class M1, class M2> std::vector<int> matroid_intersection(M1&& m1, M2&& m2);
template <class M1, class M2, class T>
std::pair<std::vector<T>, std::vector<int>> weighted_matroid_intersection(M1&& m1, M2&& m2, const std::vector<T>& w, int k = -1);
```

- `matroid_intersection(m1, m2)` は、要素数が最大の共通独立集合を、要素の番号の列で返します。
- `weighted_matroid_intersection(m1, m2, w, k)` は、要素 $x$ の重みを `w[x]` として、要素数 $j$ の共通独立集合の重みの和の最大 `best[j]` ($j = 0, 1, \dots$) と、要素数が最大で重みの和が `best.back()` の集合を返します。`k >= 0` なら要素数 `k` で打ち切ります。重みの和の最小は、符号を反転した重みを渡して求めます。

マトロイドの型は、[matroids.hpp](./matroids.md) の型か、次を持つ自作の型です。

```cpp
int size() const;                       // 台集合の大きさ n
void build(const std::vector<int>& I);  // 独立集合 I に合わせて中の構造を作り直す
bool can_add(int x);                    // x が I に無いとき、I + x が独立か
// 次の 2 つのうち、書きやすいほう (両方でもよい)
bool can_exchange(int y, int x);        // y が I にあり x が I に無いとき、I - y + x が独立か
template <class F> void circuit(int x, F f);  // I + x が独立でないとき、I + x の閉路に入る I の要素 y ごとに f(y)
```

足りないほうはライブラリが補います。

- `can_exchange` が無い型: x の閉路を `build` ごとに覚えて、y が入っているかを答えます。
- `circuit` が無い型: I の要素を全部 `can_exchange` で調べて閉路を作ります。

型によって自然に書けるほうが違います。たとえば横断マトロイドは閉路を交互路の探索 1 回で求められ、グラフ的マトロイドの双対は 2 引数の判定を $O(1)$ で書けます。

## 中身

- 重みなし: 2 引数の判定で交換グラフの層を BFS で作り、層に沿った最短の増加路を DFS で次々に探します (Cunningham の方法)。$M_2$ への問いはどれも x を固定して y を変える形にしてあるので、`can_exchange` を持たない型は $M_2$ の側に回します。判定の回数は $O(n r^{1.5})$ 回ほどです ($r$ は答えの要素数)。
- 重み付き: 閉路の列挙で交換グラフを作ります。要素の長さ (I に無い $x$ は $-w(x)$、I にある $y$ は $w(y)$) の和と辺の数の組で最小の路を、Bellman-Ford (待ち行列) で求めて 1 本ずつ増やします。`circuit` を持たない型を渡すと、増加路 1 本ごとに $O(n r)$ 回の判定が要るので遅くなります。

方式は procon-judge の self-matroid-graphic-partition、self-matroid-xor-partition、self-matroid-weighted-graphic-partition で書き比べて選びました。重み付きでは、Frank の重みの分け方をポテンシャルにした Dijkstra より、待ち行列の Bellman-Ford のほうが速くなりました。
