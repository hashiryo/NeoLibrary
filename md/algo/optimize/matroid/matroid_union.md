---
title: マトロイドの和
documentation_of: ../../neo/algo/optimize/matroid/matroid_union.hpp
---

台集合が同じ $k$ 個のマトロイド $M_0, \dots, M_{k-1}$ (要素は $0$ から $n-1$ の番号) について、互いに素な独立集合 $I_0, \dots, I_{k-1}$ を、和集合の要素数が最大になるように分けます。

```cpp
template <class M> std::vector<int> matroid_union(std::vector<M>& ms);
template <class M> std::vector<int> matroid_union(std::vector<M>& ms, const std::vector<int>& order);
```

- `matroid_union(ms)` は、要素 $x$ を入れた集合の番号 `owner[x]` (どれにも入れない要素は $-1$) を返します。入れた要素の数が和の最大の独立集合の大きさです。
- `matroid_union(ms, order)` は、`order` の順に要素を入れ、入らない要素を捨てます。和もマトロイドなので、`order` を重みの大きい順に並べると、入れた要素の重みの和が最大になります。`order` には要素を重複なく並べます (全部でなくてもかまいません)。

`ms` の型は 1 つで、同じ型の $k$ 個を渡します。各マトロイドは `build` などで中の構造を書き換えられます。型は [matroid_intersection](./matroid_intersection.md) に渡す型と同じ (`size`、`build`、`can_add` と、`can_exchange` か `circuit`) で、さらに次の `add` を持てます。

```cpp
void add(int x);  // x が I に無く、I + x が独立のとき、中の構造を作り直さずに x を足し、can_add だけを正しく保つ
```

`add` を持たない型でも動きます。持つ型は、どの集合にもそのまま入る要素を足すたびの `build` を省けるので、グラフ的マトロイドで 10 倍近く速くなりました。[GraphicMatroid、PartitionMatroid、XorMatroid](./matroids.md) が持ちます。

## 使い方

森 $k$ 本の和で、辺を $k$ 本の森に分けます。

```cpp
std::vector<GraphicMatroid> ms(k, GraphicMatroid(n, edges));
std::vector<int> owner = matroid_union(ms);
int covered = std::count_if(owner.begin(), owner.end(), [](int j) { return j >= 0; });
bool all = covered == (int)edges.size();   // 全部の辺を k 本の森に分けられるか
bool trees = covered == k * (n - 1);       // 辺素な全域木が k 本あるか (連結なグラフのとき)
```

重み付きは、重みの小さい順に渡します。辺素な $k$ 本の全域森の重みの和の最小は、次のように求まります。

```cpp
std::vector<int> ord(m);
std::iota(ord.begin(), ord.end(), 0);
std::sort(ord.begin(), ord.end(), [&](int a, int b) { return w[a] < w[b]; });
std::vector<int> owner = matroid_union(ms, ord);
```

## 中身

Edmonds のマトロイドの分割の算法です。要素を 1 つずつ入れ、どの $M_j$ にもそのまま入らなければ、交換の列を幅優先で探します。要素 $u$ を $I_j$ に入れると、$I_j + u$ の閉路の要素 $y$ が押し出されるので、$y$ をほかの $I_{j'}$ に入れ直します。どれかにそのまま入る要素に着いたら、列に沿って入れ替えます。列が見つからない要素は捨てます。

速くする工夫が 2 つあります。

- 見つからなかった探索でたどった要素の集合 $S$ は、以後の探索で飛ばします。$S$ の要素をどの $I_j$ に入れても $S$ の要素しか押し出さず、そのまま入る要素も無く、入れ替えでは $I_j \cap S$ が変わらないので、$S$ は以後のどの最短の列にも現れないためです。外すと 10 倍から 100 倍遅くなりました。
- `order` を渡さないときは、そのまま入る要素を先に全部入れてから、残りを探索で入れます。入れ替えでは各 $I_j$ の張る空間が変わらず、そのまま入る要素で広がるだけなので、一度どの $I_j$ にもそのまま入らなかった要素は以後もそのまま入らないためです。`order` を渡すときは入れる順を変えられないので使いません。

`circuit` を持たない型には、$I_j$ の要素を全部 `can_exchange` で調べて閉路を作ります。

方式は procon-judge の self-matroid-forest-union と self-matroid-xor-union で書き比べて選びました。和を「$k$ 個の写しの直和」と「同じ要素は 1 つまで」の分割マトロイドの交差に直して [matroid_intersection](./matroid_intersection.md) で解く形は、台集合が $k$ 倍になり、10 倍から 1800 倍遅くなりました。和を型にして別の算法に渡す使い方は、調べた問題 (和の 10 問) にありませんでした。
