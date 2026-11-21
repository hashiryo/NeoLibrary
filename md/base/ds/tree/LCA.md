---
title: 最小共通祖先 (O(1))
documentation_of: ../../neo/base/ds/tree/LCA.hpp
---

木の最小共通祖先を $O(1)$ で答えます。

```cpp
explicit LCA(const HeavyLightDecomposition& hld);
int operator()(int u, int v) const;  // u と v の最小共通祖先
```

木は [HeavyLightDecomposition](HeavyLightDecomposition.md) で渡します。中身は写すので、作ったあとに捨ててかまいません。同じ木の 2 頂点で使います。作るのは $O(N)$ です。

## 中身

HLD の位置は行きがけ順なので、位置 i の頂点の親の位置を $A_i$ に並べ、[RangeMin](../data_structure/RangeMin.md) を載せます。u と v の位置を $x < y$ とすると、$A_{x+1}, \ldots, A_y$ の最小の位置の頂点が最小共通祖先です。u が v の祖先なら u の位置が、そうでなければ最小共通祖先から v へ向かう子の親として、最小共通祖先の位置が範囲に現れます。範囲の頂点はどれも最小共通祖先の真の子孫なので、親の位置がそれより前に行くことはありません。表は $N = 5 \times 10^5$ で 8 MB ほどです。

HLD の `lca` は鎖を上がるので $O(\log N)$ です。procon-judge の yosupo-lca ($N, Q \le 5 \times 10^5$、木を作るところから測る) で、LCA と HLD の `lca` を比べました。Xeon 6975P-C の x64-gcc で、二分木の形のケースは 26.3 ms と 83.0 ms です (今の Library の HeavyLightDecomposition は 178.3 ms)。一直線の木のケースは 20.9 ms と 14.1 ms で、HLD の方が速くなります。記録は algo-notes の notes/rmq.md と notes/tree-queries.md にあります。
