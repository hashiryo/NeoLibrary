---
title: 木の道の積
documentation_of: ../../neo/base/ds/tree/PathProduct.hpp
---

木の頂点に値を置き、1 点の代入と、道の積を答えます。積は可換でなくてよく、u から v への道の順に掛けます。どちらも $O(\log N)$ です。

```cpp
template <class M> class PathProduct;
PathProduct(const HeavyLightDecomposition& hld, const std::vector<T>& a);  // a は頂点ごとの値
explicit PathProduct(const HeavyLightDecomposition& hld);                   // 値はすべて ti()
int size() const;
T get(int v) const;
void set(int v, T x);
T prod(int u, int v) const;  // op(a[u], ..., a[v]) (u から v への道の順)
```

構造体 M は SegmentTree と同じモノイド (`using T`、`static T ti()`、`static T op(T a, T b)`) です。木は [HeavyLightDecomposition](HeavyLightDecomposition.md) で渡します。中身は写すので、作ったあとに捨ててかまいません。`prod` は同じ木の 2 頂点で使います。

```cpp
// 1 次関数の合成 (vertex_set_path_composite)
using M = monoid::Composite<ZMod<998244353>>;
PathProduct<M> pp(HeavyLightDecomposition(n, edges), f);
pp.set(p, {c, d});
auto g = pp.prod(u, v);  // u を先に当てる
```

## 中身

HLD の鎖ごとに二分木を組みます。中国の競プロで全局平衡二叉树と呼ばれる形で、Library Checker の vertex_set_path_composite の首位 (chaihf) と同じ組み方です。HLD とセグメント木で書くと、道は $O(\log N)$ 本の鎖に分かれ、鎖ごとに $O(\log N)$ かかるので、合わせて $O(\log^2 N)$ です。

組み方は次のとおりです。

1. 軽い子を先に、重い子を最後に配る行きがけ順の位置を作る。鎖の隣どうしの位置の差は、上の頂点と、そこにぶら下がる軽い子の部分木の大きさの和になる。
2. 鎖の頂点の優先度を「位置と、鎖の次の頂点の位置の xor の最上位 bit」として、優先度の大きい方を根にした位置の順の二分木を組む。大きい部分木がぶら下がる頂点ほど根に近くなる。
3. ある頂点から鎖の二分木の根まで上がる段数は、おおよそ log (鎖全体の大きさ ÷ その頂点の重み) になる。軽い辺で上の鎖へ移ると、入った先の頂点の重みは下の鎖の全体の大きさ以上ある。そのため道に沿って段数を足すと打ち消し合って、$O(\log N)$ に収まる。

節は、自分の値と、部分木の前向きの積と後ろ向きの積、左の部分木と自分の前向きの積と後ろ向きの積を持ちます。鎖の頂点から鎖の先頭までの積は、二分木の根へ上がりながら、位置がその頂点以下の祖先の積を集めて求めます。1 点の代入は、二分木の根まで上がりながら、節ごとに op を 4 回します。

procon-judge の yosupo-vertex-set-path-composite ($N, Q \le 2 \times 10^5$) で比べました。相手は、HLD と `SegmentTree<monoid::Composite>` 2 本で書いた形です。Xeon 6975P-C の x64-gcc で、軽い辺を多く通るケースでは 106.4 ms と 153.7 ms でした (今の Library の HLD とセグメント木 2 本は 238.2 ms)。一方、一直線の木で代入が半分を占めるケースでは 104.2 ms と 54.6 ms で、セグメント木の方が速くなります。記録は algo-notes の notes/tree-queries.md にあります。
