---
title: HL 分解
documentation_of: ../../neo/graph/HeavyLightDecomposition.hpp
---

根付き木 (森でもよい) を HL 分解します。頂点を、重い子を先にした行きがけ順の位置に並べます。道は O(log N) 個の位置の区間に、部分木は 1 つの区間になります。関数の名前は今の Library の HeavyLightDecomposition と同じです。

```cpp
HeavyLightDecomposition(int n, const std::vector<std::array<int, 2>>& edges, int root = 0);  // 無向の辺の列から
explicit HeavyLightDecomposition(const std::vector<int>& par);  // 親の列から (根は par[v] < 0)
int to_seq(int v) const;     // 頂点 v の位置
int to_vertex(int i) const;  // 位置 i の頂点
std::size_t size() const;
int parent(int v) const;     // 根なら -1
int head(int v) const;       // v を含む鎖の先頭の頂点
int root(int v) const;
bool connected(int u, int v) const;
bool in_subtree(int u, int v) const;       // u が v の部分木にあるか
int subtree_size(int v) const;
std::pair<int, int> subtree(int v) const;  // v の部分木の位置の半開区間 [l, r)
int depth(int v) const;
int lca(int u, int v) const;
int la(int v, int k) const;           // v の k 個上の祖先 (0 <= k <= depth(v))
int jump(int u, int v, int k) const;  // u から v への道の k 番目 (u が 0 番目) の頂点。無ければ -1
int dist(int u, int v) const;
template <class F> void path(int u, int v, const F& f, bool edge = false) const;
```

辺の列から作るときは、root を含む木を root を根にして位置 0 から置きます。ほかの木の根は決めません。親の列から作るときは、番号の小さい根の木から順に位置を配ります。`lca`、`dist`、`jump`、`path` は同じ木の 2 頂点で使います。作るのは $O(N)$、`lca` と `path` は $O(\log N)$ です。

`path(u, v, f, edge)` は、u から v への道を位置の閉区間に分け、道の順に `f(l, r)` を呼びます。u から lca へ上る区間は `l > r` で、位置 l から r へ向かいます。lca から v へ下りる区間は `l <= r` です。`edge` が真なら lca を除きます (辺の値を子の頂点に置く使い方)。今の Library の `path` は区間の vector を返していましたが、呼ぶたびの確保をなくすため、区間ごとに関数を呼ぶ形にしました。

```cpp
// 道の和。区間和の木は位置の順に値を持つ (位置 i の頂点は to_vertex(i))
SegmentTree<monoid::Sum<long long>> seg(n, [&](int i) { return a[hld.to_vertex(i)]; });
long long s = 0;
hld.path(u, v, [&](int l, int r) { s += seg.prod(std::min(l, r), std::max(l, r) + 1); });
```

可換でない積は、位置の順の木と逆順の木の 2 本を持ち、`l > r` の区間を逆順の木で取ります。1 点の代入と道の積だけなら [PathProduct](PathProduct.md) が使えます。lca を何度も引くなら、O(1) で答える [LCA](LCA.md) があります。

## 中身

作るときは、子が親より前に来る順を作ってから、その順で部分木の大きさと重い子を決め、逆順で位置を配ります。無向の辺からは隣接を組みません。頂点ごとに次数と隣の番号の xor を持ち、次数 1 の頂点を葉から剥がします (剥がすときの xor が親)。親の列からは、$p_i < i$ (根が 0) かを分岐の無い 1 回の走査で確かめ、そうなら番号の大きい順をそのまま使います。そうでなければ、子の CSR を組んで根から BFS します。

位置ごとに、鎖の先頭の位置と先頭の親の位置を 8 byte の組で持ちます。`lca` と `path` は、鎖を 1 本上がるたびにこの組を 1 つ読むだけです。配列は 1 つの領域にまとめ、1 MB 以上なら 2 MB 境界に置いて huge page を頼みます (Linux のみ)。

procon-judge で、作り方と問い合わせの持ち方を書き比べました ($N, Q \le 5 \times 10^5$)。どの問題も木を作るところから測ります。Xeon 6975P-C の x64-gcc の最大は次のとおりです。

| 問題 | この型 | 今の Library |
| --- | --- | --- |
| yosupo-lca (`lca`) | 84.3 ms | HeavyLightDecomposition 184.5 ms |
| yosupo-vertex-add-path-sum (`path` と区間和の木) | 89.8 ms | HLD と BIT 220.5 ms |
| yosupo-vertex-add-subtree-sum (`subtree` と区間和の木) | 31.7 ms | HLD と BIT 80.7 ms |

記録は algo-notes の notes/tree-queries.md にあります。
