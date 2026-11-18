---
title: 全方位木 DP
documentation_of: ../../neo/graph/Rerooting.hpp
---

木 (森でもよい) の頂点 v ごとに、v を根にしたときの v の値を $O(N)$ で求めます。関数の受け方は今の Library の Rerooting と同じで、木は辺の列で渡します。

```cpp
template <class T> class Rerooting;
template <class U, class PE, class OP, class PV>
Rerooting(int n, const std::vector<std::array<int, 2>>& edges, const PE& put_edge, const OP& op, const U& ui, const PV& put_vertex);
const T& operator[](int v) const;           // v を根にしたときの v の値
const T* begin() const;                     // 頂点の順
const T* end() const;
const T& operator()(int root, int v) const;  // root を根にしたときの v の部分木の値
std::size_t size() const;
```

渡す関数は次の 3 つと、`op` の単位元 `ui` です。`op` は可換で結合的なものに限ります。`T` と `U` は既定のコンストラクタで作れる型にします。

| 関数 | 中身 |
| --- | --- |
| `put_edge(v, e, x) -> U` | 辺 e の向こうの部分木の値 x (型 T) を、v から見た寄与にする |
| `op(l, r) -> U` | 寄与をまとめる |
| `put_vertex(v, s) -> T` | v に集まった寄与 s から v の値を作る |

```cpp
// 頂点ごとの、木の中で最も遠い頂点までの辺の数
std::vector<std::array<int, 2>> edges(n - 1);
auto put_edge= [&](int, int, int d) { return d + 1; };
auto op= [](int l, int r) { return std::max(l, r); };
auto put_vertex= [](int, int d) { return d; };
Rerooting<int> dp(n, edges, put_edge, op, 0, put_vertex);
for(int x: dp) std::cout << x << '\n';
```

`dp(root, v)` は root と v が同じ木の頂点のときに使います。root が v の隣か内部の根なら $O(1)$ です。それ以外は、初めて呼んだときに $O(N)$ の前処理をし、そのあとは v の次数を d として $O(\log d)$ です。

## 中身

隣接の CSR を組み、根を全部先に置いて BFS した順に頂点を並べ直します。こうすると、どの頂点の子も隣り合う位置に並びます。

1. 下からの段は、位置の逆順に回ります。子の寄与を `put_edge` で作って `op` でまとめ、`put_vertex` で部分木の値を作ります。
2. 上からの段は、位置の順に回ります。頂点ごとに、親の側の寄与と、子の寄与の後ろからの積と前からの積を合わせて、子から見た親の側の値を作ります。

どちらの段も、配列を位置の順に読み書きします。`dp(root, v)` は、root が v の親か内部の根なら下からの段の値を、root が v の子ならその子から見た親の側の値を返します。それ以外は、内部の根で見た行きがけ順の番号を作っておき、root を含む v の子を兄弟の中から二分探索で探します。

procon-judge の yosupo-tree-path-composite-sum ($N \le 2 \times 10^5$、木を作るところから測る) で、作り方を 3 通り書き比べました。

| 作り方 | 最大 (Xeon 6975P-C、x64-gcc) |
| --- | --- |
| 隣接の CSR を組んで BFS する (この形) | 20.12 ms |
| 葉から剥がした順で回り、親の積を途中まで持つ | 22.68 ms |
| 葉から剥がして親を決め、子の CSR で BFS する | 25.63 ms |
| 今の Library の Rerooting | 67.04 ms |
| 和が引き算できることを使う、この問題に特化した形 | 17.88 ms |

記録は algo-notes の notes/tree-queries.md にあります。
