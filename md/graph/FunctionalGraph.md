---
title: 関数グラフ
documentation_of: ../../neo/graph/FunctionalGraph.hpp
---

頂点 x から `to[x]` へ辺が 1 本ずつ出る有向グラフ (関数グラフ) を扱います。x, f(x), f(f(x)), ... は尾をたどってから閉路を回ります。$f^k(x)$ と、x から k 歩の道を $O(\log N)$ で答えます。今の Library の Period のうち、関数を int の配列で渡す場合を写したものです。

```cpp
explicit FunctionalGraph(const std::vector<int>& to);  // to[x] = f(x) (0 <= to[x] < n)
int jump(int x, long long k) const;                   // f^k(x)
template <class F> void path(int x, long long k, const F& f) const;
int to_seq(int x) const;      // x の位置
int to_vertex(int i) const;   // 位置 i の頂点
int dist_to_cycle(int x) const;  // x から閉路に入るまでの辺の数 (閉路の上なら 0)
int cycle_length(int x) const;   // x が入る閉路の長さ
int size() const;
```

`path(x, k, f)` は、x から k 歩たどる間の k + 1 個の頂点を、位置の区間と回数の組に分けます。区間 `[l, r)` と回数 `c` ごとに、道の順で `f(l, r, c)` を呼びます。閉路を何周もする分は、閉路の区間 1 つに周回の数を付けて 1 回で渡します。今の Library の `path` が返していた 3 つの区間の列と周回の数を、1 つの関数で受ける形です。

```cpp
// x から k 歩たどる間に訪れる頂点の値の和 (同じ頂点は訪れた回数だけ足す)
FunctionalGraph g(to);
std::vector<long long> s(n + 1);
for(int v= 0; v < n; ++v) s[g.to_seq(v) + 1]= a[v];
for(int i= 0; i < n; ++i) s[i + 1]+= s[i];
long long ans= 0;
g.path(x, k, [&](int l, int r, long long c) { ans+= (s[r] - s[l]) * c; });
```

閉路に入るまでと閉路 1 周 (今の Library の `path_upto_cycle`) は、`path(x, g.dist_to_cycle(x) + g.cycle_length(x) - 1, f)` で取れます。状態が int でない漸化式 (今の Library の `Period<T>(f, {x0})`) は受けません。状態を使う側で int に直してから渡します。

## 中身

位置は、閉路ごとに閉路の頂点を f の順に 1 つの区間へ並べます。閉路へ入る木の頂点は、f(x) を親とした木の HL 分解で並べます。鎖を上る (f をたどる) と位置が 1 ずつ増える向きで、重い子を先にした行きがけ順を逆にした順です。

作るときは、まず入次数が 0 の頂点を葉から剥がし、木の頂点を、子が親より前となる順で並べます。剥がれずに残った頂点が閉路です。その順で部分木の大きさと重い子を決め、逆順で位置を配ります。位置ごとに、木の頂点なら鎖の先頭の位置とその先の頂点の位置を、閉路の頂点なら閉路の始まりの位置と長さを持ちます。`jump` と `path` は鎖を 1 本上がるたびにこの組を 1 つ読みます。$N = 2 \times 10^5$ の一様な乱択の関数で、`path` が渡す区間は多くて 12 個でした。
