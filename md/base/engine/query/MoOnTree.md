---
title: 木の上の Mo's algorithm
documentation_of: ../../neo/base/engine/query/MoOnTree.hpp
---

木のパスのクエリを、オイラーツアーの括弧列に落として処理する Mo's algorithm です。木を DFS して各頂点を 2 回（入るときと出るとき）並べた長さ $2N$ の列を作り、パス $u$ から $v$ までをその列の 1 つの区間に対応させます。区間の端点を 1 つ動かすと、その位置の頂点が入ると出るを交互に繰り返すので、`add` と `erase` ではなく `toggle(v, entering)` を渡します。

```cpp
struct MoOnTree {
  explicit MoOnTree(int n);
  void add_edge(int u, int v);
  void query(int u, int v);
  void build(int root = 0);                   // 辺とクエリを足したあとに 1 回

  bool is_on(int v) const;
  template <class Toggle, class Out>
  void run(const Toggle& toggle, const Out& out);
  template <class Toggle, class Out>
  void run_hilbert(const Toggle& toggle, const Out& out);
};
```

`toggle(v, entering)` は、`entering` が真なら頂点 `v` が区間に入るとき、偽なら出るときです。`out(id)` はクエリの答えを出すときに呼びます。パスに LCA が入らない形の区間では、`out` の直前に LCA を足して外すので、渡す側は LCA を気にしなくてかまいません。`build` は LCA の前計算（$O(N\log N)$）とクエリの区間への変換をします。木が連結であることを前提にします。

```cpp
// 木のパスの異なる色の個数
MoOnTree mo(n);
for (int i = 0; i < n - 1; ++i) {
  int u, v;
  std::cin >> u >> v;
  mo.add_edge(u, v);
}
for (int i = 0; i < q; ++i) {
  int u, v;
  std::cin >> u >> v;
  mo.query(u, v);
}
mo.build(0);
std::vector<int> cnt(V, 0), ans(q);
int distinct = 0;
auto toggle = [&](int v, bool entering) {
  if (entering) {
    if (cnt[color[v]]++ == 0) ++distinct;
  } else {
    if (--cnt[color[v]] == 0) --distinct;
  }
};
auto out = [&](int id) { ans[id] = distinct; };
mo.run(toggle, out);
```

$N$ を頂点数、$Q$ をクエリ数、$T$ を `toggle` 1 回の計算量とすると、$O((N + Q)\sqrt{N}\,T + N\log N)$ です。LCA は二分探索で 1 クエリ $O(\log N)$ です。`MoOnTreeWithUpdates` は同じ木の上に時間軸を足したもので、`query(u, v, t)` と `run(toggle, apply, undo, out)` を使います。`apply(t, on)` と `undo(t, on)` の `on` は、そのときどの頂点が区間に入っているかです。
