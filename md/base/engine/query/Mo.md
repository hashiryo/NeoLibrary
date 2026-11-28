---
title: 区間の Mo's algorithm
documentation_of: ../../neo/base/engine/query/Mo.hpp
---

区間クエリをオフラインで処理する Mo's algorithm です。半開区間 $l \le x < r$ のクエリを集め、端点の移動が少なくなる順に並べて、区間が伸び縮みするたびに `add` と `erase` を呼びます。並べ方は、$l$ のブロックと $r$ の向きを交互にするもの（奇偶ソート）と、ヒルベルト曲線の順の 2 つです。

```cpp
struct Mo {
  void query(int l, int r);                    // [l, r)
  int size() const;
  bool empty() const;

  template <class AL, class AR, class EL, class ER, class Out>
  void run(const AL& add_left, const AR& add_right, const EL& erase_left, const ER& erase_right, const Out& out) const;
  template <class A, class E, class O>
  void run(const A& add, const E& erase, const O& out) const;

  template <class AL, class AR, class EL, class ER, class Out>
  void run_hilbert(const AL& add_left, const AR& add_right, const EL& erase_left, const ER& erase_right, const Out& out) const;
  template <class A, class E, class O>
  void run_hilbert(const A& add, const E& erase, const O& out) const;
};
```

`add_left(i)` は区間が左に伸びて `i` を含むとき、`add_right(i)` は右に伸びて含むとき、`erase_left(i)` と `erase_right(i)` は縮んで外すときです。`out(id)` は `query` で `id` 番目に足したクエリの答えを出すときに呼びます。左右で `add` と `erase` が同じなら、2 引数の `run` が使えます。

```cpp
// 区間 [l, r) の異なる値の個数
std::vector<int> a(n);
Mo mo;
for (int i = 0; i < q; ++i) {
  int l, r;
  std::cin >> l >> r;
  mo.query(l, r);
}
std::vector<int> cnt(n + 1, 0), ans(q);
int distinct = 0;
auto add = [&](int i) {
  if (cnt[a[i]]++ == 0) ++distinct;
};
auto erase = [&](int i) {
  if (--cnt[a[i]] == 0) --distinct;
};
auto out = [&](int id) { ans[id] = distinct; };
mo.run(add, erase, out);
```

ブロックの幅は `max(r) / sqrt(q)` です。$N$ を数列の長さ、$Q$ をクエリ数、$T$ を `add` か `erase` 1 回の計算量とすると、`run` は $O((N + Q)\sqrt{Q}\,T)$、`run_hilbert` も同程度で、移動回数の定数が小さいことがあります。クエリの $r$ が大きい問題（$Q$ が $N$ より大きいもの）では、ヒルベルト順が効きます。`run_hilbert` は最初に全クエリの曲線上の位置を計算してから並べます。
