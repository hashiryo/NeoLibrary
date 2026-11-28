---
title: 回滚する Mo's algorithm
documentation_of: ../../neo/base/engine/query/RollbackMo.hpp
---

削除で答えを戻せない問題のための Mo's algorithm です。区間を右に伸ばすのは単調にして、左側はブロックの中で毎回足してから巻き戻します。追加だけをして、クエリごとに一時的に足したぶんを `undo()` で外します。

```cpp
struct RollbackMo {
  void query(int l, int r);                   // [l, r)
  int size() const;
  bool empty() const;

  template <class Add, class Undo, class Out>
  void run(const Add& add, const Undo& undo, const Out& out) const;
};
```

`add(i)` は位置 `i` を足します。`undo()` は直前に足したものを 1 つ外します。外す順は LIFO です。各ブロックの終わりに、そのブロックで足したものを全部 `undo()` して空に戻すので、状態は `run` の前と各ブロックの始まりで空にしておきます。

```cpp
// 区間の最頻値の出現回数。削除で最大を戻せないので、足す前の最大を控えて巻き戻す。
std::vector<int> cnt(V, 0), ans(q);
int mx = 0;
std::vector<std::pair<int, int>> hist;  // (足した位置, 足す前の mx)
auto add = [&](int i) {
  hist.emplace_back(i, mx);
  if (++cnt[a[i]] > mx) mx = cnt[a[i]];
};
auto undo = [&]() {
  auto [i, old] = hist.back();
  hist.pop_back();
  --cnt[a[i]];
  mx = old;
};
auto out = [&](int id) { ans[id] = mx; };
RollbackMo mo;
for (int i = 0; i < q; ++i) {
  int l, r;
  std::cin >> l >> r;
  mo.query(l, r);
}
mo.run(add, undo, out);
```

$N$ を列の長さ、$Q$ をクエリ数、$T$ を `add` か `undo` 1 回の計算量とすると、$O(N\sqrt{Q}\,T)$ です。ブロックの幅は $N/\sqrt{Q}$ で、右端はブロックごとに 1 回だけ伸び、左端はクエリごとにブロックの右端まで足して巻き戻します。クエリがブロックの中に収まる場合は、その区間をそのまま足して巻き戻します。
