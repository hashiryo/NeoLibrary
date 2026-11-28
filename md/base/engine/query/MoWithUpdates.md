---
title: 時間軸つきの Mo's algorithm
documentation_of: ../../neo/base/engine/query/MoWithUpdates.hpp
---

点更新と区間クエリが混ざった列を、オフラインで処理する Mo's algorithm です。クエリは $(l, r, t)$ の 3 つ組で、$t$ は適用済みの更新の数です。状態は、更新 $0$ から $t-1$ までを適用したあとの列の半開区間 $l \le x < r$ です。並べ方は $l$、$r$、$t$ の 3 本で、ブロックの幅は $N^{2/3}$ です。

```cpp
struct MoWithUpdates {
  void query(int l, int r, int t);            // [l, r), 更新 0..t-1 を適用した状態
  int size() const;
  bool empty() const;

  template <class AL, class AR, class EL, class ER, class AP, class UN, class Out>
  void run(const AL& add_left, const AR& add_right, const EL& erase_left, const ER& erase_right,
           const AP& apply, const UN& undo, const Out& out) const;
};
```

`add_left`、`add_right`、`erase_left`、`erase_right` は区間の端点を動かすときの処理です。`apply(t, l, r)` は更新 `t` を適用するとき、`undo(t, l, r)` は取り消すときに呼びます。$l$ と $r$ はそのときの区間なので、更新の位置が区間に入るかを渡す側で判定します。区間の移動が先で、時刻の移動が後です。

```cpp
// 点更新つきの、区間の異なる値の個数
std::vector<int> a(n);
std::vector<std::tuple<int, int, int>> upd;  // (位置, 古い値, 新しい値)
MoWithUpdates mo;
for (int i = 0; i < q; ++i) {
  int l, r, t;
  std::cin >> l >> r >> t;
  mo.query(l, r, t);
}
std::vector<int> cnt(V, 0), ans(q);
int distinct = 0;
auto add = [&](int i) {
  if (cnt[a[i]]++ == 0) ++distinct;
};
auto erase = [&](int i) {
  if (--cnt[a[i]] == 0) --distinct;
};
auto apply = [&](int t, int l, int r) {
  auto [p, old, nv] = upd[t];
  if (l <= p && p < r) {
    if (--cnt[old] == 0) --distinct;
    if (cnt[nv]++ == 0) ++distinct;
  }
  a[p] = nv;
};
auto undo = [&](int t, int l, int r) {
  auto [p, old, nv] = upd[t];
  if (l <= p && p < r) {
    if (--cnt[nv] == 0) --distinct;
    if (cnt[old]++ == 0) ++distinct;
  }
  a[p] = old;
};
auto out = [&](int id) { ans[id] = distinct; };
mo.run(add, add, erase, erase, apply, undo, out);
```

$N$ を列の長さ、$Q$ をクエリ数とすると、$N = Q$ のときの移動回数は $O(N^{5/3})$ です。`apply` と `undo` は更新 1 回ぶんの処理で、位置が区間に入るときだけ答えを直します。
