---
title: 区間に 1 次関数を足せる kinetic segment tree
documentation_of: ../../neo/base/ds/range/KineticSegmentTree.hpp
---

各要素が直線 $(x_i, y_i)$ の列について、区間 `[l, r)` の各 i へ `a x_i + b` を足すことと、区間の最大 (または最小) を引くことを、交互に行えます。足す `a` は 0 以上です。

```cpp
template <bool MAX = true> class KineticSegmentTree;
explicit KineticSegmentTree(const std::vector<long long>& x, const std::vector<long long>& y);
int size() const;
void add(int l, int r, long long a, long long b);  // 0 <= a。i in [l, r) に y_i += a x_i + b
void set(int i, long long x, long long y);         // (x_i, y_i) = (x, y)。y は今の時刻の値
long long prod(int l, int r);                      // MAX なら max(y_l, ..., y_{r - 1})、でなければ最小。l < r
long long all_prod();                              // 全体の最大 (最小)
long long get(int i);                              // y_i
```

`MAX = true` が最大、`false` が最小です。構築は $O(N)$ で、`add` と `prod` は償却です。入れ替わりの回数は全体で $O(N \log N)$、1 回の取り直しが $O(\log N)$ です。$x_i$, $y_i$, `nxt`, $a x_i + b$ が 64 bit に収まる範囲で使います。`a` が負のときは使えません。

## 中身

節は、部分木で勝っている直線 (傾きと切片)、その今の値、値が変わりうるまでに足せる `a` の幅 `nxt` を持ちます。`a` が `nxt` より小さい間は、節の値を `a x + b` だけずらして `nxt` から `a` を引くだけで済みます。幅に収まらない作用は子へ配って降り、2 つの子から勝っている直線を取り直します (Segment Tree Beats と同じ形)。幅は作用のたびに減るので、`a` は 0 以上でないと使えません。

区間の最大は、節の値が常に今の値なので、覆う節の値を比べるだけで取れます。端で子へ降りるときだけ、積んだ作用を配ります。点の値を取るときは、葉までの道の作用を配ってから葉の値を返します。`set` は今の時刻の値を置き換えるので、道の作用を配ってから葉を置き、上へ取り直します。

方式は procon-judge の self-kinetic-range-add-range-max / -min ($N, Q \le 10^5$) で、手書きの同じ形と比べました。手元の Apple M2 Max の local で、最大版の 10 ケースの最大が、手書きの beats で 79.5 ms、この型を使う提出で 83.1 ms でした。同じ問題には、前後両向きの幅を持つ melt、平方分割と塊ごとの凸包で持つ参照実装も置いてあります。記録は algo-notes の notes/kinetic.md にあります。
