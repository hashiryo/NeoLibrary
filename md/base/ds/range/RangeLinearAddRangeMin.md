---
title: 区間に 1 次関数を足せる区間最小
documentation_of: ../../neo/base/ds/range/RangeLinearAddRangeMin.hpp
---

整数の列 a について、区間 `[l, r)` の各 i へ `b i + c` を足すことと、区間の最小を引くことを、交互に行えます。

```cpp
class RangeLinearAddRangeMin;
explicit RangeLinearAddRangeMin(const std::vector<long long>& a);
int size() const;
void add(int l, int r, long long b, long long c);  // i in [l, r) に a[i] += b i + c。l >= r なら何もしない
long long min(int l, int r);                       // min(a[l], ..., a[r - 1])。l < r
long long get(int i) const;                        // a[i]
```

足すときの i は列の先頭からの位置です。区間の先頭からの位置で足したいときは、`c` を `c - b l` にします。計算量は、構築が $O(N \log N)$、`add` と `min` が償却で $O(\log^2 N)$、`get` が $O(\log N)$ です。途中の値の絶対値と、足した `b` の累計の絶対値に N を掛けたものが、どちらも $2^{61}$ ほどに収まる範囲で使います。

## 中身

点 (i, a_i) の下側の凸包を、節ごとに左右の子の凸包を結ぶ橋 (共通の接線の 2 点) で持つ二分木です (Overmars と van Leeuwen の形)。区間に足すのは覆う節に作用を積むだけで、1 次関数を足しても点の並びの向きは変わらないので、その節より下の橋は変わりません。端の道の節には「橋が古い」印だけを付け、橋をたどるときに古い子から探し直します。

この木に、kinetic segment tree の考え方を重ねます。要素 i を傾き i の直線と見ると、区間に `b i` を足すのは、その区間の時刻を b 進めることにあたります。節は、部分木の最小とそれを取る葉、時刻を前と後ろへどれだけ進めてもその葉が最小のままか (melt) を持ちます。作用が melt に収まる間は、節の最小をずらすだけで済みます。収まらなかった節は「最小が古い」とだけ記し、`min` で要るときに次の順で求め直します。

1. kinetic segment tree と同じく、古い子へ再帰して 2 つの子から取り直す。
2. 古い節を 32 個より多く見たら 1 をやめ、橋をたどって 1 本の道で最小の葉へ下りる。道の上の橋の傾きから、その葉が最小のままでいる幅もちょうど決まる (その葉の両隣の凸包の辺は、どれも道の上のどこかの節の橋になっている)。

kinetic segment tree だけでは、時刻を前後に行き来させる入力で計算量の保証がありません。例えば $a_i = -i$ に全体へ $+2i$ と $-2i$ を交互に足すと、どの節でも最小を取る葉が毎回入れ替わり、1 回が $O(N)$ になります。ここでは 1 節の求め直しが $O(32 + \log N)$ で済み、橋の探し直しの回数は `add` で付けた印の数で抑えられます。

方式は procon-judge の yosupo-range-linear-add-range-min ($N, Q \le 10^5$) で比べて選びました。時刻を行き来させる入力は、自作の self-segtree-range-linear-add-range-min に入れてあります。この型を使う提出は、EPYC 7763 の x64-gcc で Library Checker のケースの最大が 73.3 ms でした。同じ回で、Library Checker の首位と同じ kinetic segment tree は 91.8 ms、橋だけで持つ形は 239.0 ms です。x64-clang でも、この型が 87.3 ms、kinetic segment tree が 117.2 ms でした。自作の問題では、kinetic segment tree が時刻を行き来させるケースで 10 秒を越え、この型は最大 350.6 ms (橋だけで持つ形は 271.6 ms) でした。記録は algo-notes の notes/segment-tree-design.md にあります。
