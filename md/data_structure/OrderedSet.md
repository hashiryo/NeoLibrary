---
title: 順序付き集合 (B+ 木)
documentation_of: ../../neo/data_structure/OrderedSet.hpp
---

int の順序付き集合です。`OrderedSet<true>` は k 番目の要素と x 以下の個数も引けます。`OrderedSet<false>` はそれを持たない分、入れる、除くが少し安くなります。

```cpp
template <bool Counted = true> class OrderedSet;
OrderedSet();
explicit OrderedSet(const std::vector<int>& a);  // a の要素からなる集合。a は狭義単調増加
int size() const;
bool empty() const;
bool contains(int x) const;
bool insert(int x);                     // 入ったら true、既にあれば false
bool erase(int x);                      // 除いたら true、無ければ false
std::optional<int> prev(int x) const;   // x 以下で最大の要素
std::optional<int> next(int x) const;   // x 以上で最小の要素
// Counted = true のときだけ
std::optional<int> kth(int k) const;    // 小さい方から k 番目 (0 始まり)。範囲の外なら nullopt
int count_lt(int x) const;              // x 未満の要素の個数
int count_le(int x) const;              // x 以下の要素の個数
```

キーは int の全域を使えます。どの操作も $O(\log M)$ で、M はそれまでに入れた要素の延べの数 (初めの集合を含む) です。木の 1 段は 64 分岐です。構築は $O(N)$ です。コピーと move ができます。

## 使い分け

| 型 | 入れられる値 | メモリ | 引けるもの |
| --- | --- | --- | --- |
| OrderedSet`<true>` | int の全域 | 要素の数に比例 | 入っているか、前後の要素、k 番目、x 以下の個数 |
| OrderedSet`<false>` | int の全域 | 要素の数に比例 | 入っているか、前後の要素 |
| [FastSet](FastSet.md) | 0 以上 n 未満 (n は作るときに決める) | 値の範囲 n に比例 (n = $10^7$ で 1.25 MB) | 入っているか、前後の要素 |
| [HashMap](HashMap.md) | u64 の全域 | 要素の数に比例 | 入っているか (順序は無い) |

std::set の代わりに使います。値は int の全域を取れて、メモリは入っている要素の数で決まります。k 番目の要素や x 以下の個数が要るなら `OrderedSet<true>` にします。中央値を保つ、値の順位を答える、といった使い方です。前後の要素だけなら `OrderedSet<false>` にして、個数を直す分を払わないようにします。

値が 0 以上 n 未満 (n は $10^7$ ほどまで) に収まり、前後の要素だけが要るなら、[FastSet](FastSet.md) の方がずっと速く動きます。値を全部先に知っていれば、座標圧縮してから FastSet を使う手もあります。値がオンラインで来て範囲も広いときが、OrderedSet の出番です。入っているかだけを見るなら [HashMap](HashMap.md) で足ります。

## 中身

B+ 木です。葉はキーを 63 個まで昇順に持って前後の葉とつなぎ、内部の節点は子を 63 個まで持って、子ごとに部分木のキーの最大を持ちます。キーで下りるときは、最大の列のうち x 未満のものの数を AVX2 で 8 個ずつ比べて数え、子を選びます。入れて 64 個になった節点は半分に割り、除いて空になった節点は親から外します。小さくなった節点を隣と併せることはしません。除いても高さは増えないので、高さはそれまでに入れた延べの数で決まります。

`Counted = true` では、内部の節点が子ごとの部分木の個数と、子を 8 個ずつに分けた組の累積を持ちます。k 番目は、組の累積 8 個を 1 回比べて組を選んでから、組の中の累積をレジスタの中で作って子を選びます。x 以下の個数は、組の累積に組の中の手前の和を足します。どちらも分岐を使いません。入れるときと除くときは、通った節点ごとに子の個数 1 つと組の累積 8 個を足し引きします。

方式は procon-judge の yosupo-ordered-set で比べて選びました ($N, Q \le 5 \times 10^5$、値は $10^9$ 以下)。同じ方式の bptree_c8b は、EPYC 9V45 の x64-gcc で最大 35 ms ほどでした。`__gnu_pbds::tree` と今の Library の SegmentTree_Patricia の 14 分の 1、重み平衡木の半分です。記録は algo-notes の notes/associative-containers.md にあります。
