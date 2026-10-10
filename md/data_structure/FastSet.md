---
title: 整数の集合 (64 分木)
documentation_of: ../../neo/data_structure/FastSet.hpp
---

0 以上 n 未満の整数の集合です。入れる、除く、入っているかと、前後の要素を引けます。

```cpp
explicit FastSet(int n = 0);                                 // 空の集合
template <class F> FastSet(int n, F f);                      // f(i) が真の i からなる集合
FastSet(int n, const std::vector<unsigned long long>& bits); // bits[i / 64] の下から i % 64 番目が 1 の i からなる集合
bool contains(int i) const;
void insert(int i);
void erase(int i);
int next(int i) const;  // i 以上で最小の要素。無ければ -1
int prev(int i) const;  // i 以下で最大の要素。無ければ -1
```

`contains`、`insert`、`erase` は $0 \le i < n$ で使います。`next` と `prev` には範囲の外の i も渡せます。どの操作も $O(\log_{64} n)$ で、構築は $O(n)$ です。ビットの列から作る形は、`bits` の n 以上の位置の bit を 0 にしておきます。

## 使い分け

| 型 | 入れられる値 | メモリ | 引けるもの |
| --- | --- | --- | --- |
| FastSet | 0 以上 n 未満 (n は作るときに決める) | 値の範囲 n に比例 (n = $10^7$ で 1.25 MB) | 入っているか、前後の要素 |
| [OrderedSet](OrderedSet.md)`<false>` | int の全域 | 要素の数に比例 | 入っているか、前後の要素 |
| [OrderedSet](OrderedSet.md)`<true>` | int の全域 | 要素の数に比例 | 上に加えて、k 番目と x 以下の個数 |
| [HashMap](HashMap.md) | u64 の全域 | 要素の数に比例 | 入っているか (順序は無い) |

値が添字、位置、時刻のように 0 以上 n 未満に収まり、n が $10^7$ ほどまでなら、FastSet を使います。「i 以降でまだ使っていない番号」を探す、補グラフの BFS で残っている頂点を順に取り出す、といった使い方です。メモリは要素の数でなく値の範囲で決まるので、要素が少なくても n bit を使います。

値が大きい、まばら、負になりうるなどで、範囲の大きさの配列を取れないときは OrderedSet を使います。k 番目の要素や x 以下の個数が要るときも OrderedSet です。値を全部先に知っていれば、座標圧縮してから FastSet を使う手もあります。ただ、値がオンラインで来るならその手は使えません。前後の要素が要らず、入っているかだけを見るなら HashMap で足ります。

値が FastSet に収まるなら、FastSet の方がずっと速く動きます。1 回の操作が 64 分木の 4 段ほどをたどるだけで、各段は 64 bit の語を 1 つ読むだけだからです。

## 中身

64 分木の bitset です。段 h の i 番目の bit は、段 h - 1 の i 番目の語が 0 でないことを表し、段は 1 本の配列に根から順に並べます。前後の要素を探すときは、葉の語を読む前に 1 つ上の段の bit を見て、葉の語が 0 なら葉を読まずに上の段から探します。$n = 10^7$ で葉の段は 1.25 MB あって L2 に収まらないことが多い一方、1 つ上の段は 19.5 KB で L1 に載るので、疎な集合で葉を読みに行く分が省けます。入れるときと除くときは、上の段が変わらなくなったところで止めます。

方式は procon-judge の yosupo-predecessor-problem で比べて選びました ($n = 10^7$、$Q = 10^6$)。同じ方式の fs_td は、EPYC 9V45 の x64-gcc で最大 11 ms ほどで、今の Library の SegmentTree の 19 分の 1、std::set の 180 分の 1 でした。最大は要素がほとんど無い集合で前後を問うケースです。段の並べ方 (段ごとに別の配列、葉とその上の段を組にする形、512 分木) と huge page は、どれも 0 % から 3 % ほどの差か、遅くなるかでした。記録は algo-notes の notes/associative-containers.md にあります。
