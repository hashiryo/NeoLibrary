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

## 中身

64 分木の bitset です。段 h の i 番目の bit は、段 h - 1 の i 番目の語が 0 でないことを表し、段は 1 本の配列に根から順に並べます。前後の要素を探すときは、葉の語を読む前に 1 つ上の段の bit を見て、葉の語が 0 なら葉を読まずに上の段から探します。$n = 10^7$ で葉の段は 1.25 MB あって L2 に収まらないことが多い一方、1 つ上の段は 19.5 KB で L1 に載るので、疎な集合で葉を読みに行く分が省けます。入れるときと除くときは、上の段が変わらなくなったところで止めます。

方式は procon-judge の yosupo-predecessor-problem で比べて選びました ($n = 10^7$、$Q = 10^6$)。同じ方式の fs_td は、EPYC 9V45 の x64-gcc で最大 11 ms ほどで、今の Library の SegmentTree の 19 分の 1、std::set の 180 分の 1 でした。最大は要素がほとんど無い集合で前後を問うケースです。段の並べ方 (段ごとに別の配列、葉とその上の段を組にする形、512 分木) と huge page は、どれも 0 % から 3 % ほどの差か、遅くなるかでした。記録は algo-notes の notes/associative-containers.md にあります。
