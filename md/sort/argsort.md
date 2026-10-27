---
title: 添字の並べ替え (argsort)
documentation_of: ../../neo/sort/argsort.hpp
---

列 `a` について、`a[p[0]] <= a[p[1]] <= ...` となる添字の列 `p` を返します。等しい値は添字の小さい順に並べます (安定な並べ方)。

```cpp
template <class I = int, class T> std::vector<I> argsort(const std::vector<T>& a);
```

`I` は返す添字の型です。4 byte の整数の列は下の道で並べ、ほかの型は添字の列を `std::stable_sort` で並べます。

値と添字を 64 bit に詰めて、値の上の桁で振り分けます。塊の中は、桁より下の bit と塊の中の位置を u32 に詰めた鍵を、AVX2 のソーティングネットワークで並べます。振り分けは元の順を保つので、塊の中の位置の順は添字の順と同じです。

振り分け方は CPU で変えます。Intel 以外の CPU では、桁を塊の平均が 16 個以下になる幅 (16 bit まで) にして、振り分けを 1 回で済ませます。Intel の CPU では書き先が 65536 か所に散らばる振り分けが重いので、上の 8 bit で振り分けてから、1 段目の塊を 1 つずつ、L1 に載る小さい配列へもう 1 回振り分けます。procon-judge の self-sort-argsort-u32 (一様な u32 を $10^6$ 個) で比べると、Intel の Xeon では 2 段の形が 1.3〜1.45 倍速く出ました。AMD の EPYC では 1 回の形が 0〜4 割速く出ます。

self-sort-argsort-u32 では、EPYC 9V74 の x64-gcc で最大 7.7 ms ほどでした。添字の列を `std::stable_sort` で並べる形 (91.5 ms) の 12 倍です。Intel の Xeon での数は、まだ写し元の提出のもの (9〜15 ms ほど) しかありません。作業用に $8N$ byte の配列を使います。
