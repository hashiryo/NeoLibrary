---
title: 整数の列の並べ替え
documentation_of: ../../neo/sort/sort.hpp
---

列を昇順に並べ替えます。

```cpp
template <class T> void sort(std::vector<T>& a);
```

4 byte と 8 byte の整数の列 (`int`、`unsigned`、`long long`、`unsigned long long` など) は、長さと CPU を見て道を選びます。ほかの型は `std::sort` に任せます。作業用に、列と同じ大きさの配列を使います。

- 4 byte で 64 個以下の列は、AVX2 のソーティングネットワークで並べます。最後の 8 個をふつうに読み書きして端のベクタの lane を回すので、vpmaskmovd を使いません。
- 4 byte で $2^{15}$ 個以下の列と、Intel の CPU での長い列は、作業用の配列へ分割するクイックソートで並べます。8 個ずつ読んで軸と比べ、表を引いた vpermd で寄せてから、作業用の配列の前と後ろに書きます。64 個以下になった塊はネットワークで並べます。
- 4 byte でそれより長い列は、Intel 以外の CPU では LSD の基数ソートで並べます。要素によって違う bit だけを 11 bit 以下の桁に分け、全要素が同じ値になる桁は飛ばします。
- 8 byte の列は、上の桁で 1 回振り分けてから、塊の中を u32 の鍵のネットワークで並べます。鍵は、桁より下の bit の上 26 bit と塊の中の位置を詰めたものです。

Intel の CPU では、振り分けの散らばった書き込みが重くなります。procon-judge の self-sort-u32 (一様な u32 を $10^6$ 個) では、Xeon 8573C と 8370C で LSD がクイックソートの 1.3〜1.5 倍かかりました。AMD の EPYC では逆に、クイックソートが LSD の 1.4〜1.9 倍かかります。作業用の配列は、512 KB 以上なら Linux では huge page を頼みます。

写し元の提出の速さは、self-sort-u32 の x64-gcc で 2.5〜8.7 ms ほどで、`std::sort` の 7〜22 倍でした。8 byte の self-sort-u64 と self-sort-i64 では 5〜21 ms ほどで、`std::sort` の 3〜11 倍です。
