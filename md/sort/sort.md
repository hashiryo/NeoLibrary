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

procon-judge で測った x64-gcc の最大は次のとおりです (2026-10-27 から 11-01)。

| 問題 | 列 | CPU | この関数 | `std::sort` |
| --- | --- | --- | --- | --- |
| self-sort-u32 | 一様な u32 を $10^6$ 個 | EPYC 9V45 | 2.4 ms | 57.5 ms |
| | | Xeon 8573C | 7.5 ms | 76.4 ms |
| self-sort-small-u32 | 長さ 1 から 64 の u32 の配列 65536 本 | EPYC 7763 | 4.3 ms | 38.9 ms |
| | | Xeon 8370C | 3.9 ms | 40.6 ms |
| self-sort-large-u32 | 一様な u32 を $10^7$ 個 | EPYC 7763 | 88 ms | 778 ms |
| | | Xeon 8370C | 97 ms | 802 ms |
| self-sort-u64 | 一様な u64 を $10^6$ 個 | EPYC 7763 | 9.3 ms | 66.6 ms |
| | | Xeon 8370C | 16.4 ms | 69.7 ms |
| self-sort-i64 | 一様な i64 を $10^6$ 個 | EPYC 9V45 | 7.1 ms | 59.6 ms |
| | | Xeon 8573C | 17.8 ms | 72.0 ms |

Intel の Xeon でも、ほとんどの回で、同じ回の写し元の提出の首位と並ぶか、それより速く出ました。self-sort-i64 では、EPYC 9V45 と Xeon 8370C で写し元より 2〜4 ms 遅く出ています (9V45 では値が i64 の全域にわたるケースだけ)。理由はまだ見ていません。$10^7$ 個の列は、AMD の 9V74 と 9V45 では、上の 8 bit で振り分けてから塊ごとに LSD を回す形 (msd_lsd) が、この関数より 1.3 倍ほど速く出ています。この大きさの道はまだ詰めていません。
