---
title: 傾きが単調な Convex Hull Trick
documentation_of: ../../neo/algo/optimize/line/MonotoneCHT.hpp
---

傾きが単調な順に直線 y = a x + b を足しながら、点 x での最小 (または最大) を引きます。dp[i] = min_j (dp[j] + (h_i - h_j)^2) のように、足す直線の傾きが決まった向きに並ぶ DP のための型です。

```cpp
template <class T= long long, bool MAX= false> class MonotoneCHT;
```

足す順には約束があります。最小なら傾き a が増えない順、最大なら減らない順です (同じ傾きが続いてもかまいません)。逆の順のときは、x の符号を反転して、傾きを -a にして使います。傾きの順を問わないときや、線分を足すときは、LiChaoTree を使います。

`T` は 64 bit までの符号付き整数で、a x + b が `T` に収まる必要があります。直線を 1 本も足さないうちの問いは未定義です。

## 使い方

```cpp
MonotoneCHT<long long> c;         // 最小。傾きは足す順に増えない
c.add_line(-2 * h[0], h[0] * h[0]);
for(int i= 1; i < n; ++i) {
 dp= c.query_monotone(h[i]) + h[i] * h[i] + C;  // h は増える
 c.add_line(-2 * h[i], dp + h[i] * h[i]);
}
```

| 口 | 中身 |
| --- | --- |
| `add_line(a, b)` | 直線を足す。償却 O(1) |
| `query(x)` | 任意の x での最小 (最大)。二分探索で O(log) |
| `query_monotone(x)` | 前の問いより小さくない x での最小 (最大)。償却 O(1) |
| `reserve(n)` | 直線 n 本分の領域を取っておく |

`query` と `query_monotone` は混ぜて使えます。

## 中身

包絡に出る直線を配列に持ち、足すたびに、末尾の要らなくなった直線を除きます。要るかどうかは交点の位置の比較で、切片の差と傾きの差の積を `__int128` で比べます。直線は交点が式で出るので、判定が O(1) で済みます。`query_monotone` は、包絡の先頭から、次の直線のほうが小さくなる間だけ進みます。最大は、傾きと切片の符号を反転して最小で持ちます。

procon-judge で書き比べた核の、x64-gcc (Xeon 6975P-C) の計測区間の最大です (run 45323501621)。

| 問題 | query_monotone | query | LiChaoTreePoints | 今の Library の ConvexHullTrick |
| --- | --- | --- | --- | --- |
| atcoder-dp-z (N = 2 × 10^5) | 5.8 ms | 8.9 ms | 20.1 ms | 58.2 ms |
| yuki-409 (N ≤ 3 × 10^5) | 7.5 ms | 13.0 ms | 34.4 ms | 55.4 ms |

記録は algo-notes の notes/dp-optimization.md にあります。
