---
title: 素数の列挙
documentation_of: ../../neo/number_theory/enumerate_primes.hpp
---

$N$ 以下の素数を列挙します。

```cpp
std::vector<unsigned> enumerate_primes(unsigned N);
std::pair<unsigned, std::vector<unsigned>> enumerate_primes(unsigned N, unsigned A, unsigned B);
```

1 つ目は $N$ 以下の素数を小さい順にすべて返します。2 つ目は $N$ 以下の素数の個数 $\pi(N)$ と、小さい方から数えて $B$ 番目 (0 始まり) から $A$ 個おきの素数を返します。Library Checker の Enumerate Primes と同じ形です。$0 \le N < 2^{32}$、$0 \le B < A$ の範囲で使い、$B < A$ でなければ assert で止めます。

エラトステネスの篩を 30 の車輪のビットで持ち、32 KB の区画ごとに回します。1 バイトが 30 個の数のうち 30 と互いに素な 8 個を表します。163 以下の素数は周期の型を AND して消し、それより大きい $\sqrt{N}$ 以下の素数は倍数を車輪で進めて消します。4000 以下の素数は 210 の車輪で進め、7 の倍数を飛ばします。素数ごとに次に消す周の先頭だけを持ち、区画の中に先頭がある周を丸ごと消して、区画の後ろへのはみ出しを次の区画に持ち越します。$A$ 個おきの素数は、区画ごとに popcount の累積から取り出す語を決め、語の中の位置を pdep (BMI2) で引きます。

$N = 5 \times 10^8$、$A = 27$ の時間は、procon-judge の yosupo-enumerate-primes の x64 の CI のランナーで 22〜45 ms です。Library Checker の最速の提出の写しより 2.0〜2.5 倍速く出ました。手元の Apple M2 では、$5 \times 10^8$ の $A = 27$ が 34 ms、$10^8$ 以下の素数をすべて返すのが 14 ms、$2^{32} - 1$ までを数えるのが 0.35 s でした。メモリは区画と型の 300 KB ほどと、返す配列の分を使います。
