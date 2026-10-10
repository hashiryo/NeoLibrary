---
title: 素因数分解
documentation_of: ../../neo/number_theory/factorize.hpp
---

64 bit の符号なし整数を素因数分解します。

```cpp
constexpr std::vector<unsigned long long> factorize(unsigned long long n);
constexpr std::vector<std::pair<unsigned long long, int>> factorize_pairs(unsigned long long n);
```

`factorize(n)` は、$n$ の素因数を重複込みで昇順に並べて返します。12 なら `{2, 2, 3}` です。`factorize_pairs(n)` は、同じ結果を (素数, 指数) の組にまとめて返します。12 なら `{(2, 2), (3, 1)}` です。どちらも $n \le 1$ なら空で、$0 \le n < 2^{64}$ のすべての `n` を受けます。

97 以下の素数で試し割りしてから、残りの合成数の因数を 1 つずつ見つけて再帰します。素数かどうかは同じ number_theory の is_prime で見ます。平方数は平方根で割り、$2^{48}$ 未満は Pollard の rho 法 (Brent の方法で、2 本の列を同時に進めます)、それより大きければ楕円曲線法 (ECM) で因数を探します。

ECM では、Suyama の形の Montgomery 曲線を使います。点は x 座標だけで持ちます。stage 1 の上限は $B_1 = 150$ です。stage 2 は $B_2 \approx 3400$ までを見ます。$2^{60}$ 未満の数では、2 本の曲線を同じ手順で同時に回します。$10^{18}$ 前後で 2 つの素因数がどちらも $2^{30}$ 前後の数は、rho 法では 1 個あたり 3 万段以上かかりますが、ECM なら曲線 5 本ほどで割れます。1 個あたりは 16〜27 us ほどです (procon-judge の yosupo-factorize、EPYC 9V45 と Xeon 8573C の x64-gcc)。

constexpr なので、コンパイル時にも使えます。CI と同じ定数式の評価の上限で、$10^{18}$ 前後の 2 素数の積と、$2^{32}$ に近い 2 素数の積をコンパイル時に割れることを確かめました。要る曲線の数は数によって違うので、大きい数では上限に当たることがあります。
