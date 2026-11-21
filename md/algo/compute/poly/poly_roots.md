---
title: 多項式の根
documentation_of: ../../neo/algo/compute/poly/poly_roots.hpp
---

素数 $p$ を法とする多項式の根を、重複なしで求めます。

```cpp
template <class Z> std::vector<Z> poly_roots(std::vector<Z> f);
```

`Z` は `ZMod` で、法は素数であることが要ります。根の順は決めていません。$f$ が 0 か定数なら、空を返します。次数 $N$ の $f$ に $O(N^2 \log p)$ の時間がかかり、FFT を使いません。

## 中身

Cantor–Zassenhaus の方法です。$p$ が奇数なら、$h = x^{(p-1)/2} \bmod f$ を作り、$\gcd(f, h - 1)$ と $\gcd(f, h + 1)$ で 0 でない根を平方剰余とそうでないものに分けます。それぞれを、乱択の $a$ による $\gcd(g, (x + a)^{(p-1)/2} - 1)$ で 1 次になるまで割ります。$p = 2$ なら、0 と 1 を代入して調べます。

冪は $f$ を法とする 2 乗の繰り返しです。積は `convolve_naive` と同じ AVX2 の核で、剰余は $\mathrm{rev}(f)^{-1} \bmod x^{N-1}$ を先に作っておく Barrett の方法 (積の下半分 2 回) で求めます。

procon-judge の yosupo-polynomial-root-finding ($N \le 4000$、mod 998244353) で書き比べた naive を、法を定数にしない形で写しました。Xeon 6975P-C の x64-gcc (run 44256846647) で、この形の最大は 162.24 ms でした (写し元の naive は 201.31 ms)。
