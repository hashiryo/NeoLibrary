---
title: 拡張 gcd と逆元
documentation_of: ../../neo/number_theory/inv_gcd.hpp
---

拡張 Euclid の互除法の結果と、法 $m$ での逆元を求めます。

```cpp
template <class T> constexpr std::pair<T, T> inv_gcd(T a, T b);
template <class T> constexpr T inv_mod(T a, T m);
```

`inv_gcd(a, b)` は、$g = \gcd(a, b)$ と、$ax \equiv g \pmod b$ かつ $0 \le x < b / g$ を満たす $x$ の組 $(g, x)$ を返します。$b \ge 1$ が前提です。約束は ACL の `internal::inv_gcd` と同じです。$ax + by = g$ の $y$ が要るときは $y = (g - ax) / b$ で求まります。

`inv_mod(a, m)` は、$ax \equiv 1 \pmod m$ かつ $0 \le x < m$ を満たす $x$ を返します。$m \ge 1$ が前提で、$\gcd(a, m) = 1$ でなければ assert で止まります。

`T` は 64 bit 以下の整数型です。符号付きなら $a$ は負でもかまいません ($b$ で割った余りに直してから計算します)。符号なしの 64 bit なら、$2^{64}$ 未満のすべての値を受けます。constexpr なのでコンパイル時にも使えます。

割り算を使わない Stein の方法の拡張版で計算します。値を $2^c$ で割る代わりにもう一方の係数を $2^c$ 倍し (Kaliski の almost inverse と同じ形)、最後に $2^{-k}$ を Montgomery の還元で掛けて戻します。係数はずっと $b$ 以下に収まるので、128 bit の値を持たずに済みます。2 つの数の桁数の差が 8 bit を超えるときは、先に割り算を 1 回します。

64 bit の一様な乱数の組で 1 回 60 ns ほどで、教科書どおりの拡張 Euclid の互除法より 1.8 倍ほど速いです。値は procon-judge の self-gcd-test-inv-gcd-64 の、EPYC 9V45 の x64-gcc の記録です。
