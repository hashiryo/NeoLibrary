---
title: 最大公約数
documentation_of: ../../neo/number_theory/gcd.hpp
---

整数の最大公約数を求めます。

```cpp
template <class T> constexpr T gcd(T a, T b);
```

`T` は 64 bit 以下の整数型で、$|a|$ と $|b|$ の最大公約数を返します。`gcd(0, 0)` は 0 です。constexpr なのでコンパイル時にも使えます。

std::gcd と同じ名前です。`using namespace std` の下で呼んでも曖昧にはなりません。2 つの引数を同じ型で渡すとこの関数が選ばれます。型の違う引数を渡すと std::gcd が選ばれます。

Stein の方法 (binary GCD) で計算し、2 つの数の桁数の差が 8 bit を超えるときだけ、先に大きいほうを小さいほうで 1 回割ります。Stein の方法は大きい数と小さい数の組では 1 段に 2 bit ほどしか縮まないので、そこだけ割り算に任せます。

64 bit の一様な乱数の組で 1 回 70 ns ほど、64 bit と 8 bit の組で 17 ns ほどです。値は procon-judge の self-gcd-test-runtime-64 の、EPYC 7763 の x64-gcc の記録です。
