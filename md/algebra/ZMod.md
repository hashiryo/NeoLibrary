---
title: 剰余環 Z/MZ
documentation_of: ../../neo/algebra/ZMod.hpp
---

法 $M$ で割った余りの環 $\mathbb{Z}/M\mathbb{Z}$ の元です。$1 \le M < 2^{30}$ で、$M$ の偶奇も素数かどうかも問いません。`ZMod<0>` は、法を実行時に `set_mod` で決める型になります。

```cpp
template <unsigned M> class ZMod;
static void ZMod<0>::set_mod(unsigned m);
static constexpr unsigned mod();
template <class T> constexpr ZMod(T n);
static constexpr ZMod raw(unsigned x);
constexpr unsigned val() const;
constexpr ZMod pow(unsigned long long e) const;
constexpr ZMod inv() const;
constexpr Fixed fixed() const;
```

四則演算 (`+`、`-`、`*`、`/` と複合代入、単項の `-`)、`==`、`!=`、ストリームの入出力があります。`ZMod(n)` は 128 bit 以下の整数型を受け、負の値も $M$ で割った余りに直します。`raw(x)` は、$2M$ 未満の $x$ を、還元しないでそのまま値とします。`val()` は $M$ 未満の値を返します。`inv()` は拡張 Euclid の互除法で逆元を求め、$\gcd(\mathrm{val}(), M) = 1$ でなければ assert で止まります。

`set_mod(m)` は `ZMod<0>` だけで呼べます ($1 \le m < 2^{30}$)。法は `ZMod<0>` の全体で 1 つで、変えると、それまでに作った値と `fixed()` の前計算は使えなくなります。

## 決まった値を掛ける

同じ値を何度も掛けるときは、`fixed()` で前計算してから掛けます。前計算はループの外で 1 回だけ作ります。

```cpp
using mint= ZMod<998244353>;
mint p= ...;
auto P= p.fixed();
for(auto& x : dp) x*= P;  // x * P、P * x も書ける
```

法が奇数の定数なら Plantard の方法、それ以外 (偶数の定数と `ZMod<0>`) は Shoup の方法で前計算します。

## 中身

値は、Montgomery のような表現に直さず、$M$ で割った余りと合同な $2M$ 未満の数で持ちます。掛け算は、前計算 $\lfloor (2^{64}-1)/M \rfloor$ による Barrett の還元です。整数からの変換は、法が定数ならコンパイラに任せた割り算、`ZMod<0>` なら Barrett の還元で行います。

方式は procon-judge の self/modulo-* の問題で比べて選びました。両方の値が変わる掛け算では、Barrett、Plantard、Montgomery の差は数 % 以内でした。`fixed()` の前計算は、法が奇数の定数のとき、ベクトル化しないループの latency と throughput の両方で、前計算をしない掛け算より速くなります。

C++20 以上が要ります。
