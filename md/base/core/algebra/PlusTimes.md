---
title: 演算子つきの型の読み替え
documentation_of: ../../neo/base/core/algebra/PlusTimes.hpp
---

`int`、`double`、`ZMod` のように、`+` と `*` を演算子で書く型を、ライブラリが受ける構造 (値の型と static 関数を並べた struct) に読み替えます。`PlusTimes<X>` は `X` の (+, ×) です。

```cpp
template <class X> struct PlusTimes {
 using T= X;
 static T zero();         // X()
 static T one();          // X::one() があればそれ、無ければ X(1)
 static T add(T a, T b);  // a + b
 static T mul(T a, T b);  // a * b
 static T sub(T a, T b);  // a - b。引き算がある型だけ
};
```

構造だけを受け取る道具に、演算子つきの型をそのまま渡すときに使います。たとえば `Matrix<PlusTimes<long long>>` は、要素が `long long` で、和と積が普通の `+` と `*` の行列です。

`one()` は `X::one()` を先に見ます。`MinPlus` のように単位元を自分で持つ型ではそちらの値に、`int` や `ZMod` のように持たない型では `X(1)` になります。

引き算を持たない型でも、`add`、`mul`、`zero`、`one` だけを使う道具には渡せます。`sub` を使うのは、環を受ける道具と、`sub` から `inv` を作る [monoid](monoid.md) の `Additive` だけです。
