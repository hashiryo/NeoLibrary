---
title: 半環からモノイドを取る
documentation_of: ../../neo/base/core/algebra/monoid.hpp
---

半環 (`T`、`add`、`mul`、`zero`、`one`) から、モノイドを 1 つ取り出します。半環にはモノイドが 2 つ (和と積) 入っているので、どちらを取るかを名前で書きます。

```cpp
namespace monoid {
template <class R> struct Additive;        // 和。ti は zero()、op は add。可換で、環なら inv も持つ
template <class R> struct Multiplicative;  // 積。ti は one()、op は mul
}
```

セグメント木に半環をそのまま渡すと、`add` と `mul` のどちらを `op` と見るか決まらないので、コンパイルエラーになります。区間の max が欲しければ `Additive`、区間の和が欲しければ `Multiplicative` と書きます。

```cpp
struct MaxPlus {
 using T= long long;
 static T zero() { return -(1LL << 60); }
 static T one() { return 0; }
 static T add(T a, T b) { return std::max(a, b); }
 static T mul(T a, T b) { return a + b; }
};
SegmentTree<monoid::Additive<MaxPlus>> a(n);        // 区間の max
SegmentTree<monoid::Multiplicative<MaxPlus>> b(n);  // 区間の和
SegmentTree<PlusTimes<long long>> c(n);             // コンパイルエラー (どちらか決まらない)
```

`Additive` は可換なので `commute` を持ち、区間の反転を持つ木でも逆順の積を持ちません。`inv` は、半環が引き算 (`sub`) を持つときだけ使えます。値は `sub(zero(), a)` で、群を受ける道具 (Fenwick tree など) に渡すためのものです。

[PlusTimes](PlusTimes.md) と組み合わせると、今の `monoid::Sum<X>` に当たるモノイドを `monoid::Additive<PlusTimes<X>>` と書けます。今の `monoid::Sum` は `SegmentTree` の特殊化のキーなので、そのまま残しています。
