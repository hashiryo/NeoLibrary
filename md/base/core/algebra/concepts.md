---
title: 構造の判定
documentation_of: ../../neo/base/core/algebra/concepts.hpp
---

ライブラリが受ける構造 (値の型と static 関数を並べた struct) の判定です。セグメント木、列、重み平衡木などが、渡された struct の種類をコンパイル時に見分けるのに使います。

```cpp
template <class M> concept Monoid;       // using T、static T op(T, T)、static T ti()
template <class M> concept Action;       // using T、using E、static void mp(...)、static void cp(E&, E)
template <class M> concept Commutative;  // using commute
```

`Monoid` はモノイド (`T`、`op`、`ti`)、`Action` は作用 (`T`、`E`、`mp`、`cp`) の決まりを満たすかを見ます。決まりの中身は `neo/base/ds/range/SegmentTree.hpp` と同じです。`Commutative` は `using commute= void` があるかだけを見ます。

判定は `&M::op` のような関数のアドレスで行い、`M::op(a, a)` を呼ぶ式では書きません。名前が 1 つの関数に決まることまで見るためです。関数を呼ぶ式にすると、引数の型が変換で合う型まで受けるようになり、受ける範囲が変わります。

半環からモノイドを取る変換は [monoid](monoid.md)、演算子つきの型を構造に読み替える部品は [PlusTimes](PlusTimes.md) にあります。
