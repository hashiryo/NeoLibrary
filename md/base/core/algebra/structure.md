---
title: 構造の読み替え
documentation_of: ../../neo/base/core/algebra/structure.hpp
---

ライブラリが構造を受ける場所 (行列のテンプレート引数など) に書かれた型を、構造 (値の型と static 関数を並べた struct) に読み替えます。

```cpp
template <class X, class = void> struct StructureOf { using type= PlusTimes<X>; };
template <class X> struct StructureOf<X, std::void_t<typename X::T>> { using type= X; };
template <class X> using structure_of= typename StructureOf<X>::type;
```

`X` が値の型 `T` を持つなら、`X` を構造そのものとして使います。持たないなら、演算子つきの型 (`int`、`double`、`ZMod` など) とみなし、[PlusTimes](PlusTimes.md) (普通の + と ×) の略記として扱います。値を持つ struct (`T` が自分自身) も構造です。

```cpp
Matrix<long long> a;        // Matrix<PlusTimes<long long>>。要素は long long
Matrix<ZMod<998244353>> b;  // Matrix<PlusTimes<ZMod<998244353>>>。要素は ZMod
Matrix<MinPlus<int>> c;     // MinPlus<int> は構造。要素は int
```

構造は、値の型 `T` と、道具が要る static 関数 (`add`、`mul`、`zero`、`one`、`sub` など) を持ちます。どの関数が要るかは、道具と [concepts](concepts.md) の族で決まります。

読み替えは型だけで決まるので、`T` を持つ値の型を構造の場所に書くと、構造として読まれます。値の型に `T` という名前の別の型を置かないでください。
