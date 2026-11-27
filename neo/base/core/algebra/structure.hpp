#pragma once
#include <type_traits>
#include "neo/base/core/algebra/PlusTimes.hpp"
// ライブラリが構造を受ける場所 (Matrix や SquareMatrix のテンプレート引数) に書かれた型を、構造に読み替える。
// X が値の型 T を持つなら、X は構造 (値の型と static 関数を並べた struct) そのもの。持たないなら演算子つきの型
// (int、double、ZMod など) で、PlusTimes<X> (普通の + と ×) の略記として扱う。値を持つ struct (T が自分自身) も構造。
template <class X, class = void> struct StructureOf { using type= PlusTimes<X>; };
template <class X> struct StructureOf<X, std::void_t<typename X::T>> { using type= X; };
template <class X> using structure_of= typename StructureOf<X>::type;
