#pragma once
// 演算子つきの型を、ライブラリが受ける構造 (値の型と static 関数を並べた struct) に読み替える。
// PlusTimes<X> は X の (+, ×) で、zero() は X()、one() は X::one() があればそれ、無ければ X(1)。
// sub は引き算がある型だけ。Matrix などの半環を受ける道具に、int、double、ZMod をそのまま渡せる。
template <class X> struct PlusTimes {
 using T= X;
 static T zero() { return T(); }
 static T one() {
  if constexpr(requires { X::one(); }) return X::one();
  else return T(1);
 }
 static T add(T a, T b) { return a + b; }
 static T mul(T a, T b) { return a * b; }
 static T sub(T a, T b)
  requires requires { a - b; }
 {
  return a - b;
 }
};
