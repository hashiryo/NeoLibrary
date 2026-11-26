#pragma once
// 半環からモノイドを取り出す忘却関手。半環にはモノイドが 2 つ入っているので、どちらを取るかを名前で書く。
//   Additive:       和 (zero、add)。可換なので commute を持ち、引き算がある半環 (環) なら inv も持つ。
//   Multiplicative: 積 (one、mul)。
// 判定は neo/base/core/algebra/concepts.hpp の Monoid、Commutative がそのまま使える。
namespace monoid {
template <class R> struct Additive {
 using T= typename R::T;
 using commute= void;
 static T ti() { return R::zero(); }
 static T op(T a, T b) { return R::add(a, b); }
 static T inv(T a)
  requires requires { R::sub(a, a); }
 {
  return R::sub(R::zero(), a);
 }
};
template <class R> struct Multiplicative {
 using T= typename R::T;
 static T ti() { return R::one(); }
 static T op(T a, T b) { return R::mul(a, b); }
};
}  // namespace monoid
