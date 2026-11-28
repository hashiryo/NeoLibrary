#pragma once
#include <type_traits>
// ライブラリが受ける構造 (値の型と static 関数を並べた struct) の判定。名前は族ごとに決める。
//   Monoid:      using T と op(T, T)、ti() (演算が 1 つの族)
//   Action:      using T、using E、mp(T&, E)、cp(E&, E) (遅延セグメント木の作用)
//   ActionLen:   mp(T&, E, int) の形 (作用が区間の長さを受ける。Euler Tour Tree と部分木の遅延用)
//   Commutative: using commute (可換であることの印)
// 受ける範囲は、今の SegmentTree などの std::void_t と decltype(&M::op) の特殊化と同じにする。関数を呼ぶ式
// (requires(typename M::T a) { M::op(a, a); }) では書かない。&M::op で、名前が 1 つの関数に決まることまで見る。
template <class M>
concept Monoid= requires {
 typename M::T;
 &M::op;
 &M::ti;
};
template <class M>
concept Action= requires {
 typename M::T;
 typename M::E;
 &M::mp;
 &M::cp;
};
template <class M, bool = Action<M>>
struct ActionLenV: std::false_type {};
template <class M>
struct ActionLenV<M, true>: std::integral_constant<bool, requires(typename M::T& v, const typename M::E& e, int len) { M::mp(v, e, len); }> {};
template <class M>
concept ActionLen= ActionLenV<M>::value;
template <class M>
concept Commutative= requires { typename M::commute; };
