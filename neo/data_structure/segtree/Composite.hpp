#pragma once
#include <array>
#include "neo/data_structure/SegmentTree.hpp"
// 1 次関数の合成のモノイド monoid::Composite<Z>。値 {a, b} は x -> a x + b で、op(f, g) は f を先に、g をあとに当てた関数。
// SegmentTree に渡すと二分木の汎用の形のままで、区間の積は覆う節を分岐なしで集めてから 2 本の鎖で畳む (procon-judge の
// yosupo-point-set-range-composite で、この形が書き比べの首位と並んだ)。区間代入と合わせるなら monoid::AssignProd に渡す
// (neo/data_structure/segtree/AssignProd.hpp)。
namespace monoid {
template <class Z> struct Composite {
 using T= std::array<Z, 2>;
 static constexpr T ti() { return {Z(1), Z(0)}; }
 static constexpr T op(const T& f, const T& g) { return {g[0] * f[0], g[0] * f[1] + g[1]}; }
 // f(x)
 static constexpr Z eval(const T& f, const Z& x) { return f[0] * x + f[1]; }
};
}
