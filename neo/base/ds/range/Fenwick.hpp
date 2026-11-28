#pragma once
#include <vector>
// 1 点の更新と接頭辞の積を O(log n) で行う。M はセグメント木と同じモノイド
// (using T、static T ti()、static T op(T a, T b))。和にも最大にも使える。
// 逆元の要らないモノイド (最大、最小、gcd など) をそのまま渡せる。
template <class M> struct Fenwick {
 using T= typename M::T;
 int n;
 std::vector<T> dat;
 explicit Fenwick(int n= 0): n(n), dat(n, M::ti()) {}
 // i 番目 (0 始まり) の値を、M::op で x と合わせる。
 void add(int i, const T& x) {
  for(++i; i <= n; i+= i & -i) dat[i - 1]= M::op(dat[i - 1], x);
 }
 // [0, r) の積。
 T prod(int r) const {
  T s= M::ti();
  for(; r > 0; r-= r & -r) s= M::op(s, dat[r - 1]);
  return s;
 }
};
