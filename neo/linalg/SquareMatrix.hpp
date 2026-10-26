#pragma once
#include <array>
#include <bit>
#include <cstdint>
#include <type_traits>
#include "neo/algebra/ZMod.hpp"
// 大きさを型で決める N × N の行列。std::array に行優先で持つ (a[i * N + j] が i 行 j 列) ので、セグメント木の要素にしても
// ヒープを使わず、N が定数なので内側のループが展開される。要素の型 T は Matrix と同じく半環でよく、0 と 1 は、T に静的メンバ関数
// zero() と one() があればそれを、なければ T(0) と T(1) を使う。引き算とスカラー倍は、T がそれを持つときだけ使える。
// T が ZMod なら、積の成分ごとに N 個の積を u64 の和にためて 1 回だけ還元する。値は表現のまま (2M 未満) 掛け、和が u64 を
// 溢れうるとき (法が 2^30 に近いと N > 4) だけ M 未満に直してから掛ける。さらに N > 16 なら、17 個目からは 8 段ごとに和の
// 2^32 の位を 2^32 mod M に掛けて畳む。procon-judge の atcoder-abc246-h の zmod_fixed_acc を写したもので、積と和のたびに
// 還元する形の 1.6 倍 (x64-gcc) と 1.9 倍 (x64-clang) の速さだった。
namespace square_matrix_internal {
template <class T> constexpr T zero() {
 if constexpr(requires { T::zero(); }) return T::zero();
 else return T(0);
}
template <class T> constexpr T one() {
 if constexpr(requires { T::one(); }) return T::one();
 else return T(1);
}
template <class T> struct mod_param: std::integral_constant<int, -1> {};  // ZMod でなければ -1
template <unsigned M> struct mod_param<ZMod<M>>: std::integral_constant<int, int(M)> {};
// ZMod の値を表現のまま (2M 未満) 読む。
template <class Z> constexpr unsigned raw(Z x) { return std::bit_cast<unsigned>(x); }
// 表現のまま (2M 未満) の値どうしの積を、いくつまで u64 に足せるか。ZMod<0> は法が 2^30 未満なので 4。
template <class Z> constexpr unsigned long long raw_terms() {
 constexpr unsigned long long M= mod_param<Z>::value;
 if constexpr(M == 0) return 4;
 else return ~0ull / ((2 * M - 1) * (2 * M - 1));
}
// sum_k x[k] y[k s] を、還元を 1 回で。L は項の数、x と y は Z の配列で、y は s 個おきに読む。
template <class Z, int L, int S> constexpr Z dot(const Z* x, const Z* y) {
 using u64= unsigned long long;
 if constexpr(L <= raw_terms<Z>()) {
  u64 t= 0;
  for(int k= 0; k < L; ++k) t+= u64(raw(x[k])) * raw(y[k * S]);
  return Z(t);
 } else if constexpr(L <= 16) {
  u64 t= 0;
  for(int k= 0; k < L; ++k) t+= u64(x[k].val()) * y[k * S].val();
  return Z(t);
 } else {
  const u64 R= (u64(1) << 32) % Z::mod();
  u64 t= 0;
  for(int k= 0; k < L; ++k) {
   t+= u64(x[k].val()) * y[k * S].val();
   if(k >= 15 && (k & 7) == 7 && k + 1 < L) t= (t >> 32) * R + unsigned(t);
  }
  return Z(t);
 }
}
}
template <class T, int N> struct SquareMatrix {
 static_assert(N >= 1);
 std::array<T, N * N> a;
 constexpr SquareMatrix() { a.fill(square_matrix_internal::zero<T>()); }
 constexpr SquareMatrix(const std::array<T, N * N>& a): a(a) {}
 static constexpr SquareMatrix identity() {
  SquareMatrix r;
  for(int i= 0; i < N; ++i) r.a[i * N + i]= square_matrix_internal::one<T>();
  return r;
 }
 constexpr T* operator[](int i) { return a.data() + i * N; }
 constexpr const T* operator[](int i) const { return a.data() + i * N; }
 constexpr SquareMatrix& operator+=(const SquareMatrix& r) {
  for(int i= 0; i < N * N; ++i) a[i]+= r.a[i];
  return *this;
 }
 constexpr SquareMatrix& operator-=(const SquareMatrix& r) {
  for(int i= 0; i < N * N; ++i) a[i]-= r.a[i];
  return *this;
 }
 constexpr SquareMatrix& operator*=(const T& c) {
  for(auto& x: a) x*= c;
  return *this;
 }
 friend constexpr SquareMatrix operator+(SquareMatrix l, const SquareMatrix& r) { return l+= r; }
 friend constexpr SquareMatrix operator-(SquareMatrix l, const SquareMatrix& r) { return l-= r; }
 friend constexpr SquareMatrix operator*(SquareMatrix l, const T& c) { return l*= c; }
 friend constexpr SquareMatrix operator*(const T& c, SquareMatrix r) {
  for(auto& x: r.a) x= c * x;
  return r;
 }
 // 行列の積。ZMod なら成分ごとの内積を 1 回で還元し、ほかは i-k-j の順に足す。
 friend constexpr SquareMatrix operator*(const SquareMatrix& l, const SquareMatrix& r) {
  SquareMatrix c;
  if constexpr(square_matrix_internal::mod_param<T>::value >= 0) {
   for(int i= 0; i < N; ++i)
    for(int j= 0; j < N; ++j) c.a[i * N + j]= square_matrix_internal::dot<T, N, N>(&l.a[i * N], &r.a[j]);
  } else {
   for(int i= 0; i < N; ++i)
    for(int k= 0; k < N; ++k) {
     const T x= l.a[i * N + k];
     for(int j= 0; j < N; ++j) c.a[i * N + j]+= x * r.a[k * N + j];
    }
  }
  return c;
 }
 constexpr SquareMatrix& operator*=(const SquareMatrix& r) { return *this= *this * r; }
 // 列ベクトルとの積 A v。
 friend constexpr std::array<T, N> operator*(const SquareMatrix& l, const std::array<T, N>& v) {
  std::array<T, N> r;
  for(int i= 0; i < N; ++i) {
   if constexpr(square_matrix_internal::mod_param<T>::value >= 0) r[i]= square_matrix_internal::dot<T, N, 1>(&l.a[i * N], v.data());
   else {
    T s= square_matrix_internal::zero<T>();
    for(int j= 0; j < N; ++j) s+= l.a[i * N + j] * v[j];
    r[i]= s;
   }
  }
  return r;
 }
 // A^k。二分累乗。
 constexpr SquareMatrix pow(uint64_t k) const {
  SquareMatrix r= identity(), b= *this;
  for(; k; k>>= 1) {
   if(k & 1) r= r * b;
   if(k > 1) b= b * b;
  }
  return r;
 }
 friend constexpr bool operator==(const SquareMatrix& l, const SquareMatrix& r) { return l.a == r.a; }
 friend constexpr bool operator!=(const SquareMatrix& l, const SquareMatrix& r) { return !(l == r); }
};
