#pragma once
#include <cassert>
#include <vector>
// 素数 p を法とする階乗と、その逆元の表。Z は ZMod (法が素数で、0 <= n < p の範囲で使う)。
//   Factorial<Z> f(n);  // 0! から n! までを用意する
//   f.ensure(m);        // 必要になったら m! まで延ばす (延ばした分だけ作る)
//   f.fact(k)、f.finv(k)、f.inv(k)、f.nPr(n, r)、f.nCr(n, r)、f.nHr(n, r)
//
// 表は fac と ifac の 2 本。ifac[n] = fac[n]^{-1} を 1 回の逆元で作り、下向きの掛け算で ifac[0] まで埋める。
// 今の Library の FactorialPrecalculation は iv、fct、fiv の 3 本を持ち、fiv を 1 つずつ inv(k) を掛けて作る。
// 2 本にするとメモリが 3 分の 2 になり、fiv の作り方も掛け算だけになる (n = 10^7 で 150 MB が 80 MB)。
// n >= p ではファクトリオドが 0 になるので、nCr は Lucas の定理が要る。この表は n < p のときだけ正しい。
template <class Z> struct Factorial {
 std::vector<Z> fac_, ifac_;
 int n_= -1;
 [[gnu::always_inline]] void ensure(int n) {
  if(n <= n_) return;
  assert(0 <= n);
  assert((unsigned)n < Z::mod());
  fac_.resize(n + 1), ifac_.resize(n + 1);
  if(n_ < 0) fac_[0]= Z(1), n_= 0;
  for(int i= n_ + 1; i <= n; ++i) fac_[i]= fac_[i - 1] * Z(i);
  Z iv= fac_[n].inv();
  ifac_[n]= iv;
  for(int i= n; i > n_; --i) ifac_[i - 1]= ifac_[i] * Z(i);
  n_= n;
 }
 explicit Factorial(int n= 0) { ensure(n); }
 // 0 <= n
 [[gnu::always_inline]] Z fact(int n) { return ensure(n), fac_[n]; }
 [[gnu::always_inline]] Z finv(int n) { return ensure(n), ifac_[n]; }
 // 1 <= n
 [[gnu::always_inline]] Z inv(int n) { return ensure(n), fac_[n - 1] * ifac_[n]; }
 // nPr = n! / (n - r)!
 [[gnu::always_inline]] Z nPr(int n, int r) {
  return r < 0 || n < r ? Z(0) : (ensure(n), fac_[n] * ifac_[n - r]);
 }
 // nCr。n < p のときだけ正しい。
 [[gnu::always_inline]] Z nCr(int n, int r) {
  return r < 0 || n < r ? Z(0) : (ensure(n), fac_[n] * ifac_[r] * ifac_[n - r]);
 }
 // nHr = nCr(n + r - 1, r) で、[x^r] (1 - x)^{-n}
 [[gnu::always_inline]] Z nHr(int n, int r) {
  return !r ? Z(1) : nCr(n + r - 1, r);
 }
};
