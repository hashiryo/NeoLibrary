#pragma once
#include <cassert>
#include <utility>
#include <vector>
#include "neo/algo/compute/nt/factorize.hpp"
#include "neo/base/core/nt/is_prime.hpp"
// 素数 p (2 <= p < 2^64) を法とする乗法の位数を返す。
//   OrderFp ord(p);
//   ord(x);              // x^k ≡ 1 (mod p) となる最小の k >= 1。x ≡ 0 なら 0
//   ord.primitive_root();// 位数 p-1 の元を 1 つ (小さい方から探す)
//
// p - 1 を factorize_pairs で分解し、位数の候補 p - 1 から素因数 q を 1 つずつ落とせるか試す。
// 累乗は Montgomery (p は奇数)。今の Library の OrderFp は constexpr で、位数を素因数の積の
// 二分木でまとめて求めるが、ここでは factorize を使う実行時の形にした (constexpr は無い)。
namespace order_fp_internal {
using u64= unsigned long long;
using u128= unsigned __int128;
constexpr u64 inv64(u64 n) {
 u64 x= (3 * n) ^ 2, y= 1 - n * x;
 x*= 1 + y, y*= y;
 x*= 1 + y, y*= y;
 x*= 1 + y, y*= y;
 return x * (1 + y);
}
// 奇数の m を法とする Montgomery。値は m 未満で持ち、mul の結果も m 未満に直す。m < 2^64 で使える。
struct Mont {
 u64 m, mi, r1, r2;  // mi = m^{-1} mod 2^64、r1 = 2^64 mod m、r2 = 2^128 mod m
 static constexpr u64 rmod(u64 m) {
  u64 r= ~0ull % m + 1;
  return r == m ? 0 : r;
 }
 explicit Mont(u64 m_): m(m_), mi(inv64(m_)), r1(rmod(m_)), r2(u128(rmod(m_)) * rmod(m_) % m_) {}
 [[gnu::always_inline]] u64 mul(u64 a, u64 b) const {
  const u128 t= (u128)a * b;
  const u64 hi= u64(t >> 64), lo= u64(t);
  const u64 mhi= u64((u128)(lo * mi) * m >> 64);
  const u128 r= (u128)hi + m - mhi;
  return u64(r >= m ? r - m : r);
 }
 // a^e mod m。a は m 未満。
 [[gnu::always_inline]] u64 pow(u64 a, u64 e) const {
  u64 r= r1, b= mul(a, r2);
  for(; e; e>>= 1, b= mul(b, b))
   if(e & 1) r= mul(r, b);
  return mul(r, 1);
 }
};
}
class OrderFp {
 using u64= unsigned long long;
 u64 p_;
 std::vector<std::pair<u64, int>> fac_;
 order_fp_internal::Mont mont_;
 u64 pow_mod(u64 a, u64 e) const { return p_ == 2 ? 1 : mont_.pow(a % p_, e); }
public:
 // p は素数
 explicit OrderFp(u64 p): p_(p), fac_(factorize_pairs(p - 1)), mont_(p) { assert(2 <= p && is_prime(p)); }
 // x の位数。x ≡ 0 なら 0、x ≡ 1 なら 1。
 u64 operator()(u64 x) const {
  if(x%= p_; !x) return 0;
  if(x == 1) return 1;
  u64 ret= p_ - 1;
  for(auto [q, e]: fac_)
   for(int i= 0; i < e && ret % q == 0; ++i) {
    if(pow_mod(x, ret / q) != 1) break;
    ret/= q;
   }
  return ret;
 }
 // 位数 p - 1 の元を 1 つ返す (小さい方から探す)。p = 2 なら 1。
 u64 primitive_root() const {
  if(p_ == 2) return 1;
  for(u64 g= 2;; ++g)
   if((*this)(g) == p_ - 1) return g;
 }
};
