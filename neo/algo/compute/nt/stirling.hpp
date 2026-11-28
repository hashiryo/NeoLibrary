#pragma once
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>
#include "neo/base/core/nt/is_prime.hpp"
// 小さい素数 p (< 2^15) を法とする第 1 種・第 2 種のスターリング数を、n が大きい (<= 1e18) 範囲で返す。
//   StirlingNumber sn(p, first, second);
//   sn.nCk(n, k);  // C(n, k) mod p
//   sn.S1(n, k);   // [n k] mod p
//   sn.S2(n, k);   // {n k} mod p
//
// p 未満の n については、C(n, k) の三角形と、スターリング数の三角形を p で埋めておく。n >= p は Lucas の定理と
// 同じ形で、n と k を p 進の桁に分けて求める。今の Library の StirlingNumber と同じ算法で、三角形は
// vector<vector<uint16_t>> ではなく 1 本の配列 (行 i は i(i+1)/2 から i+1 個) に詰める。
class StirlingNumber {
 using u16= unsigned short;
 int p_;
 std::vector<u16> c_, s1_, s2_;
 std::vector<int> off_;
 u16 at(const std::vector<u16>& t, int i, int j) const { return t[off_[i] + j]; }
 void build_c() {
  c_.assign(off_[p_], 0);
  for(int i= 0; i < p_; ++i) {
   c_[off_[i]]= 1, c_[off_[i] + i]= 1;
   for(int j= 1; j < i; ++j) {
    const unsigned v= (unsigned)at(c_, i - 1, j - 1) + at(c_, i - 1, j);
    c_[off_[i] + j]= v >= unsigned(p_) ? v - p_ : v;
   }
  }
 }
 void build_s1() {
  s1_.assign(off_[p_], 0);
  for(int i= 0; i < p_; ++i) {
   s1_[off_[i] + i]= 1;
   for(int j= 1; j < i; ++j) {
    const unsigned v= (unsigned)(p_ - i + 1) * at(s1_, i - 1, j) + at(s1_, i - 1, j - 1);
    s1_[off_[i] + j]= v % p_;
   }
  }
 }
 void build_s2() {
  s2_.assign(off_[p_], 0);
  for(int i= 0; i < p_; ++i) {
   s2_[off_[i] + i]= 1;
   for(int j= 1; j < i; ++j) {
    const unsigned v= (unsigned)j * at(s2_, i - 1, j) + at(s2_, i - 1, j - 1);
    s2_[off_[i] + j]= v % p_;
   }
  }
 }
public:
 // p は素数、p < 2^15。first / second で作る三角形を選ぶ。
 explicit StirlingNumber(unsigned p, bool first= true, bool second= true): p_(p) {
  assert(2 <= p && p < (1u << 15) && is_prime(p));
  off_.resize(p_ + 1), off_[0]= 0;
  for(int i= 1; i <= p_; ++i) off_[i]= off_[i - 1] + i;
  build_c();
  if(first) build_s1();
  if(second) build_s2();
 }
 // C(n, k) mod p (Lucas)
 int nCk(unsigned long long n, unsigned long long k) const {
  if(k > n) return 0;
  int ret= 1;
  for(k= std::min(k, n - k); k; ret= (long long)ret * at(c_, int(n % p_), int(k % p_)) % p_, n/= p_, k/= p_)
   if(int(n % p_) < int(k % p_)) return 0;
  return ret;
 }
 // 第 1 種 [n k] mod p
 int S1(unsigned long long n, unsigned long long k) const {
  if(k > n) return 0;
  const unsigned long long i= n / p_;
  if(i > k) return 0;
  long long a= (k - i) / (p_ - 1);
  const u16 j= n % p_;
  u16 b= (k - i) % (p_ - 1);
  if(!b && j) b+= p_ - 1, --a;
  if(a < 0 || i < (unsigned long long)a || b > j) return 0;
  const int v= (long long)nCk(i, a) * at(s1_, j, b) % p_;
  return v && ((i + a) & 1) ? p_ - v : v;
 }
 // 第 2 種 {n k} mod p
 int S2(unsigned long long n, unsigned long long k) const {
  if(k > n) return 0;
  if(!n) return 1;
  const unsigned long long i= k / p_;
  if(n <= i) return 0;
  const unsigned long long a= (n - i - 1) / (p_ - 1);
  const u16 j= k % p_, b= (n - i) - a * (p_ - 1);
  if(j > b) return 0;
  return b == p_ - 1 && !j ? (i ? nCk(a, i - 1) : 0) : (long long)nCk(a, i) * at(s2_, b, j) % p_;
 }
};
