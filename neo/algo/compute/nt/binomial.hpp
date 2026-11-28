#pragma once
#include <cassert>
#include <cstddef>
#include <utility>
#include <vector>
#include "neo/algo/compute/nt/factorize.hpp"
#include "neo/base/core/nt/inv_gcd.hpp"
// 法 m <= 1e6 (合成数可) の二項係数を、前計算 O(m) ほどで 1 問 O(log_m n) で返す。
//   BinomialCoefficient binom(m);
//   binom.nCr(n, r);  // 0 <= n <= 1e18、0 <= r
//   binom.nHr(n, r);
//
// m = Π p^e と分け、素数べき p^e ごとに Granville の方法で C(n, r) mod p^e を求めて CRT で束ねる。
// p を取り除いた階乗の表 facts と、p 進の桁和の表 ds を持ち、n を pp = p^(e+k) で割りながら桁を見る。
// pp は pre_size 以下の最大の p の冪 (p^e の倍数) で、法が素数 1 つのときは pre_size = 2^20、それ以外は 2^14。
// 今の Library の BinomialCoefficient と同じ算法で、因数分解に neo の factorize_pairs、逆元に inv_mod を使う。
class BinomialCoefficient {
 struct PrimePower {
  int p, e, pe, pp, q, mask= 0;
  std::vector<unsigned> ppows, facts, ds;
  void init(int p_, int e_, std::size_t pre_size) {
   p= p_, e= e_, ppows.assign(e + 1, 1);
   for(int i= 1; i <= e; ++i) ppows[i]= ppows[i - 1] * p;
   for(pp= pe= ppows[e]; std::size_t(pp) * p <= pre_size;) pp*= p;
   q= pp / pe * p, facts.assign(pp, 1);
   // p の倍数を除いた 1..pp-1 の積を、p 進の桁ごとの塊で作る。
   for(int qq= 1, l= pp / p; qq < q; qq*= p, l/= p)
    for(int i= 0; i < l; ++i)
     for(int j= i * p + 1; j < i * p + p; ++j) facts[j * qq]= j;
   for(int i= 1; i < pp; ++i) facts[i]= (unsigned long long)facts[i - 1] * facts[i] % pe;
   mask= facts[pp - 1] == unsigned(pe - 1);  // p = 2 の符号
   ds.assign(q, 0);
   for(int i= 0; i < pp / pe; ++i)
    for(int j= 0, s= ds[i]; j < p; ++j) ds[i * p + j]= s + j;
  }
  // C(n, r) mod p^e。0 <= r <= n。
  int operator()(long long n, long long r) const {
   long long s= n - r;
   int num= 1, den= 1, x= 0, sg= 0;
   if(e > 1)
    for(long long n1, r1, s1; n > 0; n= n1, r= r1, s= s1) {
     n1= n / pp, r1= r / pp, s1= s / pp;
     num= (long long)num * facts[n - n1 * pp] % pp;
     den= (long long)den * facts[r - r1 * pp] % pp * facts[s - s1 * pp] % pp;
     sg+= n1 - r1 - s1, n1= n / q, r1= r / q, s1= s / q;
     x+= ds[r - r1 * q] + ds[s - s1 * q] - ds[n - n1 * q];
    }
   else
    for(long long n1, r1, s1; n > 0; n= n1, r= r1, s= s1) {
     n1= n / pp, r1= r / pp, s1= s / pp;
     const int nr= n - n1 * pp, rr= r - r1 * pp, sr= s - s1 * pp;
     num= (long long)num * facts[nr] % pp;
     den= (long long)den * facts[rr] % pp * facts[sr] % pp;
     sg+= n1 - r1 - s1, x+= ds[rr] + ds[sr] - ds[nr];
    }
   if(x >= e * (p - 1)) return 0;
   if(p > 2) x/= p - 1;
   int ret= (long long)num * inv_mod(den, pe) % pe * ppows[x] % pe;
   return (sg & mask) && ret > 0 ? pe - ret : ret;
  }
 };
 int mod_;
 std::vector<PrimePower> pp_;
 std::vector<int> iprod_;
public:
 explicit BinomialCoefficient(int m, std::size_t pre_size= 1 << 14): mod_(m) {
  assert(1 <= m && m <= 1000000);
  auto f= factorize_pairs(m);
  if(f.size() == 1) pre_size= 1 << 20;
  int prod= 1;
  for(auto [p, e]: f) {
   pp_.emplace_back(), pp_.back().init((int)p, (int)e, pre_size);
   iprod_.push_back(inv_mod(prod, pp_.back().pe));
   prod*= pp_.back().pe;
  }
 }
 // C(n, r) mod m
 int nCr(long long n, long long r) const {
  if(r < 0 || n < r) return 0;
  if(r == 0) return mod_ > 1;
  int ret= 0, prod= 1;
  for(std::size_t i= 0; i < pp_.size(); ++i) {
   const int d= pp_[i].pe;
   ret+= (long long)(pp_[i](n, r) + d - ret % d) * iprod_[i] % d * prod;
   prod*= d;
  }
  return ret;
 }
 // [x^r] (1 - x)^{-n}
 int nHr(long long n, long long r) const { return !r ? 1 : nCr(n + r - 1, r); }
};
