#pragma once
#include <cassert>
// Σ_{i=0}^{n-1} ⌊(a i + b) / m⌋ を AtCoder Library の floor_sum と同じ約束で求める (0 ≤ n < 2^32、1 ≤ m < 2^32、a と b は
// 任意。和は mod 2^64)。負の a、b と、m 以上の a、b は先に m の剰余に直し、残りは互除法で
//   a n + b < m なら終わり。Q = ⌊(a n + b) / m⌋、r = (a n + b) mod m として、n (n-1) / 2 ⌊a/m⌋ + n ⌊b/m⌋ を足し、
//   (n, m, a, b) を (Q - ⌊a/m⌋ n - ⌊b/m⌋, a mod m, m, r) にして繰り返す。
// Q は a と b を割る前の a n + b から求めるので、1 段の割り算 3 回は互いに待たない (AtCoder Library の形は a mod m と
// b mod m を作ってから割るので、割り算 2 回ぶんが 1 段の依存になる)。段の始めの n、m、a、b はどれも 32 bit に収まり、
// a n + b < 2^64。procon-judge の yosupo-sum-of-floor-of-linear の 3 回目の fs3_v0 を写したもの。
namespace floor_sum_internal {
using u32= unsigned;
using u64= unsigned long long;
// a, b < m
inline u64 floor_sum_lt(u32 n, u32 m, u32 a, u32 b) {
 u64 ans= 0;
 while(true) {
  const u64 y= (u64)a * n + b;
  if(y < m) break;
  const u64 Q= y / m;
  const u32 r= (u32)(y % m), qa= a / m, qb= b / m;
  ans+= (u64)n * (n - 1) / 2 * qa + (u64)n * qb;
  const u32 am= a - qa * m;
  const u64 nn= Q - (u64)qa * n - qb;
  if(nn == 0) break;
  n= (u32)nn, b= r, a= m, m= am;
 }
 return ans;
}
}
inline long long floor_sum(long long n, long long m, long long a, long long b) {
 using namespace floor_sum_internal;
 assert(0 <= n && n < (1ll << 32));
 assert(1 <= m && m < (1ll << 32));
 u64 ans= 0;
 const u64 N= (u64)n, M= (u64)m;
 // a、b を [0, m) に直す (負なら下へ、m 以上なら上へずらした分を足し引きする)
 u64 A, B;
 if(a < 0) {
  A= (u64)(a % m + m) % M;
  ans-= N * (N - 1) / 2 * ((A - (u64)a) / M);
 } else {
  A= (u64)a % M;
  ans+= N * (N - 1) / 2 * ((u64)a / M);
 }
 if(b < 0) {
  B= (u64)(b % m + m) % M;
  ans-= N * ((B - (u64)b) / M);
 } else {
  B= (u64)b % M;
  ans+= N * ((u64)b / M);
 }
 return (long long)(ans + floor_sum_lt((u32)N, (u32)M, (u32)A, (u32)B));
}
