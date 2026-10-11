#pragma once
#include <array>
#include <cassert>
#include <utility>
#include "neo/number_theory/factorize.hpp"
// 素数 p の原始根を 1 つ求める。n = p - 1 を factorize_pairs で分解し、n のどの素因数 q でも g^(n / q) ≠ 1 となる g を作る。
//
// q = 2 の条件は g が平方非剰余であることなので、小さい素数から Jacobi 記号で最初の平方非剰余 g_0 を探す。奇数の素因数
// q_1, ..., q_k については、g_0^(n / q_j) を二分木で一度に求める。h = g_0^(n / (q_1 ... q_k)) を作り、組を半分に分けて、
// 左の組には h を右の組の積で上げた値を、右の組には左の組の積で上げた値を渡す (2 つの累乗は同じ底の表で並べて回す)。
// g_0 がすべての q を満たせば g_0 を返す。満たさない q があれば、q ごとに splitmix64 の擬似乱数 a で a^(n / q) ≠ 1 を
// 確かめて探す (乱数なら q で落ちる確率は 1 / q。小さい素数がどれも 3 乗剰余になるような p でも詰まらない)。
// 候補 g_i が受け持つ q^e の積を A_i として、Π g_i^(n / A_i) を返す (g_i^(n / A_i) の位数は A_i ちょうどで、A_i どうしは
// 互いに素なので、積の位数は n)。
//
// 剰余の掛け算は R = 2^64 の Montgomery。p < 2^62 なら値を 2p 未満のまま持ち、それ以上なら p 未満に保つ。累乗は 4 bit の
// 窓 (指数が 2^16 未満なら 2 進) で、どちらも分岐なし。
namespace primitive_root_internal {
using u32= unsigned;
using u64= unsigned long long;
using u128= unsigned __int128;
template <bool LAZY> struct Mont {
 u64 n, ninv, r1, r2;  // r1 = 2^64 mod n、r2 = 2^128 mod n
 explicit Mont(u64 n_): n(n_) {
  u64 x= (3 * n) ^ 2, y= 1 - n * x;
  x*= 1 + y, y*= y;
  x*= 1 + y, y*= y;
  x*= 1 + y, y*= y;
  ninv= x * (1 + y);
  r1= (u64)(((u128)1 << 64) % n), r2= (u64)((u128)r1 * r1 % n);
 }
 // a b / 2^64 mod n。LAZY なら [0, 2n) (a, b < 2n、n < 2^62)、そうでなければ [0, n) (a, b < n)
 u64 mul(u64 a, u64 b) const {
  const u128 t= (u128)a * b;
  const u64 hi= (u64)(t >> 64), m= (u64)((u128)((u64)t * ninv) * n >> 64);
  if constexpr(LAZY) return hi + n - m;
  else return hi >= m ? hi - m : hi - m + n;
 }
 u64 to(u64 a) const { return mul(a, r2); }
 u64 norm(u64 a) const {
  if constexpr(LAZY) return a >= n ? a - n : a;
  else return a;
 }
 bool is_one(u64 a) const { return norm(a) == norm(r1); }
 u64 pow(u64 b, u64 e) const {
  if(e == 0) return r1;
  if(e < (1u << 16)) {
   int sh= 63 - __builtin_clzll(e);
   u64 r= b;
   for(--sh; sh >= 0; --sh) {
    r= mul(r, r);
    const u64 t= mul(r, b);
    r= (e >> sh & 1) ? t : r;
   }
   return r;
  }
  u64 tab[16];
  tab[0]= r1, tab[1]= b;
  for(int i= 2; i < 16; ++i) tab[i]= mul(tab[i - 1], b);
  int sh= (63 - __builtin_clzll(e)) & ~3;
  u64 r= tab[(e >> sh) & 15];
  for(sh-= 4; sh >= 0; sh-= 4) {
   r= mul(r, r), r= mul(r, r), r= mul(r, r), r= mul(r, r);
   r= mul(r, tab[(e >> sh) & 15]);
  }
  return r;
 }
 // (b^e1, b^e2) を、b の 4 bit の表を共有して 2 本の鎖で並べて求める
 std::pair<u64, u64> pow_pair(u64 b, u64 e1, u64 e2) const {
  const u64 e= e1 | e2;
  if(e < 16) {
   u64 x1= r1, x2= r1;
   for(int sh= 3; sh >= 0; --sh) {
    x1= mul(x1, x1), x2= mul(x2, x2);
    const u64 t1= mul(x1, b), t2= mul(x2, b);
    x1= (e1 >> sh & 1) ? t1 : x1, x2= (e2 >> sh & 1) ? t2 : x2;
   }
   return {x1, x2};
  }
  u64 tab[16];
  tab[0]= r1, tab[1]= b;
  for(int i= 2; i < 16; ++i) tab[i]= mul(tab[i - 1], b);
  int sh= (63 - __builtin_clzll(e)) & ~3;
  u64 x1= tab[(e1 >> sh) & 15], x2= tab[(e2 >> sh) & 15];
  for(sh-= 4; sh >= 0; sh-= 4) {
   for(int i= 0; i < 4; ++i) x1= mul(x1, x1), x2= mul(x2, x2);
   x1= mul(x1, tab[(e1 >> sh) & 15]), x2= mul(x2, tab[(e2 >> sh) & 15]);
  }
  return {x1, x2};
 }
};
// Jacobi 記号 (a / n)。n は奇数。
inline int jacobi(u64 a, u64 n) {
 int s= 1;
 a%= n;
 while(a) {
  const int z= __builtin_ctzll(a);
  a>>= z;
  if((z & 1) && ((n & 7) == 3 || (n & 7) == 5)) s= -s;
  if((a & 3) == 3 && (n & 3) == 3) s= -s;
  std::swap(a, n), a%= n;
 }
 return n == 1 ? s : 0;
}
// out[j] = b^(Π q / q[j]) (j < cnt)
template <class M> void tree(const M& m, u64 b, const u64* q, int cnt, u64* out) {
 if(cnt == 1) return void(out[0]= b);
 const int half= cnt / 2;
 u64 pl= 1, pr= 1;
 for(int j= 0; j < half; ++j) pl*= q[j];
 for(int j= half; j < cnt; ++j) pr*= q[j];
 const auto [xl, xr]= m.pow_pair(b, pr, pl);
 tree(m, xl, q, half, out), tree(m, xr, q + half, cnt - half, out + half);
}
template <bool LAZY> u64 solve(u64 p) {
 const Mont<LAZY> m(p);
 const u64 n= p - 1;
 u64 qs[16], qe[16];  // 奇数の素因数と、その冪 q^e
 int k= 0;
 u64 rad= 1, pow2= 1;
 for(auto [q, e]: factorize_pairs(n)) {
  u64 qq= 1;
  for(int i= 0; i < e; ++i) qq*= q;
  if(q == 2) pow2= qq;
  else qs[k]= q, qe[k]= qq, rad*= q, ++k;
 }
 // q = 2: 小さい素数 (表の先は整数) から最初の平方非剰余
 u64 g0= 2;
 while(jacobi(g0, p) != -1) g0= g0 == 2 ? 3 : g0 + 2;
 if(k == 0) return g0;
 u64 cand[17], A[17];  // 使う候補 (Montgomery 形) と、それが受け持つ冪の積
 cand[0]= m.to(g0);
 u64 out[16];
 tree(m, m.pow(cand[0], n / rad), qs, k, out);
 u32 sat= 0;
 for(int j= 0; j < k; ++j) sat|= (u32)!m.is_one(out[j]) << j;
 if(sat == (1u << k) - 1) return g0;
 A[0]= pow2;  // g_0 は q = 2 も受け持つ
 for(int j= 0; j < k; ++j)
  if(sat >> j & 1) A[0]*= qe[j];
 int nc= 1;
 u64 st= p;  // splitmix64 の状態
 for(int j= 0; j < k; ++j) {
  if(sat >> j & 1) continue;
  for(const u64 e= n / qs[j];;) {
   u64 z= (st+= 0x9e3779b97f4a7c15ull);
   z= (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
   z= (z ^ (z >> 27)) * 0x94d049bb133111ebull;
   z= (z ^ (z >> 31)) % p;
   if(z == 0) continue;
   const u64 am= m.to(z);
   if(!m.is_one(m.pow(am, e))) {
    cand[nc]= am, A[nc]= qe[j], ++nc;
    break;
   }
  }
 }
 u64 res= m.r1;
 for(int i= 0; i < nc; ++i) res= m.mul(res, m.pow(cand[i], n / A[i]));
 return m.norm(m.mul(res, 1));
}
}
// 素数 p を法とする原始根を 1 つ返す (p = 2 なら 1)。2 ≤ p < 2^64 で、p が素数であることは確かめない。
inline unsigned long long primitive_root(unsigned long long p) {
 assert(p >= 2);
 if(p == 2) return 1;
 return p < (1ull << 62) ? primitive_root_internal::solve<true>(p) : primitive_root_internal::solve<false>(p);
}
