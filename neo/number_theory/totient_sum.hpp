#pragma once
#include "neo/number_theory/mertens.hpp"
// Φ(N) = Σ_{i≤N} φ(i) = (1 + Σ_{d≤N} μ(d) ⌊N/d⌋^2) / 2 を、Mertens 関数の床の値から求める (Θ(N^{2/3}))。R = ⌊√N⌋ として、
// d ≤ R は μ(d) を直に足し、d > R は q = ⌊N/d⌋ ごとにまとめて q^2 (M(⌊N/q⌋) - M(⌊N/(q+1)⌋)) を足す。和は __int128 で
// 丸めない。⌊N/d⌋ は (N + 1/2) (1/d) の切り捨て (N < 2^51 なら割り算と一致する)。
// Φ(N)。N < 2^50。Θ(N^{2/3}) 時間、Θ(N^{1/2}) 空間
inline unsigned __int128 totient_sum(unsigned long long N) {
 using namespace mertens_internal;
 assert(N < (1ull << 50));
 if(N == 0) return 0;
 const Table T= floor_table(N);
 const u64 R= T.R;
 __int128 S= 1;
 const double Nh= (double)N + 0.5;
 for(u64 d= 1; d <= R; ++d) {
  const u64 x= (u64)(Nh * T.rt[d]);
  S+= (__int128)(T.mu[d] * (i64)1) * ((__int128)x * x);
 }
 const u64 q1= N / (R + 1);
 for(u64 q= 1; q <= q1; ++q) S+= (__int128)(q * q) * (T.W[q] - (q < q1 ? T.W[q + 1] : (i64)T.M[R]));
 return (unsigned __int128)(S / 2);
}
