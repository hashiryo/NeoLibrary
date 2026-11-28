#pragma once
#include <array>
#include <type_traits>
#include <utility>
#include "neo/algo/compute/nt/factorize.hpp"
#include "neo/base/core/nt/gcd.hpp"
#include "neo/base/core/nt/is_prime.hpp"
// a ↑↑ b (a の b 回のテトレーション) を m で割った余り。0^0 = 1、b = 0 は 1、m = 1 は 0。
//   a ↑↑ 1 = a、a ↑↑ b = a^(a ↑↑ (b-1))
//
// 指数をそのまま作れないので、指数を L(m) で割った余りに落とす。L は、2^e では e <= 3 のとき φ(2^e)、
// e >= 4 のとき λ(2^e) = 2^(e-2)、奇数の素数べきでは φ(p^e) とし、素数べきごとの値を lcm で束ねたもの。
// L(m) は log2(m) 以上なので、a と m が互いに素でなくても、指数が L 以上なら a^E mod m は E mod L で決まる。
// 再帰の結果が 0 のときは、真の指数が L の正の倍数なので L を使う。
// 今の Library の mod_tetration と同じ算法。実行時は、m < 2^32 の因数分解を 65536 未満の素数の表で
// 確保なしに試し割りする (テトレーションの再帰では法が小さくなるので、この道が大半。素数は先に
// is_prime で分ける)。表は定数評価の段数に収まらないので実行時に作り、constexpr のときは
// factorize_pairs の道を通す。
namespace tetration_internal {
using u64= unsigned long long;
using u128= unsigned __int128;
struct PrimeTable {
 static constexpr int LIM= 1 << 16, CAP= 7000;
 std::array<unsigned, CAP> p{};
 int n= 0;
};
inline const PrimeTable& primes() {
 static const PrimeTable t= [] {
  PrimeTable r{};
  std::array<bool, PrimeTable::LIM> comp{};
  for(int i= 2; i < PrimeTable::LIM; ++i)
   if(!comp[i]) {
    r.p[r.n++]= i;
    for(long long j= 1LL * i * i; j < PrimeTable::LIM; j+= i) comp[j]= true;
   }
  return r;
 }();
 return t;
}
[[gnu::always_inline]] constexpr void merge(u64& l, u64 p, int e) {
 u64 t= p - 1;
 const int ex= e - (p == 2 && e > 3 ? 1 : 0);
 for(int i= 1; i < ex; ++i) t*= p;
 l= l / gcd(l, t) * t;
}
constexpr u64 lambda_like(u64 m) {
 u64 l= 1;
 if(m < (1ull << 32) && !std::is_constant_evaluated()) {
  if(is_prime(m)) return merge(l, m, 1), l;  // 素数を 65536 まで割らずに済ませる
  const PrimeTable& P= primes();
  for(int i= 0; i < P.n; ++i) {
   const u64 p= P.p[i];
   if(p * p > m) break;
   if(m % p == 0) {
    int e= 0;
    do m/= p, ++e;
    while(m % p == 0);
    merge(l, p, e);
   }
  }
  if(m > 1) merge(l, m, 1);  // 残りは素数
  return l;
 }
 for(auto [p, e]: factorize_pairs(m)) merge(l, p, e);
 return l;
}
}
constexpr unsigned long long mod_tetration(unsigned long long a, unsigned long long b, unsigned long long m) {
 using namespace tetration_internal;
 if(m == 1) return 0;
 if(a == 0) return (b ^ 1) & 1;
 if(b == 0) return 1;
 if(b == 1) return a % m;
 const u64 l= lambda_like(m);
 u64 e= b == 2 ? a : mod_tetration(a, b - 1, l);
 if(!e) e= l;
 u64 r= 1, x= a % m;
 for(; e; e>>= 1, x= (u64)((u128)x * x % m))
  if(e & 1) r= (u64)((u128)r * x % m);
 return r;
}
