#pragma once
#include <algorithm>
#include <cmath>
#include <vector>
// i^d (0 ≤ i ≤ n) の表。0^0 = 1。T は掛け算を持つ型 (ZMod など)。
// i^d は完全乗法的なので、素数 p だけ p^d を冪で求め、合成数 m は最小の素因数 p について m^d = (m/p)^d p^d とする。
// √n 以下の素数で 2^15 個ずつの区画を篩い、合成数 m には m/p と p^d を記録する (素数は大きい順に書いて、最後に書いた最小の
// 素数を残す)。m/p は区画の前か区画の中の m より前にあるので、前から埋めればそろっている。
// procon-judge の yosupo-sum-of-exponential-times-polynomial で書き比べた se1 の表を写したもの。
template <class T> std::vector<T> pow_table(int n, unsigned long long d) {
 std::vector<T> t(n + 1);
 auto pw= [d](T b) {
  T r(1);
  for(auto e= d; e; e>>= 1, b*= b)
   if(e & 1) r*= b;
  return r;
 };
 if(d == 0) return std::fill(t.begin(), t.end(), T(1)), t;
 t[0]= T(0);
 if(n >= 1) t[1]= T(1);
 int sq= int(std::sqrt(double(n)));
 while((long long)sq * sq > n) --sq;
 while((long long)(sq + 1) * (sq + 1) <= n) ++sq;
 std::vector<char> comp(sq + 1);
 std::vector<int> ps;
 std::vector<T> pp;
 for(int i= 2; i <= sq; ++i)
  if(!comp[i]) {
   ps.push_back(i), pp.push_back(pw(T(i)));
   for(int j= i * i; j <= sq; j+= i) comp[j]= 1;
  }
 constexpr int B= 1 << 15;
 std::vector<unsigned> cof(B);
 std::vector<int> pi(B);
 for(int L= 2; L <= n; L+= B) {
  const int R= std::min(n + 1, L + B);
  std::fill(cof.begin(), cof.begin() + (R - L), 0u);
  for(int j= int(ps.size()) - 1; j >= 0; --j) {
   const int p= ps[j];
   const long long st= std::max((long long)p * p, (long long)(L + p - 1) / p * p);
   if(st >= R) continue;
   unsigned k= unsigned(st / p);
   for(int m= int(st); m < R; m+= p, ++k) cof[m - L]= k, pi[m - L]= j;
  }
  for(int m= L; m < R; ++m) {
   const unsigned c= cof[m - L];
   t[m]= c ? t[c] * pp[pi[m - L]] : pw(T(m));
  }
 }
 return t;
}
