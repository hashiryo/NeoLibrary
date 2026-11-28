#pragma once
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>
// floor(N/i) が同じになる i の区間を列挙する。N = 10 なら、q = 1, 2, 3, 5, 10 の 5 つで、
//   (1, 5, 10)、(2, 3, 5)、(3, 2, 3)、(5, 1, 2)、(10, 0, 1)
// を q の昇順で返す。tuple は (q, l, r) で、i in (l, r] のとき floor(N/i) = q。
// 区間の数は floor(2 sqrt(N)) ほど。l を 1 から進めて r = N / q で飛ばす形で、q は大きい順に出るので最後に逆順にする。
inline std::vector<std::tuple<unsigned long long, unsigned long long, unsigned long long>>
enumerate_quotients(unsigned long long N) {
 std::vector<std::tuple<unsigned long long, unsigned long long, unsigned long long>> ret;
 if(!N) return ret;
 unsigned long long s= std::sqrt((long double)N);
 while(s > 0 && s * s > N) --s;
 while((s + 1) <= N / (s + 1)) ++s;
 ret.reserve(2 * s + 1);
 for(unsigned long long l= 1;;) {
  const unsigned long long q= N / l, r= N / q;
  ret.emplace_back(q, l - 1, r);
  if(r == N) break;
  l= r + 1;
 }
 return std::reverse(ret.begin(), ret.end()), ret;
}
// F(n) = sum_{i=1}^n f(i) を、n in S = {floor(N/x) | 1 <= x <= N} の所だけ持つ配列。
// i <= K = floor(sqrt(N)) は X[i]、i > K は X[K + N/i] に置く (S の要素は K 以下と、N/j の形の 2 種類しかない)。
// ベクトルと同じく +、-、スカラー倍ができる。add(i, v) は f(i) <- f(i) + v で、F(n) を n >= i の所へ足す。
// 今の Library の CumSumQuotient の add は、i > K + N/i のときに何もしない (i が N/j の形の大きい値のとき)。
// ここでは X[K + j] (j <= N/i) と X[i..K] の両方を足す形にして、どの i でも正しく動くようにした。
template <class T> struct CumSumQuotient {
 using u64= unsigned long long;
 u64 N;
 std::size_t K;
 std::vector<T> X;
 explicit CumSumQuotient(u64 N_): N(N_) {
  u64 k= std::sqrt((long double)N_);
  while(k > 0 && k * k > N_) --k;
  while((k + 1) <= N_ / (k + 1)) ++k;
  K= k, X.assign(2 * K + 1, T());
 }
 [[gnu::always_inline]] std::size_t index(u64 n) const { return n > K ? K + N / n : n; }
 // n in S
 [[gnu::always_inline]] T& operator[](u64 n) { return X[index(n)]; }
 [[gnu::always_inline]] T operator()(u64 n) const { return X[index(n)]; }
 // f(i) <- f(i) + v。F(n) を n >= i の所へ足す。
 void add(u64 i, T v) {
  assert(1 <= i && i <= N);
  const u64 jm= std::min<u64>(N / i, K);
  for(u64 j= 1; j <= jm; ++j) X[K + j]+= v;
  for(u64 n= std::max<u64>(i, 1); n <= K; ++n) X[n]+= v;
 }
 // F(n) = sum_{i=1}^n f(i)
 [[gnu::always_inline]] T sum(u64 n) const { return X[index(n)]; }
 // F(N)。N = 1 では K + 1 が N の位置と違うので、N の位置を引く (Library は X[K + 1] を返していた)。
 [[gnu::always_inline]] T sum() const { return N ? X[index(N)] : T(); }
 CumSumQuotient& operator+=(const CumSumQuotient& r) {
  for(std::size_t i= 0; i < X.size(); ++i) X[i]+= r.X[i];
  return *this;
 }
 CumSumQuotient& operator-=(const CumSumQuotient& r) {
  for(std::size_t i= 0; i < X.size(); ++i) X[i]-= r.X[i];
  return *this;
 }
 CumSumQuotient& operator*=(T a) {
  for(std::size_t i= 0; i < X.size(); ++i) X[i]*= a;
  return *this;
 }
 CumSumQuotient operator-() const {
  CumSumQuotient ret= *this;
  for(std::size_t i= 0; i < ret.X.size(); ++i) ret.X[i]= -ret.X[i];
  return ret;
 }
 CumSumQuotient operator+(const CumSumQuotient& r) const { return CumSumQuotient(*this)+= r; }
 CumSumQuotient operator-(const CumSumQuotient& r) const { return CumSumQuotient(*this)-= r; }
 CumSumQuotient operator*(T a) const { return CumSumQuotient(*this)*= a; }
 friend CumSumQuotient operator*(T a, const CumSumQuotient& x) { return x * a; }
};
