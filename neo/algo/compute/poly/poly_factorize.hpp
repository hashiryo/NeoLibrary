#pragma once
#include <utility>
#include <vector>
#include "neo/base/core/poly/convolve_naive.hpp"
#include "neo/base/core/poly/poly_divmod.hpp"
// 素数 p を法とする多項式 f (0 でない) を、monic な既約多項式の積に分解する。返すのは (既約因子、重複度) の組で、順は決めて
// いない。f の最高次の係数は落とす (monic にしてから分解する)。Z は ZMod (法は素数)。FFT を使わない素朴な形で、次数 100
// 程度までを想定している。
// 1. 平方因子なしの分解: f' と gcd を取る Yun の方法。f' = 0 の部分は f = h(x^p) なので h に直して重複度を p 倍する。
// 2. 次数ごとの分解: 法 s での Frobenius h -> h^p を、行 x^{ip} mod s を並べた行列で O(n^2) で当て、gcd(s, x^{p^d} - x)。
// 3. 同じ次数 d の積の分解: 乱択の a について、p が奇数なら gcd(g, a^{(p^d - 1)/2} - 1)、p = 2 なら
//    gcd(g, a + a^2 + ... + a^{2^{d-1}})。
// procon-judge の yosupo-factorization-of-polynomials で書き比べた naive を、法の型を引数にしたもの。
namespace polyfactor_internal {
using u64= unsigned long long;
template <class Z> using poly= std::vector<Z>;
template <class Z> int deg(const poly<Z>& a) { return int(a.size()) - 1; }
template <class Z> void monic(poly<Z>& a) {
 if(a.empty()) return;
 const Z il= a.back().inv();
 for(auto& x: a) x*= il;
}
template <class Z> poly<Z> mul(const poly<Z>& a, const poly<Z>& b) { return convolve_naive(a, b); }
template <class Z> poly<Z> rem(const poly<Z>& a, const poly<Z>& b) { return poly_divmod(a, b).second; }
template <class Z> poly<Z> quo(const poly<Z>& a, const poly<Z>& b) { return poly_divmod(a, b).first; }
template <class Z> poly<Z> deriv(const poly<Z>& a) {
 poly<Z> r(a.size() > 1 ? a.size() - 1 : 0);
 for(size_t i= 1; i < a.size(); ++i) r[i - 1]= a[i] * Z(i);
 polydiv_internal::trim(r);
 return r;
}
template <class Z> poly<Z> powmod(poly<Z> b, u64 e, const poly<Z>& m) {
 poly<Z> r{Z(1)};
 b= rem(b, m);
 for(; e; e>>= 1, b= rem(mul(b, b), m))
  if(e & 1) r= rem(mul(r, b), m);
 return r;
}
// 法 m (n 次) での Frobenius の行列。Q[i] = x^{ip} mod m (i < n)。
template <class Z> std::vector<poly<Z>> frob_matrix(const poly<Z>& m) {
 const int n= deg(m);
 std::vector<poly<Z>> Q(n);
 Q[0]= poly<Z>{Z(1)};
 if(n == 1) return Q;
 const poly<Z> xp= powmod(poly<Z>{Z(), Z(1)}, Z::mod(), m);
 for(int i= 1; i < n; ++i) Q[i]= rem(mul(Q[i - 1], xp), m);
 return Q;
}
// h^p mod m (h は m より次数が低い)
template <class Z> poly<Z> frob(const poly<Z>& h, const std::vector<poly<Z>>& Q) {
 poly<Z> r(Q.size());
 for(size_t i= 0; i < h.size(); ++i) {
  if(h[i] == Z()) continue;
  const auto hi= h[i].fixed();
  for(size_t j= 0; j < Q[i].size(); ++j) r[j]+= Q[i][j] * hi;
 }
 polydiv_internal::trim(r);
 return r;
}
struct Rng {
 u64 x= 88172645463325252ull;
 u64 operator()() { return x^= x << 7, x^= x >> 9; }
};
template <class Z> using Out= std::vector<std::pair<poly<Z>, int>>;
template <class Z> void edf(const poly<Z>& g, int d, int mult, Out<Z>& out, Rng& rng) {
 const int n= deg(g);
 if(n == d) return out.push_back({g, mult}), void();
 const auto Q= frob_matrix(g);
 const u64 p= Z::mod();
 while(true) {
  poly<Z> a(n);
  for(auto& x: a) x= Z(rng() % p);
  polydiv_internal::trim(a);
  if(deg(a) <= 0) continue;
  poly<Z> b;
  if(p == 2) {
   poly<Z> t= a;
   b= a;
   for(int j= 1; j < d; ++j) {
    t= frob(t, Q);
    if(b.size() < t.size()) b.resize(t.size());
    for(size_t k= 0; k < t.size(); ++k) b[k]+= t[k];
   }
   polydiv_internal::trim(b);
  } else {
   poly<Z> t= a, pr= a;
   for(int j= 1; j < d; ++j) t= frob(t, Q), pr= rem(mul(pr, t), g);
   b= powmod(pr, (p - 1) / 2, g);
   if(b.empty()) b.push_back(Z());
   b[0]-= Z(1), polydiv_internal::trim(b);
  }
  const poly<Z> h= poly_gcd(g, b);
  if(0 < deg(h) && deg(h) < n) {
   edf(h, d, mult, out, rng);
   edf(quo(g, h), d, mult, out, rng);
   return;
  }
 }
}
// 平方因子のない monic な s を分解する。
template <class Z> void factor_sqfree(const poly<Z>& s, int mult, Out<Z>& out, Rng& rng) {
 if(deg(s) == 1) return out.push_back({s, mult}), void();
 const auto Q= frob_matrix(s);
 poly<Z> cur= s, h{Z(), Z(1)};
 for(int d= 1; 2 * d <= deg(cur); ++d) {
  h= rem(frob(h, Q), cur);
  poly<Z> t= h;
  if(t.size() < 2) t.resize(2);
  t[1]-= Z(1), polydiv_internal::trim(t);
  const poly<Z> g= poly_gcd(cur, t);
  if(deg(g) > 0) {
   edf(g, d, mult, out, rng);
   cur= quo(cur, g);
   h= rem(h, cur);
  }
 }
 if(deg(cur) > 0) out.push_back({cur, mult});
}
template <class Z> void sff(poly<Z> f, int mult, Out<Z>& out, Rng& rng) {
 if(deg(f) <= 0) return;
 const u64 p= Z::mod();
 const poly<Z> g= deriv(f);
 if(g.empty()) {  // f = h(x^p)
  poly<Z> h(deg(f) / p + 1);
  for(size_t i= 0; i < h.size(); ++i) h[i]= f[i * p];
  return sff(h, mult * int(p), out, rng);
 }
 poly<Z> c= poly_gcd(f, g), w= quo(f, c);
 for(int i= 1; deg(w) > 0; ++i) {
  const poly<Z> y= poly_gcd(w, c);
  poly<Z> fac= quo(w, y);
  if(deg(fac) > 0) monic(fac), factor_sqfree(fac, i * mult, out, rng);
  w= y, c= quo(c, y);
 }
 if(deg(c) > 0) {  // 残りは重複度が p の倍数の因子で、c = h(x^p)
  poly<Z> h(deg(c) / p + 1);
  for(size_t i= 0; i < h.size(); ++i) h[i]= c[i * p];
  monic(h);
  sff(h, mult * int(p), out, rng);
 }
}
}  // namespace polyfactor_internal
template <class Z> std::vector<std::pair<std::vector<Z>, int>> poly_factorize(std::vector<Z> f) {
 using namespace polyfactor_internal;
 polydiv_internal::trim(f);
 Out<Z> out;
 if(f.empty()) return out;
 monic(f);
 Rng rng;
 sff(f, 1, out, rng);
 return out;
}
