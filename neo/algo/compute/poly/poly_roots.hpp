#pragma once
#include <algorithm>
#include <vector>
#include "neo/base/core/convolve/naive.hpp"
#include "neo/base/core/poly/poly_divmod.hpp"
// 素数 p を法とする多項式 f の根を、重複なしで求める (順は決めていない)。Z は ZMod (法は素数)。f が 0 なら空を返す。
// Cantor–Zassenhaus の方法で、FFT を使わない。p が奇数なら h = x^{(p-1)/2} mod f から gcd(f, h - 1) と gcd(f, h + 1) で
// 0 でない根を平方剰余とそうでないものに分け、それぞれを乱択の a による gcd(g, (x + a)^{(p-1)/2} - 1) で 1 次まで割る。
// p = 2 なら 0 と 1 を代入して調べる。
// 冪は f を法とする剰余の 2 乗で、積は convolve_internal::conv (AVX2)、剰余は rev(f)^{-1} mod x^{N-1} を先に作っておく Barrett の
// 方法 (積の下半分 2 回)。次数 N の f で O(N^2 log p)。
// procon-judge の yosupo-polynomial-root-finding で書き比べた naive を、法を定数にしない形で写したもの。
namespace polyroots_internal {
using u32= unsigned;
using u64= unsigned long long;
using poly= std::vector<u32>;
inline void trim(poly& a) {
 while(!a.empty() && a.back() == 0) a.pop_back();
}
template <class Z> poly mul_low(const poly& a, const poly& b, int len) {
 const int n= int(a.size()), m= int(b.size());
 if(!n || !m || len <= 0) return {};
 len= std::min(len, n + m - 1);
 poly c(len);
 convolve_internal::conv<Z>(a.data(), n, b.data(), m, c.data(), len);
 return c;
}
// monic な f を法とする剰余。inv_rev = rev(f)^{-1} mod x^{N-1}。
template <class Z> struct ModF {
 poly f, inv_rev;
 int N;
 explicit ModF(const poly& f_): f(f_), N(int(f_.size()) - 1) {
  const int L= std::max(N - 1, 1);
  inv_rev.assign(L, 0), inv_rev[0]= Z(1).val();
  for(int k= 1; k < L; ++k) {
   u64 s= 0;
   Z acc(0);
   for(int j= 1, t= 0; j <= std::min(k, N); ++j) {
    s+= u64(f[N - j]) * inv_rev[k - j];
    if(++t == 15) acc+= Z(s), s= 0, t= 0;
   }
   inv_rev[k]= (-(acc + Z(s))).val();
  }
 }
 poly reduce(poly c) const {
  trim(c);
  if(int(c.size()) <= N) return c;
  const int dc= int(c.size()) - 1, L= dc - N + 1;
  poly cr(L);
  for(int i= 0; i < L; ++i) cr[i]= c[dc - i];
  poly qr= mul_low<Z>(cr, inv_rev, L);
  qr.resize(L);
  std::reverse(qr.begin(), qr.end());
  const poly qf= mul_low<Z>(qr, f, N);
  poly r(N);
  for(int k= 0; k < N; ++k) r[k]= (Z::raw(c[k]) - Z::raw(k < int(qf.size()) ? qf[k] : 0)).val();
  trim(r);
  return r;
 }
 // a (x + s) mod f
 poly mul_linear(const poly& a, u32 s) const {
  poly r(a.size() + 1);
  for(size_t i= 0; i < a.size(); ++i) r[i + 1]= (Z::raw(r[i + 1]) + Z::raw(a[i])).val(), r[i]= (Z::raw(r[i]) + Z::raw(a[i]) * Z::raw(s)).val();
  trim(r);
  if(int(r.size()) > N) {
   const Z lc= Z::raw(r[N]);
   for(int i= 0; i < N; ++i) r[i]= (Z::raw(r[i]) - lc * Z::raw(f[i])).val();
   r.resize(N), trim(r);
  }
  return r;
 }
 // (x + s)^e mod f
 poly pow_linear(u32 s, u64 e) const {
  poly r{Z(1).val()};
  for(int b= 63; b >= 0; --b) {
   if(r.size() > 1) {
    const int n= int(r.size());
    poly sq(2 * n - 1);
    convolve_internal::conv<Z>(r.data(), n, r.data(), n, sq.data(), 2 * n - 1);
    r= reduce(std::move(sq));
   } else if(r.size() == 1) r[0]= (Z::raw(r[0]) * Z::raw(r[0])).val();
   if(e >> b & 1) r= mul_linear(r, s);
  }
  return r;
 }
};
template <class Z> std::vector<Z> to_z(const poly& a) { return convolve_internal::from_u32<Z>(a); }
// h + d (d は定数)。空の vector に push_back しない形にする (gcc 15 が -Wfree-nonheap-object を誤って出すため)。
template <class Z> poly add_scalar(const poly& h, Z d) {
 poly r(std::max<std::size_t>(h.size(), 1), 0);
 std::copy(h.begin(), h.end(), r.begin());
 r[0]= (Z::raw(r[0]) + d).val();
 trim(r);
 return r;
}
template <class Z> poly gcd(const poly& a, const poly& b) { return convolve_internal::to_u32(poly_gcd(to_z<Z>(a), to_z<Z>(b))); }
template <class Z> poly divide(const poly& g, const poly& d) { return convolve_internal::to_u32(poly_divmod(to_z<Z>(g), to_z<Z>(d)).first); }
struct Rng {
 u64 x= 88172645463325252ull;
 u64 operator()() { return x^= x << 7, x^= x >> 9; }
};
template <class Z> void split(const poly& g, std::vector<Z>& roots, Rng& rng) {
 const int dg= int(g.size()) - 1;
 if(dg <= 0) return;
 if(dg == 1) return roots.push_back(-Z::raw(g[0])), void();  // x + g0 (monic)
 const ModF<Z> mf(g);
 const u32 M= Z::mod();
 while(true) {
  const u32 s= u32(rng() % M);
  const poly h= add_scalar<Z>(mf.pow_linear(s, (M - 1) / 2), -Z(1));
  const poly d= gcd<Z>(g, h);
  const int dd= int(d.size()) - 1;
  if(0 < dd && dd < dg) {
   split<Z>(d, roots, rng);
   split<Z>(divide<Z>(g, d), roots, rng);
   return;
  }
 }
}
}  // namespace polyroots_internal
template <class Z> std::vector<Z> poly_roots(std::vector<Z> f) {
 using namespace polyroots_internal;
 polydiv_internal::trim(f);
 std::vector<Z> roots;
 if(f.size() <= 1) return roots;
 int z= 0;
 while(f[z] == Z()) ++z;
 if(z) roots.push_back(Z()), f.erase(f.begin(), f.begin() + z);
 if(f.size() <= 1) return roots;
 const u32 M= Z::mod();
 if(M == 2) {  // 残りの根は 1 だけがありうる
  Z s(0);
  for(auto& x: f) s+= x;
  if(s == Z()) roots.push_back(Z(1));
  return roots;
 }
 const Z il= f.back().inv();
 for(auto& x: f) x*= il;
 poly g= convolve_internal::to_u32(f);
 Rng rng;
 const ModF<Z> mf(g);
 const poly h= mf.pow_linear(0, (M - 1) / 2);  // x^{(p-1)/2} mod f
 split<Z>(gcd<Z>(g, add_scalar<Z>(h, -Z(1))), roots, rng);
 split<Z>(gcd<Z>(g, add_scalar<Z>(h, Z(1))), roots, rng);
 return roots;
}
