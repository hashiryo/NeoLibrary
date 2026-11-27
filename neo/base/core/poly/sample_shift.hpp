#pragma once
// 標本点のずらし。y_i = f(i) (0 <= i < n、f は n-1 次以下) から f(c), f(c+1), ..., f(c+m-1) を求める。
// Lagrange の形 f(x) = P(x) Σ_i a_i/(x-i) (P(x) = Π_j (x-j)、a_i = y_i (-1)^{n-1-i}/(i! (n-1-i)!)) を、
// a と 1/(c-n+1+k) の 1 回の畳み込みと、P(c+j) の窓の積で作る。手間は O(M(n + m)) で、M は畳み込みの手間。
// c+j が標本点 (0 <= i < n) に重なるところは P(c+j) = 0 で割り算が壊れるので、f(c+j) = y_i をそのまま返す。
// v_k = c-(n-1)+k (k < n+m-1) の逆数と P の窓の積は、0 を 1 に置き換えた前積と 1 回の逆数でまとめて作る。
// 核はテンプレート引数 Conv で受け、畳み込みの共通の口 (convolve.hpp) の convolve_range に渡す。
// K は体で、ZMod なら法は素数。1/i を使うので、i < n は法より小さいこと。
#include <algorithm>
#include <cassert>
#include <vector>
#include "neo/base/core/convolve/convolve.hpp"

namespace sample_points_shift_internal {
// 0 を 1 に置き換えて、まとめて逆数を作る。
template <class K> std::vector<K> batch_inverse(const K* v, int n) {
 std::vector<K> pre(n + 1);
 pre[0]= K(1);
 for (int i= 0; i < n; ++i) pre[i + 1]= pre[i] * (v[i] == K() ? K(1) : v[i]);
 std::vector<K> res(n);
 K t= pre[n].inv();
 for (int i= n - 1; i >= 0; --i) {
  const K x= (v[i] == K() ? K(1) : v[i]);
  res[i]= t * pre[i];
  t*= x;
 }
 return res;
}
}  // namespace sample_points_shift_internal

template <class Conv= convolve_internal::Auto, class K>
std::vector<K> sample_points_shift(const std::vector<K>& y, K c, int m) {
 using sample_points_shift_internal::batch_inverse;
 const int n= int(y.size());
 std::vector<K> out(m, K());
 if (m == 0 || n == 0) return out;
 if (n == 1) {
  std::fill(out.begin(), out.end(), y[0]);
  return out;
 }
 // 階乗とその逆数。
 std::vector<K> fact(n), ifact(n);
 fact[0]= K(1);
 for (int i= 1; i < n; ++i) fact[i]= fact[i - 1] * K(i);
 ifact[n - 1]= fact[n - 1].inv();
 for (int i= n - 1; i > 0; --i) ifact[i - 1]= ifact[i] * K(i);
 std::vector<K> a(n);
 for (int i= 0; i < n; ++i) {
  const K ai= y[i] * ifact[i] * ifact[n - 1 - i];
  a[i]= ((n - 1 - i) & 1) ? -ai : ai;
 }
 // v_k = c - (n - 1) + k と、その逆数 (0 の所はあとで上書きするので 0 にしておく)。
 const int sz= n + m - 1;
 std::vector<K> v(sz);
 int z= -1;
 for (int k= 0; k < sz; ++k) {
  v[k]= c + K(k - (n - 1));
  if (v[k] == K()) z= k;
 }
 std::vector<K> b= batch_inverse(v.data(), int(v.size()));
 if (z >= 0) b[z]= K();
 // 1 回の畳み込みで Σ_i a_i/(c+j-i) を作る。
 std::vector<K> sum(m);
 convolve_range<Conv>(a.data(), n, b.data(), sz, n - 1, n - 1 + m, sum.data());
 // P(c+j) は v の窓 [j, j+n-1] の積。
 std::vector<K> pre(sz + 1);
 pre[0]= K(1);
 for (int k= 0; k < sz; ++k) pre[k + 1]= pre[k] * (v[k] == K() ? K(1) : v[k]);
 const std::vector<K> pinv= batch_inverse(pre.data(), m);
 for (int j= 0; j < m; ++j) {
  const bool covered= z >= 0 && j <= z && z <= j + n - 1;
  out[j]= covered ? K() : pre[j + n] * pinv[j] * sum[j];
 }
 // 標本点に重なった所は y をそのまま返す。
 if (z >= 0)
  for (int j= std::max(0, z - n + 1); j <= std::min(m - 1, z); ++j) out[j]= y[n - 1 - z + j];
 return out;
}
