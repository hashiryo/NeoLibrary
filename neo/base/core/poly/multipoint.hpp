#pragma once
// 多点評価と補間。
// 多点評価は、部分積木 (葉が X - x_i) と余りの木で f(x_i) を求める。余りは poly_divmod_fast (fps::inv の上) で作り、
// 手間は O(M(n) log n) (M は畳み込みの手間)。
// 補間は、根の多項式 P の導関数 P' を多点評価して P'(x_i) を作り、w_i = y_i / P'(x_i) から
// f = Σ_i w_i P/(X - x_i) を下から f_left P_right + f_right P_left でまとめる。
// 核はテンプレート引数 Conv で受け、畳み込みの共通の口 (convolve.hpp) の convolve_range に渡す。
// K は体で、ZMod なら法は素数。点は相異なること。小さい n では素直な木より lagrange::interpolate (O(n^2)) が速い。
#include <algorithm>
#include <cassert>
#include <vector>
#include "neo/base/core/convolve/convolve.hpp"
#include "neo/base/core/poly/poly_divmod_fast.hpp"

namespace multipoint_internal {
template <class Conv, class K> std::vector<K> mul(const std::vector<K>& a, const std::vector<K>& b) {
 if (a.empty() || b.empty()) return {};
 std::vector<K> c(a.size() + b.size() - 1);
 convolve_range<Conv>(a.data(), int(a.size()), b.data(), int(b.size()), 0, int(c.size()), c.data());
 return c;
}
// 部分積木。tr[sz + i] = X - x_i (i < m)、余りは 1 にして積で消す。tr[1] が根。
template <class Conv, class K> std::vector<std::vector<K>> build_tree(const std::vector<K>& x, int& sz) {
 const int m= int(x.size());
 sz= 1;
 while (sz < m) sz<<= 1;
 std::vector<std::vector<K>> tr(2 * sz);
 for (int i= 0; i < m; ++i) tr[sz + i]= {-x[i], K(1)};
 for (int i= m; i < sz; ++i) tr[sz + i]= {K(1)};
 for (int i= sz - 1; i >= 1; --i) tr[i]= mul<Conv>(tr[2 * i], tr[2 * i + 1]);
 return tr;
}
// 余りの木。rem[i] = f mod tr[i]。点の無い節は空。
template <class Conv, class K>
std::vector<std::vector<K>> remainders(const std::vector<std::vector<K>>& tr, int sz, const std::vector<K>& f) {
 std::vector<std::vector<K>> rem(2 * sz);
 rem[1]= poly_divmod_fast<Conv>(f, tr[1]).second;
 for (int i= 1; i < sz; ++i) {
  if (tr[2 * i].size() == 1) continue;  // 点の無い節 (点は前詰めなので、左が無ければ右も無い)
  rem[2 * i]= poly_divmod_fast<Conv>(rem[i], tr[2 * i]).second;
  rem[2 * i + 1]= poly_divmod_fast<Conv>(rem[i], tr[2 * i + 1]).second;
 }
 return rem;
}
}  // namespace multipoint_internal

// f の x[0..m) での値。
template <class Conv= convolve_internal::Auto, class K>
std::vector<K> multipoint_eval(const std::vector<K>& f, const std::vector<K>& x) {
 using namespace multipoint_internal;
 const int m= int(x.size());
 std::vector<K> res(m);
 if (m == 0 || f.empty()) return res;
 if (f.size() == 1) {
  std::fill(res.begin(), res.end(), f[0]);
  return res;
 }
 int sz;
 const auto tr= build_tree<Conv, K>(x, sz);
 const auto rem= remainders<Conv, K>(tr, sz, f);
 for (int i= 0; i < m; ++i) res[i]= rem[sz + i].empty() ? K() : rem[sz + i][0];
 return res;
}

// 点 (x_i, y_i) (i < n) を通る n-1 次以下の多項式の係数。点は相異なること。
template <class Conv= convolve_internal::Auto, class K>
std::vector<K> interpolate(const std::vector<K>& x, const std::vector<K>& y) {
 using namespace multipoint_internal;
 const int n= int(x.size());
 assert(int(y.size()) == n);
 if (n == 0) return {};
 int sz;
 const auto tr= build_tree<Conv, K>(x, sz);
 // P' を多点評価して P'(x_i) を作り、w_i = y_i / P'(x_i) にする。
 std::vector<K> dp(tr[1].size() - 1);
 for (int i= 1; i < int(tr[1].size()); ++i) dp[i - 1]= tr[1][i] * K(i);
 auto rem= remainders<Conv, K>(tr, sz, dp);
 for (int i= 0; i < n; ++i) rem[sz + i]= {y[i] / rem[sz + i][0]};
 for (int i= n; i < sz; ++i) rem[sz + i].clear();
 // 下から f = f_left P_right + f_right P_left。
 for (int i= sz - 1; i >= 1; --i) {
  std::vector<K> a= mul<Conv>(rem[2 * i], tr[2 * i + 1]);
  const std::vector<K> b= mul<Conv>(rem[2 * i + 1], tr[2 * i]);
  if (b.size() > a.size()) a.resize(b.size());
  for (int j= 0; j < int(b.size()); ++j) a[j]+= b[j];
  rem[i]= std::move(a);
 }
 rem[1].resize(n);
 return rem[1];
}
