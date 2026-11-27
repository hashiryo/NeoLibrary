#pragma once
// 形式的冪級数の道具。係数の型 K は std::vector<K> で受け、先頭 n 項だけを求める。
// 核はテンプレート引数 Conv で受けて、畳み込みの共通の口 (convolve.hpp) の convolve_range と convolve に渡す。
// Conv を省くと、include した核のうち K を受け持つ既定の核になる (convolve.hpp の kernel_t を参照)。
#include <algorithm>
#include <cassert>
#include <vector>
#include "neo/base/core/convolve/convolve.hpp"

namespace fps {
// f の先頭 n 項の逆数 g (f g = 1 mod x^n)。f[0] は逆元を持つこと。
// Newton 法: g <- g (2 - f g) を、長さ m から 2m ずつ倍にする。
// 1 段は、e = (f g)[m, 2m) を範囲の積で求め、g[m, 2m) = -(g e)[0, m) で直す (f g = 1 + x^m e なので g (2 - f g) = g - x^m g e)。
template <class Conv= convolve_internal::Auto, class K> std::vector<K> inv(const std::vector<K>& f, int n) {
 assert(n >= 1);
 assert(!f.empty() && f[0] != K());
 std::vector<K> g(1, f[0].inv());
 for (int m= 1; m < n;) {
  const int m2= std::min(2 * m, n), cnt= m2 - m;
  const int fm= std::min(int(f.size()), m2);  // f の m2 次以上は積の [m, m2) に入らない
  std::vector<K> e(cnt), u(cnt);
  convolve_range<Conv>(f.data(), fm, g.data(), m, m, m2, e.data());
  convolve_range<Conv>(e.data(), cnt, g.data(), m, 0, cnt, u.data());
  g.resize(m2);
  for (int i= m; i < m2; ++i) g[i]= -u[i - m];
  m= m2;
 }
 return g;
}
}  // namespace fps
