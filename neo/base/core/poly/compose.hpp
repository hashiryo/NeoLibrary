#pragma once
// 合成 f(g) mod x^n。g の定数項は 0 であること。
// Kinoshita–Li の転置の power projection で、O(M(n) log n) (M は畳み込みの手間)。
// 分母 q_t(x, y) = Π_{j < 2^t} (1 - y g(ω^j x^{1/2^t})) を y について Kronecker 置換で 1 本の配列にし、
// q_t(x, y) q_t(-x, y) の偶数次の部分を取って q_{t+1}(x, y) = q_t(√x, y) q_t(-√x, y) を作る (x の次数は (n-1)>>t まで)。
// そのあと t = T-1 から 0 まで、q_t を逆向きに掛けて窓を取る転置の段を回すと、f の係数が答えの係数に移る。
// 使うのは範囲の積 (convolve_range) だけで、核に変換が無くても動く。
// 2 変数の配列は [x の次数][y の次数] の順に並べ、積は y の次数がはみ出さない間隔 (stride) で 1 変数に潰して取る。
#include <algorithm>
#include <cassert>
#include <vector>
#include "neo/base/core/convolve/convolve.hpp"

template <class Conv= convolve_internal::Auto, class K>
std::vector<K> compose(const std::vector<K>& f, const std::vector<K>& g, int n) {
 std::vector<K> res(n, K());
 if (n <= 0) return res;
 assert(g.empty() || g[0] == K());  // g(0) = 0
 bool zero= true;
 for (int i= 0; i < n && i < int(g.size()); ++i)
  if (g[i] != K()) zero= false;
 if (zero || f.empty()) {
  res[0]= f.empty() ? K() : f[0];
  return res;
 }
 int levels= 0;
 for (int m= n - 1; m; m>>= 1) ++levels;
 // 分母 q_t。行が x の次数 (0..(n-1)>>t)、列が y の次数 (0..2^t)。q_0 = 1 - y g(x)。
 std::vector<std::vector<K>> qs(levels + 1);
 {
  std::vector<K> q(2 * n, K());
  for (int i= 0; i < n && i < int(g.size()); ++i) q[2 * i + 1]= -g[i];
  q[0]+= K(1);
  qs[0]= std::move(q);
 }
 for (int t= 0; t < levels; ++t) {
  const int mt= (n - 1) >> t, mn= mt >> 1;
  const int cols= (1 << t) + 1, stride= (1 << (t + 1)) + 1;
  const std::vector<K>& q= qs[t];
  std::vector<K> a(mt * stride + cols, K()), b(a.size(), K());
  for (int i= 0; i <= mt; ++i)
   for (int y= 0; y < cols; ++y) {
    a[i * stride + y]= q[i * cols + y];
    b[i * stride + y]= (i & 1) ? -q[i * cols + y] : q[i * cols + y];
   }
  std::vector<K> c(2 * a.size() - 1);
  convolve_range<Conv>(a.data(), int(a.size()), b.data(), int(b.size()), 0, int(c.size()), c.data());
  std::vector<K> next((mn + 1) * stride, K());
  for (int j= 0; j <= mn; ++j)
   for (int y= 0; y < stride; ++y) {
    const int pos= 2 * j * stride + y;
    if (pos < int(c.size())) next[j * stride + y]= c[pos];
   }
  qs[t + 1]= std::move(next);
 }
 // 転置の段。z は t 段で (m_t + 1) 行 × 2^t 列。
 std::vector<K> z(std::size_t(1) << levels, K());
 for (int k= 0; k < n && k < int(f.size()); ++k) z[k]= f[k];
 for (int t= levels - 1; t >= 0; --t) {
  const int mt= (n - 1) >> t, bt= mt & 1, mn= mt >> 1;
  const int w= 1 << (t + 1), cols_q= (1 << t) + 1;
  std::vector<K> r(std::size_t(mt + 1) * w, K());
  for (int j= 0; j <= mn; ++j)
   for (int u= 0; u < w; ++u) r[(2 * j + bt) * w + u]= z[j * w + u];
  const std::vector<K>& q= qs[t];
  const int len_q= mt * w + cols_q;
  std::vector<K> qt(len_q, K());
  for (int i= 0; i <= mt; ++i)
   for (int y= 0; y < cols_q; ++y) qt[i * w + y]= (i & 1) ? -q[i * cols_q + y] : q[i * cols_q + y];
  std::reverse(qt.begin(), qt.end());
  const int cols_new= 1 << t, lo= len_q - 1, hi= mt * w + (cols_new - 1) + len_q;
  std::vector<K> conv(hi - lo);
  convolve_range<Conv>(r.data(), int(r.size()), qt.data(), len_q, lo, hi, conv.data());
  std::vector<K> next(std::size_t(mt + 1) * cols_new, K());
  for (int i= 0; i <= mt; ++i)
   for (int u= 0; u < cols_new; ++u) {
    const int pos= i * w + u;
    if (pos < int(conv.size())) next[std::size_t(i) * cols_new + u]= conv[pos];
   }
  z= std::move(next);
 }
 for (int m= 0; m < n; ++m) res[m]= z[n - 1 - m];
 return res;
}
