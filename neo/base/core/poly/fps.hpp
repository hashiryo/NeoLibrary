#pragma once
// 形式的冪級数の道具。係数の型 K は std::vector<K> で受け、先頭 n 項だけを求める。
// 核はテンプレート引数 Conv で受けて、畳み込みの共通の口 (convolve.hpp) の convolve_range と convolve に渡す。
// Conv を省くと、include した核のうち K を受け持つ既定の核になる (convolve.hpp の kernel_t を参照)。
#include <algorithm>
#include <cassert>
#include <vector>
#include "neo/base/core/algebra/ZMod.hpp"
#include "neo/base/core/convolve/convolve.hpp"
#include "neo/base/core/nt/sqrt_mod.hpp"

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

// f の先頭 n 項の log (f[0] = 1)。log f = ∫ f' / f。
template <class Conv= convolve_internal::Auto, class K> std::vector<K> log(const std::vector<K>& f, int n) {
 assert(n >= 1);
 assert(!f.empty() && f[0] == K(1));
 std::vector<K> g(n);
 if (n == 1) return g;                             // log 1 = 0
 const std::vector<K> iv= inv<Conv>(f, n - 1);     // 1/f の n - 1 項
 std::vector<K> t(n - 1);                          // f' / f の n - 1 項
 const int dlen= std::min(int(f.size()), n) - 1;
 if (dlen > 0) {
  std::vector<K> df(dlen);
  for (int i= 0; i < dlen; ++i) df[i]= f[i + 1] * K(i + 1);
  convolve_range<Conv>(df.data(), dlen, iv.data(), n - 1, 0, n - 1, t.data());
 }
 for (int i= 1; i < n; ++i) g[i]= t[i - 1] / K(i);
 return g;
}

// f の先頭 n 項の exp (f[0] = 0)。g[0] = 1 から、g <- g (1 - log g + f) を倍々に。
template <class Conv= convolve_internal::Auto, class K> std::vector<K> exp(const std::vector<K>& f, int n) {
 assert(n >= 1);
 assert(f.empty() || f[0] == K());
 std::vector<K> g(1, K(1));
 for (int m= 1; m < n;) {
  const int m2= std::min(2 * m, n);
  const std::vector<K> lg= log<Conv>(g, m2);      // log g の m2 項 (g[0] = 1)
  std::vector<K> h(m2);
  h[0]= K(1);
  for (int i= 0; i < m2; ++i) h[i]+= (i < int(f.size()) ? f[i] : K()) - lg[i];
  std::vector<K> ng(m2);
  convolve_range<Conv>(g.data(), m, h.data(), m2, 0, m2, ng.data());
  g.swap(ng);
  m= m2;
 }
 return g;
}

namespace fps_internal {
// K の e 乗 (e は 64 bit)。K に pow が無い型でも使えるように、掛け算だけで作る。
template <class K> K kpow(K c, unsigned long long e) {
 K r(1), b= c;
 for (; e; e>>= 1, b*= b)
  if (e & 1) r*= b;
 return r;
}
}  // namespace fps_internal

// f の先頭 n 項の k 乗。f[0] = 0 のときは、0 でない最初の次数を v として x^{v k} を外に出し、
// 残り (先頭が 1) を exp(k log) で作る。0^0 は 1 とする。
template <class Conv= convolve_internal::Auto, class K> std::vector<K> pow(const std::vector<K>& f, unsigned long long k, int n) {
 assert(n >= 1);
 std::vector<K> g(n);
 if (k == 0) {
  g[0]= K(1);
  return g;
 }
 int v= 0;
 while (v < int(f.size()) && f[v] == K()) ++v;
 if (v == int(f.size())) return g;                          // f = 0
 if (v > 0 && k > (unsigned long long)(n - 1) / v) return g;  // x^{v k} が n 次以上
 const int sk= v * int(k);
 const K c= f[v], ci= c.inv();
 std::vector<K> u(std::min(int(f.size()), n) - v);
 for (int i= 0; i < int(u.size()); ++i) u[i]= f[v + i] * ci;  // u[0] = 1
 std::vector<K> lg= log<Conv>(u, n - sk);
 const K kk(k);
 for (auto& x: lg) x*= kk;
 const std::vector<K> t= exp<Conv>(lg, n - sk);
 const K ck= fps_internal::kpow(c, k);
 for (int i= 0; i < int(t.size()); ++i) g[sk + i]= t[i] * ck;
 return g;
}

// f の先頭 n 項の平方根。平方根が無ければ空を返す。
// 先頭が 0 のときは、0 でない最初の次数 v が偶数なら x^{v/2} を外に出す (v が奇数なら無い)。
// 残りの定数 c の平方根は、K が ZMod なら sqrt_mod で求め (無ければ空)、それ以外の型では c = 1 のときだけ扱う。
// 外に出した残りは exp(log / 2) で作る。
template <class Conv= convolve_internal::Auto, class K> std::vector<K> sqrt(const std::vector<K>& f, int n) {
 assert(n >= 1);
 int v= 0;
 while (v < int(f.size()) && f[v] == K()) ++v;
 if (v == int(f.size())) return std::vector<K>(n);  // 0 の平方根は 0
 if (v & 1) return {};                             // 位数が奇数なら無い
 const int half= v >> 1;
 if (half >= n) return std::vector<K>(n);
 const K c= f[v];
 K r0;
 if constexpr (is_zmod_v<K>) {
  const long long r= sqrt_mod(c.val(), K::mod());
  if (r < 0) return {};
  r0= K((unsigned long long)r);
 } else {
  assert(c == K(1));
  r0= K(1);
 }
 const K ci= c.inv();
 std::vector<K> u(n - half);  // u[0] = 1。f の n - half + v 次以上は 0 と見る
 for (int i= 0; i < int(u.size()) && v + i < int(f.size()); ++i) u[i]= f[v + i] * ci;
 std::vector<K> lg= log<Conv>(u, n - half);
 const K iv2= K(2).inv();
 for (auto& x: lg) x*= iv2;
 const std::vector<K> t= exp<Conv>(lg, n - half);
 std::vector<K> g(n);
 for (int i= 0; i < int(t.size()); ++i) g[half + i]= t[i] * r0;
 return g;
}
}  // namespace fps
