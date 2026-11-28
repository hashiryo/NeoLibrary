#pragma once
// 多変数の畳み込み (truncated)。dim で各次元の大きさを受け取り、n = ∏ dim 個の配列 2 本の、
// 各次元で dim_i 以上の次数を落とした積 (結果も n 個) を求める。
// ランクの変数 chi[i] = Σ_j (i / (dim_0 ... dim_{j-1})) mod k を使う。chi[i] + chi[j] ≡ chi[i+j] (mod k) で、
// 桁上がりは chi が k だけ進むことに吸収されるので、クラスごとに長さ 2n 以上の変換をして、
// (i, j) の組のクラス i + j (mod k) に足せば、truncated の積が戻る。
// 変換は単独の NTT の核 (ntt.hpp の Ntt、または dft / idft を持つ同じ形の核) の dft と idft を使う。
// K は法がコンパイル時に決まっている ZMod (ZMod<0> は使えない)。
#include <algorithm>
#include <cassert>
#include <numeric>
#include <vector>
#include "neo/base/core/algebra/ZMod.hpp"
#include "neo/base/core/convolve/ntt.hpp"

template <class Ntt, class K> std::vector<K> multivariate_convolve(const std::vector<int>& dim, const std::vector<K>& f, const std::vector<K>& g) {
 static_assert(ntt_internal::static_ntt_mod<K>::value, "多変数の畳み込みは、法がコンパイル時に決まった ntt_friendly な ZMod だけ");
 constexpr unsigned P= ntt_internal::static_ntt_mod<K>::mod;
 const int k= int(dim.size());
 int n= 1;
 for (int d: dim) n*= d;
 assert(int(f.size()) == n && int(g.size()) == n);
 if (k == 0) return {f[0] * g[0]};
 // chi[i] = Σ_j (i / (dim_0 ... dim_{j-1})) mod k
 std::vector<int> chi(n, 0);
 for (int i= n; i--;)
  for (int den= 1, j= 0; j < k; ++j) chi[i]+= i / (den*= dim[j]);
 for (int i= n; i--;) chi[i]%= k;
 int m= 1;
 while (m < 2 * n) m<<= 1;
 std::vector<K> a(size_t(k) * m), b(size_t(k) * m), c(size_t(k) * m, K());
 for (int i= 0; i < n; ++i) a[size_t(chi[i]) * m + i]= f[i], b[size_t(chi[i]) * m + i]= g[i];
 for (int t= 0; t < k; ++t) Ntt::template dft<P>(a.data() + size_t(t) * m, m);
 for (int t= 0; t < k; ++t) Ntt::template dft<P>(b.data() + size_t(t) * m, m);
 for (int i= 0; i < k; ++i)
  for (int j= 0; j < k; ++j) {
   int l= i + j;
   if (l >= k) l-= k;
   K* cl= c.data() + size_t(l) * m;
   const K* ai= a.data() + size_t(i) * m;
   const K* bj= b.data() + size_t(j) * m;
   for (int t= 0; t < m; ++t) cl[t]+= ai[t] * bj[t];
  }
 for (int t= 0; t < k; ++t) Ntt::template idft<P>(c.data() + size_t(t) * m, m);
 std::vector<K> r(n);
 for (int i= 0; i < n; ++i) r[i]= c[size_t(chi[i]) * m + i];
 return r;
}
