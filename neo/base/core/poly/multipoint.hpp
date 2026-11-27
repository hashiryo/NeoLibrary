#pragma once
// 多点評価と補間。どちらも部分積木 (葉が X - x_i) の上で、転置の形で計算する。
// 木の多項式は反転して持ち (P* = X^deg P(1/X)、定数項が 1)、内側の節は子の最高次の 1 を外して掛ける
// (親の反転は 1 + (u + w + u w) なので、積 u w は 2k - 1 項で、巡回の長さが子の 2 倍にならない)。
// 多点評価は、根で S = f mod P を作り、v = rev(S) / P* を 1 回だけ求めてから、各節のブロックを
//   左の子 = (v * P_R*)[deg P_R, deg P)、右の子 = (v * P_L*)[deg P_L, deg P)
// で子へ降りる。葉のブロックが f(x_i)。割り算は根の 1 回だけで、あとは範囲の積 (convolve_range) だけ。
// 補間は、P' を同じ降下で評価して w_i = y_i / P'(x_i) を作り、下から
//   G = G_L * P_R* + G_R * P_L*
// でまとめる (反転したまままとまるので、最後に 1 回反転する)。
// 点の無い葉は多項式 1 にして積で消す (X にすると x_i = 0 の点と根が重なる)。
// 核はテンプレート引数 Conv で受け、畳み込みの共通の口 (convolve.hpp) だけを使う。K は体で、点は相異なること。
// 小さい n では lagrange::interpolate (O(n^2)) の方が速いことがある。
#include <algorithm>
#include <cassert>
#include <vector>
#include "neo/base/core/convolve/convolve.hpp"
#include "neo/base/core/poly/fps.hpp"
#include "neo/base/core/poly/poly_divmod_fast.hpp"

namespace multipoint_internal {
template <class K> int deg(const std::vector<K>& p) { return p.empty() ? 0 : int(p.size()) - 1; }
// 部分積木 (反転した多項式)。tr[sz + i] が葉、tr[1] が根。点の無い葉は空 (多項式 1)。
template <class Conv, class K> std::vector<std::vector<K>> build_tree(const std::vector<K>& x, int& sz) {
 const int m= int(x.size());
 sz= 1;
 while (sz < m) sz<<= 1;
 std::vector<std::vector<K>> tr(2 * sz);
 for (int i= 0; i < m; ++i) tr[sz + i]= {K(1), -x[i]};
 std::vector<K> t(sz);  // 積の作業 (親の次数 - 1 項まで)
 for (int i= sz - 1; i >= 1; --i) {
  const int dl= deg(tr[2 * i]), dr= deg(tr[2 * i + 1]);
  if (dl == 0 || dr == 0) {  // 片方が 1 なら、もう片方と同じ
   if (dl + dr > 0) tr[i]= tr[dl > 0 ? 2 * i : 2 * i + 1];
   continue;
  }
  // 親の反転は 1 + (u + w + u w)。u w は dl + dr - 1 項で、巡回の長さ dl + dr に収まる。
  const K* u= tr[2 * i].data() + 1;
  const K* w= tr[2 * i + 1].data() + 1;
  const int d= dl + dr;
  tr[i].assign(d + 1, K());
  tr[i][0]= K(1);
  for (int j= 0; j < dl; ++j) tr[i][1 + j]+= u[j];
  for (int j= 0; j < dr; ++j) tr[i][1 + j]+= w[j];
  convolve_range<Conv>(u, dl, w, dr, 0, d - 1, t.data());
  for (int j= 0; j < d - 1; ++j) tr[i][2 + j]+= t[j];
 }
 return tr;
}
// v = rv / P* の先頭 d 項。rv は S の反転した係数で、S は次数 d 未満 (f mod P)。
template <class Conv, class K> std::vector<K> root_series(const std::vector<K>& pstar, const std::vector<K>& rv, int d) {
 const std::vector<K> invp= fps::inv<Conv>(pstar, d);
 std::vector<K> v(d);
 convolve_range<Conv>(rv.data(), d, invp.data(), d, 0, d, v.data());
 return v;
}
// 各節のブロックを子へ降ろす。節 o のブロックは v[start, start + deg(o)) にあり、子のブロックも同じ場所に入る。
// 点は前詰めなので、点の無い部分木のブロックは作らない。
template <class Conv, class K>
void descend(const std::vector<std::vector<K>>& tr, int sz, int m, std::vector<K>& v, std::vector<K>& t1, std::vector<K>& t2) {
 for (int cap= sz, d= 0; cap > 1; cap>>= 1, ++d)
  for (int o= 1 << d, e= 2 << d; o < e; ++o) {
   const int start= (o - (1 << d)) * cap;
   if (start >= m) continue;
   const int dl= deg(tr[2 * o]), dr= deg(tr[2 * o + 1]);
   if (dl == 0 || dr == 0) continue;  // 片方が 1 なら、ブロックはそのまま (dr == 0 だけが起きる)
   const int ds= dl + dr;
   // 左の子 = (v * P_R*)[dr, ds)、右の子 = (v * P_L*)[dl, ds)。
   convolve_range<Conv>(v.data() + start, ds, tr[2 * o + 1].data(), dr + 1, dr, ds, t1.data());
   convolve_range<Conv>(v.data() + start, ds, tr[2 * o].data(), dl + 1, dl, ds, t2.data());
   std::copy(t1.begin(), t1.begin() + dl, v.begin() + start);
   std::copy(t2.begin(), t2.begin() + dr, v.begin() + start + dl);
  }
}
// 下から G = G_L * P_R* + G_R * P_L* でまとめる。子のブロックは v[start, start + deg) にある。
template <class Conv, class K>
void combine(const std::vector<std::vector<K>>& tr, int sz, int m, std::vector<K>& v, std::vector<K>& t1, std::vector<K>& t2) {
 int d= 0;
 while ((2 << d) <= sz) ++d;
 for (--d; d >= 0; --d)
  for (int o= 1 << d, e= 2 << d; o < e; ++o) {
   const int cap= sz >> d, start= (o - (1 << d)) * cap;
   if (start >= m) continue;
   const int dl= deg(tr[2 * o]), dr= deg(tr[2 * o + 1]);
   if (dl == 0 || dr == 0) continue;
   const int ds= dl + dr;
   convolve_range<Conv>(v.data() + start, dl, tr[2 * o + 1].data(), dr + 1, 0, ds, t1.data());
   convolve_range<Conv>(v.data() + start + dl, dr, tr[2 * o].data(), dl + 1, 0, ds, t2.data());
   for (int j= 0; j < ds; ++j) v[start + j]= t1[j] + t2[j];
  }
}
}  // namespace multipoint_internal

// f の x[0..m) での値。
template <class Conv= convolve_internal::Auto, class K>
std::vector<K> multipoint_eval(const std::vector<K>& f, const std::vector<K>& x) {
 const int m= int(x.size());
 std::vector<K> res(m);
 if (m == 0 || f.empty()) return res;
 if (f.size() == 1) {
  std::fill(res.begin(), res.end(), f[0]);
  return res;
 }
 using namespace multipoint_internal;
 int sz;
 const auto tr= build_tree<Conv, K>(x, sz);
 std::vector<K> s= f;
 if (int(f.size()) > m) s= poly_divmod_fast<Conv>(f, std::vector<K>(tr[1].rbegin(), tr[1].rend())).second;
 std::vector<K> rv(m, K());
 for (int j= 0; j < int(s.size()); ++j) rv[m - 1 - j]= s[j];
 std::vector<K> v= root_series<Conv, K>(tr[1], rv, m);
 std::vector<K> t1(sz / 2), t2(sz / 2);
 descend<Conv, K>(tr, sz, m, v, t1, t2);
 std::copy(v.begin(), v.begin() + m, res.begin());
 return res;
}

// 点 (x_i, y_i) (i < n) を通る n - 1 次以下の多項式の係数。点は相異なること。
template <class Conv= convolve_internal::Auto, class K>
std::vector<K> interpolate(const std::vector<K>& x, const std::vector<K>& y) {
 const int n= int(x.size());
 assert(int(y.size()) == n);
 if (n == 0) return {};
 using namespace multipoint_internal;
 int sz;
 const auto tr= build_tree<Conv, K>(x, sz);
 // rev(P')[j] = (n - j) P*[j]。P' を多点評価して P'(x_i) を作る。
 std::vector<K> rv(n);
 for (int j= 0; j < n; ++j) rv[j]= tr[1][j] * K(n - j);
 std::vector<K> v= root_series<Conv, K>(tr[1], rv, n);
 std::vector<K> t1(sz), t2(sz);  // まとめる段では根の積が n 項になる
 descend<Conv, K>(tr, sz, n, v, t1, t2);
 for (int i= 0; i < n; ++i) v[i]= y[i] / v[i];
 combine<Conv, K>(tr, sz, n, v, t1, t2);
 std::vector<K> f(n);
 for (int j= 0; j < n; ++j) f[j]= v[n - 1 - j];
 return f;
}
