#pragma once
// 多項式の拡張 gcd (half-gcd)。f x + g y = gcd(f, g) の (d, x, y) を O(M(n) log n) で求める。
// 算法は half-gcd で、hgcd が 2x2 の多項式行列 (余りを半分の次数まで進める) を返し、cogcd がその行列を使って係数を返す。
// 今の Library の extgcd (mylib/fft/extgcd.hpp) と同じ組み立てで、商と余りを poly_divmod_fast、積を核の convolve にする。
// K は体で、多項式は末尾の 0 を取り除いた係数の vector (0 は空)。
#include <algorithm>
#include <array>
#include <cassert>
#include <utility>
#include <vector>
#include "neo/base/core/convolve/convolve.hpp"
#include "neo/base/core/poly/poly_divmod_fast.hpp"

namespace poly_extgcd_internal {
template <class K> int deg(const std::vector<K>& p) { return int(p.size()) - 1; }
template <class K> void trim(std::vector<K>& p) {
 while (!p.empty() && p.back() == K()) p.pop_back();
}
// 2x2 の多項式行列。v[0]、v[1] が 1 行目、v[2]、v[3] が 2 行目。
template <class K> struct Mat {
 std::vector<K> v[4];
};
template <class K> std::vector<K> suf(const std::vector<K>& p, int bg, int ed) {
 if (bg >= ed || bg >= int(p.size())) return {};
 ed= std::min(ed, int(p.size()));
 std::vector<K> r(p.begin() + bg, p.begin() + ed);
 trim(r);
 return r;
}
template <class Conv, class K> std::vector<K> mul(const std::vector<K>& a, const std::vector<K>& b) {
 std::vector<K> c= convolve<Conv>(a, b);
 trim(c);
 return c;
}
template <class K> std::vector<K> add(const std::vector<K>& a, const std::vector<K>& b) {
 std::vector<K> r(std::max(a.size(), b.size()));
 for (size_t i= 0; i < r.size(); ++i) r[i]= (i < a.size() ? a[i] : K()) + (i < b.size() ? b[i] : K());
 trim(r);
 return r;
}
template <class K> std::vector<K> sub(const std::vector<K>& a, const std::vector<K>& b) {
 std::vector<K> r(std::max(a.size(), b.size()));
 for (size_t i= 0; i < r.size(); ++i) r[i]= (i < a.size() ? a[i] : K()) - (i < b.size() ? b[i] : K());
 trim(r);
 return r;
}
template <class K> std::vector<K> neg(const std::vector<K>& a) {
 std::vector<K> r(a.size());
 for (size_t i= 0; i < a.size(); ++i) r[i]= -a[i];
 return r;
}
// a を b で割った商 (b は空でないこと)。
template <class Conv, class K> std::vector<K> quot(const std::vector<K>& a, const std::vector<K>& b) {
 return poly_divmod_fast<Conv>(a, b).first;
}

template <class Conv, class K> Mat<K> hgcd(const std::vector<K>& p0, const std::vector<K>& p1);

// hgcd の後半 (mylib の __COPY_HOGE)。p3 と r1 は新しい 2 つの余り、R はそこまでの行列。
template <class Conv, class K> Mat<K> hgcd_tail(Mat<K>& R, const std::vector<K>& p3, const std::vector<K>& r1, int m) {
 const int n= deg(p3), k= 2 * m - n, o= int(r1.size());
 if (o <= m) return R;
 Mat<K> S;
 if (o <= (((n + 1) * 3 + k) >> 2)) {
  auto q= quot<Conv>(suf(p3, k, n + 1), suf(r1, k, o));
  S.v[0]= R.v[2], S.v[1]= R.v[3];
  S.v[2]= sub(R.v[0], mul<Conv>(q, R.v[2]));
  S.v[3]= sub(R.v[1], mul<Conv>(q, R.v[3]));
  return S;
 }
 auto A= hgcd<Conv>(suf(p3, k, n + 1), suf(r1, k, o));
 S.v[0]= add(mul<Conv>(A.v[0], R.v[0]), mul<Conv>(A.v[1], R.v[2]));
 S.v[1]= add(mul<Conv>(A.v[0], R.v[1]), mul<Conv>(A.v[1], R.v[3]));
 S.v[2]= add(mul<Conv>(A.v[2], R.v[0]), mul<Conv>(A.v[3], R.v[2]));
 S.v[3]= add(mul<Conv>(A.v[2], R.v[1]), mul<Conv>(A.v[3], R.v[3]));
 return S;
}

// deg p0 > deg p1 のとき、M (p0, p1) の 2 行目が deg p0 / 2 未満になる行列 M を返す。
template <class Conv, class K> Mat<K> hgcd(const std::vector<K>& p0, const std::vector<K>& p1) {
 const int o= deg(p0), m= (o + 1) / 2, n= deg(p1);
 assert(o > n);
 if (n < (((o + 1) * 3 + m) >> 2)) {
  // 高次の半分で 1 回割り、その余りでもう 1 回割って、行列を作る。
  assert(n >= m);
  const auto tmp= quot<Conv>(suf(p0, m, o + 1), suf(p1, m, n + 1));
  const auto p3= sub(p0, mul<Conv>(tmp, p1));
  if (deg(p3) < 0) {
   Mat<K> R;
   R.v[0]= {K(1)}, R.v[1]= neg(tmp), R.v[2]= {}, R.v[3]= {K(1)};
   return R;
  }
  auto qr= poly_divmod_fast<Conv>(p1, p3);
  Mat<K> R;
  R.v[0]= {K(1)}, R.v[1]= neg(tmp), R.v[2]= neg(qr.first), R.v[3]= add({K(1)}, mul<Conv>(qr.first, tmp));
  return hgcd_tail<Conv, K>(R, p3, qr.second, m);
 }
 // 高次の半分で再帰し、その行列で 1 段進めてから、残りをまた半分ずつ進める。
 auto R= hgcd<Conv>(suf(p0, m, o + 1), suf(p1, m, n + 1));
 const auto p3= add(mul<Conv>(R.v[2], p0), mul<Conv>(R.v[3], p1));
 if (deg(p3) < 0) return R;
 auto qr= poly_divmod_fast<Conv>(add(mul<Conv>(R.v[0], p0), mul<Conv>(R.v[1], p1)), p3);
 Mat<K> S;
 S.v[0]= R.v[2], S.v[1]= R.v[3];
 S.v[2]= sub(R.v[0], mul<Conv>(qr.first, R.v[2]));
 S.v[3]= sub(R.v[1], mul<Conv>(qr.first, R.v[3]));
 return hgcd_tail<Conv, K>(S, p3, qr.second, m);
}

// deg p0 > deg p1 のとき、p0 x + p1 y = gcd(p0, p1) の (x, y) を返す。
template <class Conv, class K> std::pair<std::vector<K>, std::vector<K>> cogcd(const std::vector<K>& p0, const std::vector<K>& p1) {
 const int o= deg(p0), m= (o + 1) / 2, n= deg(p1), k= (o + m + 1) / 2;
 assert(o > n);
 if (n < k) {
  auto t= poly_divmod_fast<Conv>(p0, p1);
  if (deg(t.second) < 0) return {{}, {K(1)}};
  auto qr= poly_divmod_fast<Conv>(p1, t.second);
  if (deg(qr.second) < 0) return {{K(1)}, neg(t.first)};
  auto A= cogcd<Conv>(t.second, qr.second);
  auto x= sub(A.first, mul<Conv>(A.second, qr.first));
  return {x, sub(A.second, mul<Conv>(x, t.first))};
 }
 auto M= hgcd<Conv>(p0, p1);
 auto p3= add(mul<Conv>(M.v[2], p0), mul<Conv>(M.v[3], p1));
 if (deg(p3) < 0) return {M.v[0], M.v[1]};
 auto qr= poly_divmod_fast<Conv>(add(mul<Conv>(M.v[0], p0), mul<Conv>(M.v[1], p1)), p3);
 if (deg(qr.second) < 0) return {M.v[2], M.v[3]};
 auto A= cogcd<Conv>(p3, qr.second);
 auto x= sub(A.first, mul<Conv>(A.second, qr.first));
 return {add(mul<Conv>(A.second, M.v[0]), mul<Conv>(x, M.v[2])), add(mul<Conv>(A.second, M.v[1]), mul<Conv>(x, M.v[3]))};
}
}  // namespace poly_extgcd_internal

// f x + g y = gcd(f, g) を満たす (d, x, y) を求める。d は gcd を monic にしない形 (f と g の線形結合)。
// f と g は空でないこと。x と y は空 (0) のこともある。
template <class Conv= convolve_internal::Auto, class K>
std::vector<K> poly_extgcd(const std::vector<K>& f, const std::vector<K>& g, std::vector<K>& x, std::vector<K>& y) {
 using namespace poly_extgcd_internal;
 std::vector<K> a= f, b= g;
 trim(a), trim(b);
 assert(!a.empty() && !b.empty());
 if (deg(a) <= deg(b)) {
  auto qr= poly_divmod_fast<Conv>(a, b);
  if (deg(qr.second) < 0) {  // b が a を割る
   x= {}, y= {K(1)};
   return b;
  }
  auto A= cogcd<Conv>(b, qr.second);
  y= sub(A.first, mul<Conv>(A.second, qr.first));
  x= A.second;
 } else {
  auto A= cogcd<Conv>(a, b);
  x= A.first, y= A.second;
 }
 std::vector<K> d= add(mul<Conv>(a, x), mul<Conv>(b, y));
 return d;
}
