#pragma once
// 多項式係数の線形漸化式 (P-recursive) Σ_{i=0}^{d} p_i(n) a_{n+i} = 0 (n >= 0) の第 N 項 a_N。
// p[i][j] は p_i の n^j の係数で、a0 は a_0..a_{d-1}。p_d(n) != 0 (0 <= n <= N - d) を仮定する。
//
// 遷移行列は、行 r < d-1 が p_d(n) e_{r+1}、最後の行が -p_i(n) と取る。y_n = (Π_{i<n} p_d(i)) v_n と置くと
//   y_{n+1} = B(n) y_n,  a_N = (e_{d-1}^T B(N-d) ... B(0) v_0) / (Π_{n=0}^{N-d} p_d(n))
// になる。ずらしの行にも p_d を掛けるのが肝で、割り算は最後の 1 回だけになる。
//
// 算法は Bostan–Gaudry–Schost の Theorem 8 (Chudnovsky の baby-step/giant-step の値ベースの形)。
// ブロック多項式 M_t(X) = B(X) B(X+1) ... B(X+t-1) を係数では作らず、等差数列 X = jk 上の値 M_t(jk) を持ち、
//   M_{2t}(jk) = M_t(jk + t) M_t(jk)
// で倍々にする。X で t ずらすのは添字で t/k ずらすことになり、標本点のずらし (sample_points_shift) 1 回で済む。
// 最後に M_k(0), M_k(k), ..., M_k((q-1)k) を順に当てる。手間は √(e N) の準線形で、メモリは O(d^2 √(e N))。
// 定数係数 (次数が全部 0) は遷移行列が一定なので、行列の冪で進める。
// ZMod では遷移行列が法を周期に持つ (B(n + 法) = B(n)) ので、法の 4 倍より長い範囲は 1 周期分の積の冪
//   Π_{i<len} B(i) = (Π_{i<r} B(i)) (Π_{i<法} B(i))^q  (len = q·法 + r)
// に落とす (N が 10^18 でも √(法) の手間で済む)。
// 核はテンプレート引数 Conv で受け、畳み込みの共通の口 (convolve.hpp) と標本点のずらしに渡す。
// K は体で、法 (標数) は 2 e k + 1 より大きいこと (標本点のずらしが 1/j! を使う)。
#include <algorithm>
#include <cassert>
#include <vector>
#include "neo/base/core/algebra/ZMod.hpp"
#include "neo/base/core/convolve/convolve.hpp"
#include "neo/base/core/poly/sample_shift.hpp"

namespace p_recursive_internal {
template <class K> using Poly= std::vector<K>;
// d x d の多項式行列 (行が先)。
template <class K> struct PMat {
 int d;
 std::vector<Poly<K>> e;
 PMat(int d_= 0): d(d_), e(size_t(d_) * d_) {}
 Poly<K>& operator()(int i, int j) { return e[size_t(i) * d + j]; }
 const Poly<K>& operator()(int i, int j) const { return e[size_t(i) * d + j]; }
};
template <class K> K eval(const Poly<K>& f, K x) {
 K r= K();
 for (int i= int(f.size()) - 1; i >= 0; --i) r= r * x + f[i];
 return r;
}
// 長さ L の積の値を、等差数列 X = jk 上に作る。k は 2 冪で e k + 1 >= L / k。vals[e][j] は成分 e の M_L(jk)、
// vals[d*d][j] は Π_{i<k} p_d(jk+i)。j = 0..e k。
template <class Conv, class K>
void build_values(const PMat<K>& base, const std::vector<K>& pdl, int d, int edeg, unsigned long long L,
                  unsigned long long& k, std::vector<std::vector<K>>& val) {
 k= 1;
 while (k < L && (unsigned long long)edeg * k * k < L) k<<= 1;
 const int nent= d * d + 1;
 val.assign(nent, {});
 for (int e= 0; e < d * d; ++e) {
  val[e].resize(edeg + 1);
  for (int j= 0; j <= edeg; ++j) val[e][j]= eval(base.e[e], K((unsigned long long)j * k));
 }
 val[d * d].resize(edeg + 1);
 for (int j= 0; j <= edeg; ++j) val[d * d][j]= eval(pdl, K((unsigned long long)j * k));
 const K kinv= K(k).inv();
 std::vector<std::vector<K>> bs(nent), tl(nent), t2(nent), nv(nent);
 std::vector<K> A(d * d), Bv(d * d), C(d * d);
 for (unsigned long long t= 1; t < k; t<<= 1) {
  const int cnt= int(edeg * t) + 1;
  for (int e= 0; e < nent; ++e) {
   bs[e]= sample_points_shift<Conv>(val[e], K(t) * kinv, cnt);  // F(t/k + j)
   tl[e]= sample_points_shift<Conv>(val[e], K(cnt), cnt);       // F(cnt + j)
   t2[e]= sample_points_shift<Conv>(bs[e], K(cnt), cnt);        // F(t/k + cnt + j)
  }
  for (int e= 0; e < nent; ++e) nv[e].resize(2 * cnt - 1);
  for (int j= 0; j < 2 * cnt - 1; ++j) {
   for (int i= 0; i < d; ++i)
    for (int kk= 0; kk < d; ++kk) {
     A[i * d + kk]= j < cnt ? bs[i * d + kk][j] : t2[i * d + kk][j - cnt];
     Bv[i * d + kk]= j < cnt ? val[i * d + kk][j] : tl[i * d + kk][j - cnt];
    }
   for (int i= 0; i < d; ++i)
    for (int kk= 0; kk < d; ++kk) {
     K s= K();
     for (int l= 0; l < d; ++l) s+= A[i * d + l] * Bv[l * d + kk];
     C[i * d + kk]= s;
    }
   for (int i= 0; i < d * d; ++i) nv[i][j]= C[i];
   const K fd= j < cnt ? val[d * d][j] : tl[d * d][j - cnt];
   const K ft= j < cnt ? bs[d * d][j] : t2[d * d][j - cnt];
   nv[d * d][j]= fd * ft;
  }
  val.swap(nv);
 }
}
// 長さ L の積を v に掛け、D に Π_{i<L} p_d(i) を掛ける。
template <class Conv, class K>
void apply_range(const PMat<K>& base, const std::vector<K>& pdl, int d, int edeg, unsigned long long L, std::vector<K>& v, K& D) {
 if (L == 0) return;
 unsigned long long k;
 std::vector<std::vector<K>> val;
 build_values<Conv>(base, pdl, d, edeg, L, k, val);
 const unsigned long long q= L / k, r= L % k;
 std::vector<K> w(d);
 for (unsigned long long j= 0; j < q; ++j) {
  D*= val[d * d][j];
  for (int i= 0; i < d; ++i) {
   K s= K();
   for (int kk= 0; kk < d; ++kk) s+= val[i * d + kk][j] * v[kk];
   w[i]= s;
  }
  v.swap(w);
 }
 for (unsigned long long i= 0; i < r; ++i) {
  const unsigned long long n= q * k + i;
  const K nk= K(n);
  K s= K();
  for (int kk= 0; kk < d; ++kk) s+= eval(base(d - 1, kk), nk) * v[kk];
  for (int kk= 0; kk + 1 < d; ++kk) w[kk]= eval(base(kk, kk + 1), nk) * v[kk + 1];
  w[d - 1]= s;
  v.swap(w);
  D*= eval(pdl, nk);
 }
}
// 1 周期 (長さ pr) の積を d x d の行列 P として返し、Dp = Π_{i<pr} p_d(i) も返す。
template <class Conv, class K>
void period_matrix(const PMat<K>& base, const std::vector<K>& pdl, int d, int edeg, unsigned long long pr, std::vector<K>& P, K& Dp) {
 unsigned long long k;
 std::vector<std::vector<K>> val;
 build_values<Conv>(base, pdl, d, edeg, pr, k, val);
 const unsigned long long q= pr / k, r= pr % k;
 P.assign(d * d, K());
 for (int i= 0; i < d; ++i) P[i * d + i]= K(1);
 Dp= K(1);
 std::vector<K> t(d * d);
 for (unsigned long long j= 0; j < q; ++j) {
  Dp*= val[d * d][j];
  for (int i= 0; i < d; ++i)
   for (int kk= 0; kk < d; ++kk) {
    K s= K();
    for (int l= 0; l < d; ++l) s+= val[i * d + l][j] * P[l * d + kk];
    t[i * d + kk]= s;
   }
  P.swap(t);
 }
 for (unsigned long long i= 0; i < r; ++i) {
  const unsigned long long n= q * k + i;
  const K nk= K(n);
  for (int ii= 0; ii < d; ++ii)
   for (int kk= 0; kk < d; ++kk) {
    K s= K();
    if (ii + 1 == d) {
     for (int l= 0; l < d; ++l) s+= eval(base(d - 1, l), nk) * P[l * d + kk];
    } else {
     s= eval(base(ii, ii + 1), nk) * P[(ii + 1) * d + kk];  // ずらしの行にも p_d が掛かる
    }
    t[ii * d + kk]= s;
   }
  P.swap(t);
  Dp*= eval(pdl, nk);
 }
}
}  // namespace p_recursive_internal

namespace p_recursive {
template <class Conv= convolve_internal::Auto, class K>
K nth(const std::vector<std::vector<K>>& p, const std::vector<K>& a0, unsigned long long N) {
 using namespace p_recursive_internal;
 const int d= int(p.size()) - 1;
 assert(d >= 1 && int(a0.size()) >= d);
 if (N < (unsigned long long)d) return a0[N];
 const unsigned long long len= N - (unsigned long long)d + 1;
 PMat<K> base(d);
 for (int r= 0; r + 1 < d; ++r) base(r, r + 1)= p[d];
 for (int i= 0; i < d; ++i) {
  base(d - 1, i)= p[i];
  for (auto& c: base(d - 1, i)) c= -c;
 }
 int edeg= int(p[d].size()) - 1;
 for (int i= 0; i < d; ++i) edeg= std::max(edeg, int(p[i].size()) - 1);
 if (edeg == 0) {  // 定数係数: 遷移行列が一定なので、行列の冪で進める
  std::vector<K> mx(d * d), my(d * d, K()), t(d * d);
  for (int i= 0; i < d; ++i)
   for (int j= 0; j < d; ++j) mx[i * d + j]= eval(base(i, j), K());
  for (int i= 0; i < d; ++i) my[i * d + i]= K(1);
  for (unsigned long long e= len; e; e>>= 1) {
   if (e & 1) {
    for (int i= 0; i < d; ++i)
     for (int j= 0; j < d; ++j) {
      K s= K();
      for (int kk= 0; kk < d; ++kk) s+= my[i * d + kk] * mx[kk * d + j];
      t[i * d + j]= s;
     }
    my.swap(t);
   }
   for (int i= 0; i < d; ++i)
    for (int j= 0; j < d; ++j) {
     K s= K();
     for (int kk= 0; kk < d; ++kk) s+= mx[i * d + kk] * mx[kk * d + j];
     t[i * d + j]= s;
    }
   mx.swap(t);
  }
  K D(1), pwbase= p[d][0];
  for (unsigned long long x= len; x; x>>= 1, pwbase*= pwbase)
   if (x & 1) D*= pwbase;
  K r= K();
  for (int i= 0; i < d; ++i) r+= my[(d - 1) * d + i] * a0[i];
  return r / D;
 }
 unsigned long long pr= 0;
 if constexpr (is_zmod_v<K>) {
  const unsigned long long mm= K::mod();
  if (mm != 0 && len > 4 * mm) pr= mm;
 }
 std::vector<K> v(a0.begin(), a0.begin() + d);
 K D(1);
 if (pr != 0) {
  std::vector<K> P;
  K Dp;
  period_matrix<Conv>(base, p[d], d, edeg, pr, P, Dp);
  const unsigned long long q= len / pr, r= len % pr;
  std::vector<K> pw(d * d, K()), bs= P, t(d * d);
  for (int i= 0; i < d; ++i) pw[i * d + i]= K(1);
  for (unsigned long long e= q; e; e>>= 1) {
   if (e & 1) {
    for (int i= 0; i < d; ++i)
     for (int j= 0; j < d; ++j) {
      K s= K();
      for (int kk= 0; kk < d; ++kk) s+= pw[i * d + kk] * bs[kk * d + j];
      t[i * d + j]= s;
     }
    pw.swap(t);
   }
   for (int i= 0; i < d; ++i)
    for (int j= 0; j < d; ++j) {
     K s= K();
     for (int kk= 0; kk < d; ++kk) s+= bs[i * d + kk] * bs[kk * d + j];
     t[i * d + j]= s;
    }
   bs.swap(t);
  }
  K dq(1), bp= Dp;
  for (unsigned long long e= q; e; e>>= 1, bp*= bp)
   if (e & 1) dq*= bp;
  std::vector<K> w(d);
  for (int i= 0; i < d; ++i) {
   K s= K();
   for (int kk= 0; kk < d; ++kk) s+= pw[i * d + kk] * v[kk];
   w[i]= s;
  }
  v.swap(w);
  D= dq;
  apply_range<Conv>(base, p[d], d, edeg, r, v, D);
 } else {
  apply_range<Conv>(base, p[d], d, edeg, len, v, D);
 }
 return v[d - 1] / D;
}
}  // namespace p_recursive
