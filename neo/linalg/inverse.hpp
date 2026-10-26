#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>
#include "neo/algebra/ZMod.hpp"
#include "neo/linalg/Matrix.hpp"
// 正方行列の逆行列。正則でなければ std::nullopt。T は体で、ZMod なら法は素数、浮動小数点数なら絶対値が 1e-8 以下の要を 0 と
// みなす。
// ZMod は、[A | I] を前進消去すると要の行が [T | Y] (T は上三角で Y A = T) になるので、T X = Y を後退代入で解く。前進消去は、
// 行を u64 で持ち、要を 8 本ずつ貯めて残りの各行を 1 回だけ読み書きする (遅らせた更新)。積は M^2 < 2^60 なので、2^62 + 2^32
// 未満の値に 8 本足しても u64 に収まり、足したあとは 2^32 の位を 2^32 mod M に掛けて畳む。右半分の列は「行 r の単位行列の 1 の
// 位置」として付け替える。行を t 番目の要にするとき、その行の 1 を位置 t に移し、位置 t にあった行の 1 を空いた位置に移す (その
// 2 つの位置は、まだ持ち主の行にしか値がない)。すると要の行の右半分は位置 t 以下にしか値がなく、残りの行の更新は列
// [今の列, n + 要の数) で済む。前進消去の手間は 5n^3/6 から n^3/2 に減る。後退代入は下から 8 行ずつの組にして、組より下の行の
// 分を 8 行まとめて足す。n <= 16 の ZMod とほかの体は、素朴な掃き出し (浮動小数点数は絶対値が最大の要を選ぶ) にする。
// procon-judge の yosupo-inverse-matrix で書き比べた zmod_inv_relabel を写したもの。neo/linalg/adjugate.hpp もこの核を使う。
namespace inverse_internal {
using u32= unsigned;
using u64= unsigned long long;
constexpr int K= 8;
// row[j] += sum_{s<8} f[s] U[s][j] を j in [jb, w) で 4 個ずつ足してから畳む (neo/linalg/det.hpp の写し)。
[[gnu::always_inline]] inline void update_row(u64* __restrict row, const u64* const* U, const u32* f, int jb, int w, u64 R) {
 const __m256i Rv= _mm256_set1_epi64x(R), LO= _mm256_set1_epi64x(0xffffffffll);
 const __m256i f0= _mm256_set1_epi64x(f[0]), f1= _mm256_set1_epi64x(f[1]), f2= _mm256_set1_epi64x(f[2]), f3= _mm256_set1_epi64x(f[3]);
 const __m256i f4= _mm256_set1_epi64x(f[4]), f5= _mm256_set1_epi64x(f[5]), f6= _mm256_set1_epi64x(f[6]), f7= _mm256_set1_epi64x(f[7]);
 const u64 *U0= U[0], *U1= U[1], *U2= U[2], *U3= U[3], *U4= U[4], *U5= U[5], *U6= U[6], *U7= U[7];
 for(int j= jb; j < w; j+= 4) {
  __m256i acc= _mm256_loadu_si256((const __m256i*)(row + j));
  acc= _mm256_add_epi64(acc, _mm256_mul_epu32(_mm256_loadu_si256((const __m256i*)(U0 + j)), f0));
  acc= _mm256_add_epi64(acc, _mm256_mul_epu32(_mm256_loadu_si256((const __m256i*)(U1 + j)), f1));
  acc= _mm256_add_epi64(acc, _mm256_mul_epu32(_mm256_loadu_si256((const __m256i*)(U2 + j)), f2));
  acc= _mm256_add_epi64(acc, _mm256_mul_epu32(_mm256_loadu_si256((const __m256i*)(U3 + j)), f3));
  acc= _mm256_add_epi64(acc, _mm256_mul_epu32(_mm256_loadu_si256((const __m256i*)(U4 + j)), f4));
  acc= _mm256_add_epi64(acc, _mm256_mul_epu32(_mm256_loadu_si256((const __m256i*)(U5 + j)), f5));
  acc= _mm256_add_epi64(acc, _mm256_mul_epu32(_mm256_loadu_si256((const __m256i*)(U6 + j)), f6));
  acc= _mm256_add_epi64(acc, _mm256_mul_epu32(_mm256_loadu_si256((const __m256i*)(U7 + j)), f7));
  acc= _mm256_add_epi64(_mm256_mul_epu32(_mm256_srli_epi64(acc, 32), Rv), _mm256_and_si256(acc, LO));
  _mm256_storeu_si256((__m256i*)(row + j), acc);
 }
}
// acc の 8 行 (行の間隔 wx) に X の ns 行 (行の間隔 wx) を足す: acc[q][j] += sum_s c[s * 8 + q] X[s][j]
// (neo/linalg/solve.hpp の写し)。
inline void block_update(u64* __restrict acc, const u64* __restrict X, const u32* c, int ns, int wx, u64 R) {
 const __m256i Rv= _mm256_set1_epi64x(R), LO= _mm256_set1_epi64x(0xffffffffll);
 auto fold= [&](__m256i x) { return _mm256_add_epi64(_mm256_mul_epu32(_mm256_srli_epi64(x, 32), Rv), _mm256_and_si256(x, LO)); };
 for(int j= 0; j < wx; j+= 4) {
  u64* ap= acc + j;
  __m256i a0= _mm256_loadu_si256((const __m256i*)ap), a1= _mm256_loadu_si256((const __m256i*)(ap + wx));
  __m256i a2= _mm256_loadu_si256((const __m256i*)(ap + 2 * wx)), a3= _mm256_loadu_si256((const __m256i*)(ap + 3 * wx));
  __m256i a4= _mm256_loadu_si256((const __m256i*)(ap + 4 * wx)), a5= _mm256_loadu_si256((const __m256i*)(ap + 5 * wx));
  __m256i a6= _mm256_loadu_si256((const __m256i*)(ap + 6 * wx)), a7= _mm256_loadu_si256((const __m256i*)(ap + 7 * wx));
  const u64* xp= X + j;
  const u32* cp= c;
  for(int s= 0; s < ns;) {
   for(const int e= std::min(ns, s + 8); s < e; ++s, xp+= wx, cp+= 8) {
    const __m256i x= _mm256_loadu_si256((const __m256i*)xp);
    a0= _mm256_add_epi64(a0, _mm256_mul_epu32(x, _mm256_set1_epi32(cp[0])));
    a1= _mm256_add_epi64(a1, _mm256_mul_epu32(x, _mm256_set1_epi32(cp[1])));
    a2= _mm256_add_epi64(a2, _mm256_mul_epu32(x, _mm256_set1_epi32(cp[2])));
    a3= _mm256_add_epi64(a3, _mm256_mul_epu32(x, _mm256_set1_epi32(cp[3])));
    a4= _mm256_add_epi64(a4, _mm256_mul_epu32(x, _mm256_set1_epi32(cp[4])));
    a5= _mm256_add_epi64(a5, _mm256_mul_epu32(x, _mm256_set1_epi32(cp[5])));
    a6= _mm256_add_epi64(a6, _mm256_mul_epu32(x, _mm256_set1_epi32(cp[6])));
    a7= _mm256_add_epi64(a7, _mm256_mul_epu32(x, _mm256_set1_epi32(cp[7])));
   }
   a0= fold(a0), a1= fold(a1), a2= fold(a2), a3= fold(a3), a4= fold(a4), a5= fold(a5), a6= fold(a6), a7= fold(a7);
  }
  _mm256_storeu_si256((__m256i*)ap, a0), _mm256_storeu_si256((__m256i*)(ap + wx), a1);
  _mm256_storeu_si256((__m256i*)(ap + 2 * wx), a2), _mm256_storeu_si256((__m256i*)(ap + 3 * wx), a3);
  _mm256_storeu_si256((__m256i*)(ap + 4 * wx), a4), _mm256_storeu_si256((__m256i*)(ap + 5 * wx), a5);
  _mm256_storeu_si256((__m256i*)(ap + 6 * wx), a6), _mm256_storeu_si256((__m256i*)(ap + 7 * wx), a7);
 }
}
// [A | I] の前進消去 (右半分の付け替えつき)。要は列の順に探し、その列で値が 0 でない最初の行を選ぶ (まだ要にしていない行の
// 並び rem は、要を抜くたびに末尾の行で埋める)。t 番目の要の行の右半分は位置 t 以下にしか値がない。右半分の位置 q の列は、
// 単位行列の列 pos2row[q] を付け替えたもの。
template <class Z> struct Elim {
 int n, w;
 std::vector<u64> a, U;  // a は n × w (2^62 + 2^32 未満)。U は要の行を M 未満に直したもの (要にした順に w 個ずつ)
 std::vector<int> rem, prow, pcol, pos2row, row2pos;
 std::vector<u32> pval;
 explicit Elim(const Matrix<Z>& A): n(A.h), w((2 * A.h + 3) & ~3), a((size_t)n * w), rem(n), pos2row(n), row2pos(n) {
  for(int i= 0; i < n; ++i) {
   u64* r= row(i);
   for(int j= 0; j < n; ++j) r[j]= A[i][j].val();
   r[n + i]= 1;
  }
  std::iota(rem.begin(), rem.end(), 0), std::iota(pos2row.begin(), pos2row.end(), 0), std::iota(row2pos.begin(), row2pos.end(), 0);
 }
 u64* row(int r) { return a.data() + (size_t)r * w; }
 const u64* urow(int t) const { return U.data() + (size_t)t * w; }
 // 要のない列は飛ばし、その数が skip を超えたら false を返して止める。最後の組の更新も残りの行に足す。
 bool run(int skip) {
  const u32 M= Z::mod();
  const u64 R= (u64(1) << 32) % M;
  U.assign((size_t)n * w, 0);
  std::vector<u32> F((size_t)n * K), v(n);
  std::vector<u64> zero(w);
  const u64* Ub[K];
  int skipped= 0;
  for(int col= 0; col < n;) {
   int bt= 0;
   for(; bt < K && col < n; ++col) {
    const int nr= rem.size();
    int pi= -1;
    for(int i= 0; i < nr; ++i) {
     const int r= rem[i];
     const u32* fr= &F[(size_t)r * K];
     u64 x= a[(size_t)r * w + col];
     for(int s= 0; s < bt; ++s) x+= u64(fr[s]) * Ub[s][col];
     v[i]= Z(x).val();
     if(pi < 0 && v[i]) pi= i;
    }
    if(pi < 0) {
     if(++skipped > skip) return false;
     continue;
    }
    const int p= rem[pi], t= prow.size(), q= row2pos[p];
    if(q != t) {  // 行 p の 1 を位置 t に、位置 t にあった行 r2 の 1 を位置 q に移す
     const int r2= pos2row[t];
     u64 *ap= row(p), *ar= row(r2);
     ap[n + q]= 0, ap[n + t]= 1, ar[n + t]= 0, ar[n + q]= 1;
     pos2row[q]= r2, row2pos[r2]= q, pos2row[t]= p, row2pos[p]= t;
    }
    u64* up= &U[(size_t)t * w];
    const u32* fp= &F[(size_t)p * K];
    const u64* ap= row(p);
    const int je= std::min(w, (n + t + 1 + 3) & ~3);
    for(int j= col & ~3; j < je; ++j) {
     u64 x= ap[j];
     for(int s= 0; s < bt; ++s) x+= u64(fp[s]) * Ub[s][j];
     up[j]= Z(x).val();
    }
    const Z iv= Z::raw(v[pi]).inv();
    for(int i= 0; i < nr; ++i) F[(size_t)rem[i] * K + bt]= (-(Z::raw(v[i]) * iv)).val();
    prow.push_back(p), pcol.push_back(col), pval.push_back(v[pi]);
    Ub[bt++]= up;
    rem[pi]= rem.back(), rem.pop_back();
   }
   if(!bt) continue;
   for(int s= bt; s < K; ++s) Ub[s]= zero.data();
   const int jb= col & ~3, je= std::min(w, (n + (int)prow.size() + 3) & ~3);
   for(int r: rem) {
    u32* fr= &F[(size_t)r * K];
    for(int s= bt; s < K; ++s) fr[s]= 0;
    u32 any= 0;
    for(int s= 0; s < K; ++s) any|= fr[s];
    if(any) update_row(row(r), Ub, fr, jb, je, R);
   }
  }
  return true;
 }
 // 要にした行の並び (t 番目の要が行 prow[t]) の置換が奇置換か。要が n 本のときだけ使う。
 bool odd() const {
  std::vector<char> seen(n);
  int cyc= 0;
  for(int i= 0; i < n; ++i)
   if(!seen[i]) {
    ++cyc;
    for(int j= i; !seen[j]; j= prow[j]) seen[j]= 1;
   }
  return (n - cyc) & 1;
 }
};
// 前進消去のあと、要の列が作る上三角の行列 T (T[t][s] = U_t[pcol[s]]) で T X = Y を解く。Y の行 t は U_t の列 cols。
// X は要の数 × wx (wx は cols の数を 4 の倍数に切り上げたもの) で、M 未満の値を u64 で持つ。下から 8 行ずつの組にして、
// 組より下の行の分は block_update で 8 行まとめて足し、組の中は下の行から順に update_row で足す (neo/linalg/solve.hpp の写し)。
template <class Z> std::vector<u64> tri_solve(const Elim<Z>& e, const std::vector<int>& cols, int& wx) {
 const u32 M= Z::mod();
 const u64 R= (u64(1) << 32) % M;
 const int r= e.prow.size(), nc= cols.size();
 wx= (nc + 3) & ~3;
 std::vector<u64> X((size_t)r * wx), acc((size_t)K * wx), zero(wx);
 std::vector<u32> c((size_t)r * K);
 for(int t1= r; t1 > 0; t1-= K) {
  const int t0= std::max(0, t1 - K), bs= t1 - t0, ns= r - t1;
  for(int q= 0; q < K; ++q) {
   u64* aq= &acc[(size_t)q * wx];
   std::fill(aq, aq + wx, 0);
   if(q < bs) {
    const u64* ut= e.urow(t0 + q);
    for(int k= 0; k < nc; ++k) aq[k]= ut[cols[k]];
   }
  }
  for(int s= t1; s < r; ++s)
   for(int q= 0; q < K; ++q) {
    const u32 x= q < bs ? u32(e.urow(t0 + q)[e.pcol[s]]) : 0;
    c[(size_t)(s - t1) * K + q]= x ? M - x : 0;
   }
  if(ns) block_update(acc.data(), &X[(size_t)t1 * wx], c.data(), ns, wx, R);
  for(int q= bs - 1; q >= 0; --q) {
   const int t= t0 + q;
   u64* aq= &acc[(size_t)q * wx];
   const u64* Ub[K];
   u32 f[K];
   for(int i= 0; i < K; ++i) {
    const int s= t + 1 + i;
    if(s < t1) {
     const u32 x= u32(e.urow(t)[e.pcol[s]]);
     Ub[i]= &X[(size_t)s * wx], f[i]= x ? M - x : 0;
    } else Ub[i]= zero.data(), f[i]= 0;
   }
   update_row(aq, Ub, f, 0, wx, R);
   const Z iv= Z::raw(e.pval[t]).inv();
   u64* xt= &X[(size_t)t * wx];
   for(int k= 0; k < wx; ++k) xt[k]= (Z(aq[k]) * iv).val();
  }
 }
 return X;
}
template <class Z> std::optional<Matrix<Z>> zmod(const Matrix<Z>& A) {
 const int n= A.h;
 Elim<Z> e(A);
 if(!e.run(0)) return std::nullopt;
 std::vector<int> cols(n);
 std::iota(cols.begin(), cols.end(), n);
 int wx;
 const auto X= tri_solve(e, cols, wx);
 Matrix<Z> res(n, n);
 for(int t= 0; t < n; ++t) {
  Z* rt= res[t];
  const u64* xt= &X[(size_t)t * wx];
  for(int q= 0; q < n; ++q) rt[e.pos2row[q]]= Z::raw(u32(xt[q]));
 }
 return res;
}
template <class T> bool is_zero(const T& x) {
 if constexpr(std::is_floating_point_v<T>) return std::abs(x) <= T(1e-8);
 else return x == linalg_internal::zero<T>();
}
// 体の上の素朴な前進消去。[A | I] を持ち、要にした行を残りの行から引く。要のない列は飛ばし、その数が skip を超えたら止める。
template <class T> struct Field {
 int n;
 Matrix<T> a;  // n × 2n
 std::vector<int> rem, prow, pcol;
 explicit Field(const Matrix<T>& A): n(A.h), a(n, 2 * n), rem(n) {
  for(int i= 0; i < n; ++i) std::copy(A[i], A[i] + n, a[i]), a[i][n + i]= linalg_internal::one<T>();
  std::iota(rem.begin(), rem.end(), 0);
 }
 bool run(int skip) {
  const T zero= linalg_internal::zero<T>(), one= linalg_internal::one<T>();
  int skipped= 0;
  for(int c= 0; c < n; ++c) {
   int pi= -1;
   for(int i= 0; i < (int)rem.size(); ++i) {
    const T x= a[rem[i]][c];
    if constexpr(std::is_floating_point_v<T>) {
     if(pi < 0 || std::abs(x) > std::abs(a[rem[pi]][c])) pi= i;
    } else if(x != zero) {
     pi= i;
     break;
    }
   }
   if(pi < 0 || is_zero(a[rem[pi]][c])) {
    if(++skipped > skip) return false;
    continue;
   }
   const int p= rem[pi];
   rem.erase(rem.begin() + pi);
   const T* ap= a[p];
   const T iv= one / ap[c];
   for(int r: rem) {
    T* ar= a[r];
    const T f= ar[c] * iv;
    if(f == zero) continue;
    for(int j= c; j < 2 * n; ++j) ar[j]-= f * ap[j];
   }
   prow.push_back(p), pcol.push_back(c);
  }
  return true;
 }
 // T X = Y を解く。T[t][s] = a[prow[t]][pcol[s]]、Y の行 t は a[prow[t]] の列 cols。
 Matrix<T> back(const std::vector<int>& cols) const {
  const int r= prow.size(), nc= cols.size();
  Matrix<T> X(r, nc);
  for(int t= r - 1; t >= 0; --t) {
   const T* at= a[prow[t]];
   T* xt= X[t];
   for(int k= 0; k < nc; ++k) xt[k]= at[cols[k]];
   for(int s= t + 1; s < r; ++s) {
    const T f= at[pcol[s]];
    const T* xs= X[s];
    for(int k= 0; k < nc; ++k) xt[k]-= f * xs[k];
   }
   const T iv= linalg_internal::one<T>() / at[pcol[t]];
   for(int k= 0; k < nc; ++k) xt[k]*= iv;
  }
  return X;
 }
};
template <class T> std::optional<Matrix<T>> field(const Matrix<T>& A) {
 const int n= A.h;
 Field<T> f(A);
 if(!f.run(0)) return std::nullopt;
 std::vector<int> cols(n);
 std::iota(cols.begin(), cols.end(), n);
 return f.back(cols);
}
}
template <class T> std::optional<Matrix<T>> inverse(Matrix<T> A) {
 assert(A.h == A.w);
 if constexpr(linalg_internal::is_zmod_v<T>) {
  if(A.h > 16) return inverse_internal::zmod(A);
 }
 return inverse_internal::field(A);
}
