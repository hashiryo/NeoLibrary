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
// 法が奇数なら、前進消去の要を探す列、要の行、行ごとの掛ける数を AVX2 で 4 個ずつ求める (Elim::run_v、neo/linalg/solve.hpp の
// run_v と同じ形で、単位行列の列の付け替えは右半分の列だけを触る)。要の行と要の値は本当の値の 2^{-32} 倍になる (scaled) が、
// T の各行と Y の行が同じ倍率なので、T X = Y の解は変わらない。行列と要の行は 2 MB 境界の mmap に置いて huge page を頼む (hvec)。
// procon-judge の yosupo-inverse-matrix で書き比べた zmod_inv_relabel を写し、前進消去を zmod_inv_hpv の形にしたもの。
// neo/linalg/adjugate.hpp もこの核を使う。
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
// 4 個の u64 (2^63 未満) を、2^32 の位を R = 2^32 mod M に掛けて畳んでから Montgomery の還元で 2^{-32} 倍し、M 未満に直す。
struct Mont {
 __m256i R, LO, NIV, M, M32;
 explicit Mont(u32 m) {
  u32 x= m;
  for(int i= 0; i < 5; ++i) x*= 2 - m * x;
  R= _mm256_set1_epi64x((u64(1) << 32) % m), LO= _mm256_set1_epi64x(0xffffffffll), NIV= _mm256_set1_epi64x(u32(-x));
  M= _mm256_set1_epi64x(m), M32= _mm256_set1_epi32(m);
 }
 [[gnu::always_inline]] __m256i redc(__m256i x) const {
  x= _mm256_add_epi64(_mm256_mul_epu32(_mm256_srli_epi64(x, 32), R), _mm256_and_si256(x, LO));
  const __m256i t= _mm256_srli_epi64(_mm256_add_epi64(x, _mm256_mul_epu32(_mm256_mul_epu32(x, NIV), M)), 32);
  return _mm256_min_epu32(t, _mm256_sub_epi32(t, M32));  // 上の 32 bit は 0 のまま (0 - M は大きいので min で 0)
 }
};
[[gnu::always_inline]] inline __m128i pack_lo(__m256i x) { return _mm256_castsi256_si128(_mm256_permutevar8x32_epi32(x, _mm256_setr_epi32(0, 2, 4, 6, 0, 0, 0, 0))); }
// [A | I] の前進消去 (右半分の付け替えつき)。要は列の順に探し、その列で値が 0 でない最初の行を選ぶ (まだ要にしていない行の
// 並び rem は、要を抜くたびに末尾の行で埋める)。t 番目の要の行の右半分は位置 t 以下にしか値がない。右半分の位置 q の列は、
// 単位行列の列 pos2row[q] を付け替えたもの。a は後ろに 8 個の余りを持つ (run_v が行の後ろを 8 個まで読む)。
template <class Z> struct Elim {
 int n, w;
 linalg_internal::hvec<u64> a, U;  // a は n × w (2^62 + 2^32 未満)。U は要の行を M 未満に直したもの (要にした順に w 個ずつ)
 std::vector<int> rem, prow, pcol, pos2row, row2pos;
 std::vector<u32> pval;
 bool scaled= false;  // run_v で消した (要の行と要の値が 2^{-32} 倍)
 explicit Elim(const Matrix<Z>& A): n(A.h), w((2 * A.h + 3) & ~3), a((size_t)n * w + 8), rem(n), pos2row(n), row2pos(n) {
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
 // run と同じ前進消去を、要を探す列、要の行、行ごとの掛ける数を 4 個ずつ求めて行う。法は奇数。選ぶ要は run と同じ。
 // 組の初めに残りの行の次の 8 列を列ごとの配列 C (残りの行の並びの順) に 4 × 4 の転置で写し、要を探す列はそこから読む
 // (付け替えは右半分の列だけを触るので、C に写した左半分の列は変わらない)。行ごとの掛ける数 F も残りの行の並びの順に列ごとに
 // 持ち、要の行を抜くときは末尾の行の分をその位置へ移す。値は Montgomery の還元で M 未満に直し、要の行 U と要の値 pval は
 // 本当の値の 2^{-32} 倍、F は 2^32 倍で持つ。
 bool run_v(int skip) {
  const u32 M= Z::mod();
  const u64 R= (u64(1) << 32) % M;
  const Mont mt(M);
  const int n4= (n + 3) & ~3;
  scaled= true;
  U.assign((size_t)n * w, 0);
  std::vector<u64> C((size_t)K * n4), zero(w + 8);
  std::vector<u32> F((size_t)K * n4), V(n4);
  const u64* Ub[K];
  auto fill_cache= [&](int c0) {
   const int nr= rem.size();
   for(int i= 0; i < nr; i+= 4) {
    const u64* r0= row(rem[i]) + c0;
    const u64* r1= i + 1 < nr ? row(rem[i + 1]) + c0 : zero.data();
    const u64* r2= i + 2 < nr ? row(rem[i + 2]) + c0 : zero.data();
    const u64* r3= i + 3 < nr ? row(rem[i + 3]) + c0 : zero.data();
    for(int h= 0; h < 8; h+= 4) {
     const __m256i x0= _mm256_loadu_si256((const __m256i*)(r0 + h)), x1= _mm256_loadu_si256((const __m256i*)(r1 + h));
     const __m256i x2= _mm256_loadu_si256((const __m256i*)(r2 + h)), x3= _mm256_loadu_si256((const __m256i*)(r3 + h));
     const __m256i t0= _mm256_unpacklo_epi64(x0, x1), t1= _mm256_unpackhi_epi64(x0, x1);
     const __m256i t2= _mm256_unpacklo_epi64(x2, x3), t3= _mm256_unpackhi_epi64(x2, x3);
     _mm256_storeu_si256((__m256i*)&C[(size_t)(h + 0) * n4 + i], _mm256_permute2x128_si256(t0, t2, 0x20));
     _mm256_storeu_si256((__m256i*)&C[(size_t)(h + 1) * n4 + i], _mm256_permute2x128_si256(t1, t3, 0x20));
     _mm256_storeu_si256((__m256i*)&C[(size_t)(h + 2) * n4 + i], _mm256_permute2x128_si256(t0, t2, 0x31));
     _mm256_storeu_si256((__m256i*)&C[(size_t)(h + 3) * n4 + i], _mm256_permute2x128_si256(t1, t3, 0x31));
    }
   }
  };
  int skipped= 0;
  for(int col= 0; col < n;) {
   int bt= 0, c0= col;
   fill_cache(c0);
   for(; bt < K && col < n; ++col) {
    if(col == c0 + 8) fill_cache(c0= col);
    const int nr= rem.size(), nr4= (nr + 3) & ~3;
    const u64* cc= &C[(size_t)(col - c0) * n4];
    __m256i us[K];
    for(int s= 0; s < bt; ++s) us[s]= _mm256_set1_epi64x(Ub[s][col]);
    int pi= -1;
    for(int i= 0; i < nr4; i+= 4) {
     __m256i x= _mm256_loadu_si256((const __m256i*)(cc + i));
     for(int s= 0; s < bt; ++s) x= _mm256_add_epi64(x, _mm256_mul_epu32(_mm256_cvtepu32_epi64(_mm_loadu_si128((const __m128i*)&F[(size_t)s * n4 + i])), us[s]));
     const __m256i y= mt.redc(x);
     _mm_storeu_si128((__m128i*)&V[i], pack_lo(y));
     if(pi < 0) {
      const int nz= ~_mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpeq_epi64(y, _mm256_setzero_si256()))) & 15;
      if(nz) pi= i + __builtin_ctz(nz);
     }
    }
    if(pi >= nr) pi= -1;  // 余りの位置 (値は 0) は選ばない
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
    {  // 要の行 U[j] = (a[p][j] + sum_s F[s][pi] U[s][j]) 2^{-32} mod M (右半分は位置 t まで。U は 0 で埋めてある)
     __m256i fs[K];
     for(int s= 0; s < bt; ++s) fs[s]= _mm256_set1_epi64x(F[(size_t)s * n4 + pi]);
     const u64* ap= row(p);
     const int je= std::min(w, (n + t + 1 + 3) & ~3);
     for(int j= col & ~3; j < je; j+= 4) {
      __m256i x= _mm256_loadu_si256((const __m256i*)(ap + j));
      for(int s= 0; s < bt; ++s) x= _mm256_add_epi64(x, _mm256_mul_epu32(_mm256_loadu_si256((const __m256i*)(Ub[s] + j)), fs[s]));
      _mm256_storeu_si256((__m256i*)(up + j), mt.redc(x));
     }
    }
    {  // F[bt][i] = -V[i] c mod M (c = V[pi]^{-1} 2^32)。Shoup の方法で 4 個ずつ。
     const u32 c= (Z::raw(V[pi]).inv() * Z::raw(u32(R))).val(), cs= u32((u64(c) << 32) / M);
     const __m256i Cv= _mm256_set1_epi64x(c), CS= _mm256_set1_epi64x(cs);
     u32* fb= &F[(size_t)bt * n4];
     for(int i= 0; i < nr4; i+= 4) {
      const __m256i vv= _mm256_cvtepu32_epi64(_mm_loadu_si128((const __m128i*)&V[i]));
      const __m256i qq= _mm256_srli_epi64(_mm256_mul_epu32(vv, CS), 32);
      __m256i x= _mm256_sub_epi64(_mm256_mul_epu32(vv, Cv), _mm256_mul_epu32(qq, mt.M));  // [0, 2M)
      x= _mm256_min_epu32(x, _mm256_sub_epi32(x, mt.M32));
      x= _mm256_sub_epi32(mt.M32, x);
      x= _mm256_min_epu32(x, _mm256_sub_epi32(x, mt.M32));
      _mm_storeu_si128((__m128i*)(fb + i), pack_lo(x));
     }
    }
    prow.push_back(p), pcol.push_back(col), pval.push_back(V[pi]);
    Ub[bt++]= up;
    const int last= nr - 1;  // 要の行を抜き、末尾の行の分を位置 pi へ移して、末尾を 0 にする
    rem[pi]= rem[last], rem.pop_back();
    for(int s= 0; s < bt; ++s) F[(size_t)s * n4 + pi]= F[(size_t)s * n4 + last], F[(size_t)s * n4 + last]= 0;
    for(int qq= 0; qq < K; ++qq) C[(size_t)qq * n4 + pi]= C[(size_t)qq * n4 + last], C[(size_t)qq * n4 + last]= 0;
   }
   if(!bt) continue;
   for(int s= bt; s < K; ++s) Ub[s]= zero.data();
   const int jb= col & ~3, je= std::min(w, (n + (int)prow.size() + 3) & ~3), nr= rem.size();
   for(int i= 0; i < nr; ++i) {  // 最後の組も足す (余因子行列が残りの行を読む)
    u32 f[K], any= 0;
    for(int s= 0; s < K; ++s) any|= (f[s]= s < bt ? F[(size_t)s * n4 + i] : 0);
    if(any) update_row(row(rem[i]), Ub, f, jb, je, R);
   }
   for(int s= 0; s < bt; ++s) std::fill(&F[(size_t)s * n4], &F[(size_t)s * n4] + n4, 0u);
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
template <class Z> linalg_internal::hvec<u64> tri_solve(const Elim<Z>& e, const std::vector<int>& cols, int& wx) {
 const u32 M= Z::mod();
 const u64 R= (u64(1) << 32) % M;
 const int r= e.prow.size(), nc= cols.size();
 wx= (nc + 3) & ~3;
 linalg_internal::hvec<u64> X((size_t)r * wx);
 std::vector<u64> acc((size_t)K * wx), zero(wx);
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
 if(!(Z::mod() & 1 ? e.run_v(0) : e.run(0))) return std::nullopt;
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
