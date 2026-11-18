#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>
#include "neo/algebra/ZMod.hpp"
#include "neo/linalg/Matrix.hpp"
// ハフニアン haf(A)。対角が 0 の n × n の対称行列 A を辺の重みとみて、完全マッチングの重み (辺の重みの積) の和を返す。n が奇数なら
// 0、n = 0 なら 1。対角と、対角より下の成分は読まない (a[i][j] は i < j のものだけを使う)。T は可換環 (和、積、0、1)。
// 頂点 2p と 2p+1 を組 p (p < h = n/2) とする。完全マッチングに組の辺を重ねると交互の閉路に分かれるので、閉路を最大の組の小さい
// 順に 1 本ずつ作る。閉路は最大の組 w の頂点 2w から出て、マッチングの辺で別の組の頂点 k に入り、組の辺で k^1 へ移る、を繰り返し、
// 2w+1 へ戻って閉じる。c[S] を組の集合 S をちょうど覆って閉路がすべて閉じた重みの和、D_e[T] を S = {w} ∪ T で最大の組の閉路が
// いま頂点 e で開いている重みの和とすると、層 w の T ごとの
//   v_T[k] = c[T] a[2w][k] + sum_{p ∈ T} (D_{2p}[T] a[2p][k] + D_{2p+1}[T] a[2p+1][k])
// で、v_T[2w+1] が c[S]、T に入らない組 p へは D_{2p+1}[T + 2^p] = v_T[2p]、D_{2p}[T + 2^p] = v_T[2p+1] を配る。答えは c[全体]。
// 手間は O(2^{n/2} n^2)。
// ZMod で法が奇数なら、v_T を AVX2 で列の向きに 4 個ずつ求める。隣接行列は Montgomery の表現で持ち、u64 の和にためて、8 項ごとに
// 最上位の bit が立っていれば 2M 2^32 を引き、最後に 1 回還元する。D_e は頂点ごとの配列で、組 p の頂点の配列は T から bit p を
// 抜いた番号で引く (T を小さい順に見ると、どの配列も前から順に読み書きされる)。層 8 以上では、組を下の 5 個と上の残りに分け、
// 上の組の部分集合 TH ごとに、配る先の列 (下の組の頂点、2w+1、上の組で T に入らない組の頂点) と使う行 (下の組の頂点、2w、TH の
// 組の頂点) だけを並べた隣接行列で回す (全部の列を掛けると積和の半分ほどを捨てる)。TH は小さい順に回し、1 つ前の TH から出入り
// した組 (2 進の繰り上がりなので平均 2 組) の行と列だけを、末尾との入れ替えで直す (Gray 符号の順にすると D の読み書きが飛ぶ)。
// ほかの T と法が偶数の ZMod は、同じ動的計画を要素の演算子で配る形にする。procon-judge の yosupo-hafnian-of-matrix で書き比べた
// zmod_layer_h の法を実行時の奇数に広げ、下の組の数と作り直しを zmod_layer_i5 の形にしたもの。
namespace hafnian_internal {
using u32= unsigned;
using u64= unsigned long long;
// T から bit p を抜いた番号。
[[gnu::always_inline]] inline u32 drop_bit(u32 T, int p) {
#ifdef __BMI2__
 return _pext_u32(T, ~(1u << p));
#else
 const u32 low= (1u << p) - 1;
 return (T & low) | ((T >> 1) & ~low);
#endif
}
struct Consts {
 __m256i M, M2, TOP, NI;
 explicit Consts(u32 m) {
  u32 x= m;
  for(int i= 0; i < 5; ++i) x*= 2 - m * x;
  M= _mm256_set1_epi32(m), M2= _mm256_set1_epi32(2 * m), TOP= _mm256_set1_epi64x((long long)(1ull << 63)), NI= _mm256_set1_epi32(-x);
 }
};
#define HAF_REP(X) X(0) X(1) X(2) X(3) X(4) X(5) X(6) X(7) X(8) X(9)
#define HAF_KEEP(q) \
 if constexpr(NV > q) a##q= _mm256_sub_epi64(a##q, _mm256_min_epu32(k.M2, _mm256_and_si256(k.TOP, a##q)));
#define HAF_INIT(q) \
 if constexpr(NV > q) a##q= _mm256_mul_epu32(x, _mm256_load_si256((const __m256i*)(aw + 4 * q)));
#define HAF_MAC(q) \
 if constexpr(NV > q) a##q= _mm256_add_epi64(a##q, _mm256_add_epi64(_mm256_mul_epu32(x0, _mm256_load_si256((const __m256i*)(r0 + 4 * q))), _mm256_mul_epu32(x1, _mm256_load_si256((const __m256i*)(r1 + 4 * q)))));
#define HAF_RED(q) \
 if constexpr(NV > q) { \
  __m256i r= _mm256_srli_epi64(_mm256_add_epi64(a##q, _mm256_mul_epu32(_mm256_mul_epu32(a##q, k.NI), k.M)), 32); \
  r= _mm256_min_epu32(r, _mm256_sub_epi32(r, k.M2)); \
  _mm256_store_si256((__m256i*)(v + 4 * q), _mm256_min_epu32(r, _mm256_sub_epi32(r, k.M))); \
 }
// 層 w を、全部の列 (k < 2w + 2) を掛けて回す。A は行 e (e <= 2w) ごとに NV * 4 個の u64 (Montgomery の表現、残りは 0)。
template <int NV> [[gnu::noinline]] void layer(int w, const u64* __restrict A, u32* __restrict D, std::size_t stride, u32* __restrict c, const Consts& k) {
 const u32 bw= 1u << w;
 const u64* aw= A + std::size_t(2 * w) * NV * 4;
 alignas(32) u64 v[NV * 4];
 for(u32 T= 0; T < bw; ++T) {
  __m256i a0, a1, a2, a3, a4, a5, a6, a7, a8, a9;
  {
   const __m256i x= _mm256_set1_epi32(c[T]);
   HAF_REP(HAF_INIT)
  }
  int cnt= 1;
  for(u32 t= T; t; t&= t - 1) {
   const int p= __builtin_ctz(t);
   if(cnt > 6) {
    HAF_REP(HAF_KEEP)
    cnt= 0;
   }
   const u32 ix= drop_bit(T, p);
   const __m256i x0= _mm256_set1_epi32(D[std::size_t(2 * p) * stride + ix]), x1= _mm256_set1_epi32(D[std::size_t(2 * p + 1) * stride + ix]);
   const u64 *r0= A + std::size_t(2 * p) * NV * 4, *r1= r0 + NV * 4;
   HAF_REP(HAF_MAC)
   cnt+= 2;
  }
  HAF_REP(HAF_KEEP)  // 2^63 未満にしてから還元する
  HAF_REP(HAF_RED)
  c[bw + T]= u32(v[2 * w + 1]);
  for(u32 t= ~T & (bw - 1); t; t&= t - 1) {
   const int p= __builtin_ctz(t);
   const u32 ix= drop_bit(T, p);
   D[std::size_t(2 * p + 1) * stride + ix]= u32(v[2 * p]);
   D[std::size_t(2 * p) * stride + ix]= u32(v[2 * p + 1]);
  }
 }
}
// 層 w で上の組の部分集合を THbits (bit s..w-1) に決めたときの、下の組の部分集合 TL の処理。行の幅は WS (u64 の数)。
template <int NV> [[gnu::noinline]] void layer_i(int w, int s, u32 THbits, const int* thp, int nth, const int* chp, int nch, const u64* __restrict A2, std::size_t WS, u32* __restrict D, std::size_t stride, u32* __restrict c, const Consts& k) {
 const u32 nl= 1u << s;
 const u64* aw= A2 + std::size_t(2 * s) * WS;
 alignas(32) u64 v[NV * 4];
 for(u32 TL= 0; TL < nl; ++TL) {
  const u32 T= THbits | TL;
  __m256i a0, a1, a2, a3, a4, a5, a6, a7, a8, a9;
  {
   const __m256i x= _mm256_set1_epi32(c[T]);
   HAF_REP(HAF_INIT)
  }
  int cnt= 1;
  for(u32 t= TL; t; t&= t - 1) {
   const int p= __builtin_ctz(t);
   if(cnt > 6) {
    HAF_REP(HAF_KEEP)
    cnt= 0;
   }
   const u32 ix= drop_bit(T, p);
   const __m256i x0= _mm256_set1_epi32(D[std::size_t(2 * p) * stride + ix]), x1= _mm256_set1_epi32(D[std::size_t(2 * p + 1) * stride + ix]);
   const u64 *r0= A2 + std::size_t(2 * p) * WS, *r1= r0 + WS;
   HAF_REP(HAF_MAC)
   cnt+= 2;
  }
  for(int i= 0; i < nth; ++i) {
   const int p= thp[i];
   if(cnt > 6) {
    HAF_REP(HAF_KEEP)
    cnt= 0;
   }
   const u32 ix= drop_bit(T, p);
   const __m256i x0= _mm256_set1_epi32(D[std::size_t(2 * p) * stride + ix]), x1= _mm256_set1_epi32(D[std::size_t(2 * p + 1) * stride + ix]);
   const u64 *r0= A2 + std::size_t(2 * s + 1 + 2 * i) * WS, *r1= r0 + WS;
   HAF_REP(HAF_MAC)
   cnt+= 2;
  }
  HAF_REP(HAF_KEEP)
  HAF_REP(HAF_RED)
  c[(1u << w) + T]= u32(v[2 * s]);
  for(u32 t= ~TL & (nl - 1); t; t&= t - 1) {
   const int p= __builtin_ctz(t);
   const u32 ix= drop_bit(T, p);
   D[std::size_t(2 * p + 1) * stride + ix]= u32(v[2 * p]);
   D[std::size_t(2 * p) * stride + ix]= u32(v[2 * p + 1]);
  }
  for(int i= 0; i < nch; ++i) {
   const int p= chp[i];
   const u32 ix= drop_bit(T, p);
   D[std::size_t(2 * p + 1) * stride + ix]= u32(v[2 * s + 1 + 2 * i]);
   D[std::size_t(2 * p) * stride + ix]= u32(v[2 * s + 2 + 2 * i]);
  }
 }
}
#undef HAF_KEEP
#undef HAF_INIT
#undef HAF_MAC
#undef HAF_RED
// a は n × n (n = 2h、M 未満)。法は奇数。
inline u32 zmod_odd(int n, const u32* a, u32 Mod) {
 constexpr int W0= 8, SL= 5;  // W0 より下の層は全部の列を掛ける。上の層の下の組は SL 個
 const int h= n / 2;
 const Consts k(Mod);
 const std::size_t stride= std::max<std::size_t>(1, std::size_t(1) << std::max(0, h - 2));
 linalg_internal::Region D(std::size_t(std::max(1, 2 * (h - 1))) * stride), C(std::size_t(1) << h);
 u32* c= C.p;
 c[0]= 1 % Mod;
 std::vector<u64> Am((std::size_t)n * n);  // Montgomery の表現。対称にして持つ
 for(int i= 0; i < n; ++i)
  for(int j= 0; j < n; ++j) Am[(std::size_t)i * n + j]= i == j ? 0 : (u64(a[(std::size_t)std::min(i, j) * n + std::max(i, j)]) << 32) % Mod;
 std::vector<u64> buf((std::size_t)(n + 1) * 44 + 4);
 u64* Ap= reinterpret_cast<u64*>((reinterpret_cast<std::uintptr_t>(buf.data()) + 31) & ~std::uintptr_t(31));
 for(int w= 0; w < h; ++w) {
  if(w < W0) {
   const int nv= (2 * w + 2 + 3) / 4;
   for(int e= 0; e <= 2 * w; ++e)
    for(int q= 0; q < nv * 4; ++q) Ap[std::size_t(e) * nv * 4 + q]= q < 2 * w + 2 ? Am[(std::size_t)e * n + q] : 0;
   switch(nv) {
#define HAF_CASE(q) \
 case q + 1: layer<q + 1>(w, Ap, D.p, stride, c, k); break;
    HAF_REP(HAF_CASE)
#undef HAF_CASE
   }
   continue;
  }
  const int s= std::min(SL, w), nh= w - s;
  const int KM= 2 * s + 1 + 2 * nh;
  const std::size_t WS= std::size_t((KM + 3) / 4) * 4;
  int thp[32], chp[32], pth[32], pch[32], rowv[64], colv[64];
  int nth= 0, nch= 0;
  for(int e= 0; e < 2 * s; ++e) rowv[e]= e, colv[e]= e;
  rowv[2 * s]= 2 * w, colv[2 * s]= 2 * w + 1;
  for(int i= 0; i < nh; ++i) {
   const int p= s + i;
   chp[i]= p, pch[p]= i, colv[2 * s + 1 + 2 * i]= 2 * p, colv[2 * s + 2 + 2 * i]= 2 * p + 1;
  }
  nch= nh;
  auto A2= [&](int r) { return Ap + std::size_t(r) * WS; };
  for(int r= 0; r <= 2 * s; ++r)
   for(std::size_t q= 0; q < WS; ++q) A2(r)[q]= int(q) < KM ? Am[(std::size_t)rowv[r] * n + colv[q]] : 0;
  // 組 p が TH から出る: 行を末尾の組の行で埋め、列を末尾に足す
  auto leave= [&](int p) {
   const int i= pth[p], last= nth - 1;
   if(i != last) {
    const int p2= thp[last];
    for(int d= 0; d < 2; ++d) {
     std::memcpy(A2(2 * s + 1 + 2 * i + d), A2(2 * s + 1 + 2 * last + d), WS * sizeof(u64));
     rowv[2 * s + 1 + 2 * i + d]= rowv[2 * s + 1 + 2 * last + d];
    }
    thp[i]= p2, pth[p2]= i;
   }
   --nth;
   const int R= 2 * s + 1 + 2 * nth, j= nch++;
   chp[j]= p, pch[p]= j;
   for(int d= 0; d < 2; ++d) {
    const int q= 2 * s + 1 + 2 * j + d;
    colv[q]= 2 * p + d;
    for(int r= 0; r < R; ++r) A2(r)[q]= Am[(std::size_t)rowv[r] * n + colv[q]];
   }
  };
  // 組 p が TH に入る: 列を末尾の組の列で埋めて末尾を 0 にし、行を末尾に足す
  auto enter= [&](int p) {
   const int R= 2 * s + 1 + 2 * nth, i= pch[p], last= nch - 1;
   if(i != last) {
    const int p2= chp[last];
    for(int d= 0; d < 2; ++d) {
     const int qi= 2 * s + 1 + 2 * i + d, ql= 2 * s + 1 + 2 * last + d;
     for(int r= 0; r < R; ++r) A2(r)[qi]= A2(r)[ql];
     colv[qi]= colv[ql];
    }
    chp[i]= p2, pch[p2]= i;
   }
   for(int d= 0; d < 2; ++d)
    for(int r= 0; r < R; ++r) A2(r)[2 * s + 1 + 2 * last + d]= 0;
   --nch;
   const int K= 2 * s + 1 + 2 * nch, j= nth++;
   thp[j]= p, pth[p]= j;
   for(int d= 0; d < 2; ++d) {
    const int r= 2 * s + 1 + 2 * j + d;
    rowv[r]= 2 * p + d;
    u64* dst= A2(r);
    for(int q= 0; q < K; ++q) dst[q]= Am[(std::size_t)rowv[r] * n + colv[q]];
    for(std::size_t q= K; q < WS; ++q) dst[q]= 0;
   }
  };
  for(u32 TH= 0; TH < (1u << nh); ++TH) {
   if(TH) {
    const int t= __builtin_ctz(TH);
    for(int b= 0; b < t; ++b) leave(s + b);
    enter(s + t);
   }
   const int K= 2 * s + 1 + 2 * nch, nv= (K + 3) / 4;
   switch(nv) {
#define HAF_CASE(q) \
 case q + 1: layer_i<q + 1>(w, s, TH << s, thp, nth, chp, nch, Ap, WS, D.p, stride, c, k); break;
    HAF_REP(HAF_CASE)
#undef HAF_CASE
   }
  }
 }
 return c[(std::size_t(1) << h) - 1];
}
#undef HAF_REP
// 要素の演算子で配る形。dp[S][j] は組の集合 S を覆い、最大の組の閉路がいま頂点 j で開いている重みの和。
template <class T> T generic(const Matrix<T>& A) {
 const int n= A.h, h= n / 2;
 const T zero= linalg_internal::zero<T>(), one= linalg_internal::one<T>();
 auto at= [&](int i, int j) { return i < j ? A[i][j] : A[j][i]; };
 std::vector<T> dp((std::size_t(1) << h) * n, zero);
 auto D= [&](u32 S, int j) -> T& { return dp[(std::size_t)S * n + j]; };
 for(int p= 0; p < h; ++p) D(1u << p, 2 * p)= one;
 T ans= zero;
 for(int w= 0; w < h; ++w)
  for(u32 S= 1u << w; S < (2u << w); ++S)
   for(int j= 0; j < 2 * w + 1; ++j) {
    if(!(S >> (j / 2) & 1)) continue;
    const T x= D(S, j);
    for(int q= 0; q < 2 * w; ++q) {
     if(S >> (q / 2) & 1) continue;
     D(S | (1u << (q / 2)), q ^ 1)+= x * at(j, q);
    }
    const T y= x * at(j, 2 * w + 1);
    for(int q= w + 1; q < h; ++q) D(S | (1u << q), 2 * q)+= y;
    if(S + 1 == (std::size_t(1) << h)) ans+= y;
   }
 return ans;
}
}
template <class T> T hafnian(const Matrix<T>& A) {
 assert(A.h == A.w);
 const int n= A.h;
 if(n & 1) return linalg_internal::zero<T>();
 if(n == 0) return linalg_internal::one<T>();
 assert(n <= 62);  // 組の集合を u32 で持つ
 if constexpr(linalg_internal::is_zmod_v<T>) {
  if(T::mod() & 1) {
   std::vector<unsigned> a((std::size_t)n * n);
   for(int i= 0; i < n; ++i)
    for(int j= 0; j < n; ++j) a[(std::size_t)i * n + j]= A[i][j].val();
   return T::raw(hafnian_internal::zmod_odd(n, a.data(), T::mod()));
  }
 }
 return hafnian_internal::generic(A);
}
