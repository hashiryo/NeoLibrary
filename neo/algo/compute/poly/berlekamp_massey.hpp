#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <vector>
#include "neo/base/core/poly/convolve_naive.hpp"
// Berlekamp–Massey の方法。列 a_0, ..., a_{n-1} について、a_i = Σ_{j=1}^{d} c_j a_{i-j} (d ≤ i < n) を満たす最短の
// c_1, ..., c_d を O(n^2) で求める。K は体で、ZMod なら法は素数であること。
// ZMod は値を M 未満の u32 で持ち、ずれの内積は AVX2 で u64 にためて 8 回ごとに 2^32 の位を畳み、漸化式の更新
// C -= x z^m B は掛ける数が決まっているので Shoup の方法で 8 個ずつ掛ける。長さが伸びる回は、新しい C を B の置き場に後ろから
// 書いて置き場を入れ替え、写しを作らない。
// procon-judge の yosupo-find-linear-recurrence で書き比べた bm_avx2 を、法を定数にしない形で写したもの。
namespace bm_internal {
using u32= unsigned;
using u64= unsigned long long;
// Σ x[i] y[i] mod M (x、y は M 未満)
template <class Z> u32 dot(const u32* x, const u32* y, int n, __m256i R) {
 __m256i a0= _mm256_setzero_si256(), a1= _mm256_setzero_si256();
 int i= 0;
 for(int t= 0; i + 8 <= n; i+= 8) {
  const __m256i u= _mm256_loadu_si256((const __m256i*)(x + i)), v= _mm256_loadu_si256((const __m256i*)(y + i));
  a0= _mm256_add_epi64(a0, _mm256_mul_epu32(u, v));
  a1= _mm256_add_epi64(a1, _mm256_mul_epu32(_mm256_srli_epi64(u, 32), _mm256_srli_epi64(v, 32)));
  if(++t == 8) a0= poly_internal::fold(a0, R), a1= poly_internal::fold(a1, R), t= 0;
 }
 a0= _mm256_add_epi64(poly_internal::fold(a0, R), poly_internal::fold(a1, R));
 alignas(32) u64 s[4];
 _mm256_store_si256((__m256i*)s, a0);
 Z r= Z(s[0]) + Z(s[1]) + Z(s[2]) + Z(s[3]);
 for(; i < n; ++i) r+= Z::raw(x[i]) * Z::raw(y[i]);
 return r.val();
}
// a x mod M を 8 個ずつ (a < M、x < M、q = floor(x 2^32 / M))。M 未満で返す。
[[gnu::always_inline]] inline __m256i mul8(__m256i a, __m256i x, __m256i q, __m256i p) {
 const __m256i ev= _mm256_srli_epi64(_mm256_mul_epu32(a, q), 32), od= _mm256_mul_epu32(_mm256_srli_epi64(a, 32), q);
 const __m256i r= _mm256_sub_epi32(_mm256_mullo_epi32(a, x), _mm256_mullo_epi32(_mm256_blend_epi32(ev, od, 0b10101010), p));
 return _mm256_min_epu32(r, _mm256_sub_epi32(r, p));
}
[[gnu::always_inline]] inline u32 mul1(u32 a, u32 x, u32 q, u32 M) {
 const u32 r= a * x - u32((u64(a) * q) >> 32) * M;
 return r >= M ? r - M : r;
}
// d[i] = c[i] - x b[i] (i < n)。d と c は同じでよい。
inline void sub_mul(u32* d, const u32* c, const u32* b, int n, u32 x, u32 M) {
 const u32 q= u32((u64(x) << 32) / M);
 const __m256i vx= _mm256_set1_epi32(x), vq= _mm256_set1_epi32(q), p= _mm256_set1_epi32(M);
 int i= 0;
 for(; i + 8 <= n; i+= 8) {
  const __m256i r= _mm256_sub_epi32(_mm256_loadu_si256((const __m256i*)(c + i)), mul8(_mm256_loadu_si256((const __m256i*)(b + i)), vx, vq, p));
  _mm256_storeu_si256((__m256i*)(d + i), _mm256_min_epu32(r, _mm256_add_epi32(r, p)));
 }
 for(; i < n; ++i) {
  const u32 r= c[i] - mul1(b[i], x, q, M);
  d[i]= int(r) < 0 ? r + M : r;
 }
}
// e[j] = c[j] - x e[j - m] (m ≤ j < m + n) を後ろから求める。8 個の塊ごとに読みを先に済ませるので、上書きの前に読める。
inline void sub_mul_back(u32* e, const u32* c, int m, int n, u32 x, u32 M) {
 const u32 q= u32((u64(x) << 32) / M);
 const __m256i vx= _mm256_set1_epi32(x), vq= _mm256_set1_epi32(q), p= _mm256_set1_epi32(M);
 int j= m + n;
 for(; j - 8 >= m; j-= 8) {
  const __m256i bb= _mm256_loadu_si256((const __m256i*)(e + j - 8 - m));
  const __m256i r= _mm256_sub_epi32(_mm256_loadu_si256((const __m256i*)(c + j - 8)), mul8(bb, vx, vq, p));
  _mm256_storeu_si256((__m256i*)(e + j - 8), _mm256_min_epu32(r, _mm256_add_epi32(r, p)));
 }
 for(--j; j >= m; --j) {
  const u32 r= c[j] - mul1(e[j - m], x, q, M);
  e[j]= int(r) < 0 ? r + M : r;
 }
}
template <class Z> std::vector<Z> zmod(const std::vector<Z>& seq) {
 const int n= int(seq.size());
 const u32 M= Z::mod();
 std::vector<u32> ar(n);
 for(int i= 0; i < n; ++i) ar[n - 1 - i]= seq[i].val();
 // C と B は長さの先が 0 になるように保つ。
 std::vector<u32> bufC(n + 9, 0), bufB(n + 9, 0);
 u32 *C= bufC.data(), *B= bufB.data();
 C[0]= B[0]= 1 % M;
 const __m256i R= _mm256_set1_epi64x(poly_internal::r32<Z>());
 int L= 0, m= 1, lc= 1, lb= 1;
 Z binv= Z(1);
 for(int i= 0; i < n; ++i) {
  const u32 d= dot<Z>(C, ar.data() + (n - 1 - i), L + 1, R);
  if(d == 0) {
   ++m;
   continue;
  }
  const u32 x= (Z::raw(d) * binv).val();
  if(2 * L <= i) {
   const int nl= std::max(lc, lb + m);
   sub_mul_back(B, C, m, lb, x, M);
   std::copy(C, C + m, B);
   if(lb + m < lc) std::copy(C + lb + m, C + lc, B + lb + m);
   std::swap(C, B);
   lb= lc, lc= nl;
   L= i + 1 - L, binv= Z::raw(d).inv(), m= 1;
  } else {
   sub_mul(C + m, C + m, B, lb, x, M);
   lc= std::max(lc, lb + m);
   ++m;
  }
 }
 std::vector<Z> res(L);
 for(int j= 1; j <= L; ++j) res[j - 1]= -Z::raw(C[j]);
 return res;
}
template <class K> std::vector<K> generic(const std::vector<K>& a) {
 const int n= int(a.size());
 std::vector<K> C{K(1)}, B{K(1)};
 int L= 0, m= 1;
 K b= K(1);
 for(int i= 0; i < n; ++i) {
  K d= K();
  for(int j= 0; j <= L && j < int(C.size()); ++j) d+= C[j] * a[i - j];
  if(d == K()) {
   ++m;
   continue;
  }
  const K x= d / b;
  std::vector<K> T= C;
  if(C.size() < B.size() + m) C.resize(B.size() + m, K());
  for(std::size_t j= 0; j < B.size(); ++j) C[j + m]-= x * B[j];
  if(2 * L <= i) L= i + 1 - L, B.swap(T), b= d, m= 1;
  else ++m;
 }
 C.resize(L + 1, K());
 std::vector<K> res(L);
 for(int j= 1; j <= L; ++j) res[j - 1]= -C[j];
 return res;
}
}  // namespace bm_internal
template <class K> std::vector<K> berlekamp_massey(const std::vector<K>& a) {
 if constexpr(poly_internal::is_zmod_v<K>) return bm_internal::zmod(a);
 else return bm_internal::generic(a);
}
