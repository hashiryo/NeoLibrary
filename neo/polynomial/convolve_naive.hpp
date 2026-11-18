#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <type_traits>
#include <vector>
#include "neo/algebra/ZMod.hpp"
// 素朴な畳み込み c_k = Σ_{i+j=k} a_i b_j (O(nm))。T は掛け算と足し算を持つ型で、ZMod (法は定数でも実行時でもよい) なら
// AVX2 の核で求める。
// 核は出力 32 個ずつの塊ごとに、a_i を配って b の 32 個との積を u64 の 8 本にためる。b を 8 個ずつ読み、偶数番目と奇数番目を
// _mm256_mul_epu32 で別のレジスタに足す。値は M 未満なので積は M^2 < 2^60 で、8 回ごとに 2^32 の位を 2^32 mod M に掛けて
// 畳めば (2^62 未満になる) u64 に収まる。b の前後には 32 個の 0 を置き、塊の端で添字を調べない。
// 短い側が 16 以下なら、出力ごとに u64 に足すスカラーの道にする (16 M^2 < 2^64)。
// 多項式の族のほかのヘッダ (linear_recurrence、poly_roots など) も、この核 (poly_internal::conv) を使う。
namespace poly_internal {
using u32= unsigned;
using u64= unsigned long long;
template <class T> struct is_zmod: std::false_type {};
template <unsigned M> struct is_zmod<ZMod<M>>: std::true_type {};
template <class T> constexpr bool is_zmod_v= is_zmod<T>::value;
template <class Z> [[gnu::always_inline]] inline u32 r32() { return Z(u64(1) << 32).val(); }
[[gnu::always_inline]] inline __m256i fold(__m256i x, __m256i R) { return _mm256_add_epi64(_mm256_mul_epu32(_mm256_srli_epi64(x, 32), R), _mm256_and_si256(x, _mm256_set1_epi64x(0xffffffffll))); }
constexpr int PAD= 32;
// out[t] = Σ_i a[i] b[k0 + t - i] mod M (0 ≤ t < 32)。bp は b の前後に PAD 個の 0 を置いたもの。
template <class Z> void conv32(const u32* a, int n, const u32* bp, int m, int k0, u32* out, __m256i R) {
 const int ilo= std::max(0, k0 - m + 1), ihi= std::min(n - 1, k0 + 31);
 __m256i e0= _mm256_setzero_si256(), e1= e0, e2= e0, e3= e0, o0= e0, o1= e0, o2= e0, o3= e0;
 for(int i= ilo; i <= ihi;) {
  const int iend= std::min(ihi + 1, i + 8);
  for(; i < iend; ++i) {
   const __m256i x= _mm256_set1_epi32(a[i]);
   const u32* q= bp + PAD + k0 - i;
   const __m256i b0= _mm256_loadu_si256((const __m256i*)q), b1= _mm256_loadu_si256((const __m256i*)(q + 8));
   const __m256i b2= _mm256_loadu_si256((const __m256i*)(q + 16)), b3= _mm256_loadu_si256((const __m256i*)(q + 24));
   e0= _mm256_add_epi64(e0, _mm256_mul_epu32(x, b0)), o0= _mm256_add_epi64(o0, _mm256_mul_epu32(x, _mm256_srli_epi64(b0, 32)));
   e1= _mm256_add_epi64(e1, _mm256_mul_epu32(x, b1)), o1= _mm256_add_epi64(o1, _mm256_mul_epu32(x, _mm256_srli_epi64(b1, 32)));
   e2= _mm256_add_epi64(e2, _mm256_mul_epu32(x, b2)), o2= _mm256_add_epi64(o2, _mm256_mul_epu32(x, _mm256_srli_epi64(b2, 32)));
   e3= _mm256_add_epi64(e3, _mm256_mul_epu32(x, b3)), o3= _mm256_add_epi64(o3, _mm256_mul_epu32(x, _mm256_srli_epi64(b3, 32)));
  }
  e0= fold(e0, R), e1= fold(e1, R), e2= fold(e2, R), e3= fold(e3, R), o0= fold(o0, R), o1= fold(o1, R), o2= fold(o2, R), o3= fold(o3, R);
 }
 alignas(32) u64 ev[16], od[16];
 _mm256_store_si256((__m256i*)ev, e0), _mm256_store_si256((__m256i*)(ev + 4), e1), _mm256_store_si256((__m256i*)(ev + 8), e2), _mm256_store_si256((__m256i*)(ev + 12), e3);
 _mm256_store_si256((__m256i*)od, o0), _mm256_store_si256((__m256i*)(od + 4), o1), _mm256_store_si256((__m256i*)(od + 8), o2), _mm256_store_si256((__m256i*)(od + 12), o3);
 for(int j= 0; j < 16; ++j) out[2 * j]= Z(ev[j]).val(), out[2 * j + 1]= Z(od[j]).val();
}
// out[0..len) = (a b) の 0 次から len - 1 次 (a、b、out は M 未満の値、len ≤ n + m - 1)。out は a、b と重ならないこと。
template <class Z> void conv(const u32* a, int n, const u32* b, int m, u32* out, int len) {
 if(len <= 0) return;
 if(!n || !m) return std::fill(out, out + len, 0u), void();
 if(std::min(n, m) <= 16) {
  if(n > m) std::swap(a, b), std::swap(n, m);
  for(int k= 0; k < len; ++k) {
   u64 s= 0;
   const int lo= std::max(0, k - m + 1), hi= std::min(n - 1, k);
   for(int i= lo; i <= hi; ++i) s+= u64(a[i]) * b[k - i];
   out[k]= Z(s).val();
  }
  return;
 }
 static thread_local std::vector<u32> bp;
 bp.assign(m + 2 * PAD, 0);
 std::copy(b, b + m, bp.begin() + PAD);
 const __m256i R= _mm256_set1_epi64x(r32<Z>());
 int k0= 0;
 for(; k0 + 32 <= len; k0+= 32) conv32<Z>(a, n, bp.data(), m, k0, out + k0, R);
 if(k0 < len) {
  u32 t[32];
  conv32<Z>(a, n, bp.data(), m, k0, t, R);
  std::copy(t, t + (len - k0), out + k0);
 }
}
template <class Z> std::vector<u32> to_u32(const std::vector<Z>& a) {
 std::vector<u32> r(a.size());
 for(std::size_t i= 0; i < a.size(); ++i) r[i]= a[i].val();
 return r;
}
template <class Z> std::vector<Z> from_u32(const std::vector<u32>& a) {
 std::vector<Z> r(a.size());
 for(std::size_t i= 0; i < a.size(); ++i) r[i]= Z::raw(a[i]);
 return r;
}
}  // namespace poly_internal
template <class T> std::vector<T> convolve_naive(const std::vector<T>& a, const std::vector<T>& b) {
 if(a.empty() || b.empty()) return {};
 const int n= int(a.size()), m= int(b.size());
 if constexpr(poly_internal::is_zmod_v<T>) {
  const auto x= poly_internal::to_u32(a), y= poly_internal::to_u32(b);
  std::vector<poly_internal::u32> z(n + m - 1);
  poly_internal::conv<T>(x.data(), n, y.data(), m, z.data(), n + m - 1);
  return poly_internal::from_u32<T>(z);
 } else {
  std::vector<T> c(n + m - 1, T());
  for(int i= 0; i < n; ++i)
   for(int j= 0; j < m; ++j) c[i + j]+= a[i] * b[j];
  return c;
 }
}
