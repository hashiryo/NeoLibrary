#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <cassert>
#include <cstddef>
#include <type_traits>
#include <vector>
// 長さ 2^n の列の、更新付きの部分集合の和。
//   SubsetSum<T>(n, low): 1 点への加算 add(x, v) と、部分集合の和の問い sum(x) = sum_{y ⊆ x} a[y]。bit を下の low 個 (既定 n / 2) と
//     上の残りに分け、下の bit だけゼータ変換した t[h][l] = sum_{l' ⊆ l} a[h][l'] を持つ。add は 2^{low - popcount(x の下)}、sum は
//     2^{popcount(x の上)}。問いが多いなら low を大きく、加算が多いなら小さくする。上位集合の和は添字を補集合にして使う。
//     procon-judge の codechef-DANYANUM の書き比べ ds1 と ds2 の形。
//   SubsetIntervalSum<T>(n) / (a): a ⊆ c ⊆ b の c 全部への加算 add(a, b, v) と、その和の問い sum(a, b) (a ⊆ b でなければ、加算は何もせず
//     和は 0)。操作を時間の向きにバケットに分け、バケットの初めの列 v とその部分集合の和 s を持つ。バケットより前の分の問いは、自由な bit
//     (b \ a) をたどるか、a の 1 の bit を包除で外して s を使うかの少ない方 (2^{n/2} 以下)。
//     [a ⊆ c ⊆ b] = sum_{d ⊆ a} (-1)^{|d|} [c ⊆ b \ d]。バケットの中の加算 (a', b', k) は、a ∪ a' ⊆ b ∩ b' なら k 2^{|b ∩ b'| - |a ∪ a'|}
//     を問いごとに足す (unsigned int は AVX2 で 8 個ずつ)。加算は同じ分け方で「点に足す」か「部分集合全部に足す」の印に分け、作り直す
//     ときに上位集合の和で配る。作り直すのは、バケットの中の加算を問いで足した数の合計が作り直しの手間の目安 (unsigned int は 4 n 2^n、
//     ほかの型は n 2^n) を超えたとき。1 回あたり均して O(sqrt(n 2^n))。procon-judge の luogu-P6562 の書き比べ sbk4 の形。
//   T は T() が 0 で、SubsetSum は += だけ、SubsetIntervalSum は +、-、* と整数からの変換があればよい。
template <class T> struct SubsetSum {
 int n, low;
 unsigned lm;
 std::vector<T> t;
 explicit SubsetSum(int n_, int low_= -1): n(n_), low(low_ < 0 ? n_ / 2 : low_), lm((1u << low) - 1), t(std::size_t(1) << n_) { assert(0 <= low && low <= n); }
 // a から作る (a の長さは 2^n)。O(low 2^n)
 explicit SubsetSum(const std::vector<T>& a, int low_= -1): SubsetSum(__builtin_ctzll(a.size()), low_) {
  assert(a.size() == (std::size_t(1) << n));
  t= a;
  const std::size_t W= std::size_t(1) << low;
  for(std::size_t r= 0; r < t.size(); r+= W)
   for(std::size_t w= 1; w < W; w<<= 1)
    for(std::size_t p= r; p < r + W; p+= 2 * w)
     for(std::size_t k= 0; k < w; ++k) t[p + w + k]+= t[p + k];
 }
 void add(unsigned x, T v) {
  T* row= t.data() + (std::size_t(x >> low) << low);
  const unsigned lo= x & lm, fr= lm ^ lo;
  for(unsigned s= fr;; s= (s - 1) & fr) {
   row[lo | s]+= v;
   if(!s) break;
  }
 }
 T sum(unsigned x) const {
  const unsigned h= x >> low, lo= x & lm;
  T r= T();
  for(unsigned s= h;; s= (s - 1) & h) {
   r+= t[(std::size_t(s) << low) | lo];
   if(!s) break;
  }
  return r;
 }
};
namespace subset_sum_internal {
template <class T> void zeta_sub(T* f, std::size_t N) {  // f[y] <- sum_{c ⊆ y} f[c]
 for(std::size_t w= 1; w < N; w<<= 1)
  for(std::size_t p= 0; p < N; p+= 2 * w)
   for(std::size_t k= 0; k < w; ++k) f[p + w + k]+= f[p + k];
}
template <class T> void zeta_sup(T* f, std::size_t N) {  // f[x] <- sum_{c ⊇ x} f[c]
 for(std::size_t w= 1; w < N; w<<= 1)
  for(std::size_t p= 0; p < N; p+= 2 * w)
   for(std::size_t k= 0; k < w; ++k) f[p + k]+= f[p + w + k];
}
}  // namespace subset_sum_internal
template <class T> struct SubsetIntervalSum {
 explicit SubsetIntervalSum(int n_): SubsetIntervalSum(std::vector<T>(std::size_t(1) << n_)) {}
 explicit SubsetIntervalSum(const std::vector<T>& a): n(__builtin_ctzll(a.size())), full((1u << n) - 1), v(a) {
  assert(a.size() == (std::size_t(1) << n));
  const std::size_t N= a.size();
  s.assign(N, T()), mv.assign(N, T()), ms.assign(N, T());
  pw.assign(n + 1, T());
  pw[0]= T(1);
  for(int i= 1; i <= n; ++i) pw[i]= pw[i - 1] + pw[i - 1];
  budget= (long long)(SIMD ? 4 : 1) * n * (long long)N;
  build();
 }
 // a ⊆ c ⊆ b の c 全部に k を足す
 void add(unsigned a, unsigned b, T k) {
  if((a & b) != a) return;
  if(cnt % 8 == 0) ua.resize(cnt + 8, full), ub.resize(cnt + 8, 0), uk.resize(cnt + 8, T());  // 埋め草は k = 0
  ua[cnt]= a, ub[cnt]= b, uk[cnt]= k, ++cnt;
  const unsigned fr= a ^ b;
  if(__builtin_popcount(fr) <= __builtin_popcount(a)) {
   for(unsigned d= fr;; d= (d - 1) & fr) {
    mv[a | d]+= k;
    if(!d) break;
   }
  } else {
   for(unsigned d= a;; d= (d - 1) & a) {
    if(__builtin_popcount(d) & 1) ms[b ^ d]-= k;
    else ms[b ^ d]+= k;
    if(!d) break;
   }
  }
 }
 // a ⊆ c ⊆ b の c 全部の和
 T sum(unsigned a, unsigned b) {
  if((a & b) != a) return T();
  if(scan + (long long)cnt > budget) build();
  scan+= (long long)cnt;
  return base(a, b) + inblock(a, b);
 }
private:
 static constexpr bool SIMD= std::is_same_v<T, unsigned>;
 int n;
 unsigned full;
 std::vector<T> v, s, mv, ms, pw;  // v: バケットの初めの列、s: v の部分集合の和、mv と ms: バケットの中の加算の、点と部分集合全部への印
 std::vector<unsigned> ua, ub;     // バケットの中の加算 (8 個ずつに埋める)
 std::vector<T> uk;
 std::size_t cnt= 0;
 long long scan= 0, budget;
 void build() {
  using namespace subset_sum_internal;
  const std::size_t N= v.size();
  zeta_sup(ms.data(), N);
  for(std::size_t j= 0; j < N; ++j) v[j]+= mv[j] + ms[j], mv[j]= ms[j]= T();
  s= v, zeta_sub(s.data(), N);
  cnt= 0, scan= 0;
  ua.clear(), ub.clear(), uk.clear();
 }
 T base(unsigned a, unsigned b) const {
  const unsigned fr= a ^ b;
  T r= T();
  if(__builtin_popcount(fr) <= __builtin_popcount(a)) {
   for(unsigned d= fr;; d= (d - 1) & fr) {
    r+= v[a | d];
    if(!d) break;
   }
  } else {
   for(unsigned d= a;; d= (d - 1) & a) {
    if(__builtin_popcount(d) & 1) r-= s[b ^ d];
    else r+= s[b ^ d];
    if(!d) break;
   }
  }
  return r;
 }
 T inblock(unsigned a, unsigned b) const {
  if constexpr(SIMD) {
   const __m256i X= _mm256_set1_epi32(int(a)), Y= _mm256_set1_epi32(int(b)), Z= _mm256_setzero_si256();
   const __m256i lut= _mm256_setr_epi8(0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4);
   const __m256i m4= _mm256_set1_epi8(0x0F), one8= _mm256_set1_epi8(1), one16= _mm256_set1_epi16(1);
   __m256i acc= Z;
   for(std::size_t i= 0; i < cnt; i+= 8) {
    const __m256i A= _mm256_loadu_si256((const __m256i*)(ua.data() + i)), B= _mm256_loadu_si256((const __m256i*)(ub.data() + i)), K= _mm256_loadu_si256((const __m256i*)(uk.data() + i));
    const __m256i lo= _mm256_or_si256(X, A), hi= _mm256_and_si256(Y, B);
    const __m256i ok= _mm256_cmpeq_epi32(_mm256_andnot_si256(hi, lo), Z);
    const __m256i fr= _mm256_xor_si256(hi, lo);
    // fr の bit の数 (バイトごとに 4 bit ずつ表で引いて、32 bit ごとに足す)
    const __m256i c8= _mm256_add_epi8(_mm256_shuffle_epi8(lut, _mm256_and_si256(fr, m4)), _mm256_shuffle_epi8(lut, _mm256_and_si256(_mm256_srli_epi32(fr, 4), m4)));
    const __m256i c32= _mm256_madd_epi16(_mm256_maddubs_epi16(c8, one8), one16);
    acc= _mm256_add_epi32(acc, _mm256_and_si256(_mm256_sllv_epi32(K, c32), ok));
   }
   alignas(32) unsigned tmp[8];
   _mm256_store_si256((__m256i*)tmp, acc);
   return tmp[0] + tmp[1] + tmp[2] + tmp[3] + tmp[4] + tmp[5] + tmp[6] + tmp[7];
  } else {
   T r= T();
   for(std::size_t i= 0; i < cnt; ++i) {
    const unsigned lo= a | ua[i], hi= b & ub[i];
    if((lo & hi) == lo) r+= uk[i] * pw[__builtin_popcount(lo ^ hi)];
   }
   return r;
  }
 }
};
