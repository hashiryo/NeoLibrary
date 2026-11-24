#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <cstdint>
#include <new>
#include <type_traits>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
#include "neo/base/core/algebra/ZMod.hpp"
// 素朴な畳み込み c_k = Σ_{i+j=k} a_i b_j (O(nm))。T は掛け算と足し算を持つ型で、ZMod (法は定数でも実行時でもよい) なら
// AVX2 の核で求める。
// 核は出力 32 個ずつの塊ごとに、a_i を配って b の 32 個との積を u64 の 8 本にためる。b を 8 個ずつ読み、偶数番目と奇数番目を
// _mm256_mul_epu32 で別のレジスタに足す。値は M 未満なので積は M^2 < 2^60 で、8 回ごとに 2^32 の位を 2^32 mod M に掛けて
// 畳めば (2^62 未満になる) u64 に収まる。b の前後には 32 個の 0 を置き、塊の端で添字を調べない。
// 短い側が 16 以下なら、出力ごとに u64 に足すスカラーの道にする (16 M^2 < 2^64)。
// 多項式の族のほかのヘッダ (linear_recurrence、poly_roots など) も、この核 (convolve_internal::conv) を使う。
// 末尾の構造体 Naive は、同じ核を neo/base/core/convolve/convolve.hpp の口に渡す形にしたもので、積の一部の範囲だけを求められる。
namespace convolve_internal {
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
// ZMod の列をそのままの値 (2M 未満) の u32 の列として読み書きする。ZMod は u32 を 1 つだけ持つ標準レイアウトの型なので、
// 先頭の要素とポインタを相互に変換できる (neo/base/core/set_fps/ranked_xor.hpp と同じ)。核は M 未満の値を書くので、そのまま
// ZMod の値になる。
template <class Z> inline const u32* raw(const std::vector<Z>& v) {
 static_assert(std::is_standard_layout_v<Z> && sizeof(Z) == 4);
 return reinterpret_cast<const u32*>(v.data());
}
template <class Z> inline u32* raw(std::vector<Z>& v) {
 static_assert(std::is_standard_layout_v<Z> && sizeof(Z) == 4);
 return reinterpret_cast<u32*>(v.data());
}
// 2 MB 以上の作業の配列を、Linux では 2 MB 境界の mmap に置いて huge page を頼む allocator。4 KB のページに初めて触る手間
// (CI のランナーで 1 MB あたり 0.4 ms ほど) を減らす。
template <class T> struct HugeAlloc {
 using value_type= T;
 static constexpr std::size_t H= std::size_t(1) << 21, P= 4096;
 HugeAlloc()= default;
 template <class U> HugeAlloc(const HugeAlloc<U>&) {}
 T* allocate(std::size_t n) {
#ifdef __linux__
  if(n * sizeof(T) >= H) {
   const std::size_t bytes= (n * sizeof(T) + P - 1) & ~(P - 1);
   void* m= mmap(nullptr, bytes + H, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
   if(m == MAP_FAILED) throw std::bad_alloc();
   char *p= static_cast<char*>(m), *a= reinterpret_cast<char*>((reinterpret_cast<std::uintptr_t>(p) + H - 1) & ~(H - 1));
   if(a != p) munmap(p, a - p);
   munmap(a + bytes, p + H - a);
   madvise(a, bytes, MADV_HUGEPAGE);
   return reinterpret_cast<T*>(a);
  }
#endif
  return std::allocator<T>().allocate(n);
 }
 void deallocate(T* p, std::size_t n) {
#ifdef __linux__
  if(n * sizeof(T) >= H) {
   munmap(p, (n * sizeof(T) + P - 1) & ~(P - 1));
   return;
  }
#endif
  std::allocator<T>().deallocate(p, n);
 }
 template <class U> bool operator==(const HugeAlloc<U>&) const { return true; }
 template <class U> bool operator!=(const HugeAlloc<U>&) const { return false; }
};
using hvec= std::vector<u32, HugeAlloc<u32>>;
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
}  // namespace convolve_internal
template <class T> std::vector<T> convolve_naive(const std::vector<T>& a, const std::vector<T>& b) {
 if(a.empty() || b.empty()) return {};
 const int n= int(a.size()), m= int(b.size());
 if constexpr(convolve_internal::is_zmod_v<T>) {
  const auto x= convolve_internal::to_u32(a), y= convolve_internal::to_u32(b);
  std::vector<convolve_internal::u32> z(n + m - 1);
  convolve_internal::conv<T>(x.data(), n, y.data(), m, z.data(), n + m - 1);
  return convolve_internal::from_u32<T>(z);
 } else {
  std::vector<T> c(n + m - 1, T());
  for(int i= 0; i < n; ++i)
   for(int j= 0; j < m; ++j) c[i + j]+= a[i] * b[j];
  return c;
 }
}
namespace convolve_internal {
// out[k - lo] = Σ_{i+j=k} a[i] b[j] mod M (lo <= k < hi)。a と b は ZMod の列 (値は 2M 未満) で、M 未満に直してから掛ける。
// n >= 1、m >= 1、0 <= lo < hi <= n + m - 1。out は a、b と重ならないこと。
template <class Z> void conv_range(const Z* a, int n, const Z* b, int m, int lo, int hi, Z* out) {
 static_assert(std::is_standard_layout_v<Z> && sizeof(Z) == 4);
 if(std::min(n, m) <= 16) {
  if(n > m) std::swap(a, b), std::swap(n, m);
  for(int k= lo; k < hi; ++k) {
   u64 s= 0;
   for(int i= std::max(0, k - m + 1), e= std::min(n - 1, k); i <= e; ++i) s+= u64(a[i].val()) * b[k - i].val();
   out[k - lo]= Z(s);
  }
  return;
 }
 static thread_local std::vector<u32> ap, bp;
 ap.resize(n), bp.assign(m + 2 * PAD, 0);
 for(int i= 0; i < n; ++i) ap[i]= a[i].val();
 for(int j= 0; j < m; ++j) bp[PAD + j]= b[j].val();
 const __m256i R= _mm256_set1_epi64x(r32<Z>());
 u32* o= reinterpret_cast<u32*>(out);  // 核は M 未満の値を書くので、そのまま ZMod の値になる (raw と同じ)
 int k0= lo;
 for(; k0 + 32 <= hi; k0+= 32) conv32<Z>(ap.data(), n, bp.data(), m, k0, o + (k0 - lo), R);
 if(k0 < hi) {
  u32 t[32];
  conv32<Z>(ap.data(), n, bp.data(), m, k0, t, R);
  std::copy(t, t + (hi - k0), o + (k0 - lo));
 }
}
}  // namespace convolve_internal
// 素朴な積の核。neo/base/core/convolve/convolve.hpp の Conv に渡して、convolve<Naive>(a, b) や convolve_range<Naive>(a, b, lo, hi) の
// ように使う。K が ZMod なら AVX2 の核で、ほかの型 (足し算と掛け算を持ち、K() が 0 の型) は二重ループで求める。
// 既定の核にはしない (default_kernel を宣言しない)。核を include し忘れたときに、黙って O(nm) で動かないようにするため。
struct Naive {
 // out[k - lo] = Σ_{i+j=k} a[i] b[j] (lo <= k < hi)。n >= 1、m >= 1、0 <= lo < hi <= n + m - 1。out は a、b と重ならないこと。
 template <class K> static void convolve_range(const K* a, int n, const K* b, int m, int lo, int hi, K* out) {
  if constexpr(convolve_internal::is_zmod_v<K>) convolve_internal::conv_range(a, n, b, m, lo, hi, out);
  else
   for(int k= lo; k < hi; ++k) {
    K s= K();
    for(int i= std::max(0, k - m + 1), e= std::min(n - 1, k); i <= e; ++i) s+= a[i] * b[k - i];
    out[k - lo]= s;
   }
 }
};
