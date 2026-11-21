#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <bit>
#include <cstdint>
#include <memory>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 静的な列の区間最小。列を 16 個ずつの塊に分け、塊の中で左から累積した最小 pre と右から累積した最小 suf を持ち、
// 塊ごとの最小の列には sparse table を載せる。[l, r) の両端が別の塊なら suf[l]、pre[r - 1]、間の塊の sparse table の
// min、同じ塊なら塊の 16 個 (キャッシュの 1 行) を AVX2 の 2 本で読んで [l, r) の外を最大値で埋めた min を返す。
// 表は列の写しと pre と suf が N 個ずつと、塊の sparse table が N / 16 * log(N / 16) 個で、N = 5 * 10^5 で 8 MB ほど。
// 写し、pre、suf、塊の最小は塊ごとの 1 回の走査で作り、累積は AVX2 のレーンのずらしと min で求める。表は 2 MB 境界に
// 置いて huge page を頼む (Linux のみ)。procon-judge の yosupo-staticrmq で書き比べた block16_fused を写したもの。
namespace range_min_internal {
using u32= unsigned;
// u32 の n 個分の領域。1 MB 以上なら 2 MB 境界に寄せ、Linux では huge page を頼んで触る前にまとめて用意させる。
// 小さい表で 2 MB を 0 で埋めると、小さい列をたくさん作る使い方で遅くなるので、1 MB 未満は 64 byte 境界に寄せるだけ。
class Buffer {
 static constexpr std::size_t H= std::size_t(1) << 21, L= 64;
 std::unique_ptr<char[]> raw;
 u32* p= nullptr;
public:
 Buffer()= default;
 explicit Buffer(std::size_t n) {
  if(!n) return;
  const std::size_t want= n * sizeof(u32);
  const std::size_t a= want < (std::size_t(1) << 20) ? L : H, bytes= (want + a - 1) & ~(a - 1);
  raw.reset(new char[bytes + a]);
  p= reinterpret_cast<u32*>((reinterpret_cast<std::uintptr_t>(raw.get()) + a - 1) & ~std::uintptr_t(a - 1));
#ifdef __linux__
  if(a == H) {
   madvise(p, bytes, MADV_HUGEPAGE);
   madvise(p, bytes, 23);  // MADV_POPULATE_WRITE (Linux 5.14+)。古い kernel では何もしない
  }
#endif
 }
 u32* get() const { return p; }
};
// 8 個の左からの累積の最小。128 bit の半分ごとに 1 個と 2 個ずらして min を取り、半分の境を越える分は端を配る。
inline __m256i prefix8(__m256i x, __m256i ones) {
 x= _mm256_min_epu32(x, _mm256_alignr_epi8(x, ones, 12));
 x= _mm256_min_epu32(x, _mm256_alignr_epi8(x, ones, 8));
 return _mm256_min_epu32(x, _mm256_permute2x128_si256(_mm256_shuffle_epi32(x, 0xff), ones, 0x02));
}
// 8 個の右からの累積の最小。
inline __m256i suffix8(__m256i x, __m256i ones) {
 x= _mm256_min_epu32(x, _mm256_alignr_epi8(ones, x, 4));
 x= _mm256_min_epu32(x, _mm256_alignr_epi8(ones, x, 8));
 return _mm256_min_epu32(x, _mm256_permute2x128_si256(_mm256_shuffle_epi32(x, 0x00), ones, 0x31));
}
}
class RangeMin {
 using u32= unsigned;
 static constexpr int LB= 4, B= 1 << LB;
 int n= 0, nb= 0, lg= 0;
 range_min_internal::Buffer buf;
 u32 *v= nullptr, *pre= nullptr, *suf= nullptr, *st= nullptr;
 // 塊 b の 16 個 (x0, x1) から、写し、pre、suf、塊の最小を書く。
 void block(int b, __m256i x0, __m256i x1) {
  using namespace range_min_internal;
  const __m256i ones= _mm256_set1_epi32(-1), seven= _mm256_set1_epi32(7), zero= _mm256_setzero_si256();
  const std::size_t o= std::size_t(b) << LB;
  _mm256_store_si256((__m256i*)(v + o), x0), _mm256_store_si256((__m256i*)(v + o + 8), x1);
  const __m256i p0= prefix8(x0, ones), p1= _mm256_min_epu32(prefix8(x1, ones), _mm256_permutevar8x32_epi32(p0, seven));
  const __m256i s1= suffix8(x1, ones), s0= _mm256_min_epu32(suffix8(x0, ones), _mm256_permutevar8x32_epi32(s1, zero));
  _mm256_store_si256((__m256i*)(pre + o), p0), _mm256_store_si256((__m256i*)(pre + o + 8), p1);
  _mm256_store_si256((__m256i*)(suf + o), s0), _mm256_store_si256((__m256i*)(suf + o + 8), s1);
  st[b]= u32(_mm256_cvtsi256_si32(s0));
 }
 // 塊の [x, y) の最小。x < y。
 u32 blocks(int x, int y) const {
  const int k= std::bit_width(unsigned(y - x)) - 1;
  const u32* row= st + std::size_t(k) * nb;
  return std::min(row[x], row[y - (1 << k)]);
 }
 // 塊 b の lane [lo, hi] の最小。
 u32 inblock(int b, int lo, int hi) const {
  const __m256i* p= (const __m256i*)(v + (std::size_t(b) << LB));
  const __m256i vlo= _mm256_set1_epi32(lo), vhi= _mm256_set1_epi32(hi), eight= _mm256_set1_epi32(8);
  __m256i idx= _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7), m= _mm256_set1_epi32(-1);
  for(int t= 0; t < B / 8; ++t) {
   const __m256i out= _mm256_or_si256(_mm256_cmpgt_epi32(vlo, idx), _mm256_cmpgt_epi32(idx, vhi));
   m= _mm256_min_epu32(m, _mm256_or_si256(_mm256_load_si256(p + t), out));
   idx= _mm256_add_epi32(idx, eight);
  }
  __m128i h= _mm_min_epu32(_mm256_castsi256_si128(m), _mm256_extracti128_si256(m, 1));
  h= _mm_min_epu32(h, _mm_shuffle_epi32(h, 0x4e));
  h= _mm_min_epu32(h, _mm_shuffle_epi32(h, 0xb1));
  return u32(_mm_cvtsi128_si32(h));
 }
public:
 RangeMin()= default;
 explicit RangeMin(const std::vector<u32>& a): n(a.size()), nb((n + B - 1) >> LB), lg(nb ? std::bit_width(unsigned(nb)) - 1 : 0), buf((std::size_t(nb) << LB) * 3 + std::size_t(lg + 1) * nb) {
  if(!n) return;
  const std::size_t m= std::size_t(nb) << LB;
  v= buf.get(), pre= v + m, suf= pre + m, st= suf + m;
  const u32* src= a.data();
  const int full= n >> LB;
  for(int b= 0; b < full; ++b) block(b, _mm256_loadu_si256((const __m256i*)(src + (std::size_t(b) << LB))), _mm256_loadu_si256((const __m256i*)(src + (std::size_t(b) << LB) + 8)));
  if(full < nb) {
   alignas(32) u32 t[B];
   std::fill(t, t + B, ~0u);
   std::copy(src + (std::size_t(full) << LB), src + n, t);
   block(full, _mm256_load_si256((const __m256i*)t), _mm256_load_si256((const __m256i*)(t + 8)));
  }
  for(int k= 0; k < lg; ++k) {
   const u32* s= st + std::size_t(k) * nb;
   u32* d= st + std::size_t(k + 1) * nb;
   const int h= 1 << k, len= nb - 2 * h + 1;
   int i= 0;
   for(; i + 8 <= len; i+= 8) _mm256_storeu_si256((__m256i*)(d + i), _mm256_min_epu32(_mm256_loadu_si256((const __m256i*)(s + i)), _mm256_loadu_si256((const __m256i*)(s + i + h))));
   for(; i < len; ++i) d[i]= std::min(s[i], s[i + h]);
  }
 }
 int size() const { return n; }
 // min(a[l], ..., a[r - 1])。0 <= l < r <= size()。
 u32 min(int l, int r) const {
  const int j= r - 1, bl= l >> LB, bj= j >> LB;
  if(bl == bj) return inblock(bl, l & (B - 1), j & (B - 1));
  u32 x= std::min(suf[l], pre[j]);
  if(bj - bl > 1) x= std::min(x, blocks(bl + 1, bj));
  return x;
 }
};
