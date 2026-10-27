#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <new>
#include <numeric>
#include <type_traits>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 列 a について、a[p[0]] <= a[p[1]] <= ... となる添字の列 p を返す。等しい値は添字の小さい順 (安定な並べ方)。
// 4 byte の整数は、値と添字を u64 に詰めて値の上の桁で振り分け、塊の中を u32 の鍵 (桁より下の bit を上に、塊の中の位置を下の 6 bit
// に詰めたもの) の AVX2 のネットワークで並べ、鍵の下の 6 bit で詰めた u64 を引いて添字を書く。振り分けは元の順を保つので、塊の中の
// 位置の順は添字の順と同じ。桁より下の bit が 26 bit 以下になるように桁を 6 bit 以上にする。AMD などでは、桁を塊の平均が 16 個以下に
// なる幅 (16 bit まで) にして振り分けを 1 回で済ませる。Intel の CPU では書き先が 65536 か所に散らばる振り分けが重い (10^6 個で
// 1.3〜1.45 倍) ので、上の 8 bit で振り分けてから、1 段目の塊ごとに次の桁で L1 に載る小さい配列へ振り分けて並べる。64 個を超えた
// 塊は std::sort。符号付きの値は最上位 bit を反転した鍵で並べる。作業用の配列と返す列は、2^19 byte 以上なら Linux では huge page を
// 頼む。ほかの型は、添字の列を std::stable_sort で並べる。procon-judge の self-sort-argsort-u32 の pack_msd16_net32_ai と
// pack_msd8x2_net32 を写したもの。
namespace argsort_internal {
using u32= unsigned;
using u64= unsigned long long;
// n 個の T の作業用の領域。0 では埋めない。2^19 byte 以上なら、Linux では 2 MB 境界の mmap に置いて huge page を頼み、触る前に
// まとめて用意させる (MADV_POPULATE_WRITE)。
template <class T> struct Buf {
 T* p;
 std::size_t mapped= 0;
 explicit Buf(std::size_t n) {
#ifdef __linux__
  constexpr std::size_t H= std::size_t(1) << 21;
  if(n * sizeof(T) >= (std::size_t(1) << 19)) {
   const std::size_t bytes= (n * sizeof(T) + H - 1) & ~(H - 1);
   void* m= mmap(nullptr, bytes + H, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
   if(m != MAP_FAILED) {
    char *q= static_cast<char*>(m), *a= reinterpret_cast<char*>((reinterpret_cast<std::uintptr_t>(q) + H - 1) & ~(H - 1));
    if(a != q) munmap(q, a - q);
    munmap(a + bytes, q + H - a);
    madvise(a, bytes, MADV_HUGEPAGE);
    madvise(a, bytes, 23);  // MADV_POPULATE_WRITE (Linux 5.14 から。古い kernel では何もしない)
    p= reinterpret_cast<T*>(a), mapped= bytes;
    return;
   }
  }
#endif
  p= static_cast<T*>(::operator new(n * sizeof(T)));
 }
 ~Buf() {
#ifdef __linux__
  if(mapped) {
   munmap(p, mapped);
   return;
  }
#endif
  ::operator delete(p);
 }
 Buf(const Buf&)= delete;
 Buf& operator=(const Buf&)= delete;
};
// Intel の CPU か。
inline bool is_intel() {
#if defined(__GNUC__) && defined(__x86_64__)
 static const bool v= [] {
  __builtin_cpu_init();
  return bool(__builtin_cpu_is("intel"));
 }();
 return v;
#else
 return false;
#endif
}
// 長さ n の返す列。2^19 byte 以上なら、確保だけしてから 2 MB 境界に揃った中の部分に huge page を頼み、ページの境界に揃った中の部分を
// まとめて用意させてから 0 で埋める (返す vector は std::allocator で確保するしかないので、この順にする)。
template <class I> std::vector<I> make_result(std::size_t n) {
 std::vector<I> p;
#ifdef __linux__
 if(n * sizeof(I) >= (std::size_t(1) << 19)) {
  p.reserve(n);
  constexpr std::uintptr_t H= std::uintptr_t(1) << 21, P= 4096;
  const std::uintptr_t s= reinterpret_cast<std::uintptr_t>(p.data()), e= s + n * sizeof(I);
  const std::uintptr_t hs= (s + H - 1) & ~(H - 1), he= e & ~(H - 1);
  if(hs < he) madvise(reinterpret_cast<void*>(hs), he - hs, MADV_HUGEPAGE);
  const std::uintptr_t ps= (s + P - 1) & ~(P - 1), pe= e & ~(P - 1);
  if(ps < pe) madvise(reinterpret_cast<void*>(ps), pe - ps, 23);  // MADV_POPULATE_WRITE
 }
#endif
 p.resize(n);
 return p;
}
// u32 の bitonic のソーティングネットワーク。stepD<MASK> は距離 D の lane どうしを比べ、MASK の bit が立つ lane に大きい方を置く。
namespace net {
[[gnu::always_inline]] inline __m256i vmin(__m256i a, __m256i b) { return _mm256_min_epu32(a, b); }
[[gnu::always_inline]] inline __m256i vmax(__m256i a, __m256i b) { return _mm256_max_epu32(a, b); }
template <int MASK> [[gnu::always_inline]] inline __m256i step1(__m256i v) {
 const __m256i q= _mm256_shuffle_epi32(v, 0xB1);
 return _mm256_blend_epi32(vmin(v, q), vmax(v, q), MASK);
}
template <int MASK> [[gnu::always_inline]] inline __m256i step2(__m256i v) {
 const __m256i q= _mm256_shuffle_epi32(v, 0x4E);
 return _mm256_blend_epi32(vmin(v, q), vmax(v, q), MASK);
}
template <int MASK> [[gnu::always_inline]] inline __m256i step4(__m256i v) {
 const __m256i q= _mm256_permute4x64_epi64(v, 0x4E);
 return _mm256_blend_epi32(vmin(v, q), vmax(v, q), MASK);
}
// bitonic な 8 個を昇順にする。
[[gnu::always_inline]] inline __m256i merge8(__m256i v) { return step1<0xAA>(step2<0xCC>(step4<0xF0>(v))); }
// 任意の 8 個を昇順にする。2 個ずつ昇順と降順、4 個ずつ昇順と降順に並べて bitonic にしてから merge8。
[[gnu::always_inline]] inline __m256i sort8(__m256i v) { return merge8(step1<0x5A>(step2<0x3C>(step1<0x66>(v)))); }
[[gnu::always_inline]] inline __m256i reverse8(__m256i v) { return _mm256_permutevar8x32_epi32(v, _mm256_setr_epi32(7, 6, 5, 4, 3, 2, 1, 0)); }
// a, b を並べた 16 個が bitonic なら昇順にする。
[[gnu::always_inline]] inline void merge16(__m256i& a, __m256i& b) {
 const __m256i l= vmin(a, b), h= vmax(a, b);
 a= merge8(l), b= merge8(h);
}
[[gnu::always_inline]] inline void sort16(__m256i& a, __m256i& b) {
 a= sort8(a), b= reverse8(sort8(b));
 merge16(a, b);
}
[[gnu::always_inline]] inline void sort32(__m256i& a, __m256i& b, __m256i& c, __m256i& d) {
 sort16(a, b), sort16(c, d);
 // 後ろの 16 個を逆順にして (reverse8(d), reverse8(c)) とつなぐと 32 個が bitonic になる。
 const __m256i rc= reverse8(d), rd= reverse8(c);
 __m256i l0= vmin(a, rc), h0= vmax(a, rc), l1= vmin(b, rd), h1= vmax(b, rd);
 merge16(l0, l1), merge16(h0, h1);
 a= l0, b= l1, c= h0, d= h1;
}
// a, b, c, d を並べた 32 個が bitonic なら昇順にする。
[[gnu::always_inline]] inline void merge32(__m256i& a, __m256i& b, __m256i& c, __m256i& d) {
 __m256i l0= vmin(a, c), h0= vmax(a, c), l1= vmin(b, d), h1= vmax(b, d);
 merge16(l0, l1), merge16(h0, h1);
 a= l0, b= l1, c= h0, d= h1;
}
// 8 本のベクタの 64 個を昇順にする。前後の 32 個を並べ、後ろを逆順にしてつないだ bitonic な 64 個を merge する。
[[gnu::always_inline]] inline void sort64(__m256i* v) {
 sort32(v[0], v[1], v[2], v[3]), sort32(v[4], v[5], v[6], v[7]);
 const __m256i r4= reverse8(v[7]), r5= reverse8(v[6]), r6= reverse8(v[5]), r7= reverse8(v[4]);
 __m256i l0= vmin(v[0], r4), h0= vmax(v[0], r4), l1= vmin(v[1], r5), h1= vmax(v[1], r5);
 __m256i l2= vmin(v[2], r6), h2= vmax(v[2], r6), l3= vmin(v[3], r7), h3= vmax(v[3], r7);
 merge32(l0, l1, l2, l3), merge32(h0, h1, h2, h3);
 v[0]= l0, v[1]= l1, v[2]= l2, v[3]= l3, v[4]= h0, v[5]= h1, v[6]= h2, v[7]= h3;
}
}  // namespace net
// 塊 s[0, m) の 8t 番目からの 8 個の鍵を作る。詰めた u64 を 4 個ずつ 2 回読み、上位 32 bit (値の鍵) を元の順に集めてから、鍵の bit
// [lo, lo + r) を上に、塊の中の位置を下の 6 bit に詰める。m 個より後ろの lane は最大値で埋める (読むのは s[0, 8t + 8) で、塊の後ろの
// 要素を読んでも鍵には使わない)。
[[gnu::always_inline]] inline __m256i make_keys(const u64* s, int m, int t, __m128i lo, __m256i rmask) {
 const int base= 8 * t;
 const __m256 v0= _mm256_castsi256_ps(_mm256_loadu_si256((const __m256i*)(s + base)));
 const __m256 v1= _mm256_castsi256_ps(_mm256_loadu_si256((const __m256i*)(s + base + 4)));
 // shuffle_ps で 128 bit の中ごとに奇数番目の 32 bit を集めると (x0 x1 x4 x5 | x2 x3 x6 x7) になるので、64 bit 単位で並べ直す。
 const __m256i x= _mm256_permute4x64_epi64(_mm256_castps_si256(_mm256_shuffle_ps(v0, v1, 0xDD)), 0xD8);
 const __m256i iota= _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
 const __m256i k= _mm256_or_si256(_mm256_slli_epi32(_mm256_and_si256(_mm256_srl_epi32(x, lo), rmask), 6), _mm256_add_epi32(iota, _mm256_set1_epi32(base)));
 if(m - base >= 8) return k;
 return _mm256_blendv_epi8(_mm256_set1_epi32(-1), k, _mm256_cmpgt_epi32(_mm256_set1_epi32(m - base), iota));
}
// 塊 s[0, m) (2 <= m <= 64) を値の順に並べ、添字を d[0, m) に書く。
template <class I> [[gnu::always_inline]] inline void sort_bucket(const u64* s, int m, __m128i lo, __m256i rmask, I* d) {
 alignas(32) u32 k[64];
 const __m256i pad= _mm256_set1_epi32(-1);
 if(m <= 8) {
  _mm256_store_si256((__m256i*)k, net::sort8(make_keys(s, m, 0, lo, rmask)));
 } else if(m <= 16) {
  __m256i x= make_keys(s, m, 0, lo, rmask), y= make_keys(s, m, 1, lo, rmask);
  net::sort16(x, y);
  _mm256_store_si256((__m256i*)k, x), _mm256_store_si256((__m256i*)(k + 8), y);
 } else if(m <= 32) {
  __m256i v[4];
  for(int t= 0; t < 4; ++t) v[t]= 8 * t < m ? make_keys(s, m, t, lo, rmask) : pad;
  net::sort32(v[0], v[1], v[2], v[3]);
  for(int t= 0; t < 4; ++t) _mm256_store_si256((__m256i*)(k + 8 * t), v[t]);
 } else {
  __m256i v[8];
  for(int t= 0; t < 8; ++t) v[t]= 8 * t < m ? make_keys(s, m, t, lo, rmask) : pad;
  net::sort64(v);
  for(int t= 0; t < 8; ++t) _mm256_store_si256((__m256i*)(k + 8 * t), v[t]);
 }
 for(int t= 0; t < m; ++t) d[t]= I(u32(s[k[t] & 63]));
}
// Intel の道の 2 段目: 1 段目の塊 s[0, m) (m > 64、鍵の bit [lo, lo + r) より上はそろっている) を並べ、添字を d[0, m) に書く。
// t は m + 64 個の作業用の領域。桁は塊の平均が 16 個以下になる幅 (11 bit まで) にする。
template <class I> inline void sort_level2(const u64* s, std::size_t m, int lo, int r, u64* t, I* d) {
 int db= 1;
 while(db < 11 && (m >> db) > 16) ++db;
 db= std::min(db, r);
 const int sh= lo + r - db, r2= r - db;
 const u32 mk= (1u << db) - 1;
 u32 pos[2048];
 std::fill(pos, pos + mk + 1, 0u);
 for(std::size_t j= 0; j < m; ++j) ++pos[u32(s[j] >> 32) >> sh & mk];
 for(u32 c= 0, acc= 0; c <= mk; ++c) {
  const u32 x= pos[c];
  pos[c]= acc, acc+= x;
 }
 for(std::size_t j= 0; j < m; ++j) {
  const u64 y= s[j];
  t[pos[u32(y >> 32) >> sh & mk]++]= y;
 }
 // 振り分けたあとの pos[c] は小さい塊 c の終わり。
 const __m128i lov= _mm_cvtsi32_si128(lo);
 const __m256i rmask= _mm256_set1_epi32(int((1u << r2) - 1));
 for(u32 c= 0, off= 0; c <= mk; ++c) {
  const u32 end= pos[c], k= end - off;
  if(k <= 1 || !r2) {  // 1 個以下か、小さい塊の中が全部同じ値
   for(u32 j= off; j < end; ++j) d[j]= I(u32(t[j]));
  } else if(k <= 64) {
   sort_bucket(t + off, int(k), lov, rmask, d + off);
  } else {
   if(!std::is_sorted(t + off, t + end)) std::sort(t + off, t + end);
   for(u32 j= off; j < end; ++j) d[j]= I(u32(t[j]));
  }
  off= end;
 }
}
// 値の鍵の bit [lo, hi] が要素によって違いうる (hi > lo か hi == lo)。桁を塊の平均が 16 個以下になる幅 (6 bit から 16 bit) にして
// 1 回振り分け、添字を q[0, n) に書く。6 bit 以上にするのは、桁より下の bit を 26 bit 以下にして、塊の中の位置 (6 bit) と u32 に
// 詰めるため。
template <class I, class T> inline void by_one(const T* a, std::size_t n, int hi, int lo, I* q) {
 constexpr u32 F= std::is_signed_v<T> ? 0x80000000u : 0;
 int db= 6;
 while(db < 16 && (n >> db) > 16) ++db;
 db= std::min(db, hi - lo + 1);
 const int sh= hi - db + 1, r= sh - lo;
 const u32 mk= (1u << db) - 1;
 std::vector<u32> pos(std::size_t(1) << db);
 for(std::size_t i= 0; i < n; ++i) ++pos[(u32(a[i]) ^ F) >> sh & mk];
 for(std::size_t d= 0, s= 0; d < pos.size(); ++d) {
  const u32 t= pos[d];
  pos[d]= u32(s), s+= t;
 }
 if(!r) {  // 桁より下の bit が全要素で同じ: 振り分け 1 回で添字を直接書く
  for(std::size_t i= 0; i < n; ++i) q[pos[(u32(a[i]) ^ F) >> sh & mk]++]= I(i);
  return;
 }
 Buf<u64> bb(n + 64);  // 塊を 8 個ずつ読むので、最後の塊の後ろを 63 個まで読む
 u64* b= bb.p;
 for(std::size_t i= 0; i < n; ++i) {
  const u32 x= u32(a[i]) ^ F;
  b[pos[x >> sh & mk]++]= u64(x) << 32 | i;
 }
 // 振り分けたあとの pos[d] は塊 d の終わり。
 const __m128i lov= _mm_cvtsi32_si128(lo);
 const __m256i rmask= _mm256_set1_epi32(int((1u << r) - 1));
 for(std::size_t d= 0, off= 0; d < pos.size(); ++d) {
  const std::size_t end= pos[d], m= end - off;
  if(m <= 1) {
   if(m) q[off]= I(u32(b[off]));
  } else if(m <= 64) {
   sort_bucket(b + off, int(m), lov, rmask, q + off);
  } else {
   if(!std::is_sorted(b + off, b + end)) std::sort(b + off, b + end);
   for(std::size_t t= off; t < end; ++t) q[t]= I(u32(b[t]));
  }
  off= end;
 }
}
// Intel の道。上の 8 bit で振り分け、64 個を超えた 1 段目の塊は sort_level2 で並べる。
template <class I, class T> inline void by_two(const T* a, std::size_t n, int hi, int lo, I* q) {
 constexpr u32 F= std::is_signed_v<T> ? 0x80000000u : 0;
 const int db= std::min(8, hi - lo + 1), sh= hi - db + 1, r= sh - lo;  // r は 1 段目の桁より下で、要素によって違いうる bit の数 (24 bit 以下)
 const u32 mk= (1u << db) - 1;
 u32 pos[256]= {};
 for(std::size_t i= 0; i < n; ++i) ++pos[(u32(a[i]) ^ F) >> sh & mk];
 u32 maxm= 0;
 for(u32 c= 0, acc= 0; c <= mk; ++c) {
  const u32 x= pos[c];
  pos[c]= acc, acc+= x, maxm= std::max(maxm, x);
 }
 if(!r) {  // 桁より下の bit が全要素で同じ: 振り分け 1 回で添字を直接書く
  for(std::size_t i= 0; i < n; ++i) q[pos[(u32(a[i]) ^ F) >> sh & mk]++]= I(i);
  return;
 }
 Buf<u64> bb(n + 64);  // 塊を 8 個ずつ読むので、最後の塊の後ろを 63 個まで読む
 u64* b= bb.p;
 for(std::size_t i= 0; i < n; ++i) {
  const u32 x= u32(a[i]) ^ F;
  b[pos[x >> sh & mk]++]= u64(x) << 32 | i;
 }
 // 振り分けたあとの pos[c] は塊 c の終わり。
 Buf<u64> tmp(std::size_t(maxm) + 64);
 const __m128i lov= _mm_cvtsi32_si128(lo);
 const __m256i rmask= _mm256_set1_epi32(int((1u << r) - 1));
 for(u32 c= 0, off= 0; c <= mk; ++c) {
  const u32 end= pos[c], m= end - off;
  if(m <= 1) {
   if(m) q[off]= I(u32(b[off]));
  } else if(m <= 64) {
   sort_bucket(b + off, int(m), lov, rmask, q + off);
  } else {
   sort_level2(b + off, m, lo, r, tmp.p, q + off);
  }
  off= end;
 }
}
// 4 byte の列。two は 8 bit ずつ 2 段で振り分けるか (Intel の道)。
template <class I, class T> inline std::vector<I> argsort4(const T* a, std::size_t n, bool two) {
 constexpr u32 F= std::is_signed_v<T> ? 0x80000000u : 0;
 if(n <= 64) {
  std::vector<I> p(n);
  u64 b[64];
  for(std::size_t i= 0; i < n; ++i) b[i]= u64(u32(a[i]) ^ F) << 32 | i;
  std::sort(b, b + n);
  for(std::size_t i= 0; i < n; ++i) p[i]= I(u32(b[i]));
  return p;
 }
 std::vector<I> p= make_result<I>(n);
 I* q= p.data();
 u32 o= 0, e= ~0u;
 for(std::size_t i= 0; i < n; ++i) {
  const u32 k= u32(a[i]) ^ F;
  o|= k, e&= k;
 }
 const u32 diff= o ^ e;
 if(!diff) {  // 全要素が同じ値
  for(std::size_t i= 0; i < n; ++i) q[i]= I(i);
  return p;
 }
 const int hi= 31 - __builtin_clz(diff), lo= __builtin_ctz(diff);
 if(two) by_two(a, n, hi, lo, q);
 else by_one(a, n, hi, lo, q);
 return p;
}
}  // namespace argsort_internal
// a[p[0]] <= a[p[1]] <= ... となる添字の列 p を返す。等しい値は添字の小さい順。I は添字の型。
template <class I= int, class T> std::vector<I> argsort(const std::vector<T>& a) {
 if constexpr(std::is_integral_v<T> && !std::is_same_v<T, bool> && sizeof(T) == 4) {
  return argsort_internal::argsort4<I>(a.data(), a.size(), argsort_internal::is_intel());
 } else {
  std::vector<I> p(a.size());
  std::iota(p.begin(), p.end(), I(0));
  std::stable_sort(p.begin(), p.end(), [&](I i, I j) { return a[i] < a[j]; });
  return p;
 }
}
