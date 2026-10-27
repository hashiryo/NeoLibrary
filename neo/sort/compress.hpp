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
#include <type_traits>
#include <utility>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 整数の列 a の各要素を順位 (a の中でその要素より小さい値の種類数) に置き換え、重複を除いて昇順に並べた値の列 xs を返す。返したあとは
// xs[a[i]] が元の a[i] になる。値と添字を値の上の桁で振り分け、塊の中を u32 の鍵 (桁より下の bit の上の 26 bit までを上に、塊の中の
// 位置を下の 6 bit に詰めたもの) の AVX2 のネットワークで並べたら、並べた順に順位を数えて a の添字の位置に書く。塊は値の順に並び、
// 違う塊の値は必ず違うので、順位は塊をまたいで数え続ければよい。
//   4 byte: 値と添字を u64 に詰めて振り分ける。桁は要素によって違う bit の上から、塊の平均が 16 個以下になる幅 (6 bit から 16 bit)。
//     10^6 個ほどの列は、間をあけて取った 1024 個から桁の位置を見込み、全要素の bit or と bit and を取る走査でヒストグラムも作る。
//   8 byte: 値と添字を別の配列の同じ位置に振り分ける (合わせて u64 に収まらないため)。桁は、要素によって違う bit の幅と最大値と
//     最小値の差の幅の狭い方の上から選ぶ。26 bit が等しい隣どうしがあれば、そこだけ挿入ソートで値の順に直す。
// どちらも、桁より下の bit が全要素で同じなら振り分けずに塊の番号から順位を引く表で書き換え、値の種類が 1024 個以下ならハッシュ表で
// 順位を引く。64 個を超えた塊は std::sort。符号付きの値は最上位 bit を反転した鍵で並べる。作業用の配列と返す xs は、2^19 byte 以上なら
// Linux では huge page を頼む。ほかの型は、写しを並べて重複を除き、各要素を二分探索で引く。procon-judge の self-sort-compress-u32 の
// pack_msd16_net32_fused_all と、self-sort-compress-i64 の msd16_net32_fused_ph を写したもの。
namespace compress_internal {
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
// 返す xs を、長さ 0、容量 n で作る。2^19 byte 以上なら、確保だけしてから 2 MB 境界に揃った中の部分に huge page を頼み、ページの境界に
// 揃った中の部分をまとめて用意させる。
template <class T> std::vector<T> alloc_xs(std::size_t n) {
 std::vector<T> xs;
 xs.reserve(n);
#ifdef __linux__
 if(n * sizeof(T) >= (std::size_t(1) << 19)) {
  constexpr std::uintptr_t H= std::uintptr_t(1) << 21, P= 4096;
  const std::uintptr_t s= reinterpret_cast<std::uintptr_t>(xs.data()), e= s + n * sizeof(T);
  const std::uintptr_t hs= (s + H - 1) & ~(H - 1), he= e & ~(H - 1);
  if(hs < he) madvise(reinterpret_cast<void*>(hs), he - hs, MADV_HUGEPAGE);
  const std::uintptr_t ps= (s + P - 1) & ~(P - 1), pe= e & ~(P - 1);
  if(ps < pe) madvise(reinterpret_cast<void*>(ps), pe - ps, 23);  // MADV_POPULATE_WRITE
 }
#endif
 return xs;
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
// 4 byte の値。値の鍵 (符号付きなら最上位 bit を反転したもの) を上位 32 bit、添字を下位 32 bit に詰めた u64 を振り分ける。
namespace c4 {
template <class T> constexpr u32 F= std::is_signed_v<T> ? 0x80000000u : 0;
// 塊 s[0, m) の 8t 番目からの 8 個の鍵を作る。詰めた u64 を 4 個ずつ 2 回読み、上位 32 bit を元の順に集めてから、値の鍵の bit
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
// 塊 s[0, m) (2 <= m <= 64) の鍵を並べて k[0, m) に書く。
[[gnu::always_inline]] inline void sort_keys(const u64* s, int m, __m128i lo, __m256i rmask, u32* k) {
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
}
// 値の順に並んだ詰めた u64 の列 s[0, n) (n >= 1) を先頭から走査する。値が変わるたびに順位を 1 増やして xs に足し、a の添字の位置に
// 順位を書く。
template <class T> inline std::vector<T> scan(const u64* s, std::size_t n, T* a) {
 std::vector<T> xs;
 xs.reserve(n);
 u32 prev= u32(s[0] >> 32), r= 0;
 xs.push_back(T(prev ^ F<T>));
 for(std::size_t i= 0; i < n; ++i) {
  const u64 x= s[i];
  const u32 v= u32(x >> 32);
  if(v != prev) xs.push_back(T(v ^ F<T>)), prev= v, ++r;
  a[u32(x)]= T(r);
 }
 return xs;
}
// 塊 s[0, m) (2 <= m <= 64) を鍵のネットワークで並べ、並べた順に順位を数えて a の添字の位置に書き、新しい値を xs に足す。rk は
// それまでの値の種類数 (次に付ける順位)。塊の中では桁とそれより上の bit、lo より下の bit がそろっているので、鍵の上の 26 bit 以下
// (値の鍵の bit [lo, sh)) が等しいことと値が等しいことは同じ。
template <class T> [[gnu::always_inline]] inline void rank_bucket(const u64* s, int m, __m128i lo, __m256i rmask, T* a, std::vector<T>& xs, u32& rk) {
 alignas(32) u32 k[64];
 sort_keys(s, m, lo, rmask, k);
 u32 prev= ~0u;  // 鍵の上の bit は 26 bit 以下なので、塊の最初の要素は必ず新しい値になる
 for(int t= 0; t < m; ++t) {
  const u32 v= k[t] >> 6;
  const u64 x= s[k[t] & 63];
  if(v != prev) xs.push_back(T(u32(x >> 32) ^ F<T>)), prev= v, ++rk;
  a[u32(x)]= T(rk - 1);
 }
}
// 値の種類が 1024 個以下なら、各要素をハッシュ表で引いた順位に書き換えて xs を作り、true を返す。種類が 1024 個を超えた時点で、
// a に手を付けずに false を返す。表は 4096 個の u64 (上位 32 bit が値、下位 32 bit が埋まった印の 1、あとで順位) の線形探査で、
// 消すことはないので、値 x を引くと、x の場所に着くまでに空の場所を通らない。ほとんどの入力ではすぐにやめる道なので、展開させない。
template <class T> [[gnu::noinline]] inline bool few_values(T* a, std::size_t n, std::vector<T>& xs) {
 constexpr u32 LIM= 1024, B= 12, S= 1u << B;
 std::vector<u64> tab(S);
 u64* t= tab.data();
 const auto slot= [](u32 x) { return (x * 0x9E3779B1u) >> (32 - B); };
 u32 k= 0;
 for(std::size_t i= 0; i < n; ++i) {
  const u32 x= u32(a[i]);
  const u64 want= u64(x) << 32 | 1;
  for(u32 h= slot(x);; h= (h + 1) & (S - 1)) {
   const u64 y= t[h];
   if(y == want) break;
   if(!y) {
    if(++k > LIM) return false;
    t[h]= want;
    break;
   }
  }
 }
 xs.reserve(k);
 for(u32 h= 0; h < S; ++h)
  if(t[h]) xs.push_back(T(u32(t[h] >> 32)));
 std::sort(xs.begin(), xs.end());
 for(u32 r= 0; r < k; ++r) {
  const u32 x= u32(xs[r]);
  u32 h= slot(x);
  while(u32(t[h] >> 32) != x) h= (h + 1) & (S - 1);
  t[h]= u64(x) << 32 | r;
 }
 for(std::size_t i= 0; i < n; ++i) {
  const u32 x= u32(a[i]);
  u32 h= slot(x);
  while(u32(t[h] >> 32) != x) h= (h + 1) & (S - 1);
  a[i]= T(u32(t[h]));
 }
 return true;
}
template <class T> inline std::vector<T> compress(T* a, std::size_t n) {
 constexpr u32 FL= F<T>;
 if(!n) return {};
 if(n <= 64) {
  u64 b[64];
  for(std::size_t i= 0; i < n; ++i) b[i]= u64(u32(a[i]) ^ FL) << 32 | i;
  std::sort(b, b + n);
  return scan(b, n, a);
 }
 // 桁の幅が 16 bit になる大きさなら、間をあけて取った 1024 個の bit or と bit and から桁の位置を見込み、全要素の bit or と bit and
 // を取る走査で、見込んだ桁のヒストグラムも作る。取らなかった要素にもっと上の bit があると見込みが外れるので、そのときは作り直す。
 int ssh= -1;  // 見込んだ桁の位置 (-1 は見込まない)
 if((n >> 15) > 16) {
  u32 so= 0, se= ~0u;
  const std::size_t step= n / 1024;
  for(std::size_t i= 0; i < n; i+= step) {
   const u32 x= u32(a[i]) ^ FL;
   so|= x, se&= x;
  }
  const u32 sd= so ^ se;
  if(sd && 32 - __builtin_clz(sd) - __builtin_ctz(sd) >= 16) ssh= 31 - __builtin_clz(sd) - 15;
 }
 u32 o= 0, e= ~0u;
 std::vector<u32> pos;
 if(ssh >= 0) {
  pos.assign(std::size_t(1) << 16, 0);
  u32* c= pos.data();
  for(std::size_t i= 0; i < n; ++i) {
   const u32 x= u32(a[i]) ^ FL;
   o|= x, e&= x, ++c[x >> ssh & 0xFFFF];
  }
 } else {
  for(std::size_t i= 0; i < n; ++i) {
   const u32 x= u32(a[i]) ^ FL;
   o|= x, e&= x;
  }
 }
 const u32 diff= o ^ e;
 if(!diff) {  // 全要素が同じ値
  std::vector<T> xs(1, a[0]);
  std::fill(a, a + n, T(0));
  return xs;
 }
 if(n >= (std::size_t(1) << 14)) {
  std::vector<T> xs;
  if(few_values(a, n, xs)) return xs;
 }
 const int hi= 31 - __builtin_clz(diff), lo= __builtin_ctz(diff), w= hi - lo + 1;
 // 6 bit 以上にするのは、桁より下の bit を 26 bit 以下にして、塊の中の位置 (6 bit) と u32 に詰めるため。
 int db= 6;
 while(db < 16 && (n >> db) > 16) ++db;
 db= std::min(db, w);
 const int sh= hi - db + 1, r= sh - lo;
 const u32 mk= (1u << db) - 1;
 if(sh != ssh || db != 16) {  // 見込まなかったか、見込みが外れた
  pos.assign(std::size_t(1) << db, 0);
  for(std::size_t i= 0; i < n; ++i) ++pos[(u32(a[i]) ^ FL) >> sh & mk];
 }
 if(!r) {
  // 桁が値を決めるので、空でない塊ごとに値が 1 つある。値の鍵は、桁の外の bit (全要素で同じ) に塊の番号を入れたもの。
  std::vector<T> xs;
  xs.reserve(std::min(n, pos.size()));
  const u32 base= e & ~(mk << sh);
  u32 rk= 0;
  for(std::size_t d= 0; d < pos.size(); ++d) {
   const u32 c= pos[d];
   pos[d]= rk;
   if(c) xs.push_back(T((base | u32(d) << sh) ^ FL)), ++rk;
  }
  for(std::size_t i= 0; i < n; ++i) a[i]= T(pos[(u32(a[i]) ^ FL) >> sh & mk]);
  return xs;
 }
 for(std::size_t d= 0, s= 0; d < pos.size(); ++d) {
  const u32 t= pos[d];
  pos[d]= u32(s), s+= t;
 }
 Buf<u64> bb(n + 64);  // 塊を 8 個ずつ読むので、最後の塊の後ろを 63 個まで読む
 u64* b= bb.p;
 for(std::size_t i= 0; i < n; ++i) {
  const u32 x= u32(a[i]) ^ FL;
  b[pos[x >> sh & mk]++]= u64(x) << 32 | i;
 }
 // 振り分けたあとの pos[d] は塊 d の終わり。
 const __m128i lov= _mm_cvtsi32_si128(lo);
 const __m256i rmask= _mm256_set1_epi32(int((1u << r) - 1));
 std::vector<T> xs= alloc_xs<T>(n);
 u32 rk= 0;
 for(std::size_t d= 0, off= 0; d < pos.size(); ++d) {
  const std::size_t end= pos[d], m= end - off;
  if(m == 1) {
   const u64 x= b[off];
   xs.push_back(T(u32(x >> 32) ^ FL));
   a[u32(x)]= T(rk++);
  } else if(m >= 2 && m <= 64) {
   rank_bucket(b + off, int(m), lov, rmask, a, xs, rk);
  } else if(m > 64) {
   if(!std::is_sorted(b + off, b + end)) std::sort(b + off, b + end);
   u32 prev= u32(b[off] >> 32);
   xs.push_back(T(prev ^ FL)), ++rk;
   for(std::size_t t= off; t < end; ++t) {
    const u64 x= b[t];
    const u32 v= u32(x >> 32);
    if(v != prev) xs.push_back(T(v ^ FL)), prev= v, ++rk;
    a[u32(x)]= T(rk - 1);
   }
  }
  off= end;
 }
 return xs;
}
}  // namespace c4
// 8 byte の値。値と添字を別の配列の同じ位置に振り分ける。
namespace c8 {
template <class T> constexpr u64 F= std::is_signed_v<T> ? u64(1) << 63 : 0;
// 桁を取り出すときの鍵。鍵の符号なしの順が値の順になるようにする。
template <class T> [[gnu::always_inline]] inline u64 key(T x) { return u64(x) ^ F<T>; }
// 塊の鍵の作り方。鍵 (key(x) - base) の bit [ks, ks + rk) を上に、塊の中の位置を下の 6 bit に詰める (rk <= 26)。
struct KeyParams {
 __m256i flip, base, mask;  // flip は符号付きのとき最上位 bit、mask は下の rk bit (64 bit の lane ごと)
 __m128i ks;
};
// 塊 s[0, m) の 8t 番目からの 8 個の鍵を作る。値を 4 個ずつ 2 回読み、64 bit のまま鍵の bit を取り出してから、下の 32 bit を元の順に
// 集める。m 個より後ろの lane は最大値で埋める (読むのは s[0, 8t + 8) で、塊の後ろの要素を読んでも鍵には使わない)。
template <class T> [[gnu::always_inline]] inline __m256i make_keys(const T* s, int m, int t, const KeyParams& kp) {
 const int at= 8 * t;
 const auto part= [&](const T* p) {
  const __m256i v= _mm256_sub_epi64(_mm256_xor_si256(_mm256_loadu_si256((const __m256i*)p), kp.flip), kp.base);
  return _mm256_castsi256_ps(_mm256_and_si256(_mm256_srl_epi64(v, kp.ks), kp.mask));
 };
 // shuffle_ps で 128 bit の中ごとに偶数番目の 32 bit を集めると (y0 y1 y4 y5 | y2 y3 y6 y7) になるので、64 bit 単位で並べ直す。
 const __m256i y= _mm256_permute4x64_epi64(_mm256_castps_si256(_mm256_shuffle_ps(part(s + at), part(s + at + 4), 0x88)), 0xD8);
 const __m256i iota= _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
 const __m256i k= _mm256_or_si256(_mm256_slli_epi32(y, 6), _mm256_add_epi32(iota, _mm256_set1_epi32(at)));
 if(m - at >= 8) return k;
 return _mm256_blendv_epi8(_mm256_set1_epi32(-1), k, _mm256_cmpgt_epi32(_mm256_set1_epi32(m - at), iota));
}
// 並べた鍵のベクタ v[0, nv) の先頭 m 個に、上の 26 bit が等しい隣どうしがあるか。1 lane ずらした鍵と比べる。
[[gnu::always_inline]] inline bool has_tie(const __m256i* v, int nv, int m) {
 const __m256i up= _mm256_setr_epi32(7, 0, 1, 2, 3, 4, 5, 6), last= _mm256_set1_epi32(7);
 __m256i carry= _mm256_set1_epi32(-1);  // 先頭の要素の前には、どの鍵とも等しくならない値を置く
 for(int t= 0; t < nv; ++t) {
  const int rem= m - 8 * t;
  if(rem <= 0) break;
  const __m256i p= _mm256_srli_epi32(v[t], 6);
  const __m256i q= _mm256_blend_epi32(_mm256_permutevar8x32_epi32(p, up), carry, 0x01);
  int bits= _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpeq_epi32(p, q)));
  if(rem < 8) bits&= (1 << rem) - 1;
  if(bits) return true;
  carry= _mm256_permutevar8x32_epi32(p, last);
 }
 return false;
}
// 並べた値 d と添字 di の鍵 k[0, m) の上の 26 bit が等しい並びごとに、挿入ソートで値の順に直す (添字も一緒に動かす)。
template <class T> inline void fix_ties(const u32* k, int m, T* d, u32* di) {
 for(int i= 0; i < m;) {
  int j= i + 1;
  while(j < m && (k[j] ^ k[i]) < 64) ++j;
  for(int x= i + 1; x < j; ++x) {
   const T y= d[x];
   const u32 yi= di[x];
   int z= x;
   for(; z > i && y < d[z - 1]; --z) d[z]= d[z - 1], di[z]= di[z - 1];
   d[z]= y, di[z]= yi;
  }
  i= j;
 }
}
// 値の順に並んだ値 d と添字 di の [0, m) (m >= 1) を走査し、値が変わるたびに順位を 1 増やして xs に足し、a の添字の位置に順位を
// 書く。rk はそれまでの値の種類数 (次に付ける順位)。最初の要素は、前の塊と値が必ず違うので新しい値になる。
template <class T, class D, class I> [[gnu::always_inline]] inline void rank_run(const D& d, const I& di, std::size_t m, T* a, std::vector<T>& xs, u32& rk) {
 T prev= d(0);
 xs.push_back(prev), ++rk;
 a[di(0)]= T(rk - 1);
 for(std::size_t i= 1; i < m; ++i) {
  const T v= d(i);
  if(v != prev) xs.push_back(v), prev= v, ++rk;
  a[di(i)]= T(rk - 1);
 }
}
// 塊 (値 s、添字 si) の [0, m) (2 <= m <= 64) を値の順に並べ、順位を数えて a に書く。fix は鍵が桁より下の bit を全部は持たないとき。
template <class T> [[gnu::always_inline]] inline void rank_bucket(const T* s, const u32* si, int m, const KeyParams& kp, bool fix, T* a, std::vector<T>& xs, u32& rk) {
 alignas(32) u32 k[64];
 __m256i v[8];
 const __m256i pad= _mm256_set1_epi32(-1);
 int nv;
 if(m <= 8) {
  nv= 1, v[0]= net::sort8(make_keys(s, m, 0, kp));
 } else if(m <= 16) {
  nv= 2, v[0]= make_keys(s, m, 0, kp), v[1]= make_keys(s, m, 1, kp);
  net::sort16(v[0], v[1]);
 } else if(m <= 32) {
  nv= 4;
  for(int t= 0; t < 4; ++t) v[t]= 8 * t < m ? make_keys(s, m, t, kp) : pad;
  net::sort32(v[0], v[1], v[2], v[3]);
 } else {
  nv= 8;
  for(int t= 0; t < 8; ++t) v[t]= 8 * t < m ? make_keys(s, m, t, kp) : pad;
  net::sort64(v);
 }
 for(int t= 0; t < nv; ++t) _mm256_store_si256((__m256i*)(k + 8 * t), v[t]);
 T d[64];
 u32 di[64];
 for(int i= 0; i < m; ++i) d[i]= s[k[i] & 63], di[i]= si[k[i] & 63];
 if(fix && has_tie(v, nv, m)) fix_ties(k, m, d, di);
 rank_run([&](std::size_t i) { return d[i]; }, [&](std::size_t i) { return di[i]; }, std::size_t(m), a, xs, rk);
}
// 値と添字の組 p[0, m) を std::sort で並べてから走査する (64 個以下の配列と、64 個を超える塊)。
template <class T> inline void rank_pairs(std::pair<T, u32>* p, std::size_t m, T* a, std::vector<T>& xs, u32& rk) {
 std::sort(p, p + m);
 rank_run([&](std::size_t i) { return p[i].first; }, [&](std::size_t i) { return p[i].second; }, m, a, xs, rk);
}
// 値の種類が 1024 個以下なら、各要素をハッシュ表で引いた順位に書き換えて xs を作り、true を返す。種類が 1024 個を超えた時点で、
// a に手を付けずに false を返す。表は 4096 か所の線形探査で、値 kv と印 tag (0 は空、埋まったら 1、あとで順位 + 1) を別の配列に
// 持つ。消すことはないので、値 x を引くと、x の場所に着くまでに空の場所を通らない。ほとんどの入力ではすぐにやめる道なので、展開させない。
template <class T> [[gnu::noinline]] inline bool few_values(T* a, std::size_t n, std::vector<T>& xs) {
 constexpr u32 LIM= 1024, B= 12, S= 1u << B;
 std::vector<u64> kv(S);
 std::vector<u32> tag(S);
 const auto slot= [](u64 x) { return u32((x * 0x9E3779B97F4A7C15ull) >> (64 - B)); };
 u32 k= 0;
 for(std::size_t i= 0; i < n; ++i) {
  const u64 x= u64(a[i]);
  for(u32 h= slot(x);; h= (h + 1) & (S - 1)) {
   if(!tag[h]) {
    if(++k > LIM) return false;
    kv[h]= x, tag[h]= 1;
    break;
   }
   if(kv[h] == x) break;
  }
 }
 xs.reserve(k);
 for(u32 h= 0; h < S; ++h)
  if(tag[h]) xs.push_back(T(kv[h]));
 std::sort(xs.begin(), xs.end());
 for(u32 r= 0; r < k; ++r) {
  const u64 x= u64(xs[r]);
  u32 h= slot(x);
  while(kv[h] != x || !tag[h]) h= (h + 1) & (S - 1);
  tag[h]= r + 1;
 }
 for(std::size_t i= 0; i < n; ++i) {
  const u64 x= u64(a[i]);
  u32 h= slot(x);
  while(kv[h] != x || !tag[h]) h= (h + 1) & (S - 1);
  a[i]= T(tag[h] - 1);
 }
 return true;
}
template <class T> inline std::vector<T> compress(T* a, std::size_t n) {
 std::vector<T> xs;
 u32 rk= 0;
 if(n <= 64) {
  std::pair<T, u32> p[64];
  for(std::size_t i= 0; i < n; ++i) p[i]= {a[i], u32(i)};
  xs.reserve(n);
  if(n) rank_pairs(p, n, a, xs, rk);
  return xs;
 }
 u64 o= 0, e= ~0ull, mn= ~0ull, mx= 0;
 for(std::size_t i= 0; i < n; ++i) {
  const u64 k= key(a[i]);
  o|= k, e&= k, mn= std::min(mn, k), mx= std::max(mx, k);
 }
 const u64 diff= o ^ e;
 if(!diff) {  // 全要素が同じ値
  xs.assign(1, a[0]);
  std::fill(a, a + n, T(0));
  return xs;
 }
 if(n >= (std::size_t(1) << 14) && few_values(a, n, xs)) return xs;
 // 鍵から base を引いた値の bit [lo, hi] の上から db bit を桁にする。要素によって違う bit の幅と、最大値と最小値の差の幅の狭い方を使う。
 int hi= 63 - __builtin_clzll(diff), lo= __builtin_ctzll(diff);
 u64 base= 0;
 bool by_range= false;
 if(const int wm= 64 - __builtin_clzll(mx - mn); wm < hi - lo + 1) lo= 0, hi= wm - 1, base= mn, by_range= true;
 int db= 1;
 while(db < 16 && (n >> db) > 16) ++db;
 db= std::min(db, hi - lo + 1);
 const int sh= hi - db + 1, r= sh - lo, kr= std::min(r, 26);
 const u64 mk= (u64(1) << db) - 1;
 std::vector<u32> pos(std::size_t(1) << db);
 for(std::size_t i= 0; i < n; ++i) ++pos[(key(a[i]) - base) >> sh & mk];
 if(!r) {
  // 桁が値を決めるので、空でない塊ごとに値が 1 つある。最大値と最小値の差で桁を選んだなら鍵は base に桁の値を足したもの、そうでなければ
  // 桁の外の bit (全要素で同じ) に桁の値を入れたもの。
  xs.reserve(std::min(n, pos.size()));
  for(std::size_t d= 0; d < pos.size(); ++d) {
   const u32 c= pos[d];
   pos[d]= rk;
   if(c) {
    const u64 kd= by_range ? base + (u64(d) << sh) : (e & ~(mk << sh)) | (u64(d) << sh);
    xs.push_back(T(kd ^ F<T>)), ++rk;
   }
  }
  for(std::size_t i= 0; i < n; ++i) a[i]= T(pos[(key(a[i]) - base) >> sh & mk]);
  return xs;
 }
 for(std::size_t d= 0, s= 0; d < pos.size(); ++d) {
  const u32 t= pos[d];
  pos[d]= u32(s), s+= t;
 }
 Buf<T> bb(n + 64);  // 塊を 8 個ずつ読むので、最後の塊の後ろを 63 個まで読む
 Buf<u32> ib(n);
 T* b= bb.p;
 u32* bi= ib.p;
 for(std::size_t i= 0; i < n; ++i) {
  const T x= a[i];
  const u32 p= pos[(key(x) - base) >> sh & mk]++;
  b[p]= x, bi[p]= u32(i);
 }
 // 振り分けたあとの pos[d] は塊 d の終わり。
 const KeyParams kp{_mm256_set1_epi64x((long long)F<T>), _mm256_set1_epi64x((long long)base), _mm256_set1_epi64x((long long)((u64(1) << kr) - 1)), _mm_cvtsi32_si128(sh - kr)};
 const bool fix= r > kr;
 xs= alloc_xs<T>(n);
 for(std::size_t d= 0, off= 0; d < pos.size(); ++d) {
  const std::size_t end= pos[d], m= end - off;
  if(m == 1) {
   xs.push_back(b[off]);
   a[bi[off]]= T(rk++);
  } else if(m >= 2 && m <= 64) {
   rank_bucket(b + off, bi + off, int(m), kp, fix, a, xs, rk);
  } else if(m > 64) {
   std::vector<std::pair<T, u32>> p(m);
   for(std::size_t i= 0; i < m; ++i) p[i]= {b[off + i], bi[off + i]};
   rank_pairs(p.data(), m, a, xs, rk);
  }
  off= end;
 }
 return xs;
}
}  // namespace c8
}  // namespace compress_internal
// a の各要素を順位 (a の中でその要素より小さい値の種類数) に置き換え、重複を除いて昇順に並べた値の列 xs を返す。xs[a[i]] が元の a[i]。
template <class T> std::vector<T> compress(std::vector<T>& a) {
 static_assert(std::is_integral_v<T>, "compress は整数の列を受ける");
 if constexpr(!std::is_same_v<T, bool> && sizeof(T) == 4) {
  return compress_internal::c4::compress(a.data(), a.size());
 } else if constexpr(sizeof(T) == 8) {
  return compress_internal::c8::compress(a.data(), a.size());
 } else {
  std::vector<T> xs(a);
  std::sort(xs.begin(), xs.end());
  xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
  for(T& x: a) x= T(std::lower_bound(xs.begin(), xs.end(), x) - xs.begin());
  return xs;
 }
}
