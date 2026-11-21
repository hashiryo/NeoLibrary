#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <type_traits>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 整数の列を昇順に並べ替える。4 byte と 8 byte の整数は、長さと CPU で道を選ぶ。ほかの型は std::sort に任せる。
//   4 byte で 64 個以下: AVX2 の bitonic のソーティングネットワーク (8 / 16 / 32 / 64 個)。最後の 8 個をふつうに読み書きして端の
//     ベクタの lane を回し、vpmaskmovd を使わない。8 個未満は挿入ソート。
//   4 byte で 2^15 個以下か、Intel の CPU: 作業用の配列へ分割する AVX2 のクイックソート。8 個ずつ読んで軸と比べ、movemask で表を
//     引いた vpermd で軸以下を前、軸より大きいものを後ろに寄せ、作業用の配列の前と後ろの両方に書く。64 個以下の塊はネットワーク。
//   4 byte でそれより長い (AMD など): LSD の基数ソート。要素によって違う bit だけを 11 bit 以下の桁に等分し、全要素が同じ値になる
//     桁は飛ばす。
//   8 byte: 上の桁で 1 回振り分けてから、塊 (平均 16 個以下) の中を、桁より下の bit の上の 26 bit と塊の中の位置 (6 bit) を詰めた
//     u32 の鍵のネットワークで並べる。桁は、要素によって違う bit の幅と最大値と最小値の差の幅の狭い方の上から選ぶ。64 個を超えた
//     塊は 64 bit のクイックソート。
// Intel の CPU では振り分けの散らばった書き込みが重く、4 byte の 10^6 個で LSD がクイックソートの 1.3〜1.5 倍かかる (AMD では逆に
// クイックソートが LSD の 1.4〜1.9 倍かかる)。作業用の配列は 2^19 byte 以上なら、Linux では 2 MB 境界の mmap に置いて huge page を
// 頼み、触る前にまとめて用意させる。procon-judge の self-sort-u32 の lsd_adaptive_hp_cfu と qsort_avx2_oop_all、self-sort-small-u32
// の net_avx2_ov、self-sort-u64 と self-sort-i64 の msd16_net32_ai を写したもの。
namespace sort_internal {
using u32= unsigned;
using u64= unsigned long long;
// n 個の T の作業用の領域。0 では埋めない。2^19 byte 以上なら、Linux では 2 MB 境界の mmap に置いて huge page を頼み、触る前に
// まとめて用意させる (MADV_POPULATE_WRITE)。振り分けは書き先が散らばるので、4 KB のページでは TLB を外しやすい。
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
// Intel の CPU か。4 byte の長い列は、Intel ならクイックソート、ほかは LSD で並べる。
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
// 4 byte の値の bitonic のソーティングネットワーク。S は符号付きか。stepD<MASK> は距離 D の lane どうしを比べ、MASK の bit が立つ
// lane に大きい方を置く。
namespace net4 {
template <bool S> [[gnu::always_inline]] inline __m256i vmin(__m256i a, __m256i b) {
 if constexpr(S) return _mm256_min_epi32(a, b);
 else return _mm256_min_epu32(a, b);
}
template <bool S> [[gnu::always_inline]] inline __m256i vmax(__m256i a, __m256i b) {
 if constexpr(S) return _mm256_max_epi32(a, b);
 else return _mm256_max_epu32(a, b);
}
template <bool S, int MASK> [[gnu::always_inline]] inline __m256i step1(__m256i v) {
 const __m256i q= _mm256_shuffle_epi32(v, 0xB1);
 return _mm256_blend_epi32(vmin<S>(v, q), vmax<S>(v, q), MASK);
}
template <bool S, int MASK> [[gnu::always_inline]] inline __m256i step2(__m256i v) {
 const __m256i q= _mm256_shuffle_epi32(v, 0x4E);
 return _mm256_blend_epi32(vmin<S>(v, q), vmax<S>(v, q), MASK);
}
template <bool S, int MASK> [[gnu::always_inline]] inline __m256i step4(__m256i v) {
 const __m256i q= _mm256_permute4x64_epi64(v, 0x4E);
 return _mm256_blend_epi32(vmin<S>(v, q), vmax<S>(v, q), MASK);
}
// bitonic な 8 個を昇順にする。
template <bool S> [[gnu::always_inline]] inline __m256i merge8(__m256i v) { return step1<S, 0xAA>(step2<S, 0xCC>(step4<S, 0xF0>(v))); }
// 任意の 8 個を昇順にする。2 個ずつ昇順と降順、4 個ずつ昇順と降順に並べて bitonic にしてから merge8。
template <bool S> [[gnu::always_inline]] inline __m256i sort8(__m256i v) { return merge8<S>(step1<S, 0x5A>(step2<S, 0x3C>(step1<S, 0x66>(v)))); }
[[gnu::always_inline]] inline __m256i reverse8(__m256i v) { return _mm256_permutevar8x32_epi32(v, _mm256_setr_epi32(7, 6, 5, 4, 3, 2, 1, 0)); }
// a, b を並べた 16 個が bitonic なら昇順にする。
template <bool S> [[gnu::always_inline]] inline void merge16(__m256i& a, __m256i& b) {
 const __m256i l= vmin<S>(a, b), h= vmax<S>(a, b);
 a= merge8<S>(l), b= merge8<S>(h);
}
template <bool S> [[gnu::always_inline]] inline void sort16(__m256i& a, __m256i& b) {
 a= sort8<S>(a), b= reverse8(sort8<S>(b));
 merge16<S>(a, b);
}
template <bool S> [[gnu::always_inline]] inline void sort32(__m256i& a, __m256i& b, __m256i& c, __m256i& d) {
 sort16<S>(a, b), sort16<S>(c, d);
 // 後ろの 16 個を逆順にして (reverse8(d), reverse8(c)) とつなぐと 32 個が bitonic になる。
 const __m256i rc= reverse8(d), rd= reverse8(c);
 __m256i l0= vmin<S>(a, rc), h0= vmax<S>(a, rc), l1= vmin<S>(b, rd), h1= vmax<S>(b, rd);
 merge16<S>(l0, l1), merge16<S>(h0, h1);
 a= l0, b= l1, c= h0, d= h1;
}
// a, b, c, d を並べた 32 個が bitonic なら昇順にする。
template <bool S> [[gnu::always_inline]] inline void merge32(__m256i& a, __m256i& b, __m256i& c, __m256i& d) {
 __m256i l0= vmin<S>(a, c), h0= vmax<S>(a, c), l1= vmin<S>(b, d), h1= vmax<S>(b, d);
 merge16<S>(l0, l1), merge16<S>(h0, h1);
 a= l0, b= l1, c= h0, d= h1;
}
// 8 本のベクタの 64 個を昇順にする。前後の 32 個を並べ、後ろを逆順にしてつないだ bitonic な 64 個を merge する。
template <bool S> [[gnu::always_inline]] inline void sort64(__m256i* v) {
 sort32<S>(v[0], v[1], v[2], v[3]), sort32<S>(v[4], v[5], v[6], v[7]);
 const __m256i r4= reverse8(v[7]), r5= reverse8(v[6]), r6= reverse8(v[5]), r7= reverse8(v[4]);
 __m256i l0= vmin<S>(v[0], r4), h0= vmax<S>(v[0], r4), l1= vmin<S>(v[1], r5), h1= vmax<S>(v[1], r5);
 __m256i l2= vmin<S>(v[2], r6), h2= vmax<S>(v[2], r6), l3= vmin<S>(v[3], r7), h3= vmax<S>(v[3], r7);
 merge32<S>(l0, l1, l2, l3), merge32<S>(h0, h1, h2, h3);
 v[0]= l0, v[1]= l1, v[2]= l2, v[3]= l3, v[4]= h0, v[5]= h1, v[6]= h2, v[7]= h3;
}
}  // namespace net4
// 4 byte の 64 個以下の列。
namespace small4 {
template <class T> inline void insertion(T* a, int n) {
 for(int i= 1; i < n; ++i) {
  const T x= a[i];
  int j= i;
  for(; j > 0 && a[j - 1] > x; --j) a[j]= a[j - 1];
  a[j]= x;
 }
}
// a[0, n) (8 <= n <= 64、8 (NV / 2) < n <= 8 NV) を NV 本のベクタのネットワークで並べる。q は 8 個そろったベクタの本数、c (1 から 8) は
// その次の端のベクタに入る個数。最後の 8 個 a[n - 8, n) はいつでもふつうに読み書きできる。読むときは最後の 8 個を読み、端のベクタに入る
// c 個を下の lane へ回してから、残りの lane を最大値で埋める。書くときは、並べ終えた最後の 8 個を、端の 1 つ前のベクタの上の 8 - c 個と
// 端のベクタの下の c 個から回して作り、a[n - 8, n) に書く (前のベクタと重なる所は同じ値を書き直すだけ)。
template <class T, int NV> [[gnu::always_inline]] inline void sort_ov(T* a, int n) {
 constexpr bool S= std::is_signed_v<T>;
 const int q= (n - 1) >> 3, c= n - 8 * q;
 const __m256i pad= _mm256_set1_epi32(int(std::numeric_limits<T>::max())), iota= _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
 const __m256i in_idx= _mm256_add_epi32(iota, _mm256_set1_epi32(8 - c)), out_idx= _mm256_add_epi32(iota, _mm256_set1_epi32(c));
 const __m256i in_mask= _mm256_cmpgt_epi32(_mm256_set1_epi32(c), iota), out_mask= _mm256_cmpgt_epi32(iota, _mm256_set1_epi32(7 - c));
 __m256i v[NV];
#pragma GCC unroll 8
 for(int t= 0; t < NV; ++t) {
  if(t < q) v[t]= _mm256_loadu_si256((const __m256i*)(a + 8 * t));
  else if(t == q) v[t]= _mm256_blendv_epi8(pad, _mm256_permutevar8x32_epi32(_mm256_loadu_si256((const __m256i*)(a + n - 8)), in_idx), in_mask);
  else v[t]= pad;
 }
 if constexpr(NV == 1) v[0]= net4::sort8<S>(v[0]);
 else if constexpr(NV == 2) net4::sort16<S>(v[0], v[1]);
 else if constexpr(NV == 4) net4::sort32<S>(v[0], v[1], v[2], v[3]);
 else net4::sort64<S>(v);
#pragma GCC unroll 8
 for(int t= 0; t < NV; ++t) {
  if(t < q) {
   _mm256_storeu_si256((__m256i*)(a + 8 * t), v[t]);
  } else if(t == q) {
   if constexpr(NV == 1) _mm256_storeu_si256((__m256i*)a, v[0]);
   else _mm256_storeu_si256((__m256i*)(a + n - 8), t ? _mm256_blendv_epi8(_mm256_permutevar8x32_epi32(v[t ? t - 1 : 0], out_idx), _mm256_permutevar8x32_epi32(v[t], out_idx), out_mask) : v[0]);
  }
 }
}
template <class T> [[gnu::always_inline]] inline void sort(T* a, int n) {
 if(n < 8) return insertion(a, n);
 if(n == 8) return sort_ov<T, 1>(a, 8);
 if(n <= 16) return sort_ov<T, 2>(a, n);
 if(n <= 32) return sort_ov<T, 4>(a, n);
 sort_ov<T, 8>(a, n);
}
}  // namespace small4
// 4 byte の列の、作業用の配列へ分割するクイックソート。分割のたびに a と作業用の配列の役を入れ替え、並べ終えた塊が作業用の配列の側に
// あれば a へ写す。読む位置は先頭から順に進むだけなので、その場で分けるときの左右の選び分けも、両端を先に取る手間も要らない。軸は
// 8192 個以上なら 64 個の標本の中央値、128 個以上なら 9 点、それ未満は 3 点の中央値。全要素が軸以下なら軸は最大値なので、軸と等しい
// 要素を後ろに分けて前だけを続ける (重複の多い入力で止まらないため)。再帰が深くなりすぎたら std::sort に任せる。
namespace quick4 {
// 表の m 番目: bit i が立っている lane (軸より大きい) を後ろに、立っていない lane を前に集める vpermd の添字を 3 bit ずつ詰めたもの。
struct PermTable {
 u32 t[256];
 constexpr PermTable(): t() {
  for(int m= 0; m < 256; ++m) {
   u32 packed= 0;
   int k= 0;
   for(int i= 0; i < 8; ++i)
    if(!(m >> i & 1)) packed|= u32(i) << (3 * k++);
   for(int i= 0; i < 8; ++i)
    if(m >> i & 1) packed|= u32(i) << (3 * k++);
   t[m]= packed;
  }
 }
};
inline constexpr PermTable PERM{};
// x[0, n) を、p 以下を y の前から、p より大きいものを y の後ろから詰めて書き、前に書いた個数を返す。前の書き込み [l, l + 8) と後ろの
// 書き込み [r - 8, r) が重なると後ろの書き込みが前に置いた値を上書きするので、16 個以上残っている間だけベクタで進め、残りは 1 個ずつ置く。
// 符号なしの値は、最上位 bit を反転して符号付きで比べる (AVX2 に符号なしの比較が無いため)。
template <class T> inline std::size_t partition(const T* x, T* y, std::size_t n, T p) {
 constexpr u32 F= std::is_signed_v<T> ? 0 : 0x80000000u;
 const __m256i sign= _mm256_set1_epi32(int(F));
 const __m256i pv= _mm256_set1_epi32(int(u32(p) ^ F));
 const __m256i shifts= _mm256_setr_epi32(0, 3, 6, 9, 12, 15, 18, 21);
 std::size_t l= 0, r= n, i= 0;
 for(; i + 16 <= n; i+= 8) {
  const __m256i v= _mm256_loadu_si256((const __m256i*)(x + i));
  const int m= _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpgt_epi32(_mm256_xor_si256(v, sign), pv)));
  const __m256i w= _mm256_permutevar8x32_epi32(v, _mm256_srlv_epi32(_mm256_set1_epi32(int(PERM.t[m])), shifts));
  const int nr= __builtin_popcount(m);
  _mm256_storeu_si256((__m256i*)(y + l), w);
  _mm256_storeu_si256((__m256i*)(y + r - 8), w);
  l+= 8 - nr, r-= nr;
 }
 for(; i < n; ++i) {
  if(x[i] <= p) y[l++]= x[i];
  else y[--r]= x[i];
 }
 return l;
}
// 端の半端なベクタは vpmaskmovd で読み書きし、配列の外には触れない。rem はそのベクタから先に残っている個数。読まない lane は最大値。
template <class T> inline __m256i load_pad(const T* a, int rem) {
 if(rem >= 8) return _mm256_loadu_si256((const __m256i*)a);
 const __m256i m= _mm256_cmpgt_epi32(_mm256_set1_epi32(rem), _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7));
 return _mm256_blendv_epi8(_mm256_set1_epi32(int(std::numeric_limits<T>::max())), _mm256_maskload_epi32((const int*)a, m), m);
}
template <class T> inline void store_part(T* a, int rem, __m256i v) {
 if(rem >= 8) return _mm256_storeu_si256((__m256i*)a, v);
 const __m256i m= _mm256_cmpgt_epi32(_mm256_set1_epi32(rem), _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7));
 _mm256_maskstore_epi32((int*)a, m, v);
}
// a[0, n) (n <= 64) を、8 個ずつのベクタに最大値で埋めて読み、8 / 16 / 32 / 64 個のネットワークで並べて書き戻す。
template <class T> inline void small_sort(T* a, std::size_t n) {
 constexpr bool S= std::is_signed_v<T>;
 if(n <= 1) return;
 const int r= int(n);
 if(n <= 8) {
  store_part(a, r, net4::sort8<S>(load_pad(a, r)));
 } else if(n <= 16) {
  __m256i x= load_pad(a, 8), y= load_pad(a + 8, r - 8);
  net4::sort16<S>(x, y);
  store_part(a, 8, x), store_part(a + 8, r - 8, y);
 } else if(n <= 32) {
  __m256i x= load_pad(a, 8), y= load_pad(a + 8, 8), z= load_pad(a + 16, r - 16), w= load_pad(a + 24, r - 24);
  net4::sort32<S>(x, y, z, w);
  store_part(a, 8, x), store_part(a + 8, 8, y), store_part(a + 16, r - 16, z), store_part(a + 24, r - 24, w);
 } else {
  __m256i v[8];
  for(int j= 0; j < 8; ++j) v[j]= load_pad(a + 8 * j, r - 8 * j);
  net4::sort64<S>(v);
  for(int j= 0; j < 8; ++j) store_part(a + 8 * j, r - 8 * j, v[j]);
 }
}
template <class T> inline T median3(T a, T b, T c) { return std::max(std::min(a, b), std::min(std::max(a, b), c)); }
template <class T> inline T choose_pivot(const T* a, std::size_t n) {
 if(n >= 8192) {
  alignas(32) T s[64];
  const std::size_t step= n / 64;
  for(int i= 0; i < 64; ++i) s[i]= a[i * step + step / 2];
  __m256i v[8];
  for(int j= 0; j < 8; ++j) v[j]= _mm256_load_si256((const __m256i*)(s + 8 * j));
  net4::sort64<std::is_signed_v<T>>(v);
  return T(_mm256_extract_epi32(v[4], 0));
 }
 if(n < 128) return median3(a[0], a[n / 2], a[n - 1]);
 const std::size_t s= n / 8;
 return median3(median3(a[0], a[s], a[2 * s]), median3(a[3 * s], a[4 * s], a[5 * s]), median3(a[6 * s], a[7 * s], a[n - 1]));
}
// x[0, n) を並べる。y は同じ長さの作業用の領域で、分割のたびに x と y を入れ替える。x_is_a は x が元の配列の側かで、並べ終えた値は
// 元の配列の側に置く (元の配列の側でなければ最後に y へ写す)。
template <class T> inline void rec(T* x, T* y, std::size_t n, int depth, bool x_is_a) {
 while(n > 64) {
  if(depth-- == 0) {
   std::sort(x, x + n);
   if(!x_is_a) std::memcpy(y, x, n * sizeof(T));
   return;
  }
  const T p= choose_pivot(x, n);
  std::size_t k= partition(x, y, n, p);
  std::swap(x, y), x_is_a= !x_is_a;
  if(k == n) {
   // 全要素が p 以下なので p は最大値。p と等しい要素 (並べ終わり) を後ろに分け、前だけを続ける。
   if(p == std::numeric_limits<T>::min()) {  // 全要素が最小値
    if(!x_is_a) std::memcpy(y, x, n * sizeof(T));
    return;
   }
   k= partition(x, y, n, T(p - 1));
   std::swap(x, y), x_is_a= !x_is_a;
   if(!x_is_a) std::memcpy(y + k, x + k, (n - k) * sizeof(T));
   n= k;
   continue;
  }
  // 小さい側を再帰で、大きい側をループで並べる (スタックの深さを log n に抑える)。
  if(k < n - k) rec(x, y, k, depth, x_is_a), x+= k, y+= k, n-= k;
  else rec(x + k, y + k, n - k, depth, x_is_a), n= k;
 }
 small_sort(x, n);
 if(!x_is_a) std::memcpy(y, x, n * sizeof(T));
}
template <class T> inline void sort(T* a, std::size_t n) {
 Buf<T> buf(n);
 int log= 0;
 for(std::size_t t= n; t>>= 1;) ++log;
 rec(a, buf.p, n, 2 * log, true);
}
}  // namespace quick4
// 4 byte の列の LSD の基数ソート。最初の走査で全要素の bit or と bit and を取り、要素によって違う bit (or ^ and) の最下位から最上位まで
// だけを 11 bit 以下の桁に等分する。値が 30 bit に収まる入力なら 10 bit × 3 回、下位 16 bit がすべて 0 なら 8 bit × 2 回になる。さらに
// ヒストグラムを見て、全要素が同じ値になる桁は飛ばす。振り分けの回数が奇数なら最初に作業用の配列へ写し、最後が元の配列で終わるように
// する。写しはヒストグラムを作る走査でまとめて行う。符号付きの値は最上位 bit を反転した鍵で桁を取る。振り分けでは、数え上げの表の
// ポインタを空の asm に通し、表への書き込みを書き先への書き込みより先に置かせる (GCC は表が書き先と重ならないと分かると順を入れ替え、
// EPYC 7763 では 1 割 5 分ほど遅かった)。振り分けとヒストグラムのループは GCC にも 4 回ずつ展開させる。
namespace radix4 {
template <class T> inline void sort(T* a, std::size_t n) {
 constexpr u32 F= std::is_signed_v<T> ? 0x80000000u : 0;
 u32 o= 0, e= ~0u;
 for(std::size_t i= 0; i < n; ++i) {
  const u32 k= u32(a[i]) ^ F;
  o|= k, e&= k;
 }
 const u32 diff= o ^ e;
 if(!diff) return;  // 全要素が同じ値
 const int lo= __builtin_ctz(diff), hi= 31 - __builtin_clz(diff);
 const int w= hi - lo + 1, p= (w + 10) / 11;
 int sh[3], bits[3];
 for(int d= 0, s= lo; d < p; ++d) {
  bits[d]= w / p + (d < w % p);
  sh[d]= s, s+= bits[d];
 }
 u32 cnt[3][2048];
 for(int d= 0; d < p; ++d) std::memset(cnt[d], 0, sizeof(u32) << bits[d]);
 Buf<T> buf(n);
 T* b= buf.p;
 const bool copied= p & 1;
 if(p == 1) {
  const u32 m0= (1u << bits[0]) - 1;
#pragma GCC unroll 4
  for(std::size_t i= 0; i < n; ++i) {
   const T x= a[i];
   b[i]= x;
   ++cnt[0][(u32(x) ^ F) >> sh[0] & m0];
  }
 } else if(p == 2) {
  const u32 m0= (1u << bits[0]) - 1, m1= (1u << bits[1]) - 1;
#pragma GCC unroll 4
  for(std::size_t i= 0; i < n; ++i) {
   const u32 x= u32(a[i]) ^ F;
   ++cnt[0][x >> sh[0] & m0], ++cnt[1][x >> sh[1] & m1];
  }
 } else {
  const u32 m0= (1u << bits[0]) - 1, m1= (1u << bits[1]) - 1, m2= (1u << bits[2]) - 1;
#pragma GCC unroll 4
  for(std::size_t i= 0; i < n; ++i) {
   const T x= a[i];
   b[i]= x;
   const u32 k= u32(x) ^ F;
   ++cnt[0][k >> sh[0] & m0], ++cnt[1][k >> sh[1] & m1], ++cnt[2][k >> sh[2] & m2];
  }
 }
 // 全要素が同じ値になる桁を除き、残りの桁のヒストグラムを先頭からの位置に直す。
 int use[3], k= 0;
 for(int d= 0; d < p; ++d) {
  const u32 sz= 1u << bits[d];
  bool trivial= false;
  for(u32 i= 0; i < sz; ++i)
   if(cnt[d][i]) {
    trivial= cnt[d][i] == n;
    break;
   }
  if(trivial) continue;
  u32 s= 0;
  for(u32 i= 0; i < sz; ++i) {
   const u32 t= cnt[d][i];
   cnt[d][i]= s, s+= t;
  }
  use[k++]= d;
 }
 // 奇数回なら写しから、偶数回なら元の配列から始めると、最後が元の配列で終わる。
 T *src= a, *dst= b;
 if(k & 1) {
  if(!copied) std::memcpy(b, a, n * sizeof(T));
  src= b, dst= a;
 }
 for(int j= 0; j < k; ++j) {
  const int d= use[j], s= sh[d];
  const u32 m= (1u << bits[d]) - 1;
  u32* c= cnt[d];
  asm("" : "+r"(c));  // c が dst と重なるかもしれないと思わせ、下の 2 つの書き込みの順を保たせる
#pragma GCC unroll 4
  for(std::size_t i= 0; i < n; ++i) {
   const T x= src[i];
   const u32 t= (u32(x) ^ F) >> s & m, pos= c[t];
   c[t]= pos + 1;
   dst[pos]= x;
  }
  std::swap(src, dst);
 }
}
}  // namespace radix4
// 65 個以上の 4 byte の列。
template <class T> inline void long4(T* a, std::size_t n) {
 if(n <= (std::size_t(1) << 15) || is_intel()) quick4::sort(a, n);
 else radix4::sort(a, n);
}
// 8 byte の列。
namespace msd8 {
template <class T> [[gnu::always_inline]] inline u64 key(T x) { return std::is_signed_v<T> ? u64(x) ^ (u64(1) << 63) : u64(x); }
// 表の m 番目: bit j が立っている lane (軸より大きい) を後ろに、立っていない lane を前に集める vpermd の添字 (64 bit の lane j は
// 32 bit の lane 2j と 2j + 1)。
struct PermTable {
 u32 t[16][8];
 constexpr PermTable(): t() {
  for(int m= 0; m < 16; ++m) {
   int k= 0;
   for(int pass= 0; pass < 2; ++pass)
    for(int j= 0; j < 4; ++j)
     if((m >> j & 1) == pass) t[m][2 * k]= u32(2 * j), t[m][2 * k + 1]= u32(2 * j + 1), ++k;
  }
 }
};
inline constexpr PermTable PERM{};
// 符号付きの比較に直すときに足す値 (符号なしなら最上位 bit を反転する)。
template <class T> inline __m256i bias() { return std::is_signed_v<T> ? _mm256_setzero_si256() : _mm256_set1_epi64x((long long)(1ull << 63)); }
// a > b の lane が全 bit 1。
template <class T> inline __m256i gt(__m256i a, __m256i b) {
 if constexpr(std::is_signed_v<T>) return _mm256_cmpgt_epi64(a, b);
 const __m256i s= bias<T>();
 return _mm256_cmpgt_epi64(_mm256_xor_si256(a, s), _mm256_xor_si256(b, s));
}
// x[0, n) を、p 以下を y の前から、p より大きいものを y の後ろから詰めて書き、前に書いた個数を返す。
template <class T> inline std::size_t partition(const T* x, T* y, std::size_t n, T p) {
 const __m256i s= bias<T>();
 const __m256i pv= _mm256_xor_si256(_mm256_set1_epi64x((long long)p), s);
 std::size_t l= 0, r= n, i= 0;
 for(; i + 8 <= n; i+= 4) {
  const __m256i v= _mm256_loadu_si256((const __m256i*)(x + i));
  const int m= _mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpgt_epi64(_mm256_xor_si256(v, s), pv)));
  const __m256i w= _mm256_permutevar8x32_epi32(v, _mm256_loadu_si256((const __m256i*)PERM.t[m]));
  const int nr= __builtin_popcount(m);
  _mm256_storeu_si256((__m256i*)(y + l), w);
  _mm256_storeu_si256((__m256i*)(y + r - 4), w);
  l+= 4 - nr, r-= nr;
 }
 for(; i < n; ++i) {
  if(x[i] <= p) y[l++]= x[i];
  else y[--r]= x[i];
 }
 return l;
}
// 64 bit の値の bitonic のソーティングネットワーク (AVX2 に 64 bit の min と max が無いので、比較と blend で組む)。a と b の lane ごとの
// 小さい方を a に、大きい方を b に置く。
template <class T> [[gnu::always_inline]] inline void minmax(__m256i& a, __m256i& b) {
 const __m256i g= gt<T>(a, b);
 const __m256i lo= _mm256_blendv_epi8(a, b, g), hi= _mm256_blendv_epi8(b, a, g);
 a= lo, b= hi;
}
// 距離 1 と 2 の lane どうしを比べ、MASK (32 bit の lane 単位) の bit が立つ lane に大きい方を置く。
template <class T, int MASK> [[gnu::always_inline]] inline __m256i step1(__m256i v) {
 __m256i q= _mm256_permute4x64_epi64(v, 0xB1), lo= v, hi= q;
 minmax<T>(lo, hi);
 return _mm256_blend_epi32(lo, hi, MASK);
}
template <class T, int MASK> [[gnu::always_inline]] inline __m256i step2(__m256i v) {
 __m256i q= _mm256_permute4x64_epi64(v, 0x4E), lo= v, hi= q;
 minmax<T>(lo, hi);
 return _mm256_blend_epi32(lo, hi, MASK);
}
// bitonic な 4 個を昇順にする。
template <class T> [[gnu::always_inline]] inline __m256i merge4(__m256i v) { return step1<T, 0xCC>(step2<T, 0xF0>(v)); }
// 任意の 4 個を昇順にする。2 個ずつ昇順と降順に並べて bitonic にしてから merge4。
template <class T> [[gnu::always_inline]] inline __m256i sort4(__m256i v) { return merge4<T>(step1<T, 0x3C>(v)); }
[[gnu::always_inline]] inline __m256i reverse4(__m256i v) { return _mm256_permute4x64_epi64(v, 0x1B); }
// a, b を並べた 8 個が bitonic なら昇順にする。
template <class T> [[gnu::always_inline]] inline void merge8(__m256i& a, __m256i& b) {
 minmax<T>(a, b);
 a= merge4<T>(a), b= merge4<T>(b);
}
template <class T> [[gnu::always_inline]] inline void sort8(__m256i& a, __m256i& b) {
 a= sort4<T>(a), b= reverse4(sort4<T>(b));
 merge8<T>(a, b);
}
// a, b, c, d を並べた 16 個が bitonic なら昇順にする。
template <class T> [[gnu::always_inline]] inline void merge16(__m256i& a, __m256i& b, __m256i& c, __m256i& d) {
 minmax<T>(a, c), minmax<T>(b, d);
 merge8<T>(a, b), merge8<T>(c, d);
}
template <class T> [[gnu::always_inline]] inline void sort16(__m256i& a, __m256i& b, __m256i& c, __m256i& d) {
 sort8<T>(a, b), sort8<T>(c, d);
 // 後ろの 8 個を逆順にして (reverse4(d), reverse4(c)) とつなぐと 16 個が bitonic になる。
 __m256i rc= reverse4(d), rd= reverse4(c);
 minmax<T>(a, rc), minmax<T>(b, rd);
 merge8<T>(a, b), merge8<T>(rc, rd);
 c= rc, d= rd;
}
// 8 本のベクタの 32 個を昇順にする。
template <class T> [[gnu::always_inline]] inline void sort32(__m256i* v) {
 sort16<T>(v[0], v[1], v[2], v[3]), sort16<T>(v[4], v[5], v[6], v[7]);
 __m256i r4= reverse4(v[7]), r5= reverse4(v[6]), r6= reverse4(v[5]), r7= reverse4(v[4]);
 minmax<T>(v[0], r4), minmax<T>(v[1], r5), minmax<T>(v[2], r6), minmax<T>(v[3], r7);
 merge16<T>(v[0], v[1], v[2], v[3]), merge16<T>(r4, r5, r6, r7);
 v[4]= r4, v[5]= r5, v[6]= r6, v[7]= r7;
}
// a[0, n) (n <= 32) を、4 個ずつのベクタに最大値で埋めて読み、4 / 8 / 16 / 32 個のネットワークで並べて書き戻す。端の半端な
// ベクタは vpmaskmovq で読み書きし、配列の外には触れない。rem はそのベクタから先に残っている個数。
template <class T> inline __m256i load_pad(const T* a, long long rem) {
 if(rem >= 4) return _mm256_loadu_si256((const __m256i*)a);
 const __m256i m= _mm256_cmpgt_epi64(_mm256_set1_epi64x(rem), _mm256_setr_epi64x(0, 1, 2, 3));
 const __m256i pad= _mm256_set1_epi64x((long long)std::numeric_limits<T>::max());
 return _mm256_blendv_epi8(pad, _mm256_maskload_epi64((const long long*)a, m), m);
}
template <class T> inline void store_part(T* a, long long rem, __m256i v) {
 if(rem >= 4) return _mm256_storeu_si256((__m256i*)a, v);
 const __m256i m= _mm256_cmpgt_epi64(_mm256_set1_epi64x(rem), _mm256_setr_epi64x(0, 1, 2, 3));
 _mm256_maskstore_epi64((long long*)a, m, v);
}
template <class T> inline void small_sort(T* a, std::size_t n) {
 if(n <= 1) return;
 const long long r= (long long)n;
 if(n <= 4) {
  store_part(a, r, sort4<T>(load_pad(a, r)));
 } else if(n <= 8) {
  __m256i x= load_pad(a, 4), y= load_pad(a + 4, r - 4);
  sort8<T>(x, y);
  store_part(a, 4, x), store_part(a + 4, r - 4, y);
 } else if(n <= 16) {
  __m256i x= load_pad(a, 4), y= load_pad(a + 4, 4), z= load_pad(a + 8, r - 8), w= load_pad(a + 12, r - 12);
  sort16<T>(x, y, z, w);
  store_part(a, 4, x), store_part(a + 4, 4, y), store_part(a + 8, r - 8, z), store_part(a + 12, r - 12, w);
 } else {
  __m256i v[8];
  for(int j= 0; j < 8; ++j) v[j]= load_pad(a + 4 * j, r - 4 * j);
  sort32<T>(v);
  for(int j= 0; j < 8; ++j) store_part(a + 4 * j, r - 4 * j, v[j]);
 }
}
template <class T> inline T median3(T a, T b, T c) { return std::max(std::min(a, b), std::min(std::max(a, b), c)); }
// 8192 個以上の塊は、等間隔に取った 64 個の標本の中央値 (std::nth_element) を軸にする。128 個以上は 9 点、それ未満は 3 点の中央値。
template <class T> inline T choose_pivot(const T* a, std::size_t n) {
 if(n >= 8192) {
  T s[64];
  const std::size_t step= n / 64;
  for(int i= 0; i < 64; ++i) s[i]= a[i * step + step / 2];
  std::nth_element(s, s + 32, s + 64);
  return s[32];
 }
 if(n < 128) return median3(a[0], a[n / 2], a[n - 1]);
 const std::size_t s= n / 8;
 return median3(median3(a[0], a[s], a[2 * s]), median3(a[3 * s], a[4 * s], a[5 * s]), median3(a[6 * s], a[7 * s], a[n - 1]));
}
// x[0, n) を並べる。y は同じ長さの作業用の領域で、分割のたびに x と y を入れ替える (quick4::rec と同じ形)。
template <class T> inline void rec(T* x, T* y, std::size_t n, int depth, bool x_is_a) {
 while(n > 32) {
  if(depth-- == 0) {
   std::sort(x, x + n);
   if(!x_is_a) std::memcpy(y, x, n * sizeof(T));
   return;
  }
  const T p= choose_pivot(x, n);
  std::size_t k= partition(x, y, n, p);
  std::swap(x, y), x_is_a= !x_is_a;
  if(k == n) {
   if(p == std::numeric_limits<T>::min()) {
    if(!x_is_a) std::memcpy(y, x, n * sizeof(T));
    return;
   }
   k= partition(x, y, n, T(p - 1));
   std::swap(x, y), x_is_a= !x_is_a;
   if(!x_is_a) std::memcpy(y + k, x + k, (n - k) * sizeof(T));
   n= k;
   continue;
  }
  if(k < n - k) rec(x, y, k, depth, x_is_a), x+= k, y+= k, n-= k;
  else rec(x + k, y + k, n - k, depth, x_is_a), n= k;
 }
 small_sort(x, n);
 if(!x_is_a) std::memcpy(y, x, n * sizeof(T));
}
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
// 並べた鍵 k[0, m) の上の 26 bit が等しい並びごとに、d の同じ範囲を挿入ソートで値の順に直す。
template <class T> inline void fix_ties(const u32* k, int m, T* d) {
 for(int i= 0; i < m;) {
  int j= i + 1;
  while(j < m && (k[j] ^ k[i]) < 64) ++j;
  for(int x= i + 1; x < j; ++x) {
   const T y= d[x];
   int z= x;
   for(; z > i && y < d[z - 1]; --z) d[z]= d[z - 1];
   d[z]= y;
  }
  i= j;
 }
}
// 塊 s[0, m) (2 <= m <= 64) を並べて d[0, m) に書く。fix は鍵が桁より下の bit を全部は持たないとき。
template <class T> [[gnu::always_inline]] inline void sort_bucket(const T* s, int m, const KeyParams& kp, bool fix, T* d) {
 alignas(32) u32 k[64];
 __m256i v[8];
 const __m256i pad= _mm256_set1_epi32(-1);
 int nv;
 if(m <= 8) {
  nv= 1, v[0]= net4::sort8<false>(make_keys(s, m, 0, kp));
 } else if(m <= 16) {
  nv= 2, v[0]= make_keys(s, m, 0, kp), v[1]= make_keys(s, m, 1, kp);
  net4::sort16<false>(v[0], v[1]);
 } else if(m <= 32) {
  nv= 4;
  for(int t= 0; t < 4; ++t) v[t]= 8 * t < m ? make_keys(s, m, t, kp) : pad;
  net4::sort32<false>(v[0], v[1], v[2], v[3]);
 } else {
  nv= 8;
  for(int t= 0; t < 8; ++t) v[t]= 8 * t < m ? make_keys(s, m, t, kp) : pad;
  net4::sort64<false>(v);
 }
 for(int t= 0; t < nv; ++t) _mm256_store_si256((__m256i*)(k + 8 * t), v[t]);
 for(int i= 0; i < m; ++i) d[i]= s[k[i] & 63];
 if(fix && has_tie(v, nv, m)) fix_ties(k, m, d);
}
template <class T> inline void sort(T* a, std::size_t n) {
 if(n <= 32) return small_sort(a, n);
 u64 o= 0, e= ~0ull, mn= ~0ull, mx= 0;
 for(std::size_t i= 0; i < n; ++i) {
  const u64 k= key(a[i]);
  o|= k, e&= k, mn= std::min(mn, k), mx= std::max(mx, k);
 }
 const u64 diff= o ^ e;
 if(!diff) return;  // 全要素が同じ値
 // 鍵から base を引いた値の bit [lo, hi] の上から db bit を桁にする。要素によって違う bit の幅と、最大値と最小値の差の幅の狭い方を使う。
 int hi= 63 - __builtin_clzll(diff), lo= __builtin_ctzll(diff);
 u64 base= 0;
 if(const int wm= 64 - __builtin_clzll(mx - mn); wm < hi - lo + 1) lo= 0, hi= wm - 1, base= mn;
 int db= 1;
 while(db < 16 && (n >> db) > 16) ++db;
 db= std::min(db, hi - lo + 1);
 const int sh= hi - db + 1, r= sh - lo, rk= std::min(r, 26);
 const u64 mk= (u64(1) << db) - 1;
 std::vector<u32> pos(std::size_t(1) << db);
 for(std::size_t i= 0; i < n; ++i) ++pos[(key(a[i]) - base) >> sh & mk];
 for(std::size_t d= 0, s= 0; d < pos.size(); ++d) {
  const u32 t= pos[d];
  pos[d]= u32(s), s+= t;
 }
 Buf<T> buf(n + 64);  // 塊を 8 個ずつ読むので、最後の塊の後ろを 63 個まで読む
 T* b= buf.p;
 for(std::size_t i= 0; i < n; ++i) {
  const T x= a[i];
  b[pos[(key(x) - base) >> sh & mk]++]= x;
 }
 if(!r) {  // 桁より下の bit が全要素で同じなので、振り分けた順で並んでいる
  std::memcpy(a, b, n * sizeof(T));
  return;
 }
 // 振り分けたあとの pos[d] は塊 d の終わり。
 const KeyParams kp{std::is_signed_v<T> ? _mm256_set1_epi64x((long long)(1ull << 63)) : _mm256_setzero_si256(), _mm256_set1_epi64x((long long)base),
                    _mm256_set1_epi64x((long long)((u64(1) << rk) - 1)), _mm_cvtsi32_si128(sh - rk)};
 const bool fix= r > rk;
 for(std::size_t d= 0, off= 0; d < pos.size(); ++d) {
  const std::size_t end= pos[d], m= end - off;
  if(m <= 1) {
   if(m) a[off]= b[off];
  } else if(m <= 64) {
   sort_bucket(b + off, int(m), kp, fix, a + off);
  } else {
   int log= 0;
   for(std::size_t t= m; t>>= 1;) ++log;
   rec(b + off, a + off, m, 2 * log, false);
  }
  off= end;
 }
}
}  // namespace msd8
}  // namespace sort_internal
// a を昇順に並べ替える。4 byte と 8 byte の整数の列は上の道で並べ、ほかの型は std::sort に任せる。
template <class T> [[gnu::always_inline]] inline void sort(std::vector<T>& a) {
 if constexpr(std::is_integral_v<T> && !std::is_same_v<T, bool> && sizeof(T) == 4) {
  const std::size_t n= a.size();
  if(n <= 64) sort_internal::small4::sort(a.data(), int(n));
  else sort_internal::long4(a.data(), n);
 } else if constexpr(std::is_integral_v<T> && sizeof(T) == 8) {
  sort_internal::msd8::sort(a.data(), a.size());
 } else {
  std::sort(a.begin(), a.end());
 }
}
