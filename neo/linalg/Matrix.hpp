#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <bit>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <new>
#include <type_traits>
#include <utility>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
#include "neo/algebra/ZMod.hpp"
// 行列。大きさを実行時に決め、行優先の 1 本の配列 a に持つ (a[i * w + j] が i 行 j 列)。a は 2 MB 以上なら、Linux では
// 2 MB 境界の mmap に置いて huge page を頼む (linalg_internal::HugeAlloc)。要素の型 T は和 (+)、積 (*)、0、1 を
// 持てばよい (半環)。0 と 1 は、T に静的メンバ関数 zero() と one() があればそれを、なければ T(0) と T(1) を使う
// (neo/algebra/semiring.hpp の MinPlus などは前者)。引き算とスカラー倍は、T がそれを持つときだけ使える。
// T が ZMod なら、行列の積、ベクトルとの積、累乗は linalg_internal::zmat の核を使う。積は値を M 未満に直して u64 の和にため、
// 最後に 1 回だけ還元する。小さい行列はスカラーで、ほかは AVX2 の 4 行 × 8 列のタイルで掛け、1 辺を 2 の冪に広げても
// 無駄の少ない大きい積 (法が奇数) は Strassen の再帰 (Winograd の形) にする。累乗は、法が素数で辺と指数が大きければ Frobenius 標準形を経由し、
// それ以外は二分累乗。
namespace linalg_internal {
// vector の allocator。2^19 byte 以上なら、Linux では 2 MB 境界の mmap (2 MB の倍数に切り上げる) に置いて MADV_HUGEPAGE を
// 頼む。頼まない領域を 4 KB のページに置く判定機では、大きな行列を初めて触るときのページ 1 枚ごとの費用 (1 枚 1.7 us ほど) が
// 重いため (n = 500 の u64 の行列で 0.85 ms)。
template <class T> struct HugeAlloc {
 using value_type= T;
 static constexpr std::size_t H= std::size_t(1) << 21, L= std::size_t(1) << 19;
 HugeAlloc()= default;
 template <class U> HugeAlloc(const HugeAlloc<U>&) {}
 T* allocate(std::size_t n) {
#ifdef __linux__
  if(n * sizeof(T) >= L) {
   const std::size_t bytes= (n * sizeof(T) + H - 1) & ~(H - 1);
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
  if(n * sizeof(T) >= L) {
   munmap(p, (n * sizeof(T) + H - 1) & ~(H - 1));
   return;
  }
#endif
  std::allocator<T>().deallocate(p, n);
 }
 template <class U> bool operator==(const HugeAlloc<U>&) const { return true; }
 template <class U> bool operator!=(const HugeAlloc<U>&) const { return false; }
};
// 大きな作業の配列 (掃き出しの行列など) に使う vector。
template <class T> using hvec= std::vector<T, HugeAlloc<T>>;
// n 個の u32 の作業の領域。0 では埋めない。Linux で 2^19 byte 以上なら、2 MB 境界の mmap に置いて huge page を頼み、先に
// ページを用意させる (MADV_POPULATE_WRITE、Linux 5.14 から。古い kernel では何もしない)。
struct Region {
 unsigned* p= nullptr;
 std::size_t mapped= 0;
 explicit Region(std::size_t n) {
#ifdef __linux__
  constexpr std::size_t H= std::size_t(1) << 21;
  if(n * sizeof(unsigned) >= (std::size_t(1) << 19)) {
   const std::size_t bytes= (n * sizeof(unsigned) + H - 1) & ~(H - 1);
   void* m= mmap(nullptr, bytes + H, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
   if(m != MAP_FAILED) {
    char *q= static_cast<char*>(m), *a= reinterpret_cast<char*>((reinterpret_cast<std::uintptr_t>(q) + H - 1) & ~(H - 1));
    if(a != q) munmap(q, a - q);
    munmap(a + bytes, q + H - a);
    madvise(a, bytes, MADV_HUGEPAGE);
    madvise(a, bytes, 23);
    p= reinterpret_cast<unsigned*>(a), mapped= bytes;
    return;
   }
  }
#endif
  p= static_cast<unsigned*>(::operator new(n * sizeof(unsigned)));
 }
 ~Region() {
#ifdef __linux__
  if(mapped) {
   munmap(p, mapped);
   return;
  }
#endif
  ::operator delete(p);
 }
 Region(const Region&)= delete;
 Region& operator=(const Region&)= delete;
};
template <class T> constexpr T zero() {
 if constexpr(requires { T::zero(); }) return T::zero();
 else return T(0);
}
template <class T> constexpr T one() {
 if constexpr(requires { T::one(); }) return T::one();
 else return T(1);
}
// 行列の積の速い道の 1 行分 c[j] += T::fast_mul(x, r[j]) (j < n)。c と r が重ならないことを伝えないと、x64 の GCC 15 は
// 64 bit の値の min のループをベクトル化しない。
template <class T> [[gnu::always_inline]] inline void fast_axpy(T* __restrict c, const T* __restrict r, T x, int n) {
 for(int j= 0; j < n; ++j) c[j]+= T::fast_mul(x, r[j]);
}
#ifdef __x86_64__
// MinPlus<int> と MaxPlus<int> の速い道の積 c = a b を、4 行 × 16 列の和を 8 本の ymm に持つタイルで求める。a は h4 行 × m 列
// (h4 は h を 4 の倍数に切り上げたもので、広げた行は 0 の値 z)、b は m 行 × w16 列 (w16 は w を 16 の倍数に切り上げたもので、広げた
// 列は z)、c は h 行 × w 列 (行の間隔 w) に書く。MIN なら min、そうでなければ max で足す。和 x + y は範囲の判定を通った値どうしなので
// int を溢れない。b は 16 列の帯ごとに k の向きに読む (辺 300 なら帯は 19 KB で L1 に残る)。k のループは 2 回ずつ展開し、k ごとに
// 8 本の和を空の asm に通す。通さないと x64 の GCC 15 は和を毎回別のレジスタに求めてから戻す vmovdqa を 8 本付け、1 本はスタックにも
// 書いていた (procon-judge の cf-147-b の maxplus_tile)。
template <bool MIN> [[gnu::noinline]] inline void semiring_tile(const int* __restrict a, const int* __restrict b, int* __restrict c, int h, int m, int w, int w16, int z) {
 const __m256i Z= _mm256_set1_epi32(z);
 auto op= [](__m256i x, __m256i y) {
  if constexpr(MIN) return _mm256_min_epi32(x, y);
  else return _mm256_max_epi32(x, y);
 };
 for(int j0= 0; j0 < w16; j0+= 16)
  for(int i0= 0; i0 < h; i0+= 4) {
   const int *x0= a + std::size_t(i0) * m, *x1= x0 + m, *x2= x1 + m, *x3= x2 + m;
   const int* yb= b + j0;
   __m256i c00= Z, c01= Z, c10= Z, c11= Z, c20= Z, c21= Z, c30= Z, c31= Z;
#pragma GCC unroll 2
   for(int k= 0; k < m; ++k) {
    const __m256i y0= _mm256_loadu_si256((const __m256i*)(yb + std::size_t(k) * w16)), y1= _mm256_loadu_si256((const __m256i*)(yb + std::size_t(k) * w16 + 8));
    __m256i t= _mm256_set1_epi32(x0[k]);
    c00= op(c00, _mm256_add_epi32(t, y0)), c01= op(c01, _mm256_add_epi32(t, y1));
    t= _mm256_set1_epi32(x1[k]);
    c10= op(c10, _mm256_add_epi32(t, y0)), c11= op(c11, _mm256_add_epi32(t, y1));
    t= _mm256_set1_epi32(x2[k]);
    c20= op(c20, _mm256_add_epi32(t, y0)), c21= op(c21, _mm256_add_epi32(t, y1));
    t= _mm256_set1_epi32(x3[k]);
    c30= op(c30, _mm256_add_epi32(t, y0)), c31= op(c31, _mm256_add_epi32(t, y1));
    asm("" : "+x"(c00), "+x"(c01), "+x"(c10), "+x"(c11), "+x"(c20), "+x"(c21), "+x"(c30), "+x"(c31));
   }
   if(i0 + 4 <= h && j0 + 16 <= w) {
    int* cr= c + std::size_t(i0) * w + j0;
    _mm256_storeu_si256((__m256i*)cr, c00), _mm256_storeu_si256((__m256i*)(cr + 8), c01), cr+= w;
    _mm256_storeu_si256((__m256i*)cr, c10), _mm256_storeu_si256((__m256i*)(cr + 8), c11), cr+= w;
    _mm256_storeu_si256((__m256i*)cr, c20), _mm256_storeu_si256((__m256i*)(cr + 8), c21), cr+= w;
    _mm256_storeu_si256((__m256i*)cr, c30), _mm256_storeu_si256((__m256i*)(cr + 8), c31);
   } else {  // 端のタイルは、写してから入る分だけ書く
    alignas(32) int tmp[64];
    _mm256_store_si256((__m256i*)tmp, c00), _mm256_store_si256((__m256i*)(tmp + 8), c01);
    _mm256_store_si256((__m256i*)(tmp + 16), c10), _mm256_store_si256((__m256i*)(tmp + 24), c11);
    _mm256_store_si256((__m256i*)(tmp + 32), c20), _mm256_store_si256((__m256i*)(tmp + 40), c21);
    _mm256_store_si256((__m256i*)(tmp + 48), c30), _mm256_store_si256((__m256i*)(tmp + 56), c31);
    for(int r= 0; r < 4 && i0 + r < h; ++r)
     for(int q= 0; q < 16 && j0 + q < w; ++q) c[std::size_t(i0 + r) * w + j0 + q]= tmp[r * 16 + q];
   }
  }
}
#endif
template <class T> struct is_zmod: std::false_type {};
template <unsigned M> struct is_zmod<ZMod<M>>: std::true_type {};
template <class T> constexpr bool is_zmod_v= is_zmod<T>::value;
// ZMod の行列の核。行列は行優先の配列で受ける。
namespace zmat {
using u32= unsigned;
using u64= unsigned long long;
using u128= unsigned __int128;
template <class Z> struct mod_param;
template <unsigned M> struct mod_param<ZMod<M>>: std::integral_constant<unsigned, M> {};
// 2^32 未満の n が素数か。底 2、7、61 の Miller-Rabin (4759123141 未満で正しい)。
constexpr bool is_prime32(u32 n) {
 if(n < 64) return 0x28208a20a08a28acull >> n & 1;
 if(n % 2 == 0 || n % 3 == 0 || n % 5 == 0 || n % 7 == 0) return false;
 u32 d= n - 1;
 int s= 0;
 while(!(d & 1)) d>>= 1, ++s;
 for(u64 b: {2ull, 7ull, 61ull}) {
  u64 x= 1, base= b % n;
  if(!base) continue;
  for(u32 e= d; e; e>>= 1, base= base * base % n)
   if(e & 1) x= x * base % n;
  bool ok= x == 1 || x == n - 1;
  for(int r= 1; r < s && !ok; ++r) x= x * x % n, ok= x == n - 1;
  if(!ok) return false;
 }
 return true;
}
// 法が素数か。ZMod<0> は、最後に調べた法と結果を覚えておく。
template <class Z> bool prime_mod() {
 if constexpr(mod_param<Z>::value != 0) {
  constexpr bool r= is_prime32(mod_param<Z>::value);
  return r;
 } else {
  static u32 m= 0;
  static bool r= false;
  if(m != Z::mod()) m= Z::mod(), r= is_prime32(m);
  return r;
 }
}
// 法 M (1 <= M < 2^30) の定数。r32 = 2^32 mod M。M が奇数なら niv = -M^{-1} mod 2^32、r2 = 2^64 mod M。
struct Mod {
 u32 m, r32, niv= 0, r2= 0;
 bool odd;
 explicit Mod(u32 m): m(m), r32(u32((u64(1) << 32) % m)), odd(m & 1) {
  if(odd) {
   u32 x= m;
   for(int i= 0; i < 5; ++i) x*= 2 - m * x;
   niv= -x, r2= u32(u64(r32) * r32 % m);
  }
 }
 // x 2^{-32} mod M を M 未満で (x < 2^32 M、M は奇数)。
 u32 redc(u64 x) const {
  const u32 r= u32((x + u64(u32(x) * niv) * m) >> 32);
  return r >= m ? r - m : r;
 }
 u32 to_mont(u32 x) const { return redc(u64(x) * r2); }  // x 2^32 mod M (x < M)
};
template <class I> inline u32 vin(const I& x) {
 if constexpr(std::is_same_v<I, u32>) return x;
 else return x.val();
}
// vin の 8 個版。p から I を 8 個読み、M 未満に直す。I が Z なら値は 2M 未満で、M 以上なら M を引く (min(x, x - M))。I は Z か
// u32 で、どちらも 4 byte の値 1 つなので、8 個をそのまま 32 byte で読む。m は vin_mod<I>() (u32 なら使わない)。x64 の GCC 15 は
// 1 個ずつの vin のループをベクトル化しなかった (Z の配列と書き込み先が重なりうるため)。
template <class I> inline __m256i vin_mod() {
 if constexpr(std::is_same_v<I, u32>) return _mm256_setzero_si256();
 else return _mm256_set1_epi32(I::mod());
}
template <class I> [[gnu::always_inline]] inline __m256i vin8(const I* p, __m256i m) {
 static_assert(sizeof(I) == sizeof(u32));
 const __m256i x= _mm256_loadu_si256((const __m256i*)p);
 if constexpr(std::is_same_v<I, u32>) return (void)m, x;
 else return _mm256_min_epu32(x, _mm256_sub_epi32(x, m));
}
template <class O> inline O vout(u32 x) {
 if constexpr(std::is_same_v<O, u32>) return x;
 else return O::raw(x);
}
// u64 の x の剰余を O で。
template <class Z, class O> inline O vout64(u64 x) {
 if constexpr(std::is_same_v<O, u32>) return Z(x).val();
 else return Z(x);
}
// x が 2^63 以上なら 2M 2^32 を引く。M の倍数なので剰余は変わらず、8 個の積 (8 M^2 < 2M 2^32) を足した値も 2^63 未満に戻る。
inline __m256i keep(__m256i x, __m256i m2, __m256i top) { return _mm256_sub_epi64(x, _mm256_min_epu32(m2, _mm256_and_si256(top, x))); }
// 2^63 未満の x を 2^32 の位で畳む: (x >> 32) (2^32 mod M) + (x mod 2^32) < 2^31 (M + 1)。
inline __m256i fold(__m256i x, __m256i r32, __m256i lo) { return _mm256_add_epi64(_mm256_mul_epu32(_mm256_srli_epi64(x, 32), r32), _mm256_and_si256(x, lo)); }
// 偶数番目の列の和 e と奇数番目の列の和 o を Montgomery の還元で 2^32 で割り、8 個の u32 に並べる。e、o が 2^31 (M + 1) 未満なら
// 結果は 1.5M + 1 未満、M^2 未満なら 1.25M 未満。
inline __m256i redc8(__m256i e, __m256i o, __m256i niv, __m256i m) {
 const __m256i ke= _mm256_mul_epu32(_mm256_mul_epu32(e, niv), m), ko= _mm256_mul_epu32(_mm256_mul_epu32(o, niv), m);
 return _mm256_blend_epi32(_mm256_srli_epi64(_mm256_add_epi64(e, ke), 32), _mm256_add_epi64(o, ko), 0xaa);
}
inline __m256i shrk(__m256i x, __m256i m) { return _mm256_min_epu32(x, _mm256_sub_epi32(x, m)); }
// p から読んだ 8 個の u32 z の奇数番目を、64 bit の各々の下の 32 bit に置く。x86 では p + 1 から読み直す (vpmuludq は下の
// 32 bit しか見ないので、ずらす命令が要らない。p の後ろを 4 byte 読む)。SIMDe の道では読み直すほうが遅いので、32 bit ずらす。
[[gnu::always_inline]] inline __m256i odd_lanes(const u32* p, __m256i z) {
#ifdef __x86_64__
 (void)z;
 return _mm256_loadu_si256((const __m256i*)(p + 1));
#else
 (void)p;
 return _mm256_srli_epi64(z, 32);
#endif
}
inline __m256i dilt(__m256i x, __m256i m) { return _mm256_min_epu32(x, _mm256_add_epi32(x, m)); }
// C (行 n、列 p、行の間隔 ldc) = A B を 4 行 × 8 列のタイルで求める。a は A の値 (M 未満) を、行を 4 の倍数に、列を 8 の倍数 m8 に
// 0 で広げた行優先の配列。bp は B を 8 列ずつの帯 [帯][k][8] (k < m8) に詰め、列を 8 の倍数に 0 で広げたもの。B の行の 8 個を
// 1 回読み、偶数番目はそのまま vpmuludq に渡す。奇数番目は x86 では 4 byte ずらした番地から読み直して渡す (odd_lanes。Intel の
// 機種で 3% から 5% 速く、bp の後ろを 4 byte 読む)。u64 の和 8 本 (4 行 × 偶奇) に足し、8 段ごとに keep で 2^63 未満に保つ。最後に、MODE 0 (偶数の法) は和を 1 つずつ Barrett で還元する。奇数の法は畳んでから Montgomery の
// 還元で 2^32 で割り、MODE 1 はさらに 2^64 mod M を Montgomery の掛け算で掛けて正しい値に戻す。MODE 2 は 2^{-32} 倍のまま
// 返す (Montgomery の表現どうしの積)。値は M 未満で書く。O は u32 か Z。
template <int MODE, class Z, class O> void tile(const u32* a, const u32* bp, int n, int m8, int p, O* c, std::size_t ldc, const Mod& md) {
 const __m256i M1= _mm256_set1_epi32(md.m), M2= _mm256_set1_epi32(2 * md.m), TOP= _mm256_set1_epi64x((long long)(1ull << 63));
 const __m256i NIV= _mm256_set1_epi32(md.niv), R32= _mm256_set1_epi64x(md.r32), R2= _mm256_set1_epi64x(md.r2), LO= _mm256_set1_epi64x(0xffffffffll);
 alignas(32) u32 t[8];
 alignas(32) u64 te[4], to[4];
 for(int jb= 0; jb < p; jb+= 8) {
  const u32* band= bp + std::size_t(jb) * m8;
  const int lim= std::min(8, p - jb);
  for(int i= 0; i < n; i+= 4) {
   const u32 *a0= a + std::size_t(i) * m8, *a1= a0 + m8, *a2= a1 + m8, *a3= a2 + m8;
   __m256i e0= _mm256_setzero_si256(), o0= e0, e1= e0, o1= e0, e2= e0, o2= e0, e3= e0, o3= e0;
   for(int k= 0; k < m8; k+= 8) {
    for(int s= k; s < k + 8; ++s) {
     const __m256i z= _mm256_loadu_si256((const __m256i*)(band + s * 8)), zo= odd_lanes(band + s * 8, z);
     __m256i x= _mm256_set1_epi32(a0[s]);
     e0= _mm256_add_epi64(e0, _mm256_mul_epu32(x, z)), o0= _mm256_add_epi64(o0, _mm256_mul_epu32(x, zo));
     x= _mm256_set1_epi32(a1[s]);
     e1= _mm256_add_epi64(e1, _mm256_mul_epu32(x, z)), o1= _mm256_add_epi64(o1, _mm256_mul_epu32(x, zo));
     x= _mm256_set1_epi32(a2[s]);
     e2= _mm256_add_epi64(e2, _mm256_mul_epu32(x, z)), o2= _mm256_add_epi64(o2, _mm256_mul_epu32(x, zo));
     x= _mm256_set1_epi32(a3[s]);
     e3= _mm256_add_epi64(e3, _mm256_mul_epu32(x, z)), o3= _mm256_add_epi64(o3, _mm256_mul_epu32(x, zo));
    }
    e0= keep(e0, M2, TOP), o0= keep(o0, M2, TOP), e1= keep(e1, M2, TOP), o1= keep(o1, M2, TOP);
    e2= keep(e2, M2, TOP), o2= keep(o2, M2, TOP), e3= keep(e3, M2, TOP), o3= keep(o3, M2, TOP);
   }
   auto put= [&](__m256i e, __m256i o, int r) {
    if(i + r >= n) return;
    O* dst= c + std::size_t(i + r) * ldc + jb;
    if constexpr(MODE == 0) {
     _mm256_store_si256((__m256i*)te, e), _mm256_store_si256((__m256i*)to, o);
     for(int j= 0; j < lim; ++j) dst[j]= vout64<Z, O>(j & 1 ? to[j >> 1] : te[j >> 1]);
    } else {
     __m256i v= shrk(redc8(fold(e, R32, LO), fold(o, R32, LO), NIV, M1), M1);
     if constexpr(MODE == 1) v= shrk(redc8(_mm256_mul_epu32(v, R2), _mm256_mul_epu32(_mm256_srli_epi64(v, 32), R2), NIV, M1), M1);
     if(lim == 8) _mm256_storeu_si256((__m256i*)dst, v);
     else {
      _mm256_store_si256((__m256i*)t, v);
      for(int j= 0; j < lim; ++j) dst[j]= vout<O>(t[j]);
     }
    }
   };
   put(e0, o0, 0), put(e1, o1, 1), put(e2, o2, 2), put(e3, o3, 3);
  }
 }
}
// A (行 h、列 w、行の間隔 ld) の値を M 未満に直し、行の間隔 w8 の a に写す。a の残り (広げた分) は 0 にしておく。
template <class I> void load_rows(const I* A, int h, int w, std::size_t ld, u32* a, int w8) {
 const __m256i m= vin_mod<I>();
 for(int i= 0; i < h; ++i) {
  const I* src= A + std::size_t(i) * ld;
  u32* dst= a + std::size_t(i) * w8;
  int j= 0;
  for(; j + 8 <= w; j+= 8) _mm256_storeu_si256((__m256i*)(dst + j), vin8(src + j, m));
  for(; j < w; ++j) dst[j]= vin(src[j]);
 }
}
// B (行 h、列 w、行の間隔 ld) の値を M 未満に直し、8 列ずつの帯 [帯][k][8] (k < h8) に詰める。bp の残りは 0 にしておく。
template <class I> void pack_bands(const I* B, int h, int w, std::size_t ld, u32* bp, int h8) {
 const __m256i m= vin_mod<I>();
 const int w8= w & ~7;
 for(int k= 0; k < h; ++k) {
  const I* src= B + std::size_t(k) * ld;
  for(int j= 0; j < w8; j+= 8) _mm256_storeu_si256((__m256i*)(bp + std::size_t(j) * h8 + k * 8), vin8(src + j, m));
  for(int j= w8; j < w; ++j) bp[std::size_t(w8) * h8 + k * 8 + (j & 7)]= vin(src[j]);
 }
}
// C = A B をタイルの核で。A、B は I (u32 なら M 未満) の行優先、C は O の行優先。
template <class Z, class I, class O> void tile_mul(const I* A, const I* B, O* C, int n, int m, int p, const Mod& md) {
 const int n4= (n + 3) & ~3, m8= (m + 7) & ~7, p8= (p + 7) & ~7;
 std::vector<u32> a(std::size_t(n4) * m8), bp(std::size_t(p8) * m8 + 8);  // bp の後ろは、奇数番目の列を 4 byte ずらして読む分
 load_rows(A, n, m, m, a.data(), m8), pack_bands(B, m, p, p, bp.data(), m8);
 if(md.odd) tile<1, Z>(a.data(), bp.data(), n, m8, p, C, p, md);
 else tile<0, Z>(a.data(), bp.data(), n, m8, p, C, p, md);
}
// 小さい行列の積。i-k-j の順に C の行を u64 の和で持つ。値は M 未満に直すので、16 個の積までは u64 に収まる。そのあとは
// 8 段ごとに 2^32 の位を 2^32 mod M に掛けて畳む (2^62 + 2^32 未満になり、8 個の積を足しても収まる)。作業の配列は、小さければ
// スタックに置く。
template <class Z> void small_mul(const Z* A, const Z* B, Z* C, int n, int m, int p) {
 constexpr int SB_= 1024, SA_= 128;
 u32 bs[SB_];
 u64 as[SA_];
 std::vector<u32> bv;
 std::vector<u64> av;
 u32* b= bs;
 u64* s= as;
 if(std::size_t(m) * p > SB_) bv.resize(std::size_t(m) * p), b= bv.data();
 if(p > SA_) av.resize(p), s= av.data();
 for(std::size_t i= 0, e= std::size_t(m) * p; i < e; ++i) b[i]= B[i].val();
 const u64 R= m > 16 ? (u64(1) << 32) % Z::mod() : 0;
 for(int i= 0; i < n; ++i) {
  std::fill(s, s + p, 0);
  for(int k= 0; k < m; ++k) {
   const u64 x= A[std::size_t(i) * m + k].val();
   const u32* y= b + std::size_t(k) * p;
   for(int j= 0; j < p; ++j) s[j]+= x * y[j];
   if((k & 7) == 7 && k >= 15 && k + 1 < m)
    for(int j= 0; j < p; ++j) s[j]= (s[j] >> 32) * R + u32(s[j]);
  }
  for(int j= 0; j < p; ++j) C[std::size_t(i) * p + j]= Z(s[j]);
 }
}
// Strassen の再帰。1 辺 SB のブロックを Morton 順に並べ、ブロックの積 c = a b 2^{-32} を 4 行 × 8 列のタイル (tile と同じく、
// 奇数番目の列は 4 byte ずらして読む。b の後ろは同じ領域の続きなので、はみ出しても読める) で求める。
// 再帰の中は足し算と引き算だけなので、どのブロックも同じ 2^{-32} 倍を持ち、最後に 2^64 mod M を掛けて戻す。値は M 未満。
// 7 回の積は Winograd の形で、1 段の足し算と引き算は 15 回 (並べ方は Boyer、Dumas、Pernet、Zhou (2009) の表 1 の、作業の
// 領域 X、Y の 2 つで C = A B を求める順)。作業の配列は 1 つの領域にまとめて huge page を頼む (Region)。
// procon-judge の yosupo-matrix-product の zmod_winograd_hpo (_strassen_hp.hpp) の法を実行時の奇数に広げたもの。
// SB は 64。32 にすると再帰が 1 段深くなり、積和が 7/8 になる代わりに足し算と引き算が 1 段分増える (procon-judge の
// yosupo-matrix-product で、1024 次の積を lib-neo-core32 と lib-neo-core で比べている)。
template <int SB> [[gnu::noinline]] void s_kernel(const u32* __restrict a, const u32* __restrict b, u32* __restrict c, const Mod& md) {
 const __m256i M1= _mm256_set1_epi32(md.m), M2= _mm256_set1_epi32(2 * md.m), TOP= _mm256_set1_epi64x((long long)(1ull << 63));
 const __m256i NIV= _mm256_set1_epi32(md.niv), R32= _mm256_set1_epi64x(md.r32), LO= _mm256_set1_epi64x(0xffffffffll);
 for(int i= 0; i < SB; i+= 4) {
  const u32 *a0= a + i * SB, *a1= a0 + SB, *a2= a1 + SB, *a3= a2 + SB;
  for(int j= 0; j < SB; j+= 8) {
   __m256i e0= _mm256_setzero_si256(), o0= e0, e1= e0, o1= e0, e2= e0, o2= e0, e3= e0, o3= e0;
   for(int k= 0; k < SB; k+= 8) {
    for(int w= k; w < k + 8; ++w) {
     const __m256i z= _mm256_loadu_si256((const __m256i*)(b + w * SB + j)), zo= odd_lanes(b + w * SB + j, z);
     __m256i x= _mm256_set1_epi32(a0[w]);
     e0= _mm256_add_epi64(e0, _mm256_mul_epu32(x, z)), o0= _mm256_add_epi64(o0, _mm256_mul_epu32(x, zo));
     x= _mm256_set1_epi32(a1[w]);
     e1= _mm256_add_epi64(e1, _mm256_mul_epu32(x, z)), o1= _mm256_add_epi64(o1, _mm256_mul_epu32(x, zo));
     x= _mm256_set1_epi32(a2[w]);
     e2= _mm256_add_epi64(e2, _mm256_mul_epu32(x, z)), o2= _mm256_add_epi64(o2, _mm256_mul_epu32(x, zo));
     x= _mm256_set1_epi32(a3[w]);
     e3= _mm256_add_epi64(e3, _mm256_mul_epu32(x, z)), o3= _mm256_add_epi64(o3, _mm256_mul_epu32(x, zo));
    }
    e0= keep(e0, M2, TOP), o0= keep(o0, M2, TOP), e1= keep(e1, M2, TOP), o1= keep(o1, M2, TOP);
    e2= keep(e2, M2, TOP), o2= keep(o2, M2, TOP), e3= keep(e3, M2, TOP), o3= keep(o3, M2, TOP);
   }
   _mm256_storeu_si256((__m256i*)(c + i * SB + j), shrk(redc8(fold(e0, R32, LO), fold(o0, R32, LO), NIV, M1), M1));
   _mm256_storeu_si256((__m256i*)(c + (i + 1) * SB + j), shrk(redc8(fold(e1, R32, LO), fold(o1, R32, LO), NIV, M1), M1));
   _mm256_storeu_si256((__m256i*)(c + (i + 2) * SB + j), shrk(redc8(fold(e2, R32, LO), fold(o2, R32, LO), NIV, M1), M1));
   _mm256_storeu_si256((__m256i*)(c + (i + 3) * SB + j), shrk(redc8(fold(e3, R32, LO), fold(o3, R32, LO), NIV, M1), M1));
  }
 }
}
inline void s_add(u32* x, const u32* y, const u32* z, std::size_t d, __m256i m) {
 for(std::size_t i= 0; i < d; i+= 8) _mm256_storeu_si256((__m256i*)(x + i), shrk(_mm256_add_epi32(_mm256_loadu_si256((const __m256i*)(y + i)), _mm256_loadu_si256((const __m256i*)(z + i))), m));
}
inline void s_sub(u32* x, const u32* y, const u32* z, std::size_t d, __m256i m) {
 for(std::size_t i= 0; i < d; i+= 8) _mm256_storeu_si256((__m256i*)(x + i), dilt(_mm256_sub_epi32(_mm256_loadu_si256((const __m256i*)(y + i)), _mm256_loadu_si256((const __m256i*)(z + i))), m));
}
// 1 辺 n のブロック列 (Morton 順) の積 c = a b 2^{-32}。S1 = A21 + A22、S2 = S1 - A11、S3 = A11 - A21、S4 = A12 - S2、
// T1 = B12 - B11、T2 = B22 - T1、T3 = B22 - B12、T4 = T2 - B21、P1 = A11 B11、P2 = A12 B21、P3 = S4 B22、P4 = A22 T4、
// P5 = S1 T1、P6 = S2 T2、P7 = S3 T3、U2 = P1 + P6、U3 = U2 + P7、U4 = U2 + P5 として C11 = P1 + P2、C12 = U4 + P3、
// C21 = U3 - P4、C22 = U3 + P5。X、Y は 1 辺 n / 2 の作業の領域で、深い段はその後ろを使う (どちらも n^2 / 3 ほど)。
template <int SB> void s_mul(const u32* a, const u32* b, u32* c, std::size_t n, u32* X, u32* Y, const Mod& md) {
 if(n == SB) return s_kernel<SB>(a, b, c, md);
 const __m256i m= _mm256_set1_epi32(md.m);
 const std::size_t h= n / 2, d= h * h;
 const u32 *a11= a, *a12= a + d, *a21= a + 2 * d, *a22= a + 3 * d;
 const u32 *b11= b, *b12= b + d, *b21= b + 2 * d, *b22= b + 3 * d;
 u32 *c11= c, *c12= c + d, *c21= c + 2 * d, *c22= c + 3 * d, *X2= X + d, *Y2= Y + d;
 s_sub(X, a11, a21, d, m), s_sub(Y, b22, b12, d, m), s_mul<SB>(X, Y, c21, h, X2, Y2, md);  // C21 = P7
 s_add(X, a21, a22, d, m), s_sub(Y, b12, b11, d, m), s_mul<SB>(X, Y, c22, h, X2, Y2, md);  // C22 = P5
 s_sub(X, X, a11, d, m), s_sub(Y, b22, Y, d, m), s_mul<SB>(X, Y, c12, h, X2, Y2, md);  // X = S2、Y = T2、C12 = P6
 s_sub(X, a12, X, d, m), s_mul<SB>(X, b22, c11, h, X2, Y2, md);  // X = S4、C11 = P3
 s_mul<SB>(a11, b11, X, h, X2, Y2, md);  // X = P1
 s_add(c12, X, c12, d, m), s_add(c21, c12, c21, d, m), s_add(c12, c12, c22, d, m), s_add(c22, c21, c22, d, m), s_add(c12, c12, c11, d, m);
 s_sub(Y, Y, b21, d, m), s_mul<SB>(a22, Y, c11, h, X2, Y2, md), s_sub(c21, c21, c11, d, m);  // Y = T4、C11 = P4、C21 = U6
 s_mul<SB>(a12, b21, c11, h, X2, Y2, md), s_add(c11, X, c11, d, m);  // C11 = P2、U1
}
// 行 rows、列 cols の行列 src (行の間隔 ld) の値を M 未満に直し、1 辺 n のブロック列 (Morton 順) に写す。範囲の外は 0。
template <int SB, class Z> void s_place(const Z* src, int rows, int cols, std::size_t ld, std::size_t n, u32* dst) {
 if(n == SB) {
  const __m256i m= vin_mod<Z>();
  for(int i= 0; i < SB; ++i) {
   u32* d= dst + i * SB;
   const int w= i < rows ? std::clamp(cols, 0, SB) : 0;
   if(w > 0) {
    const Z* r= src + i * ld;
    int j= 0;
    for(; j + 8 <= w; j+= 8) _mm256_storeu_si256((__m256i*)(d + j), vin8(r + j, m));
    for(; j < w; ++j) d[j]= r[j].val();
   }
   std::fill(d + w, d + SB, 0u);
  }
  return;
 }
 const std::size_t h= n / 2, d= h * h;
 s_place<SB>(src, rows, cols, ld, h, dst);
 s_place<SB>(src + h, rows, int(cols - h), ld, h, dst + d);
 s_place<SB>(src + h * ld, int(rows - h), cols, ld, h, dst + 2 * d);
 s_place<SB>(src + h * ld + h, int(rows - h), int(cols - h), ld, h, dst + 3 * d);
}
// ブロック列 src (1 辺 n、値は a b 2^{-32}) に 2^64 mod M を Montgomery の掛け算で掛け、行列 dst (行の間隔 ld) の範囲内に書く。
template <int SB, class Z> void s_unplace(const u32* src, int rows, int cols, std::size_t ld, std::size_t n, Z* dst, const Mod& md) {
 if(rows <= 0 || cols <= 0) return;
 if(n == SB) {
  const __m256i r2= _mm256_set1_epi64x(md.r2), m= _mm256_set1_epi32(md.m), niv= _mm256_set1_epi32(md.niv);
  alignas(32) u32 t[8];
  for(int i= 0; i < std::min(rows, SB); ++i)
   for(int j= 0; j < std::min(cols, SB); j+= 8) {
    const __m256i x= _mm256_loadu_si256((const __m256i*)(src + i * SB + j));
    const __m256i y= shrk(redc8(_mm256_mul_epu32(x, r2), _mm256_mul_epu32(_mm256_srli_epi64(x, 32), r2), niv, m), m);
    if(j + 8 <= cols) _mm256_storeu_si256((__m256i*)(dst + i * ld + j), y);
    else {
     _mm256_store_si256((__m256i*)t, y);
     for(int q= 0; j + q < cols; ++q) dst[i * ld + j + q]= Z::raw(t[q]);
    }
   }
  return;
 }
 const std::size_t h= n / 2, d= h * h;
 s_unplace<SB>(src, rows, cols, ld, h, dst, md);
 s_unplace<SB>(src + d, rows, int(cols - h), ld, h, dst + h, md);
 s_unplace<SB>(src + 2 * d, int(rows - h), cols, ld, h, dst + h * ld, md);
 s_unplace<SB>(src + 3 * d, int(rows - h), int(cols - h), ld, h, dst + h * ld + h, md);
}
// Strassen の 1 辺 (2 の冪に広げたもの)。タイルの核より積和が少ないと見込めるときだけ 0 より大きい値を返す。
// 積和の数は 7^d 64^3 (d はブロックの 1 辺を 64 としたときの再帰の段数) で、足し算と並べ替えの分を 1.25 倍と見て、タイルの
// n p m (広げた分を含む) と比べる。
inline std::size_t strassen_side(int n, int m, int p) {
 constexpr int B= 64;
 std::size_t N= B, c= std::size_t(B) * B * B;
 while(N < std::size_t(std::max({n, m, p}))) N*= 2, c*= 7;
 const double tile= double((n + 3) & ~3) * ((m + 7) & ~7) * ((p + 7) & ~7);
 return N >= 4 * B && 1.25 * double(c) < tile ? N : 0;
}
template <int SB, class Z> void strassen_mul(const Z* a, const Z* b, Z* c, int n, int m, int p, std::size_t N, const Mod& md) {
 const std::size_t NN= N * N, W= NN / 3 + SB * SB;
 Region R(3 * NN + 2 * W);  // A、B、C、X、Y
 u32 *A= R.p, *B= A + NN, *C= B + NN, *X= C + NN, *Y= X + W;
 s_place<SB>(a, n, m, m, N, A), s_place<SB>(b, m, p, p, N, B);
 s_mul<SB>(A, B, C, N, X, Y, md);
 s_unplace<SB>(C, n, p, p, N, c, md);
}
// C = A B (A は n × m、B は m × p、C は n × p)。
template <class Z> void mul(const Z* A, const Z* B, Z* C, int n, int m, int p) {
 if(n == 0 || p == 0) return;
 // タイルの核は行を 4 の倍数に、列と k を 8 の倍数に広げるので、広げた分が 2.5 倍を超えるならスカラーで掛ける。
 const double vol= double(n) * m * p, padded= double((n + 3) & ~3) * ((m + 7) & ~7) * ((p + 7) & ~7);
 if(vol < 128 || padded > 2.5 * vol) return small_mul(A, B, C, n, m, p);
 const Mod md(Z::mod());
 if(md.odd)
  if(const std::size_t N= strassen_side(n, m, p)) return strassen_mul<64>(A, B, C, n, m, p, N, md);
 tile_mul<Z>(A, B, C, n, m, p, md);
}
// 列ベクトルとの積 r = A v。値を M 未満に直して、行ごとの内積を 16 個ずつ u64 にためる。
template <class Z> void mul_vec(const Z* A, const Z* v, Z* r, int n, int m) {
 std::vector<u32> x(m);
 for(int j= 0; j < m; ++j) x[j]= v[j].val();
 for(int i= 0; i < n; ++i) {
  const Z* ai= A + std::size_t(i) * m;
  Z t;
  for(int j0= 0; j0 < m; j0+= 16) {
   u64 s= 0;
   for(int j= j0, je= std::min(m, j0 + 16); j < je; ++j) s+= u64(ai[j].val()) * x[j];
   t+= Z(s);
  }
  r[i]= t;
 }
}
// 以下は Frobenius 標準形で A^K を求める道。procon-judge の yosupo-pow-of-matrix の zmod_frob (_frob.hpp) の法を、静的な法と
// ZMod<0> の素数に広げたもの。剰余は Z の変換 (静的な法なら定数の割り算、ZMod<0> なら Barrett) で取る。
// y[i] += c x[i] (i < len)。x[i] < M、c < M。
inline void axpy(u64* __restrict y, const u32* __restrict x, u32 c, int len) {
 const __m256i C= _mm256_set1_epi64x(c);
 int i= 0;
 for(; i + 4 <= len; i+= 4) {
  const __m256i v= _mm256_cvtepu32_epi64(_mm_loadu_si128((const __m128i*)(x + i)));
  _mm256_storeu_si256((__m256i*)(y + i), _mm256_add_epi64(_mm256_loadu_si256((const __m256i*)(y + i)), _mm256_mul_epu32(v, C)));
 }
 for(; i < len; ++i) y[i]+= u64(c) * x[i];
}
inline void fold_all(u64* y, int len, u64 r32) {
 const __m256i R= _mm256_set1_epi64x(r32), LO= _mm256_set1_epi64x(0xffffffffll);
 int i= 0;
 for(; i + 4 <= len; i+= 4) {
  const __m256i v= _mm256_loadu_si256((const __m256i*)(y + i));
  _mm256_storeu_si256((__m256i*)(y + i), _mm256_add_epi64(_mm256_mul_epu32(_mm256_srli_epi64(v, 32), R), _mm256_and_si256(v, LO)));
 }
 for(; i < len; ++i) y[i]= (y[i] >> 32) * r32 + u32(y[i]);
}
template <class Z> struct Frob {
 static u32 rd(u64 x) { return Z(x).val(); }
 static u64 fold1(u64 x, u64 r32) { return (x >> 32) * r32 + u32(x); }
 // sum_i x[i] y[i] mod M (x[i], y[i] < M)。最後は 4 本の和を 2 回畳んで (各 2^60 + 2^32 未満) u64 のまま足し、剰余を 1 回だけ
 // 取る。Z の足し算 (2M 以上なら引く) と val() で 4 本をまとめると、x64 の GCC はその比べを分岐にし、値が乱数に近いので予測を外す
 // (Krylov の鎖を作るあいだに n^2 回呼ぶので、yosupo-pow-of-matrix の EPYC 7763 で合計が 7% 増えていた)。
 static u32 dot(const u32* x, const u32* y, int n, u64 r32) {
  const __m256i Rv= _mm256_set1_epi64x(r32), LO= _mm256_set1_epi64x(0xffffffffll);
  __m256i acc= _mm256_setzero_si256();
  int i= 0, t= 0;
  for(; i + 4 <= n; i+= 4) {
   acc= _mm256_add_epi64(acc, _mm256_mul_epu32(_mm256_cvtepu32_epi64(_mm_loadu_si128((const __m128i*)(x + i))), _mm256_cvtepu32_epi64(_mm_loadu_si128((const __m128i*)(y + i)))));
   if(++t == 8) acc= _mm256_add_epi64(_mm256_mul_epu32(_mm256_srli_epi64(acc, 32), Rv), _mm256_and_si256(acc, LO)), t= 0;
  }
  acc= _mm256_add_epi64(_mm256_mul_epu32(_mm256_srli_epi64(acc, 32), Rv), _mm256_and_si256(acc, LO));
  acc= _mm256_add_epi64(_mm256_mul_epu32(_mm256_srli_epi64(acc, 32), Rv), _mm256_and_si256(acc, LO));
  alignas(32) u64 s[4];
  _mm256_store_si256((__m256i*)s, acc);
  u64 q= s[0] + s[1] + s[2] + s[3];
  for(; i < n; ++i) q+= u64(x[i]) * y[i];
  return rd(q);
 }
 // t[0..len) (畳んだ u64) を、f (monic、次数 d) で割った余りに直して r[0..d) (M 未満) に書く。nf[j] = -f[j] mod M。
 static void reduce(std::vector<u64>& t, int len, const std::vector<u32>& nf, int d, std::vector<u32>& r, u64 r32) {
  int cnt= 0;
  for(int k= len - 1; k >= d; --k) {
   const u32 c= rd(t[k]);
   if(c) {
    u64* q= &t[k - d];
    for(int j= 0; j < d; ++j) q[j]+= u64(c) * nf[j];
   }
   if(++cnt == 8) {
    for(int j= std::max(0, k - d); j < k; ++j) t[j]= fold1(t[j], r32);
    cnt= 0;
   }
  }
  for(int j= 0; j < d; ++j) r[j]= rd(t[j]);
 }
 // x^K mod f (K >= 1、f は monic で次数 d)。上の bit から 2 乗と x 倍を重ねる。
 static std::vector<u32> xpow(u64 K, const std::vector<u32>& f, int d, u64 r32) {
  const u32 M= Z::mod();
  std::vector<u32> nf(d), r(d);
  for(int j= 0; j < d; ++j) nf[j]= f[j] ? M - f[j] : 0;
  std::vector<u64> t(2 * d);
  r[0]= 1 % M;
  for(int b= 63 - std::countl_zero(K); b >= 0; --b) {
   std::fill(t.begin(), t.end(), 0);
   int cnt= 0;
   for(int i= 0; i < d; ++i) {
    if(!r[i]) continue;
    const u64 a= r[i];
    u64* q= &t[i];
    for(int j= 0; j < d; ++j) q[j]+= a * r[j];
    if(++cnt == 8) {
     for(auto& v: t) v= fold1(v, r32);
     cnt= 0;
    }
   }
   for(auto& v: t) v= fold1(v, r32);
   reduce(t, 2 * d - 1, nf, d, r, r32);
   if(K >> b & 1) {  // r = x r mod f
    const u64 c= r[d - 1];
    for(int j= d - 1; j >= 1; --j) r[j]= rd(r[j - 1] + c * nf[j]);
    r[0]= rd(c * nf[0]);
   }
  }
  return r;
 }
 static u64 rng() {
  static u64 s= u64(std::chrono::steady_clock::now().time_since_epoch().count());
  u64 z= (s+= 0x9e3779b97f4a7c15ull);
  z= (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull, z= (z ^ (z >> 27)) * 0x94d049bb133111ebull;
  return z ^ (z >> 31);
 }
 // 乱数のベクトル v から Krylov の鎖 v, Av, A^2 v, ... を作り、それまでに集めたベクトルの張る空間に入ったところで止める。
 // 止まったときの関係 f(A) v = w (w は前の鎖の張る空間の元) から鎖の多項式 f が分かる。w が 0 でなければ、前の鎖ごとに w の
 // 成分 g_c を f で割り (乱数がうまく選べていれば割り切れる)、u = sum_c (g_c / f)(A) v_c を引いた v - u から鎖を作り直す。
 // 全体を鎖に分けると、鎖を並べた基底 T で A = T F T^{-1} (F は同伴行列を並べたブロック対角) になり、A^K = T F^K T^{-1}。
 // F^K のブロックの列 j は x^{K+j} mod f。鎖のベクトルは、要の位置で正規化した掃き出しの行 [e | α] (e = sum_j α_j t_j) として
 // 持ち、新しいベクトルを入れるたびに前の行で消す。T^{-1} の列 piv_k は、e の部分を後ろから単位ベクトルに直したときの α の
 // 部分。やり直しても割り切れなければ false を返す。a、out は n × n (M 未満)。
 static bool pow(int n, u64 K, const u32* a, u32* out) {
  const u32 M= Z::mod();
  const u64 R32= (u64(1) << 32) % M;
  // 行は [e (n 個) | α (n + 1 個)]。α の最後は、空間が埋まったあとの関係の分。
  const int W= 2 * n + 1;
  std::vector<u32> Y(std::size_t(n) * W), T(std::size_t(n) * n);  // T の鎖のベクトル t_j は T[j n .. j n + n)
  std::vector<int> piv(n), len(n);
  std::vector<u64> y(W);
  int r= 0;
  // x (= t_j) を入れる。新しい行になれば true。前の行の張る空間に入れば false で、y[n..n+j] に関係 (y[n+j] = 1) が残る。
  auto insert= [&](const u32* x, int j) -> bool {
   for(int i= 0; i < n; ++i) y[i]= x[i];
   std::fill(y.begin() + n, y.begin() + n + j + 1, 0);
   y[n + j]= 1;
   int cnt= 0;
   for(int k= 0; k < r; ++k) {
    const u32 c= rd(y[piv[k]]);
    if(!c) continue;
    axpy(y.data(), &Y[std::size_t(k) * W], M - c, len[k]);
    if(++cnt == 8) fold_all(y.data(), n + j + 1, R32), cnt= 0;
   }
   int p= -1;
   for(int i= 0; i < n + j + 1; ++i) {
    y[i]= rd(y[i]);
    if(p < 0 && i < n && y[i]) p= i;
   }
   if(p < 0) return false;
   const u64 iv= Z::raw(u32(y[p])).inv().val();
   u32* row= &Y[std::size_t(r) * W];
   for(int i= 0; i < n + j + 1; ++i) row[i]= rd(y[i] * iv);
   std::fill(row + n + j + 1, row + W, 0);
   piv[r]= p, len[r]= n + j + 1, ++r;
   return true;
  };
  std::vector<u32> cur(n), nxt(n);
  // v から鎖を作る (s 番目から)。前の空間に入ったときの鎖の長さ d を返し、y に関係を残す。
  auto gen= [&](const std::vector<u32>& v, int s) -> int {
   r= s;
   std::copy(v.begin(), v.end(), cur.begin());
   for(int j= s;; ++j) {
    if(!insert(cur.data(), j)) return j - s;
    std::copy(cur.begin(), cur.end(), &T[std::size_t(j) * n]);
    for(int i= 0; i < n; ++i) nxt[i]= dot(a + std::size_t(i) * n, cur.data(), n, R32);
    std::swap(cur, nxt);
   }
  };
  struct Block {
   int s, d;
   std::vector<u32> f;  // f[0..d]、f[d] = 1
  };
  std::vector<Block> blocks;
  std::vector<u32> v(n), w(n);
  std::vector<u64> acc(n);
  // s 番目からの鎖を 1 本作る。乱数のベクトルを TB 本試して鎖の最も長いものを選ぶ。前の鎖が A の最小多項式の次数の大きい順に
  // 選べていれば、割り算は割り切れる。法が小さいと乱数のベクトルが最大の次数を外しやすいので、TB を増やす。
  const int TB= M < 1024 ? 4 : 1;
  auto build= [&](int s) -> bool {
   int d= -1;
   bool last= false;
   for(int t= 0; t < TB && d < n - s; ++t) {
    for(auto& x: w) x= rd(rng());
    const int e= gen(w, s);
    if((last= e > d)) d= e, std::swap(v, w);
   }
   if(!last) d= gen(v, s);
   if(d == 0) return false;
   bool clean= true;
   for(int j= 0; j < s && clean; ++j) clean= y[n + j] == 0;
   if(!clean) {  // 前の鎖ごとに g_c = -(関係の成分) を f で割り、u = sum_c h_c(A) v_c を v から引く
    std::vector<u32> f(d + 1);
    for(int i= 0; i <= d; ++i) f[i]= u32(y[n + s + i]);
    std::fill(acc.begin(), acc.end(), 0);
    int cnt= 0;
    for(const auto& B: blocks) {
     std::vector<u64> g(B.d);
     for(int i= 0; i < B.d; ++i) g[i]= y[n + B.s + i] ? M - y[n + B.s + i] : 0;
     for(int deg= B.d - 1; deg >= d; --deg) {  // g -= q x^{deg - d} f
      const u32 q= rd(g[deg]);
      g[deg]= 0;
      if(!q) continue;
      for(int i= 0; i < d; ++i) g[deg - d + i]= rd(g[deg - d + i] + u64(q) * (f[i] ? M - f[i] : 0));
      axpy(acc.data(), &T[std::size_t(B.s + deg - d) * n], q, n);
      if(++cnt == 8) fold_all(acc.data(), n, R32), cnt= 0;
     }
     for(int i= 0; i < std::min(d, B.d); ++i)
      if(rd(g[i])) return false;
    }
    for(int i= 0; i < n; ++i) {
     const u32 t= rd(acc[i]);
     v[i]= v[i] >= t ? v[i] - t : v[i] + M - t;
    }
    d= gen(v, s);
    if(d == 0) return false;
    for(int j= 0; j < s; ++j)
     if(y[n + j]) return false;
   }
   std::vector<u32> f(d + 1);
   for(int i= 0; i <= d; ++i) f[i]= u32(y[n + s + i]);
   blocks.push_back({s, d, f});
   return true;
  };
  // 割り切れなければ、前の鎖の選び方が悪かったとみて最初からやり直す (法が小さければ 8 回、大きければ 2 回まで)。
  for(int g= 0;; ++g) {
   blocks.clear(), r= 0;
   bool ok= true;
   while(r < n && (ok= build(r)));
   if(ok) break;
   if(g + 1 == (M < 1024 ? 8 : 2)) return false;
  }
  // T^{-1}: α'_k = α_k - sum_{kk > k} e_k[piv_kk] α'_kk を k の大きい順に。T^{-1} の列 piv_k が α'_k。
  std::vector<u32> Ap(std::size_t(n) * n);
  for(int k= n - 1; k >= 0; --k) {
   const u32* row= &Y[std::size_t(k) * W];
   for(int j= 0; j < n; ++j) acc[j]= row[n + j];
   int cnt= 0;
   for(int kk= k + 1; kk < n; ++kk) {
    const u32 c= row[piv[kk]];
    if(!c) continue;
    axpy(acc.data(), &Ap[std::size_t(kk) * n], M - c, n);
    if(++cnt == 8) fold_all(acc.data(), n, R32), cnt= 0;
   }
   for(int j= 0; j < n; ++j) Ap[std::size_t(k) * n + j]= rd(acc[j]);
  }
  std::vector<u32> Tm(std::size_t(n) * n), Fk(std::size_t(n) * n), Ti(std::size_t(n) * n), P(std::size_t(n) * n);
  for(int j= 0; j < n; ++j)
   for(int i= 0; i < n; ++i) Tm[std::size_t(i) * n + j]= T[std::size_t(j) * n + i];
  for(int k= 0; k < n; ++k)
   for(int j= 0; j < n; ++j) Ti[std::size_t(j) * n + piv[k]]= Ap[std::size_t(k) * n + j];
  for(const auto& B: blocks) {
   std::vector<u32> col= xpow(K, B.f, B.d, R32), nf(B.d);
   for(int i= 0; i < B.d; ++i) nf[i]= B.f[i] ? M - B.f[i] : 0;
   for(int jj= 0; jj < B.d; ++jj) {
    for(int ii= 0; ii < B.d; ++ii) Fk[std::size_t(B.s + ii) * n + B.s + jj]= col[ii];
    const u64 c= col[B.d - 1];  // col = x col mod f
    for(int ii= B.d - 1; ii >= 1; --ii) col[ii]= rd(col[ii - 1] + c * nf[ii]);
    col[0]= rd(c * nf[0]);
   }
  }
  const Mod md(M);
  tile_mul<Z>(Tm.data(), Fk.data(), P.data(), n, n, n, md);
  tile_mul<Z>(P.data(), Ti.data(), out, n, n, n, md);
  return true;
 }
};
// Frobenius 標準形の道を通すか。手間を積 (n^3) の回数で見て、二分累乗は c = (K の bit の数 - 1) + (1 の bit の数 - 1)、
// Frobenius 標準形は 14 + 200 / n + 12 b / n (b は K の bit の数) とする。手元の M2 (SIMDe) で測った比を、x64 では
// スカラーの部分が 1.7 倍ほど重いと見て直したもの。
inline bool frob_worth(int n, u64 K) {
 const int b= std::bit_width(K), c= b - 1 + std::popcount(K) - 1;
 return n >= 16 && double(c) * n > 14.0 * n + 200 + 12.0 * b;
}
// R = A^K (n × n)。
template <class Z> void pow(const Z* A, int n, u64 K, Z* R) {
 if(n == 0) return;
 if(K == 0) {
  for(int i= 0; i < n; ++i)
   for(int j= 0; j < n; ++j) R[std::size_t(i) * n + j]= Z(i == j);
  return;
 }
 const std::size_t nn= std::size_t(n) * n;
 if(n < 8) {
  std::vector<Z> b(A, A + nn), r(A, A + nn), t(nn);
  for(int i= 63 - std::countl_zero(K) - 1; i >= 0; --i) {
   mul(r.data(), r.data(), t.data(), n, n, n), std::swap(r, t);
   if(K >> i & 1) mul(r.data(), b.data(), t.data(), n, n, n), std::swap(r, t);
  }
  std::copy(r.begin(), r.end(), R);
  return;
 }
 std::vector<u32> a(nn);
 for(std::size_t i= 0; i < nn; ++i) a[i]= A[i].val();
 if(K >= 2 && prime_mod<Z>() && frob_worth(n, K)) {
  std::vector<u32> o(nn);
  if(Frob<Z>::pow(n, K, a.data(), o.data())) {
   for(std::size_t i= 0; i < nn; ++i) R[i]= Z::raw(o[i]);
   return;
  }
 }
 // 二分累乗。上の bit から 2 乗と A 倍を重ねる。法が奇数なら、値を Montgomery の表現 (x 2^32 mod M) で持って掛け続ける。
 const Mod md(Z::mod());
 const int w= (n + 7) & ~7;
 const std::size_t ww= std::size_t(w) * w;
 std::vector<u32> r(ww), t(ww), bp(ww + 8), ap(ww + 8);
 for(int i= 0; i < n; ++i)
  for(int j= 0; j < n; ++j) r[std::size_t(i) * w + j]= md.odd ? md.to_mont(a[std::size_t(i) * n + j]) : a[std::size_t(i) * n + j];
 pack_bands(r.data(), n, n, w, ap.data(), w);
 for(int i= 63 - std::countl_zero(K) - 1; i >= 0; --i) {
  pack_bands(r.data(), n, n, w, bp.data(), w);
  if(md.odd) tile<2, Z>(r.data(), bp.data(), n, w, n, t.data(), w, md);
  else tile<0, Z>(r.data(), bp.data(), n, w, n, t.data(), w, md);
  std::swap(r, t);
  if(K >> i & 1) {
   if(md.odd) tile<2, Z>(r.data(), ap.data(), n, w, n, t.data(), w, md);
   else tile<0, Z>(r.data(), ap.data(), n, w, n, t.data(), w, md);
   std::swap(r, t);
  }
 }
 for(int i= 0; i < n; ++i)
  for(int j= 0; j < n; ++j) R[std::size_t(i) * n + j]= Z::raw(md.odd ? md.redc(r[std::size_t(i) * w + j]) : r[std::size_t(i) * w + j]);
}
}
}
template <class T> struct Matrix {
 int h= 0, w= 0;
 std::vector<T, linalg_internal::HugeAlloc<T>> a;
 Matrix()= default;
 Matrix(int h, int w): h(h), w(w), a((size_t)h * w, linalg_internal::zero<T>()) {}
 Matrix(int h, int w, const T& v): h(h), w(w), a((size_t)h * w, v) {}
 Matrix(const std::vector<std::vector<T>>& v): h(v.size()), w(v.empty() ? 0 : v[0].size()), a((size_t)h * w) {
  for(int i= 0; i < h; ++i) {
   assert((int)v[i].size() == w);
   for(int j= 0; j < w; ++j) a[(size_t)i * w + j]= v[i][j];
  }
 }
 static Matrix identity(int n) {
  Matrix r(n, n);
  for(int i= 0; i < n; ++i) r.a[(size_t)i * n + i]= linalg_internal::one<T>();
  return r;
 }
 int height() const { return h; }
 int width() const { return w; }
 T* operator[](int i) { return a.data() + (size_t)i * w; }
 const T* operator[](int i) const { return a.data() + (size_t)i * w; }
 Matrix& operator+=(const Matrix& r) {
  assert(h == r.h && w == r.w);
  for(size_t i= 0; i < a.size(); ++i) a[i]+= r.a[i];
  return *this;
 }
 Matrix& operator-=(const Matrix& r) {
  assert(h == r.h && w == r.w);
  for(size_t i= 0; i < a.size(); ++i) a[i]-= r.a[i];
  return *this;
 }
 Matrix& operator*=(const T& c) {
  for(auto& x: a) x*= c;
  return *this;
 }
 friend Matrix operator+(Matrix l, const Matrix& r) { return l+= r; }
 friend Matrix operator-(Matrix l, const Matrix& r) { return l-= r; }
 friend Matrix operator*(Matrix l, const T& c) { return l*= c; }
 friend Matrix operator*(const T& c, Matrix r) {
  for(auto& x: r.a) x= c * x;
  return r;
 }
 // 行列の積。T が ZMod なら linalg_internal::zmat::mul、ほかは i-k-j の順に足す。
 friend Matrix operator*(const Matrix& l, const Matrix& r) {
  assert(l.w == r.h);
  Matrix c(l.h, r.w);
  if constexpr(linalg_internal::is_zmod_v<T>) linalg_internal::zmat::mul(l.a.data(), r.a.data(), c.a.data(), l.h, l.w, r.w);
  else {
   // T が速い道 (neo/algebra/semiring.hpp の MinPlus と MaxPlus の fast_one_sided、fast_bounded、fast_mul、fast_fix) を持ち、
   // 両方の行列の値が同じ範囲にあれば、積を fast_mul (INF の比べのない和) にして足し、行ごとに fast_fix で直す (範囲の中では
   // 厳密な積と同じ値)。片側の範囲に入れば、fast_bounded は調べない。
   if constexpr(requires {
                 T::fast_mul(l.a[0], r.a[0]);
                 T::fast_one_sided(l.a.data(), size_t(0));
                 T::fast_bounded(l.a.data(), size_t(0));
                 T::fast_fix(l.a[0], false);
                }) {
    const bool one= T::fast_one_sided(l.a.data(), l.a.size()) && T::fast_one_sided(r.a.data(), r.a.size());
    if(one || (T::fast_bounded(l.a.data(), l.a.size()) && T::fast_bounded(r.a.data(), r.a.size()))) {
#ifdef __x86_64__
     // MinPlus<int> と MaxPlus<int> で、どの辺も 16 以上なら、行と列を広げた int の配列に写してタイルで掛け、答えの行列に直接書く。
     // 写す配列は呼び出しをまたいで使い回す。積ごとに確保して答えも別の配列に書いた 115be19 は、cf-147-b の EPYC 7763 で写し元の
     // 1.65 倍かかっていた (4 KB のページを初めて触る分が、積 1 回あたり 0.6 ms ほどと見ている)。
     if constexpr(requires { T::zero().x; } && sizeof(T) == sizeof(int) && std::is_same_v<std::remove_cvref_t<decltype(T::zero().x)>, int>) {
      if(l.h >= 16 && l.w >= 16 && r.w >= 16) {
       constexpr int z= T::zero().x;
       static linalg_internal::hvec<int> a, b;
       const int h= l.h, h4= (h + 3) & ~3, m= l.w, w= r.w, w16= (w + 15) & ~15;
       if(a.size() < std::size_t(h4) * m) a.resize(std::size_t(h4) * m);
       if(b.size() < std::size_t(m) * w16) b.resize(std::size_t(m) * w16);
       for(int i= 0; i < h; ++i)
        for(int k= 0; k < m; ++k) a[std::size_t(i) * m + k]= l[i][k].x;
       std::fill(a.begin() + std::size_t(h) * m, a.begin() + std::size_t(h4) * m, z);
       for(int k= 0; k < m; ++k) {
        int* bk= b.data() + std::size_t(k) * w16;
        for(int j= 0; j < w; ++j) bk[j]= r[k][j].x;
        std::fill(bk + w, bk + w16, z);
       }
       linalg_internal::semiring_tile<(z > 0)>(a.data(), b.data(), reinterpret_cast<int*>(c.a.data()), h, m, w, w16, z);
       if constexpr(requires { T::fast_fix_all(c.a.data(), c.a.size(), true); }) T::fast_fix_all(c.a.data(), c.a.size(), !one);
       else
        for(auto& v: c.a) v= T::fast_fix(v, !one);
       return c;
      }
     }
#endif
     for(int i= 0; i < l.h; ++i) {
      T* ci= c[i];
      for(int k= 0; k < l.w; ++k) linalg_internal::fast_axpy(ci, r[k], l[i][k], r.w);
      if constexpr(requires { T::fast_fix_all(ci, std::size_t(r.w), true); }) T::fast_fix_all(ci, r.w, !one);
      else
       for(int j= 0; j < r.w; ++j) ci[j]= T::fast_fix(ci[j], !one);
     }
     return c;
    }
   }
   for(int i= 0; i < l.h; ++i) {
    T* ci= c[i];
    for(int k= 0; k < l.w; ++k) {
     const T x= l[i][k];
     const T* rk= r[k];
     for(int j= 0; j < r.w; ++j) ci[j]+= x * rk[j];
    }
   }
  }
  return c;
 }
 Matrix& operator*=(const Matrix& r) { return *this= *this * r; }
 // 列ベクトルとの積 A v。
 friend std::vector<T> operator*(const Matrix& l, const std::vector<T>& v) {
  assert(l.w == (int)v.size());
  std::vector<T> r(l.h, linalg_internal::zero<T>());
  if constexpr(linalg_internal::is_zmod_v<T>) linalg_internal::zmat::mul_vec(l.a.data(), v.data(), r.data(), l.h, l.w);
  else {
   for(int i= 0; i < l.h; ++i) {
    const T* li= l[i];
    T s= linalg_internal::zero<T>();
    for(int j= 0; j < l.w; ++j) s+= li[j] * v[j];
    r[i]= s;
   }
  }
  return r;
 }
 // A^k。T が ZMod なら linalg_internal::zmat::pow、ほかは二分累乗。
 Matrix pow(uint64_t k) const {
  assert(h == w);
  if constexpr(linalg_internal::is_zmod_v<T>) {
   Matrix r(h, h);
   linalg_internal::zmat::pow(a.data(), h, k, r.a.data());
   return r;
  } else {
   Matrix r= identity(h), b= *this;
   for(; k; k>>= 1) {
    if(k & 1) r= r * b;
    if(k > 1) b= b * b;
   }
   return r;
  }
 }
 Matrix transpose() const {
  Matrix r(w, h);
  for(int i= 0; i < h; ++i)
   for(int j= 0; j < w; ++j) r.a[(size_t)j * h + i]= a[(size_t)i * w + j];
  return r;
 }
 friend bool operator==(const Matrix& l, const Matrix& r) { return l.h == r.h && l.w == r.w && l.a == r.a; }
 friend bool operator!=(const Matrix& l, const Matrix& r) { return !(l == r); }
};
