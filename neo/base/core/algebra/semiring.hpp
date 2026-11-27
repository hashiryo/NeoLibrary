#pragma once
#if defined(__x86_64__) && defined(__AVX2__)
#include <immintrin.h>
#endif
#include <cstddef>
#include <limits>
#include <type_traits>
// 半環の構造 (値の型 T と static 関数を並べた struct)。Matrix や SquareMatrix のテンプレート引数にそのまま書き、
// 要素の型は T になる (neo/base/core/algebra/structure.hpp の読み替え)。
//   MinPlus<X>: 和は min、積は + (最短路)。0 は INF、1 は 0。
//   MaxPlus<X>: 和は max、積は + (最長路)。0 は -INF、1 は 0。
//   MaxMin<X>: 和は max、積は min (道の幅の最大)。0 は X の最小、1 は X の最大。
//   OrAnd: 和は or、積は and (到達)。
// 整数の INF は X の最大の半分で、積は片方が INF (MaxPlus は -INF) ならそのまま INF にし、和が INF を超えれば INF に丸める。
// 負の重みがあっても INF は INF のまま。浮動小数点数の INF は infinity。
// 整数の MinPlus と MaxPlus は、行列の積の速い道を持つ。積を INF の比べのない和 x + y (fast_mul) にして和 (min か max) で足し、
// 最後に成分ごとに fast_fix で INF に直す。掛ける 2 つの行列の値がどれも fast_one_sided の範囲 (MinPlus は 0 以上 INF 以下、
// MaxPlus は -INF 以上 0 以下) か、どれも fast_bounded の範囲 (±INF か、絶対値が INF / 4 以下) にあれば、厳密な積と同じ値になる。
namespace semiring_internal {
template <class T> constexpr T inf() {
 if constexpr(std::is_floating_point_v<T>) return std::numeric_limits<T>::infinity();
 else return std::numeric_limits<T>::max() / 2;
}
// 4 byte か 8 byte の整数の値を持つ p と q の n 個のうち、先頭から 32 byte ずつの塊の値 (NEG なら符号を反転したもの) の or が
// mask の bit を持たないかを、x86 で AVX2 が使えれば vptest で調べ、調べた個数を返す (mask の bit があれば -1)。残りは呼び出し元
// がスカラーで調べる。塊は先頭から読むので、直前の積が行列を先頭から 32 byte ずつ書いていれば、書き込みからの受け渡しが効く
// (最後の塊を末尾から重ねて読む形は、3 × 3 の行列では 4 byte ずれて 2 つの書き込みにまたがり、procon-judge の cf-1609-e で
// スカラーの or より遅かった)。使えなければ 0 を返す。
template <class T, bool NEG> inline long long chunks_within(const void* p, const void* q, std::size_t n, T mask) {
#if defined(__x86_64__) && defined(__AVX2__)
 if constexpr(sizeof(T) == 4 || sizeof(T) == 8) {
  constexpr std::size_t S= 32 / sizeof(T);
  if(n < S) return 0;
  const char *a= static_cast<const char*>(p), *b= static_cast<const char*>(q);
  auto ld= [](const char* r) {
   const __m256i v= _mm256_loadu_si256((const __m256i*)r);
   if constexpr(!NEG) return v;
   else if constexpr(sizeof(T) == 4) return _mm256_sub_epi32(_mm256_setzero_si256(), v);
   else return _mm256_sub_epi64(_mm256_setzero_si256(), v);
  };
  const __m256i m= sizeof(T) == 4 ? _mm256_set1_epi32(int(mask)) : _mm256_set1_epi64x((long long)mask);
  __m256i o= _mm256_or_si256(ld(a), ld(b));
  std::size_t i= S;
  for(; i + S <= n; i+= S) o= _mm256_or_si256(o, _mm256_or_si256(ld(a + i * sizeof(T)), ld(b + i * sizeof(T))));
  return _mm256_testz_si256(o, m) ? (long long)i : -1;
 }
#endif
 (void)p, (void)q, (void)n, (void)mask;
 return 0;
}
// 4 byte か 8 byte の整数の値を持つ p の n 個のうち、先頭から 32 byte ずつの塊の値がどれも z (INF か -INF) か、絶対値が L 以下かを、
// x86 で AVX2 が使えれば 8 個 (4 個) ずつ調べ、調べた個数を返す (外れがあれば -1)。残りは呼び出し元がスカラーで調べる。使えなければ
// 0 を返す。スカラーのループは x64 の GCC 15 がベクトル化せず、1 個に 11 命令かかっていた。
template <class T> inline long long chunks_bounded(const void* p, std::size_t n, T z, T L) {
#if defined(__x86_64__) && defined(__AVX2__)
 if constexpr(sizeof(T) == 4 || sizeof(T) == 8) {
  constexpr std::size_t S= 32 / sizeof(T);
  const char* a= static_cast<const char*>(p);
  std::size_t i= 0;
  if constexpr(sizeof(T) == 4) {
   const __m256i Z= _mm256_set1_epi32(int(z)), Lv= _mm256_set1_epi32(int(L)), L2= _mm256_set1_epi32(int(2 * L));
   __m256i in= _mm256_set1_epi32(-1);
   for(; i + S <= n; i+= S) {
    const __m256i v= _mm256_loadu_si256((const __m256i*)(a + i * 4)), t= _mm256_add_epi32(v, Lv);
    in= _mm256_and_si256(in, _mm256_or_si256(_mm256_cmpeq_epi32(v, Z), _mm256_cmpeq_epi32(_mm256_max_epu32(t, L2), L2)));  // t <= 2L (符号なし)
   }
   return _mm256_testc_si256(in, _mm256_set1_epi32(-1)) ? (long long)i : -1;
  } else {
   const __m256i Z= _mm256_set1_epi64x((long long)z), Lv= _mm256_set1_epi64x((long long)L), SB= _mm256_set1_epi64x((long long)(1ull << 63));
   const __m256i L2= _mm256_xor_si256(_mm256_set1_epi64x((long long)(2 * L)), SB);
   __m256i out= _mm256_setzero_si256();
   for(; i + S <= n; i+= S) {
    const __m256i v= _mm256_loadu_si256((const __m256i*)(a + i * 8)), t= _mm256_xor_si256(_mm256_add_epi64(v, Lv), SB);
    out= _mm256_or_si256(out, _mm256_andnot_si256(_mm256_cmpeq_epi64(v, Z), _mm256_cmpgt_epi64(t, L2)));  // t > 2L (符号なし)
   }
   return _mm256_testz_si256(out, out) ? (long long)i : -1;
  }
 }
#endif
 (void)p, (void)n, (void)z, (void)L;
 return 0;
}
// 4 byte か 8 byte の符号付き整数の値を持つ p の n 個のうち、先頭から 32 byte ずつの塊の値を、MIN なら t 以上を、そうでなければ
// t 以下を z に直し、直した個数を返す (x86 で AVX2 が使えるときだけ。使えなければ 0)。残りは呼び出し元が直す。
template <class T, bool MIN> inline std::size_t chunks_fix(void* p, std::size_t n, T t, T z) {
#if defined(__x86_64__) && defined(__AVX2__)
 if constexpr(sizeof(T) == 4 || sizeof(T) == 8) {
  constexpr std::size_t S= 32 / sizeof(T);
  char* a= static_cast<char*>(p);
  const T u= MIN ? T(t - 1) : T(t + 1);  // MIN: v > t - 1、MAX: t + 1 > v
  const __m256i Tv= sizeof(T) == 4 ? _mm256_set1_epi32(int(u)) : _mm256_set1_epi64x((long long)u);
  const __m256i Z= sizeof(T) == 4 ? _mm256_set1_epi32(int(z)) : _mm256_set1_epi64x((long long)z);
  std::size_t i= 0;
  for(; i + S <= n; i+= S) {
   __m256i* q= (__m256i*)(a + i * sizeof(T));
   const __m256i v= _mm256_loadu_si256(q);
   __m256i m;
   if constexpr(sizeof(T) == 4) m= MIN ? _mm256_cmpgt_epi32(v, Tv) : _mm256_cmpgt_epi32(Tv, v);
   else m= MIN ? _mm256_cmpgt_epi64(v, Tv) : _mm256_cmpgt_epi64(Tv, v);
   _mm256_storeu_si256(q, _mm256_blendv_epi8(v, Z, m));
  }
  return i;
 }
#endif
 (void)p, (void)n, (void)t, (void)z;
 return 0;
}
}
template <class X> struct MinPlus {
 using T= X;
 static constexpr T INF= semiring_internal::inf<X>();
 static constexpr T zero() { return INF; }
 static constexpr T one() { return T(0); }
 static constexpr T add(T a, T b) { return a < b ? a : b; }
 // 積。max(x, y, x + y) >= INF なら INF。
 static constexpr T mul(T a, T b) {
  if constexpr(std::is_floating_point_v<T>) return T(a + b);
  else {
   using U= std::make_unsigned_t<T>;
   const T s= T(U(a) + U(b)), t= a > b ? a : b;
   const U m= U(0) - U((t > s ? t : s) >= INF);
   return T((U(s) & ~m) | (U(INF) & m));
  }
 }
 // 行列の積の速い道。fast_one_sided は p の値がすべて 0 以上 INF 以下か。このとき積は min(x + y, INF) と同じで (INF に 0 以上を
 // 足しても INF 以上)、min は単調なので、和 x + y の min を取ってから INF 以上を INF にしてよい。INF は 2^k - 1 の形なので、
 // 範囲に入るかは値の or が INF の外の bit を持たないかで分かる。
 static constexpr bool fast_one_sided(const T* p, std::size_t n) { return fast_one_sided(p, p, n); }
 // p と q の n 個ずつを 1 本のループで調べる (大きさを型で決める行列の積で、判定を 1 回の or の畳み込みにする)。
 static constexpr bool fast_one_sided(const T* p, const T* q, std::size_t n) {
  if constexpr(std::is_floating_point_v<T>) return false;
  else {
   using U= std::make_unsigned_t<T>;
   long long i= 0;
   if(!std::is_constant_evaluated() && (i= semiring_internal::chunks_within<T, false>(p, q, n, T(~U(INF)))) < 0) return false;
   U o= 0;
   for(std::size_t j= i; j < n; ++j) o|= U(p[j]) | U(q[j]);
   return !(o & ~U(INF));
  }
 }
 // fast_bounded は p の値がすべて INF か、絶対値が L = INF / 4 以下か。このとき有限の値どうしの和は絶対値が INF / 2 未満で、
 // INF を含む和は INF - L 以上 (2 INF 以下なので溢れない) なので、和の min を取ってから INF / 2 以上を INF にすれば、厳密な積の
 // min と同じになる。
 static constexpr bool fast_bounded(const T* p, std::size_t n) {
  if constexpr(std::is_floating_point_v<T> || !std::is_signed_v<T>) return false;
  else {
   using U= std::make_unsigned_t<T>;
   constexpr U L= U(INF) >> 2;
   long long i= 0;
   if(!std::is_constant_evaluated() && (i= semiring_internal::chunks_bounded<T>(p, n, INF, T(L))) < 0) return false;
   bool ok= true;
   for(std::size_t j= i; j < n; ++j) ok&= (U(p[j]) == U(INF)) | (U(U(p[j]) + L) <= 2 * L);
   return ok;
  }
 }
 static constexpr T fast_mul(T a, T b) { return T(a + b); }
 // 速い道で足した和を直す。bounded は fast_bounded の範囲で掛けたか。
 static constexpr T fast_fix(T v, bool bounded) {
  const T t= bounded ? T(INF >> 1) : INF;
  return v >= t ? INF : v;
 }
 // p の n 個に fast_fix を当てる (行列の積の速い道の最後。x86 で AVX2 が使えれば 32 byte ずつ)。
 static void fast_fix_all(T* p, std::size_t n, bool bounded) {
  std::size_t i= 0;
  if constexpr(!std::is_floating_point_v<T> && std::is_signed_v<T>) i= semiring_internal::chunks_fix<T, true>(p, n, bounded ? T(INF >> 1) : INF, INF);
  for(; i < n; ++i) p[i]= fast_fix(p[i], bounded);
 }
 static constexpr int fast_tile= 0;  // 要素が int のときの行列の積の AVX2 のタイル (0 は min、1 は max)
};
template <class X> struct MaxPlus {
 using T= X;
 static constexpr T INF= semiring_internal::inf<X>();
 static constexpr T zero() { return T(-INF); }
 static constexpr T one() { return T(0); }
 static constexpr T add(T a, T b) { return a > b ? a : b; }
 // 積。min(x, y, x + y) <= -INF なら -INF。
 static constexpr T mul(T a, T b) {
  if constexpr(std::is_floating_point_v<T>) return T(a + b);
  else {
   using U= std::make_unsigned_t<T>;
   const T s= T(U(a) + U(b)), t= a < b ? a : b;
   const U m= U(0) - U((t < s ? t : s) <= -INF);
   return T((U(s) & ~m) | (U(-INF) & m));
  }
 }
 // 行列の積の速い道。MinPlus の符号を逆にしたもので、fast_one_sided は p の値がすべて -INF 以上 0 以下か、fast_bounded は
 // すべて -INF か、絶対値が INF / 4 以下か。
 static constexpr bool fast_one_sided(const T* p, std::size_t n) { return fast_one_sided(p, p, n); }
 static constexpr bool fast_one_sided(const T* p, const T* q, std::size_t n) {
  if constexpr(std::is_floating_point_v<T> || !std::is_signed_v<T>) return false;
  else {
   using U= std::make_unsigned_t<T>;
   long long i= 0;
   if(!std::is_constant_evaluated() && (i= semiring_internal::chunks_within<T, true>(p, q, n, T(~U(INF)))) < 0) return false;
   U o= 0;
   for(std::size_t j= i; j < n; ++j) o|= (U(0) - U(p[j])) | (U(0) - U(q[j]));
   return !(o & ~U(INF));
  }
 }
 static constexpr bool fast_bounded(const T* p, std::size_t n) {
  if constexpr(std::is_floating_point_v<T> || !std::is_signed_v<T>) return false;
  else {
   using U= std::make_unsigned_t<T>;
   constexpr U L= U(INF) >> 2;
   long long i= 0;
   if(!std::is_constant_evaluated() && (i= semiring_internal::chunks_bounded<T>(p, n, T(-INF), T(L))) < 0) return false;
   bool ok= true;
   for(std::size_t j= i; j < n; ++j) ok&= (U(p[j]) == U(-INF)) | (U(U(p[j]) + L) <= 2 * L);
   return ok;
  }
 }
 static constexpr T fast_mul(T a, T b) { return T(a + b); }
 static constexpr T fast_fix(T v, bool bounded) {
  const T t= bounded ? T(INF >> 1) : INF;
  return v <= -t ? T(-INF) : v;
 }
 static void fast_fix_all(T* p, std::size_t n, bool bounded) {
  std::size_t i= 0;
  if constexpr(!std::is_floating_point_v<T> && std::is_signed_v<T>) i= semiring_internal::chunks_fix<T, false>(p, n, bounded ? T(-(INF >> 1)) : T(-INF), T(-INF));
  for(; i < n; ++i) p[i]= fast_fix(p[i], bounded);
 }
 static constexpr int fast_tile= 1;  // 要素が int のときの行列の積の AVX2 のタイル (0 は min、1 は max)
};
template <class X> struct MaxMin {
 using T= X;
 static constexpr T zero() { return std::numeric_limits<X>::lowest(); }
 static constexpr T one() { return std::numeric_limits<X>::max(); }
 static constexpr T add(T a, T b) { return a > b ? a : b; }
 static constexpr T mul(T a, T b) { return a < b ? a : b; }
};
struct OrAnd {
 using T= bool;
 static constexpr T zero() { return false; }
 static constexpr T one() { return true; }
 static constexpr T add(T a, T b) { return a | b; }
 static constexpr T mul(T a, T b) { return a & b; }
};
