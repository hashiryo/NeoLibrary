#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <cmath>
#include <type_traits>
#include <vector>
#include "neo/base/core/algebra/ZMod.hpp"
// i^d (0 ≤ i ≤ n) の表。0^0 = 1。T は掛け算を持つ型 (ZMod など) で、整数から作れること。
// i^d は完全乗法的なので、素数 p だけ p^d を冪で求め、合成数 m は最小の素因数 p について m^d = (m/p)^d p^d とする。
// √n 以下の素数で 2^15 個ずつの区画を篩い、合成数 m には m/p と p の番号を記録する (素数は大きい順に書いて、最後に書いた最小の
// 素数を残す)。m/p は区画の前か区画の中の m より前にあるので、前から埋めればそろっている。
// ZMod (法が奇数) では、2 番目以降の区画の m/p はどれも前の区画にあるので、区画の素数を集めて、32 個ずつ AVX2 の Montgomery の
// 積で冪を求める (8 個ずつの 4 本を並べて待ちを隠す)。値は答えの vector に u32 の列として直接書く。
// procon-judge の yosupo-sum-of-exponential-times-polynomial で書き比べた se1 と se2 の表を写したもの。
namespace pow_table_internal {
using u32= unsigned;
using u64= unsigned long long;
// 型が ZMod かの判定は neo/base/core/algebra/ZMod.hpp。使う場所の名前は変えない。
using ::is_zmod;
// 8 個の x (M 未満) の d 乗を 4 組 (32 個) まとめて求め、M 未満で out に書く。
template <class Z> void pow32(const u32* x, u32* out, u64 d, __m256i mv, __m256i mi, u32 R1, u32 R2) {
 auto mont8= [&](__m256i a, __m256i b) {
  const __m256i te= _mm256_mul_epu32(a, b), to= _mm256_mul_epu32(_mm256_srli_epi64(a, 32), _mm256_srli_epi64(b, 32));
  const __m256i ue= _mm256_add_epi64(te, _mm256_mul_epu32(_mm256_mul_epu32(te, mi), mv));
  const __m256i uo= _mm256_add_epi64(to, _mm256_mul_epu32(_mm256_mul_epu32(to, mi), mv));
  return _mm256_blend_epi32(_mm256_srli_epi64(ue, 32), uo, 0b10101010);
 };
 __m256i b[4], r[4];
 const __m256i r2= _mm256_set1_epi32(R2), one= _mm256_set1_epi32(R1);
 for(int q= 0; q < 4; ++q) b[q]= mont8(_mm256_loadu_si256((const __m256i*)(x + 8 * q)), r2), r[q]= one;
 int top= 63;
 while(top >= 0 && !(d >> top & 1)) --top;
 for(int i= top; i >= 0; --i) {
  for(int q= 0; q < 4; ++q) r[q]= mont8(r[q], r[q]);
  if(d >> i & 1)
   for(int q= 0; q < 4; ++q) r[q]= mont8(r[q], b[q]);
 }
 for(int q= 0; q < 4; ++q) {
  const __m256i v= mont8(r[q], _mm256_set1_epi32(1));
  _mm256_storeu_si256((__m256i*)(out + 8 * q), _mm256_min_epu32(v, _mm256_sub_epi32(v, mv)));
 }
}
template <class T> std::vector<T> generic(int n, unsigned long long d) {
 std::vector<T> t(n + 1);
 auto pw= [d](T b) {
  T r(1);
  for(auto e= d; e; e>>= 1, b*= b)
   if(e & 1) r*= b;
  return r;
 };
 if(d == 0) return std::fill(t.begin(), t.end(), T(1)), t;
 t[0]= T(0);
 if(n >= 1) t[1]= T(1);
 int sq= int(std::sqrt(double(n)));
 while((long long)sq * sq > n) --sq;
 while((long long)(sq + 1) * (sq + 1) <= n) ++sq;
 std::vector<char> comp(sq + 1);
 std::vector<int> ps;
 std::vector<T> pp;
 for(int i= 2; i <= sq; ++i)
  if(!comp[i]) {
   ps.push_back(i), pp.push_back(pw(T(i)));
   for(int j= i * i; j <= sq; j+= i) comp[j]= 1;
  }
 constexpr int B= 1 << 15;
 std::vector<unsigned> cof(B);
 std::vector<int> pi(B);
 for(int L= 2; L <= n; L+= B) {
  const int R= std::min(n + 1, L + B);
  std::fill(cof.begin(), cof.begin() + (R - L), 0u);
  for(int j= int(ps.size()) - 1; j >= 0; --j) {
   const int p= ps[j];
   const long long st= std::max((long long)p * p, (long long)(L + p - 1) / p * p);
   if(st >= R) continue;
   unsigned k= unsigned(st / p);
   for(int m= int(st); m < R; m+= p, ++k) cof[m - L]= k, pi[m - L]= j;
  }
  for(int m= L; m < R; ++m) {
   const unsigned c= cof[m - L];
   t[m]= c ? t[c] * pp[pi[m - L]] : pw(T(m));
  }
 }
 return t;
}
template <class Z> std::vector<Z> zmod(int n, u64 d) {
 std::vector<Z> res(n + 1);
 static_assert(std::is_standard_layout_v<Z> && sizeof(Z) == 4);
 u32* t= reinterpret_cast<u32*>(res.data());  // ZMod の値 (M 未満の u32) として直接書く
 const u32 M= Z::mod();
 u32 minv= M;
 for(int i= 0; i < 5; ++i) minv*= 2 - M * minv;
 const __m256i mv= _mm256_set1_epi32(M), mi= _mm256_set1_epi32(-minv);
 const u32 R1= Z(u64(1) << 32).val(), R2= (Z::raw(R1) * Z::raw(R1)).val();
 if(d == 0) return std::fill(res.begin(), res.end(), Z(1)), res;
 t[0]= 0;
 if(n >= 1) t[1]= Z(1).val();
 int sq= int(std::sqrt(double(n)));
 while((long long)sq * sq > n) --sq;
 while((long long)(sq + 1) * (sq + 1) <= n) ++sq;
 std::vector<char> comp(sq + 1);
 std::vector<int> ps;
 std::vector<u32> pp;
 for(int i= 2; i <= sq; ++i)
  if(!comp[i]) {
   ps.push_back(i), pp.push_back(Z(i).pow(d).val());
   for(int j= i * i; j <= sq; j+= i) comp[j]= 1;
  }
 constexpr int B= 1 << 15;
 std::vector<u32> cof(B), ppv(B), pr(B + 32), prv(B + 32);
 for(int L= 2; L <= n; L+= B) {
  const int R= std::min(n + 1, L + B);
  std::fill(cof.begin(), cof.begin() + (R - L), 0u);
  for(int j= int(ps.size()) - 1; j >= 0; --j) {
   const int p= ps[j];
   const long long st= std::max((long long)p * p, (long long)(L + p - 1) / p * p);
   if(st >= R) continue;
   u32 k= u32(st / p);
   const u32 v= pp[j];
   for(int m= int(st); m < R; m+= p, ++k) cof[m - L]= k, ppv[m - L]= v;
  }
  if(L == 2) {  // 最初の区画は m/p が区画の中にもあるので、素数もその場で求める
   for(int m= L; m < R; ++m) {
    const u32 c= cof[m - L];
    t[m]= c ? (Z::raw(t[c]) * Z::raw(ppv[m - L])).val() : Z(m).pow(d).val();
   }
   continue;
  }
  int np= 0;
  for(int m= L; m < R; ++m) {
   const u32 c= cof[m - L];
   if(c) t[m]= (Z::raw(t[c]) * Z::raw(ppv[m - L])).val();
   else pr[np++]= Z(m).val();
  }
  for(int i= np; i % 32; ++i) pr[i]= 1;
  for(int i= 0; i < np; i+= 32) pow32<Z>(pr.data() + i, prv.data() + i, d, mv, mi, R1, R2);
  for(int i= 0, m= L; i < np; ++m)
   if(!cof[m - L]) t[m]= prv[i++];
 }
 return res;
}
}  // namespace pow_table_internal
template <class T> std::vector<T> pow_table(int n, unsigned long long d) {
 if constexpr(pow_table_internal::is_zmod<T>::value) {
  if(T::mod() & 1) return pow_table_internal::zmod<T>(n, d);
 }
 return pow_table_internal::generic<T>(n, d);
}
