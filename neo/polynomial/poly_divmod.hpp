#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>
#include "neo/polynomial/convolve_naive.hpp"
// 多項式の商と余り、gcd (素朴な O(nm))。K は体で、ZMod なら法は素数であること。多項式は係数の vector (k 次の係数が [k]) で、
// 末尾の 0 は取り除いて返す (0 は空の vector)。
// ZMod は、割る式を monic にしてから、商の係数ごとに a[i - db + j] += (-c) b[j] を、掛ける数が決まっているので Shoup の方法で
// 8 個ずつ足す。
namespace polydiv_internal {
using u32= unsigned;
using u64= unsigned long long;
template <class K> void trim(std::vector<K>& a) {
 while(!a.empty() && a.back() == K()) a.pop_back();
}
// d[j] += x b[j] (j < n、値は M 未満、x < M)。
inline void add_mul(u32* d, const u32* b, int n, u32 x, u32 M) {
 const u32 q= u32((u64(x) << 32) / M);
 const __m256i vx= _mm256_set1_epi32(x), vq= _mm256_set1_epi32(q), p= _mm256_set1_epi32(M);
 int j= 0;
 for(; j + 8 <= n; j+= 8) {
  const __m256i bb= _mm256_loadu_si256((const __m256i*)(b + j));
  const __m256i ev= _mm256_srli_epi64(_mm256_mul_epu32(bb, vq), 32), od= _mm256_mul_epu32(_mm256_srli_epi64(bb, 32), vq);
  __m256i r= _mm256_sub_epi32(_mm256_mullo_epi32(bb, vx), _mm256_mullo_epi32(_mm256_blend_epi32(ev, od, 0b10101010), p));
  r= _mm256_min_epu32(r, _mm256_sub_epi32(r, p));
  r= _mm256_add_epi32(r, _mm256_loadu_si256((const __m256i*)(d + j)));
  _mm256_storeu_si256((__m256i*)(d + j), _mm256_min_epu32(r, _mm256_sub_epi32(r, p)));
 }
 for(; j < n; ++j) {
  u32 t= b[j] * x - u32((u64(b[j]) * q) >> 32) * M;
  t= t >= M ? t - M : t;
  t+= d[j];
  d[j]= t >= M ? t - M : t;
 }
}
// a を monic な b (次数 db ≥ 1) で割る。a[0..db) に余りが残り、q が空でなければ商の係数を置く (q[i] は x^i の係数)。
template <class Z> void divmod_monic(std::vector<u32>& a, const std::vector<u32>& b, std::vector<u32>* q) {
 const u32 M= Z::mod();
 const int db= int(b.size()) - 1;
 if(q) q->assign(std::max(0, int(a.size()) - db), 0);
 for(int i= int(a.size()) - 1; i >= db; --i) {
  const u32 c= a[i];
  if(q) (*q)[i - db]= c;
  if(!c) continue;
  add_mul(a.data() + i - db, b.data(), db, M - c, M);
  a[i]= 0;
 }
 a.resize(std::min<std::size_t>(a.size(), db));
 trim(a);
}
}  // namespace polydiv_internal
// a = q b + r (deg r < deg b) の (q, r)。b は 0 でないこと。
template <class K> std::pair<std::vector<K>, std::vector<K>> poly_divmod(std::vector<K> a, std::vector<K> b) {
 using namespace polydiv_internal;
 trim(a), trim(b);
 assert(!b.empty());
 if(a.size() < b.size()) return {std::vector<K>{}, a};
 const K il= K(1) / b.back();
 for(auto& x: b) x*= il;
 if constexpr(poly_internal::is_zmod_v<K>) {
  if(b.size() == 1) {
   for(auto& x: a) x*= il;
   return {a, std::vector<K>{}};
  }
  auto x= poly_internal::to_u32(a);
  const auto y= poly_internal::to_u32(b);
  std::vector<u32> q;
  divmod_monic<K>(x, y, &q);
  auto Q= poly_internal::from_u32<K>(q);
  for(auto& v: Q) v*= il;
  return {Q, poly_internal::from_u32<K>(x)};
 } else {
  const int db= int(b.size()) - 1;
  std::vector<K> q(a.size() - db, K());
  for(int i= int(a.size()) - 1; i >= db; --i) {
   const K c= a[i];
   q[i - db]= c * il;
   if(c == K()) continue;
   for(int j= 0; j < db; ++j) a[i - db + j]-= c * b[j];
   a[i]= K();
  }
  a.resize(db), trim(a);
  return {q, a};
 }
}
// monic な gcd (両方 0 なら空)。
template <class K> std::vector<K> poly_gcd(std::vector<K> a, std::vector<K> b) {
 using namespace polydiv_internal;
 trim(a), trim(b);
 auto monic= [](std::vector<K>& v) {
  const K il= K(1) / v.back();
  for(auto& x: v) x*= il;
 };
 if constexpr(poly_internal::is_zmod_v<K>) {
  while(!b.empty()) {
   monic(b);
   if(b.size() == 1) return b;
   auto x= poly_internal::to_u32(a);
   divmod_monic<K>(x, poly_internal::to_u32(b), nullptr);
   a= poly_internal::from_u32<K>(x);
   std::swap(a, b);
  }
 } else {
  while(!b.empty()) a= poly_divmod(a, b).second, std::swap(a, b);
 }
 if(!a.empty()) monic(a);
 return a;
}
