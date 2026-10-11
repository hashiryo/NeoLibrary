#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iterator>
#include <string>
#include <vector>
// 列の各中心の最長の回文の長さを Manacher で求める。中心 c (0 以上 2n - 1 未満) の回文を [lo, hi) とすると lo + hi = c + 1 なので、
// 回文は右端だけで決まり、長さは 2 hi - c - 1。2n - 1 個の中心を 1 本のループで回し、右端が最も右にある回文 (箱) の中心 C と右端 H
// を持つ。c / 2 + 1 < H なら、C について対称な中心 c' = 2C - c の回文の右端 h' から、hm = c - C + h' までは回文だと分かり、hm < H
// ならそれが答え、そうでなければ H から伸ばす。文字列の版は、両隣が一致しない中心 (文字 s[i] の中心なら長さ 1、s[i] と s[i + 1] の間
// なら 0) が箱の中でも自明なのを使い、出力の偶数番目を 1、奇数番目を 0 で埋めてから、両隣が一致する中心だけを中心の順に調べる。
// 両隣の一致は 32 文字ずつ AVX2 で比べた 2 つの印を 1 bit おきに組み合わせて拾い、64 個の中心がすべて一致したら印を見ずに続けて
// 回す。伸ばすときは次の 1 文字を比べてから、右側の 8 byte と、左側の 8 byte を byte の順で反転したものの排他的論理和で一致を数える。
// procon-judge の yosupo-enumerate-palindromes で書き比べた pal_skip を写したもの。
namespace palindrome_lengths_internal {
inline uint64_t load(const unsigned char* p) {
 uint64_t x;
 std::memcpy(&x, p, 8);
 return x;
}
// 32 bit の x の各 bit を、64 bit の偶数番目の bit に広げる。
inline uint64_t spread(uint32_t x) {
 uint64_t v= x;
 v= (v | v << 16) & 0x0000FFFF0000FFFFull;
 v= (v | v << 8) & 0x00FF00FF00FF00FFull;
 v= (v | v << 4) & 0x0F0F0F0F0F0F0F0Full;
 v= (v | v << 2) & 0x3333333333333333ull;
 v= (v | v << 1) & 0x5555555555555555ull;
 return v;
}
// 1, 0, 1, 0, ... を出すイテレータ (出力の配列を、自明な答えで 1 回だけ書いて作る)。
struct Alt {
 using iterator_category= std::forward_iterator_tag;
 using value_type= int;
 using difference_type= std::ptrdiff_t;
 using pointer= const int*;
 using reference= int;
 std::ptrdiff_t c;
 int operator*() const { return int(~c & 1); }
 Alt& operator++() { return ++c, *this; }
 Alt operator++(int) {
  Alt r= *this;
  ++c;
  return r;
 }
 bool operator==(const Alt& o) const { return c == o.c; }
 bool operator!=(const Alt& o) const { return c != o.c; }
 std::ptrdiff_t operator-(const Alt& o) const { return c - o.c; }
};
// t[0, n) の各中心の長さを、自明な答えで埋めた len[0, 2n - 1) に書く。
inline void bytes(const unsigned char* t, int n, int* len) {
 int C= 0, H= 0;
 auto at= [&](int c) __attribute__((always_inline)) {
  int h= c / 2 + 1;  // 最も短い回文の右端 (c が偶数なら 1 文字、奇数なら空)
  if(h < H) {
   const int c2= 2 * C - c, hm= c - C + ((len[c2] + c2 + 1) >> 1);
   if(hm < H) {
    len[c]= 2 * hm - c - 1;
    return;
   }
   h= H;
  }
  int lo= c + 1 - h;
  if(lo > 0 && h < n && t[lo - 1] == t[h]) {
   --lo, ++h;
   for(;;) {
    if(lo < 8 || h + 8 > n) {
     while(lo > 0 && h < n && t[lo - 1] == t[h]) --lo, ++h;
     break;
    }
    if(const uint64_t x= load(t + h) ^ __builtin_bswap64(load(t + lo - 8))) {
     const int m= std::countr_zero(x) >> 3;
     lo-= m, h+= m;
     break;
    }
    lo-= 8, h+= 8;
   }
  }
  len[c]= h - lo;
  C= c, H= h;
 };
 // 文字の中心 c = 2i (1 <= i <= n - 2) は t[i - 1] = t[i + 1] のとき、間の中心 c = 2i + 1 (0 <= i <= n - 2) は t[i] = t[i + 1] のとき
 // に調べる。
 int i= 1;
 if(n >= 2 && t[0] == t[1]) at(1);
 for(; i + 33 <= n; i+= 32) {
  const __m256i a= _mm256_loadu_si256(reinterpret_cast<const __m256i*>(t + i - 1)), b= _mm256_loadu_si256(reinterpret_cast<const __m256i*>(t + i)), d= _mm256_loadu_si256(reinterpret_cast<const __m256i*>(t + i + 1));
  const uint32_t odd= uint32_t(_mm256_movemask_epi8(_mm256_cmpeq_epi8(a, d))), even= uint32_t(_mm256_movemask_epi8(_mm256_cmpeq_epi8(b, d)));
  uint64_t m= spread(odd) | spread(even) << 1;
  if(m == ~0ull) {
   for(int c= 2 * i; c < 2 * i + 64; ++c) at(c);
   continue;
  }
  for(; m; m&= m - 1) at(2 * i + std::countr_zero(m));
 }
 for(; i + 1 < n; ++i) {
  if(t[i - 1] == t[i + 1]) at(2 * i);
  if(t[i] == t[i + 1]) at(2 * i + 1);
 }
}
}
// 文字列 s の、2n - 1 個の中心ごとの最長の回文の長さを返す。c 番目は位置 c / 2 を中心とする回文で、c が偶数なら文字 s[c / 2] が中心
// の奇数の長さ、奇数なら s[c / 2] と s[c / 2 + 1] の間が中心の偶数の長さ。文字は unsigned char として比べる。
inline std::vector<int> palindrome_lengths(const std::string& s) {
 using namespace palindrome_lengths_internal;
 const int n= s.size();
 std::vector<int> len(Alt{0}, Alt{n > 0 ? 2 * n - 1 : 0});
 bytes(reinterpret_cast<const unsigned char*>(s.data()), n, len.data());
 return len;
}
// 比べられる値の列 s の、2n - 1 個の中心ごとの最長の回文の長さを返す。値は == で比べる。
template <class T> std::vector<int> palindrome_lengths(const std::vector<T>& s) {
 const int n= s.size();
 std::vector<int> len(n > 0 ? 2 * n - 1 : 0);
 for(int c= 0, C= 0, H= 0; c < 2 * n - 1; ++c) {
  int h= c / 2 + 1;
  if(h < H) {
   const int c2= 2 * C - c, hm= c - C + ((len[c2] + c2 + 1) >> 1);
   if(hm < H) {
    len[c]= 2 * hm - c - 1;
    continue;
   }
   h= H;
  }
  int lo= c + 1 - h;
  while(lo > 0 && h < n && s[lo - 1] == s[h]) --lo, ++h;
  len[c]= h - lo;
  C= c, H= h;
 }
 return len;
}
