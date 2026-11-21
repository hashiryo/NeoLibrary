#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <bit>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
// Z 配列。z[i] は s と s[i, n) の最長共通接頭辞の長さで、z[0] = n。右端が最も右にある一致の区間 [l, r) (箱) を持ち、i が箱の中で
// z[i - l] < r - i なら z[i] = z[i - l] として比べず、そうでなければ r - i (箱の外なら 0) から伸ばして箱を [i, i + z[i]) に替える。
// 文字列の版は、i > 0 で z[i] が 1 以上になるのは s[i] = s[0] のときだけ (箱の中でも s[i] = s[i - l] なので変わらない) なのを使い、
// 0 で埋めた出力に、s[0] と同じ文字の位置だけを書く。その位置は 32 byte ずつ AVX2 で比べた印から拾い、32 位置がすべて一致したら
// 印を見ずに続けて回す。伸ばすときは次の 1 文字を比べてから、8 byte の排他的論理和の下位から 0 の byte を数える。procon-judge の
// yosupo-zalgorithm で書き比べた z_skip を写したもの。
namespace z_algorithm_internal {
inline uint64_t load(const unsigned char* p) {
 uint64_t x;
 std::memcpy(&x, p, 8);
 return x;
}
// t[0, n) の Z 配列を、0 で埋めた z[0, n) に書く。
inline void z_bytes(const unsigned char* t, int n, int* z) {
 if(n == 0) return;
 z[0]= n;
 int l= 0, r= 0;
 // s[i] = s[0] の位置 i (> 0) の z[i] を求める。
 auto at= [&](int i) __attribute__((always_inline)) {
  int k;
  if(i < r) {
   k= z[i - l];
   if(k < r - i) {
    z[i]= k;
    return;
   }
   k= r - i;
  } else k= 1;  // t[i] = t[0] は分かっている
  if(i + k == n || t[k] != t[i + k]) {  // 8 byte を読む前に、次の 1 文字を比べる
   z[i]= k, l= i, r= i + k;
   return;
  }
  for(;;) {
   if(i + k + 8 > n) {
    while(i + k < n && t[k] == t[i + k]) ++k;
    break;
   }
   if(const uint64_t x= load(t + k) ^ load(t + i + k)) {
    k+= std::countr_zero(x) >> 3;
    break;
   }
   k+= 8;
  }
  z[i]= k, l= i, r= i + k;
 };
 const __m256i c0= _mm256_set1_epi8(char(t[0]));
 int i= 1;
 for(; i + 32 <= n; i+= 32) {
  unsigned m= unsigned(_mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(t + i)), c0)));
  if(m == ~0u) {
   for(int j= i; j < i + 32; ++j) at(j);
   continue;
  }
  for(; m; m&= m - 1) at(i + std::countr_zero(m));
 }
 for(; i < n; ++i)
  if(t[i] == t[0]) at(i);
}
}
// 文字列 s の Z 配列を返す。文字は unsigned char として比べる。
inline std::vector<int> z_algorithm(const std::string& s) {
 std::vector<int> z(s.size(), 0);
 z_algorithm_internal::z_bytes(reinterpret_cast<const unsigned char*>(s.data()), int(s.size()), z.data());
 return z;
}
// 比べられる値の列 s の Z 配列を返す。値は == で比べる。
template <class T> std::vector<int> z_algorithm(const std::vector<T>& s) {
 const int n= s.size();
 std::vector<int> z(n, 0);
 if(n == 0) return z;
 z[0]= n;
 for(int i= 1, l= 0, r= 0; i < n; ++i) {
  int k= 0;
  if(i < r) {
   k= z[i - l];
   if(k < r - i) {
    z[i]= k;
    continue;
   }
   k= r - i;
  }
  while(i + k < n && s[k] == s[i + k]) ++k;
  z[i]= k, l= i, r= i + k;
 }
 return z;
}
