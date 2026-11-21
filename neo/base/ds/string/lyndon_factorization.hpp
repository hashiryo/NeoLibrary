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
// Lyndon 分解を Duval の方法 (1983) で求める。i から始まる部分を読み、j を進めながら周期の頭 k と比べ、s[k] = s[j] なら k も進め、
// s[k] < s[j] なら k を i に戻し、s[k] > s[j] で止まったら周期 j - k の Lyndon 語を i から切り出す。文字列の版は 3 つの手を足す。
// k が i に戻っている間は s[j] > s[i] の位置で j が進むだけなので、s[j] <= s[i] になる最初の j を 32 byte ずつ AVX2 で探して飛ぶ。
// s[k] = s[j] の続く長さは、1 文字目を比べてから 8 byte の排他的論理和の下位から 0 の byte を数えて進む。切れ目の列は n + 1 個分を
// 先に取り、同じ周期の因子が続く分は等差の列としてまとめて書く。procon-judge の yosupo-lyndon-factorization で書き比べた duval_skip
// を写したもの。
namespace lyndon_factorization_internal {
inline uint64_t load(const unsigned char* p) {
 uint64_t x;
 std::memcpy(&x, p, 8);
 return x;
}
// j 以上で t[j] <= c となる最初の j (無ければ n)。
inline int skip_greater(const unsigned char* t, int j, int n, unsigned char c) {
 const __m256i v= _mm256_set1_epi8(char(c));
 for(; j + 32 <= n; j+= 32) {
  const __m256i x= _mm256_loadu_si256(reinterpret_cast<const __m256i*>(t + j));
  if(const unsigned m= unsigned(_mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_max_epu8(x, v), v)))) return j + std::countr_zero(m);
 }
 while(j < n && t[j] > c) ++j;
 return j;
}
// 周期 p の因子を c 個、i から切り出して cut の末尾に足す。
inline void emit(std::vector<int>& cut, int i, int p, int c) {
 const std::size_t o= cut.size();
 cut.resize(o + c);
 int* w= cut.data() + o;
 for(int q= 0; q < c; ++q) w[q]= i + (q + 1) * p;
}
}
// 文字列 s の Lyndon 分解の切れ目 0 = a_0 < a_1 < ... < a_k = n を返す。s[a_{t - 1}, a_t) が t 番目の Lyndon 語で、辞書順で減って
// いかない。文字は unsigned char として比べる。O(n)。
inline std::vector<int> lyndon_factorization(const std::string& s) {
 using namespace lyndon_factorization_internal;
 const int n= s.size();
 const unsigned char* t= reinterpret_cast<const unsigned char*>(s.data());
 std::vector<int> cut;
 cut.reserve(n + 1);
 cut.push_back(0);
 for(int i= 0; i < n;) {
  int j= i + 1, k= i;
  while(j < n) {
   if(k == i) {
    j= skip_greater(t, j, n, t[i]);
    if(j == n) break;
   }
   if(t[k] != t[j]) {
    if(t[k] > t[j]) break;
    k= i, ++j;  // s[k] < s[j]: s[i, j + 1) は Lyndon 語
    continue;
   }
   ++k, ++j;
   while(j + 8 <= n) {
    if(const uint64_t x= load(t + k) ^ load(t + j)) {
     const int m= std::countr_zero(x) >> 3;
     k+= m, j+= m;
     break;
    }
    k+= 8, j+= 8;
   }
   while(j < n && t[k] == t[j]) ++k, ++j;
   if(j == n || t[k] > t[j]) break;
   k= i, ++j;
  }
  const int p= j - k, c= (k - i) / p + 1;
  emit(cut, i, p, c);
  i+= c * p;
 }
 return cut;
}
// 比べられる値の列 s の Lyndon 分解の切れ目を返す。値は < で比べる。O(n)。
template <class T> std::vector<int> lyndon_factorization(const std::vector<T>& s) {
 const int n= s.size();
 std::vector<int> cut{0};
 for(int i= 0; i < n;) {
  int j= i + 1, k= i;
  for(; j < n && !(s[j] < s[k]); ++j) k= s[k] < s[j] ? i : k + 1;
  const int p= j - k, c= (k - i) / p + 1;
  lyndon_factorization_internal::emit(cut, i, p, c);
  i+= c * p;
 }
 return cut;
}
