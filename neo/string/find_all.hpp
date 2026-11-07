#pragma once
#include <cstring>
#include <string>
#include <vector>
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
// テキストの中のパターンの出現位置をすべて、昇順に返す。パターンが 8 文字以下なら、テキストの 32 か所の位置について、パターンの全部の
// 文字を AVX2 で比べて and を取り、そのまま答えにする。9 文字以上なら、先頭の 4 文字と末尾の 1 文字で候補を絞り (W. Muła の
// SIMD-friendly algorithms for substring searching の generic な形を広げたもの)、候補の位置で残りを 32 byte ずつ比べる。候補が密な
// 入力 (同じ文字だけのテキストと長いパターンなど) では、比べた 32 byte の塊の数が n / 2 + m を超えた所から KMP に切り替えるので、
// 全体で O(n + m)。KMP の一致の判定は分岐で書く (k += (t[i] == p[k]) と書くと、次に読む p[k] の位置が前の読み込みを待つ鎖になる)。
// procon-judge の self-string-pattern-search で書き比べた simd_kmp5 を写したもの。
namespace find_all_internal {
// テキスト t の位置 i から先を KMP で走査し、出現位置を pos に足す。
inline void kmp_from(const unsigned char* T, int n, const unsigned char* P, int m, int i, std::vector<int>& pos) {
 std::vector<int> pi(m);
 pi[0]= 0;
 for(int j= 1, k= 0; j < m; ++j) {
  if(P[j] == P[k]) {
   ++k;
  } else if(k > 0) {
   do k= pi[k - 1];
   while(k > 0 && P[j] != P[k]);
   if(P[j] == P[k]) ++k;
  }
  pi[j]= k;
 }
 for(int k= 0; i < n; ++i) {
  if(T[i] == P[k]) {
   if(++k == m) pos.push_back(i - m + 1), k= pi[m - 1];
  } else if(k > 0) {
   do k= pi[k - 1];
   while(k > 0 && T[i] != P[k]);
   if(T[i] == P[k] && ++k == m) pos.push_back(i - m + 1), k= pi[m - 1];
  }
 }
}
// パターンが M 文字 (M <= 8) のとき、32 か所ずつ全部の文字を比べる。SIMD で調べ終えた次の位置を返す。
template <int M> inline int exact(const char* T, const char* P, int n, std::vector<int>& pos) {
 __m256i f[M];
 for(int j= 0; j < M; ++j) f[j]= _mm256_set1_epi8(P[j]);
 int i= 0;
 for(; i + M + 31 <= n; i+= 32) {
  __m256i acc= _mm256_cmpeq_epi8(_mm256_loadu_si256((const __m256i*)(T + i)), f[0]);
  for(int j= 1; j < M; ++j) acc= _mm256_and_si256(acc, _mm256_cmpeq_epi8(_mm256_loadu_si256((const __m256i*)(T + i + j)), f[j]));
  for(unsigned mask= (unsigned)_mm256_movemask_epi8(acc); mask; mask&= mask - 1) pos.push_back(i + __builtin_ctz(mask));
 }
 return i;
}
}
inline std::vector<int> find_all(const std::string& text, const std::string& pattern) {
 using namespace find_all_internal;
 const int n= text.size(), m= pattern.size();
 std::vector<int> pos;
 if(m > n) return pos;
 pos.reserve(n - m + 1);
 if(m == 0) {
  for(int i= 0; i <= n; ++i) pos.push_back(i);
  return pos;
 }
 const char* T= text.data();
 std::string pb= pattern;
 pb.append(64, '\0');
 const char* P= pb.data();
 int i= 0;
 switch(m) {
  case 1: i= exact<1>(T, P, n, pos); break;
  case 2: i= exact<2>(T, P, n, pos); break;
  case 3: i= exact<3>(T, P, n, pos); break;
  case 4: i= exact<4>(T, P, n, pos); break;
  case 5: i= exact<5>(T, P, n, pos); break;
  case 6: i= exact<6>(T, P, n, pos); break;
  case 7: i= exact<7>(T, P, n, pos); break;
  case 8: i= exact<8>(T, P, n, pos); break;
  default: {
   const __m256i f0= _mm256_set1_epi8(P[0]), f1= _mm256_set1_epi8(P[1]), f2= _mm256_set1_epi8(P[2]), f3= _mm256_set1_epi8(P[3]), l= _mm256_set1_epi8(P[m - 1]);
   const int body= m - 5, nc= (body + 31) / 32;
   const unsigned last= body % 32 ? (1u << (body % 32)) - 1 : 0xffffffffu;
   long long budget= n / 2 + m;
   for(; i + m + 63 <= n; i+= 32) {
    __m256i acc= _mm256_and_si256(_mm256_cmpeq_epi8(_mm256_loadu_si256((const __m256i*)(T + i)), f0), _mm256_cmpeq_epi8(_mm256_loadu_si256((const __m256i*)(T + i + m - 1)), l));
    acc= _mm256_and_si256(acc, _mm256_cmpeq_epi8(_mm256_loadu_si256((const __m256i*)(T + i + 1)), f1));
    acc= _mm256_and_si256(acc, _mm256_cmpeq_epi8(_mm256_loadu_si256((const __m256i*)(T + i + 2)), f2));
    acc= _mm256_and_si256(acc, _mm256_cmpeq_epi8(_mm256_loadu_si256((const __m256i*)(T + i + 3)), f3));
    for(unsigned mask= (unsigned)_mm256_movemask_epi8(acc); mask; mask&= mask - 1) {
     const int k= __builtin_ctz(mask);
     const char *x= T + i + k + 4, *y= P + 4;
     int c= 0;
     for(; c < nc; ++c) {
      const __m256i u= _mm256_loadu_si256((const __m256i*)(x + 32 * c)), v= _mm256_loadu_si256((const __m256i*)(y + 32 * c));
      const unsigned eq= (unsigned)_mm256_movemask_epi8(_mm256_cmpeq_epi8(u, v)) | ~(c + 1 == nc ? last : 0xffffffffu);
      if(eq != 0xffffffffu) break;
     }
     budget-= c + 1;
     if(c == nc) pos.push_back(i + k);
    }
    if(budget < 0) {
     kmp_from(reinterpret_cast<const unsigned char*>(T), n, reinterpret_cast<const unsigned char*>(P), m, i + 32, pos);
     return pos;
    }
   }
  }
 }
 for(; i + m <= n; ++i)
  if(std::memcmp(T + i, P, m) == 0) pos.push_back(i);
 return pos;
}
