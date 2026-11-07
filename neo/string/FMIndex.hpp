#pragma once
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include "neo/string/suffix_array.hpp"
// 文字列の中のパターンの出現回数を答える索引。接尾辞配列から、末尾に最小の番兵を足した文字列の BWT を作り、FM-index (Ferragina と
// Manzini、2000) で答える。文字は使われているものに 1 から番号を振り直し、文字ごとに BWT の位置の bit の列と、64 個ごとのその文字の
// 累積の数を 16 byte の組で持つ (rank は表を 1 回読んで popcount)。末尾の短い接尾辞は足りない文字を 0 として、3 文字の組ごとに
// それより小さい接尾辞の数 P も持つ (文字が k 種類なら (k + 1)^3 個、(k + 1)^3 <= 4n のときだけ)。count は、3 文字以下なら P だけで
// 答え、32 文字以下なら最後の 3 文字の範囲から FM-index で後ろへ絞り、33 文字以上なら最初の 3 文字の範囲の中を接尾辞配列の二分探索
// (8 byte ずつの比べ) で求める。長いパターンでは、二分探索の log n 回の比べは前から続けて読むだけなので、FM-index の 1 文字ごとの
// 飛び飛びの読み込みより速い。procon-judge の self-string-count-occurrences で書き比べた fm_sa3 を写したもの。
namespace fm_index_internal {
using u32= std::uint32_t;
using u64= std::uint64_t;
struct Ent {
 u64 bits;
 u32 cnt, pad;
};
}
class FMIndex {
 using u32= fm_index_internal::u32;
 using u64= fm_index_internal::u64;
 int n= 0, sigma= 0, K= 1, words= 0;
 unsigned char code[256]= {};  // 0 は使われていない文字
 int C[258]= {};               // C[x] は番号が x 未満の文字の数に 1 (番兵の行) を足したもの
 std::vector<int> sa, P;
 std::string t;
 std::vector<fm_index_internal::Ent> occ;
 // BWT の行 [0, i) の中の番号 c (0 始まり) の数
 [[gnu::always_inline]] int rank(int c, int i) const {
  const fm_index_internal::Ent& e= occ[std::size_t(c) * words + (i >> 6)];
  return int(e.cnt) + __builtin_popcountll(e.bits & ((u64(1) << (i & 63)) - 1));
 }
 // s[i, ...) の先頭 m 文字のうち from 文字目からと、p の from 文字目からの大小。p は末尾に 0 を 8 個足したもの。
 int cmp(int i, int from, int m, const char* p) const {
  const char* x= t.data() + i;
  for(int j= from; j < m; j+= 8) {
   u64 a, b;
   std::memcpy(&a, x + j, 8), std::memcpy(&b, p + j, 8);
   if(a != b) {
    const int d= j + (__builtin_ctzll(a ^ b) >> 3);
    if(d >= m) return 0;
    return (unsigned char)x[d] < (unsigned char)p[d] ? -1 : 1;
   }
  }
  return 0;
 }
public:
 FMIndex()= default;
 explicit FMIndex(const std::string& s): n(s.size()), sa(suffix_array(s)), t(s) {
  const unsigned char* u= reinterpret_cast<const unsigned char*>(s.data());
  int freq[256]= {};
  for(int i= 0; i < n; ++i) ++freq[u[i]];
  C[1]= 1;
  for(int c= 0; c < 256; ++c)
   if(freq[c]) code[c]= (unsigned char)++sigma, C[sigma + 1]= C[sigma] + freq[c];
  K= sigma + 1;
  if(std::size_t(K) * K * K <= 4 * std::size_t(n)) {
   P.assign(std::size_t(K) * K * K + 1, 0);
   for(int i= 0; i < n; ++i) {
    const int a= code[u[i]], b= i + 1 < n ? code[u[i + 1]] : 0, c= i + 2 < n ? code[u[i + 2]] : 0;
    ++P[(std::size_t(a) * K + b) * K + c + 1];
   }
   for(std::size_t k= 1; k < P.size(); ++k) P[k]+= P[k - 1];
  }
  // BWT の行 r は、行 0 が番兵の接尾辞で、行 r >= 1 が sa[r - 1]。
  words= (n + 1) / 64 + 1;
  occ.assign(std::size_t(sigma) * words, fm_index_internal::Ent{0, 0, 0});
  auto put= [&](int r, int i) {
   if(i > 0) occ[std::size_t(code[u[i - 1]] - 1) * words + (r >> 6)].bits|= u64(1) << (r & 63);
  };
  put(0, n);
  for(int r= 1; r <= n; ++r) put(r, sa[r - 1]);
  for(int c= 0; c < sigma; ++c) {
   fm_index_internal::Ent* e= occ.data() + std::size_t(c) * words;
   u32 a= 0;
   for(int w= 0; w < words; ++w) e[w].cnt= a, a+= u32(__builtin_popcountll(e[w].bits));
  }
  t.append(8, '\0');
 }
 int size() const { return n; }
 // s[i, i + |p|) = p となる i の数。p が空なら n + 1。
 int count(const std::string& p) const {
  const int m= p.size();
  if(m == 0) return n + 1;
  const unsigned char* q= reinterpret_cast<const unsigned char*>(p.data());
  if(P.empty() || m <= 32) {
   int lo, hi, j;
   if(!P.empty() && m >= 3) {
    const int a= code[q[m - 3]], b= code[q[m - 2]], c= code[q[m - 1]];
    if(!a || !b || !c) return 0;
    const std::size_t k= (std::size_t(a) * K + b) * K + c;
    lo= P[k] + 1, hi= P[k + 1] + 1, j= m - 4;
   } else if(!P.empty()) {
    const int a= code[q[0]], b= m > 1 ? code[q[1]] : 0;
    if(!a || (m > 1 && !b)) return 0;
    const std::size_t k= (std::size_t(a) * K + b) * K;
    return P[k + (m == 2 ? K : std::size_t(K) * K)] - P[k];
   } else {
    const int c= code[q[m - 1]];
    if(!c) return 0;
    lo= C[c], hi= C[c + 1], j= m - 2;
   }
   for(; j >= 0 && lo < hi; --j) {
    const int x= code[q[j]];
    if(!x) return 0;
    lo= C[x] + rank(x - 1, lo), hi= C[x] + rank(x - 1, hi);
   }
   return lo < hi ? hi - lo : 0;
  }
  if(m > n) return 0;
  const int a= code[q[0]], b= code[q[1]], c= code[q[2]];
  if(!a || !b || !c) return 0;
  const std::size_t k= (std::size_t(a) * K + b) * K + c;
  std::string pb= p;
  pb.append(8, '\0');
  const int *lo= sa.data() + P[k], *hi= sa.data() + P[k + 1];
  const int* x= std::partition_point(lo, hi, [&](int i) { return cmp(i, 3, m, pb.data()) < 0; });
  const int* y= std::partition_point(x, hi, [&](int i) { return cmp(i, 3, m, pb.data()) == 0; });
  return int(y - x);
 }
};
