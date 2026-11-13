#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <cassert>
#include <limits>
#include <type_traits>
#include <vector>
// 割り当て問題。n 行 m 列の行列 a (行優先) で、各行に違う列を 1 つずつ割り当てる (n > m なら各列に違う行を 1 つずつ) とき、
// 割り当てた成分の和を最小 (maximize なら最大) にする。返すのは和、行ごとの列 (割り当てられない行は -1)、行と列の双対
// u, v。最小化なら u_i + v_j ≤ a_ij で、割り当てた組は等号になる。最大化なら不等号の向きが逆になる。n < m なら割り当てない
// 列の v は 0、n > m なら割り当てない行の u は 0。
// 中身は Hungarian 法で、行を 1 つずつ足して最短の増加路を Dijkstra で探す。列の側の走査は 1 回で、整数の費用は AVX2 で
// 4 列ずつ回す。浮動小数点数の費用は同じ形をスカラーで回す。正方行列のときだけ、列のポテンシャルを列ごとの最小で始める
// (n < m で使うと、割り当てない列のポテンシャルがそろわず正しくない)。縮約費用が同じ列のうちでは、空いている列を選ぶ。
// 行が列より多ければ転置して解く。整数の費用は、絶対値の最大を A として (n + m) A < 2^56 であること (走査の鍵に、
// 2 倍した縮約費用とポテンシャルの和を i64 で持つため)。
// 方式は procon-judge の yosupo-assignment で比べて選んだ (hungarian_avx2_cr)。記録は algo-notes の notes/network-flow.md。
template <class Cost> struct Assignment {
 Cost cost;
 std::vector<int> col;    // 行 i の列。割り当てられない行は -1
 std::vector<Cost> u, v;  // 行と列の双対
};
namespace assignment_internal {
using i64= long long;
// n ≤ m。a は n 行 m 列 (行優先)。p[j] は列 j (1 から m) の行 (1 から n、空なら 0)。u は 1 から n、v は 1 から m。
inline void solve_i64(int n, int m, const i64* a, std::vector<int>& p, std::vector<i64>& u, std::vector<i64>& v) {
 constexpr i64 OFF= i64(1) << 60, USED= OFF / 2, INF= OFF / 2;
 std::vector<i64> veff(m + 1), minv(m + 1), tj(m + 1), way(m + 1);
 std::vector<int> usedcols;
 usedcols.reserve(m + 1);
 p.assign(m + 1, 0), u.assign(n + 1, 0), v.assign(m + 1, 0);
 if(n == m) {
  for(int j= 1; j <= m; ++j) v[j]= a[j - 1];
  for(int r= 1; r < n; ++r)
   for(int j= 1; j <= m; ++j) v[j]= std::min(v[j], a[(size_t)r * m + j - 1]);
 }
 for(int i= 1; i <= n; ++i) {
  p[0]= i;
  int j0= 0;
  i64 T= 0;
  // 鍵は「縮約費用 × 2 + 割り当て済みなら 1」。割り当て済みかどうかは 1 行を足す間は変わらないので veff に入れておく。
  for(int j= 0; j <= m; ++j) veff[j]= 2 * v[j] - (p[j] != 0), minv[j]= INF;
  usedcols.clear();
  do {
   usedcols.push_back(j0);
   tj[j0]= T, veff[j0]= -OFF, minv[j0]= USED;
   const int i0= p[j0];
   const i64* row= a + (size_t)(i0 - 1) * m - 1;
   const i64 c= 2 * (T - u[i0]);
   const __m256i vc= _mm256_set1_epi64x(c), vj0= _mm256_set1_epi64x(j0), four= _mm256_set1_epi64x(4);
   __m256i vbest= _mm256_set1_epi64x(INF), vbidx= _mm256_setzero_si256(), vidx= _mm256_setr_epi64x(1, 2, 3, 4);
   int j= 1;
   for(; j + 3 <= m; j+= 4) {
    const __m256i r= _mm256_loadu_si256((const __m256i*)(row + j));
    const __m256i ve= _mm256_loadu_si256((const __m256i*)(veff.data() + j));
    __m256i mv= _mm256_loadu_si256((const __m256i*)(minv.data() + j));
    const __m256i cur= _mm256_add_epi64(_mm256_sub_epi64(_mm256_add_epi64(r, r), ve), vc);
    const __m256i lt= _mm256_cmpgt_epi64(mv, cur);
    mv= _mm256_blendv_epi8(mv, cur, lt);
    _mm256_storeu_si256((__m256i*)(minv.data() + j), mv);
    const __m256i w= _mm256_loadu_si256((const __m256i*)(way.data() + j));
    _mm256_storeu_si256((__m256i*)(way.data() + j), _mm256_blendv_epi8(w, vj0, lt));
    const __m256i lb= _mm256_cmpgt_epi64(vbest, mv);
    vbest= _mm256_blendv_epi8(vbest, mv, lb);
    vbidx= _mm256_blendv_epi8(vbidx, vidx, lb);
    vidx= _mm256_add_epi64(vidx, four);
   }
   alignas(32) i64 bb[4], bi[4];
   _mm256_store_si256((__m256i*)bb, vbest), _mm256_store_si256((__m256i*)bi, vbidx);
   i64 best= bb[0];
   int j1= (int)bi[0];
   for(int k= 1; k < 4; ++k)
    if(bb[k] < best) best= bb[k], j1= (int)bi[k];
   for(; j <= m; ++j) {
    const i64 cur= 2 * row[j] - veff[j] + c;
    if(cur < minv[j]) minv[j]= cur, way[j]= j0;
    if(minv[j] < best) best= minv[j], j1= j;
   }
   T= best >> 1;
   j0= j1;
  } while(p[j0] != 0);
  for(int jj: usedcols) {
   const i64 d= T - tj[jj];
   u[p[jj]]+= d, v[jj]-= d;
  }
  do {
   const int j1= (int)way[j0];
   p[j0]= p[j1], j0= j1;
  } while(j0);
 }
}
// 浮動小数点数の費用。solve_i64 と同じ形をスカラーで回す。縮約費用が同じ列のうちでは空いている列を選ぶ。
template <class T> void solve_float(int n, int m, const T* a, std::vector<int>& p, std::vector<T>& u, std::vector<T>& v) {
 const T INF= std::numeric_limits<T>::infinity();
 std::vector<T> minv(m + 1), tj(m + 1);
 std::vector<int> way(m + 1), usedcols;
 std::vector<char> used(m + 1);
 usedcols.reserve(m + 1);
 p.assign(m + 1, 0), u.assign(n + 1, 0), v.assign(m + 1, 0);
 if(n == m) {
  for(int j= 1; j <= m; ++j) v[j]= a[j - 1];
  for(int r= 1; r < n; ++r)
   for(int j= 1; j <= m; ++j) v[j]= std::min(v[j], a[(size_t)r * m + j - 1]);
 }
 for(int i= 1; i <= n; ++i) {
  p[0]= i;
  int j0= 0;
  T t= 0;
  std::fill(minv.begin(), minv.end(), INF), std::fill(used.begin(), used.end(), 0);
  usedcols.clear();
  do {
   usedcols.push_back(j0);
   tj[j0]= t, used[j0]= 1;
   const int i0= p[j0];
   const T* row= a + (size_t)(i0 - 1) * m - 1;
   const T c= t - u[i0];
   T best= INF;
   int j1= -1;
   bool bfree= false;
   for(int j= 1; j <= m; ++j) {
    if(used[j]) continue;
    const T cur= row[j] - v[j] + c;
    if(cur < minv[j]) minv[j]= cur, way[j]= j0;
    const bool fr= p[j] == 0;
    if(minv[j] < best || (minv[j] == best && fr && !bfree)) best= minv[j], j1= j, bfree= fr;
   }
   t= best;
   j0= j1;
  } while(p[j0] != 0);
  for(int jj: usedcols) {
   const T d= t - tj[jj];
   u[p[jj]]+= d, v[jj]-= d;
  }
  do {
   const int j1= way[j0];
   p[j0]= p[j1], j0= j1;
  } while(j0);
 }
}
}
// n 行 m 列の行列 a (行優先、大きさ n * m) の割り当て問題を解く。
template <class Cost> Assignment<Cost> assignment(int n, int m, const std::vector<Cost>& a, bool maximize= false) {
 static_assert(std::is_arithmetic_v<Cost>);
 using namespace assignment_internal;
 assert((long long)a.size() == (long long)n * m);
 const bool tr= n > m;
 const int N= tr ? m : n, M= tr ? n : m;
 Assignment<Cost> ret{0, std::vector<int>(n, -1), std::vector<Cost>(n), std::vector<Cost>(m)};
 std::vector<int> p;
 auto at= [&](int i, int j) { return tr ? a[(size_t)j * m + i] : a[(size_t)i * m + j]; };  // 解く向きの (i, j)
 auto finish= [&](const auto& u, const auto& v) {
  // p[j] は解く向きの列 j の行 (1 始まり)。元の向きに戻し、最大化なら双対の符号を戻す。
  const Cost sg= maximize ? -1 : 1;
  for(int j= 1; j <= M; ++j)
   if(p[j]) {
    const int i= p[j] - 1, jj= j - 1;
    if(tr) ret.col[jj]= i;
    else ret.col[i]= jj;
   }
  for(int i= 0; i < N; ++i) (tr ? ret.v : ret.u)[i]= sg * Cost(u[i + 1]);
  for(int j= 0; j < M; ++j) (tr ? ret.u : ret.v)[j]= sg * Cost(v[j + 1]);
  for(int i= 0; i < n; ++i)
   if(ret.col[i] >= 0) ret.cost+= a[(size_t)i * m + ret.col[i]];
 };
 if constexpr(std::is_integral_v<Cost>) {
  std::vector<i64> u, v, b;
  auto copy= [&] {
   b.resize((size_t)N * M);
   for(int i= 0; i < N; ++i)
    for(int j= 0; j < M; ++j) b[(size_t)i * M + j]= maximize ? -i64(at(i, j)) : i64(at(i, j));
   solve_i64(N, M, b.data(), p, u, v);
  };
  if constexpr(std::is_same_v<Cost, i64>) {
   if(!tr && !maximize) solve_i64(N, M, a.data(), p, u, v);
   else copy();
  } else copy();
  finish(u, v);
 } else {
  std::vector<Cost> u, v;
  if(!tr && !maximize) solve_float(N, M, a.data(), p, u, v);
  else {
   std::vector<Cost> b((size_t)N * M);
   for(int i= 0; i < N; ++i)
    for(int j= 0; j < M; ++j) b[(size_t)i * M + j]= maximize ? -at(i, j) : at(i, j);
   solve_float(N, M, b.data(), p, u, v);
  }
  finish(u, v);
 }
 return ret;
}
