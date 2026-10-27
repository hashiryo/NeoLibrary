#pragma once
#include "neo/data_structure/WaveletMatrix.hpp"
// 静的な整数の列の、重み付きの wavelet matrix。区間で値が x 未満のものの個数と重みの和、[lo, hi) の個数と和、
// 小さい方 (大きい方) から k 個の重みの和、k 番目を引く。重みを渡さなければ値そのものを重みにする。
// 1 段で 1 bit ずつ分ける 2 分の形で、段ごとに、分けたあとの並びの重みの累積を持つ。和は、値の上限の bit が 1 の段
// だけ累積を読む (分岐で飛ばす)。4 分の形より、和の配列の読み出しが少ない。bit 列は 64 bit の語の列と語ごとの累積の
// 列。値は WaveletMatrix と同じく順位に詰め、値から順位への変換も同じ表で引く。2 MB 以上の配列は、Linux では
// 2 MB 境界の mmap に置いて MADV_HUGEPAGE を頼む。procon-judge の yosupo-static-range-sum-with-upper-bound などで
// 書き比べた wm_sep_brtbl と、atcoder-abc281-e の wm_sep を写したもの。
template <class T= int, class W= long long> class WaveletMatrixSum {
 using u32= unsigned;
 using u64= unsigned long long;
 int n= 0, lg= 1, nw= 0;
 u32 sigma= 0;
 wavelet_matrix_internal::Values<T> val;
 wavelet_matrix_internal::Vec<u64> bw;  // 段 h の語 j は bw[h * nw + j]
 wavelet_matrix_internal::Vec<u32> bc;  // 語より前の 1 の数
 std::vector<int> zero;                 // 段 h の 0 の数
 wavelet_matrix_internal::Vec<W> sums;  // 段 h で分けたあとの並びの累積は sums[h (n + 1) + i]、元の並びの累積は sums[lg (n + 1) + i]
 [[gnu::always_inline]] int rank1(int h, int i) const {
  const std::size_t o= std::size_t(h) * nw + (i >> 6);
  return bc[o] + __builtin_popcountll(bw[o] & ((u64(1) << (i & 63)) - 1));
 }
 [[gnu::always_inline]] const W* S(int h) const { return sums.data() + std::size_t(h) * (n + 1); }
 void build(const std::vector<T>& a, const std::vector<W>& w) {
  std::vector<u32> cur= val.build(a), nxt(n);
  sigma= val.vals.size();
  while((u32(1) << lg) < sigma) ++lg;
  nw= n / 64 + 1;
  bw.assign(std::size_t(nw) * lg, 0), bc.assign(std::size_t(nw) * lg, 0), zero.assign(lg, 0);
  sums.assign(std::size_t(lg + 1) * (n + 1), W());
  std::vector<W> cw(w), nx(n);
  W* top= sums.data() + std::size_t(lg) * (n + 1);
  for(int i= 0; i < n; ++i) top[i + 1]= top[i] + cw[i];
  for(int h= lg - 1; h >= 0; --h) {
   u64* B= bw.data() + std::size_t(h) * nw;
   u32* C= bc.data() + std::size_t(h) * nw;
   for(int j= 0; j * 64 < n; ++j) {
    const int m= std::min(64, n - j * 64);
    u64 x= 0;
    for(int t= 0; t < m; ++t) x|= u64(cur[j * 64 + t] >> h & 1) << t;
    B[j]= x;
   }
   u32 ones= 0;
   for(int j= 0; j < nw; ++j) C[j]= ones, ones+= __builtin_popcountll(B[j]);
   const int z= n - ones;
   zero[h]= z;
   int i0= 0, i1= z;
   for(int i= 0; i < n; ++i) {
    const u32 bit= cur[i] >> h & 1;
    const int j= bit ? i1 : i0;
    nxt[j]= cur[i], nx[j]= cw[i];
    i1+= bit, i0+= bit ^ 1;
   }
   cur.swap(nxt), cw.swap(nx);
   W* Sh= sums.data() + std::size_t(h) * (n + 1);
   for(int i= 0; i < n; ++i) Sh[i + 1]= Sh[i] + cw[i];
  }
 }
 // [l, r) の順位 x 未満の個数と和。x は 0 以上 sigma 以下。
 [[gnu::always_inline]] std::pair<int, W> count_sum_rank(int l, int r, u32 x) const {
  if(x >= sigma) return {r - l, S(lg)[r] - S(lg)[l]};
  int c= 0;
  W s= W();
  for(int h= lg - 1; h >= 0; --h) {
   const int l1= rank1(h, l), r1= rank1(h, r), l0= l - l1, r0= r - r1;
   if(x >> h & 1) {
    c+= r0 - l0, s+= S(h)[r0] - S(h)[l0];
    l= zero[h] + l1, r= zero[h] + r1;
   } else l= l0, r= r0;
  }
  return {c, s};
 }
public:
 WaveletMatrixSum(): WaveletMatrixSum(std::vector<T>()) {}
 // 値そのものを重みにする。
 explicit WaveletMatrixSum(const std::vector<T>& a): n(a.size()) { build(a, std::vector<W>(a.begin(), a.end())); }
 // a[i] の重みを w[i] にする。
 WaveletMatrixSum(const std::vector<T>& a, const std::vector<W>& w): n(a.size()) { build(a, w); }
 int size() const { return n; }
 // [l, r) の重みの和。
 [[gnu::always_inline]] W sum(int l, int r) const { return S(lg)[r] - S(lg)[l]; }
 // [l, r) のうち値が x 未満のものの個数と重みの和。
 [[gnu::always_inline]] std::pair<int, W> count_sum_lt(int l, int r, T x) const { return count_sum_rank(l, r, val.template rank<false>(x)); }
 // [l, r) のうち値が x 以下のものの個数と重みの和。
 [[gnu::always_inline]] std::pair<int, W> count_sum_le(int l, int r, T x) const { return count_sum_rank(l, r, val.template rank<true>(x)); }
 // [l, r) のうち値が [lo, hi) に入るものの個数と重みの和。lo >= hi なら {0, 0}。
 std::pair<int, W> count_sum(int l, int r, T lo, T hi) const {
  if(!(lo < hi)) return {0, W()};
  const auto [c1, s1]= count_sum_lt(l, r, hi);
  const auto [c0, s0]= count_sum_lt(l, r, lo);
  return {c1 - c0, s1 - s0};
 }
 // [l, r) の値の小さい方から k 個の重みの和 (同じ値は位置の順)。0 <= k <= r - l。
 [[gnu::always_inline]] W sum_smallest(int l, int r, int k) const {
  W s= W();
  for(int h= lg - 1; h >= 0; --h) {
   const int l1= rank1(h, l), r1= rank1(h, r), l0= l - l1, r0= r - r1, z= r0 - l0;
   if(k >= z) {
    s+= S(h)[r0] - S(h)[l0], k-= z;
    l= zero[h] + l1, r= zero[h] + r1;
   } else l= l0, r= r0;
  }
  return s + S(0)[l + k] - S(0)[l];
 }
 // [l, r) の値の大きい方から k 個の重みの和。0 <= k <= r - l。
 [[gnu::always_inline]] W sum_largest(int l, int r, int k) const { return sum(l, r) - sum_smallest(l, r, r - l - k); }
 // [l, r) の小さい方から k 番目 (0 始まり)。0 <= k < r - l。
 [[gnu::always_inline]] T kth_smallest(int l, int r, int k) const {
  u32 res= 0;
  for(int h= lg - 1; h >= 0; --h) {
   const int l1= rank1(h, l), r1= rank1(h, r), z= (r - l) - (r1 - l1);
   if(k >= z) k-= z, res|= u32(1) << h, l= zero[h] + l1, r= zero[h] + r1;
   else l-= l1, r-= r1;
  }
  return val.vals[res];
 }
 // [l, r) の大きい方から k 番目 (0 始まり)。0 <= k < r - l。
 [[gnu::always_inline]] T kth_largest(int l, int r, int k) const { return kth_smallest(l, r, r - l - 1 - k); }
};
