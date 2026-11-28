#pragma once
#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <utility>
#include <vector>
// 区間クエリをオフラインで処理する Mo's algorithm。クエリ [l, r) を集め、端点の移動が少なくなる順に並べて
// add と erase を呼ぶ。並べ方は、l のブロックと r の向きを交互にするもの (run) と、ヒルベルト曲線の順 (run_hilbert)。
//   Mo mo;
//   mo.query(l, r);                    // [l, r)
//   mo.run(add_left, add_right, erase_left, erase_right, out);  // 左右で add と erase が違うとき
//   mo.run(add, erase, out);           // 左右が同じとき
//   mo.run_hilbert(add, erase, out);   // ヒルベルト順
// add_left(i) は区間が左に伸びて i を含むとき、erase_right(i) は右が縮んで i を外すとき。out(id) は id 番目のクエリ。
// 計算量は O((N + Q) sqrt(Q) T) で、T は add か erase 1 回の計算量。ヒルベルト順も同じ程度で、定数が小さいことが多い。
namespace mo_internal {
// ヒルベルト曲線の順序。x と y は [0, 2^lg) の整数で、lg は座標の bit 数。曲線に沿って並べたときの位置を返す。
// 標準の xy2d (Hilbert curve) をそのまま使う。s は今見ている bit、n は座標の上限 2^lg。
[[gnu::always_inline]] inline unsigned long long hilbert_order(unsigned x, unsigned y, int lg) {
 if(lg <= 0) return 0;
 const unsigned n= 1u << (lg - 1);
 unsigned long long d= 0;
 for(unsigned s= n; s; s>>= 1) {
  unsigned rx= (x & s) != 0, ry= (y & s) != 0;
  d+= (unsigned long long)s * s * ((3 * rx) ^ ry);
  if(ry == 0) {
   if(rx) x= (n << 1) - 1 - x, y= (n << 1) - 1 - y;
   std::swap(x, y);
  }
 }
 return d;
}
}  // namespace mo_internal

struct Mo {
 std::vector<int> ql_, qr_;
 int n_= 0;
 // 半開区間 [l, r) のクエリを足す。n_ は右端の最大 (数列の長さ) になる。
 void query(int l, int r) {
  ql_.push_back(l), qr_.push_back(r);
  if(n_ < r) n_= r;
 }
 int size() const { return (int)ql_.size(); }
 bool empty() const { return ql_.empty(); }
 // l のブロックを左から、同じブロックの中は r を昇順と降順の交互に並べる (奇偶ソート)。ブロックの幅は n / sqrt(q)。
 template <class AL, class AR, class EL, class ER, class Out>
 void run(const AL& add_left, const AR& add_right, const EL& erase_left, const ER& erase_right, const Out& out) const {
  const int q= (int)ql_.size();
  if(q == 0) return;
  int bs= n_ / std::max(1, (int)std::sqrt((double)q));
  if(bs < 1) bs= 1;
  std::vector<int> ord(q);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b) {
   const int ab= ql_[a] / bs, bb= ql_[b] / bs;
   if(ab != bb) return ab < bb;
   return (ab & 1) ? qr_[a] > qr_[b] : qr_[a] < qr_[b];
  });
  int l= 0, r= 0;
  for(int id: ord) {
   const int ql= ql_[id], qr= qr_[id];
   while(l > ql) add_left(--l);
   while(r < qr) add_right(r++);
   while(l < ql) erase_left(l++);
   while(r > qr) erase_right(--r);
   out(id);
  }
 }
 template <class A, class E, class O>
 void run(const A& add, const E& erase, const O& out) const {
  run(add, add, erase, erase, out);
 }
 // ヒルベルト曲線の順。クエリ (l, r) を 2 次元の点と見て、曲線に沿って並べる。lg は max(r) の bit 数。
 template <class AL, class AR, class EL, class ER, class Out>
 void run_hilbert(const AL& add_left, const AR& add_right, const EL& erase_left, const ER& erase_right, const Out& out) const {
  const int q= (int)ql_.size();
  if(q == 0) return;
  const int lg= (int)std::bit_width((unsigned)n_);
  std::vector<unsigned long long> key(q);
  for(int i= 0; i < q; ++i) key[i]= mo_internal::hilbert_order((unsigned)ql_[i], (unsigned)qr_[i], lg);
  std::vector<int> ord(q);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b) { return key[a] < key[b]; });
  int l= 0, r= 0;
  for(int id: ord) {
   const int ql= ql_[id], qr= qr_[id];
   while(l > ql) add_left(--l);
   while(r < qr) add_right(r++);
   while(l < ql) erase_left(l++);
   while(r > qr) erase_right(--r);
   out(id);
  }
 }
 template <class A, class E, class O>
 void run_hilbert(const A& add, const E& erase, const O& out) const {
  run_hilbert(add, add, erase, erase, out);
 }
};
