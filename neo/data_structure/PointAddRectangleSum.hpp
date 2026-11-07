#pragma once
#include <utility>
#include <vector>
#include "neo/data_structure/WaveletMatrix.hpp"
#include "neo/sort/argsort.hpp"
// 座標が先に決まった 2 次元の点の集合について、点の重みに足すことと、長方形 [xl, xr) × [yl, yr) に入る点の重みの和を引くことを
// 交互にできる。点を x の順に並べ (neo/sort/argsort.hpp)、並べた順の y の順位の列から 1 段で 1 bit ずつ分ける 2 分の wavelet matrix
// を作り、段ごとに、分けたあとの並びの重みを BIT で持つ (元の並びの重みの BIT も 1 本持つ)。足すときは、点の位置を段ごとに rank で
// 追って各段の BIT に足す。和は、y の 2 つの境の bit が同じ段を 1 回で下り (両方の和が打ち消し合うので BIT を引かない)、違う段から
// 下は 2 本を同じループで下りて、境の bit が 1 の段だけ 0 の側へ行く範囲を BIT で引く。x の範囲と y の境の順位は、値の範囲を最大
// 65536 個の区間に分けた表と cmov の二分探索で引く。procon-judge の yosupo-point-add-rectangle-sum で書き比べた wm_bit の形。
template <class T= int, class W= long long> class PointAddRectangleSum {
 using u32= unsigned;
 using u64= unsigned long long;
 int n= 0, lg= 1, nw= 0;
 u32 sigma= 0;
 wavelet_matrix_internal::Lookup<T> xs;  // 点の x の昇順 (重複を含む) と、x 未満の点の数を引く表
 wavelet_matrix_internal::Lookup<T> ys;  // y の異なる値の昇順と、y 未満の値の数を引く表
 std::vector<int> pos;                   // 点 i の x の順の位置
 std::vector<u32> yid;                   // x の順に並べた点の y の順位
 wavelet_matrix_internal::Vec<u64> bw;   // 段 h の語 j は bw[h * nw + j]
 wavelet_matrix_internal::Vec<u32> bc;   // 語より前の 1 の数
 std::vector<int> zero;                  // 段 h の 0 の数
 wavelet_matrix_internal::Vec<W> bit;    // 段 h の BIT は bit[h (n + 1) + 1 ..]。段 lg は元の並び
 [[gnu::always_inline]] int rank1(int h, int i) const {
  const std::size_t o= std::size_t(h) * nw + (i >> 6);
  return bc[o] + __builtin_popcountll(bw[o] & ((u64(1) << (i & 63)) - 1));
 }
 // 段 h の並びの先頭 i 個の重みの和。
 [[gnu::always_inline]] W pre(int h, int i) const {
  const W* b= bit.data() + std::size_t(h) * (n + 1);
  W s= W();
  for(; i > 0; i&= i - 1) s+= b[i];
  return s;
 }
 [[gnu::always_inline]] void upd(int h, int i, W v) {
  W* b= bit.data() + std::size_t(h) * (n + 1);
  for(++i; i <= n; i+= i & -i) b[i]+= v;
 }
 void make(int h, const std::vector<W>& w) {
  W* b= bit.data() + std::size_t(h) * (n + 1);
  for(int i= 1; i <= n; ++i) b[i]= w[i - 1];
  for(int i= 1; i <= n; ++i)
   if(const int j= i + (i & -i); j <= n) b[j]+= b[i];
 }
 // x の順の位置 [l, r) のうち y の順位が x 未満の点の重みの和。
 [[gnu::always_inline]] W sum_rank(int l, int r, u32 x) const {
  if(x >= sigma) return pre(lg, r) - pre(lg, l);
  W s= W();
  for(int h= lg - 1; h >= 0; --h) {
   const int l1= rank1(h, l), r1= rank1(h, r), l0= l - l1, r0= r - r1;
   if(x >> h & 1) {
    s+= pre(h, r0) - pre(h, l0);
    l= zero[h] + l1, r= zero[h] + r1;
   } else l= l0, r= r0;
  }
  return s;
 }
 // x の順の位置 [l, r) のうち y の順位が [a, b) に入る点の重みの和 (a < b)。
 [[gnu::always_inline]] W sum_range_rank(int l, int r, u32 a, u32 b) const {
  if(b >= sigma) return pre(lg, r) - pre(lg, l) - sum_rank(l, r, a);
  for(int h= lg - 1; h >= 0; --h) {
   const int l1= rank1(h, l), r1= rank1(h, r);
   if((a ^ b) >> h & 1) {
    const int l0= l - l1, r0= r - r1;
    W s= pre(h, r0) - pre(h, l0);
    int lb= zero[h] + l1, rb= zero[h] + r1, la= l0, ra= r0;
    for(int g= h - 1; g >= 0; --g) {
     const int p1= rank1(g, lb), q1= rank1(g, rb), p0= lb - p1, q0= rb - q1;
     const int u1= rank1(g, la), v1= rank1(g, ra), u0= la - u1, v0= ra - v1;
     if(b >> g & 1) s+= pre(g, q0) - pre(g, p0), lb= zero[g] + p1, rb= zero[g] + q1;
     else lb= p0, rb= q0;
     if(a >> g & 1) s-= pre(g, v0) - pre(g, u0), la= zero[g] + u1, ra= zero[g] + v1;
     else la= u0, ra= v0;
    }
    return s;
   }
   if(b >> h & 1) l= zero[h] + l1, r= zero[h] + r1;
   else l-= l1, r-= r1;
  }
  return W();
 }
public:
 PointAddRectangleSum(): PointAddRectangleSum(std::vector<T>(), std::vector<T>(), std::vector<W>()) {}
 // 重み w[i] の点 (x[i], y[i]) の集合から作る。あとで重みを足す点も、ここで (重み 0 などで) 入れておく。
 PointAddRectangleSum(const std::vector<T>& x, const std::vector<T>& y, const std::vector<W>& w): n(x.size()) {
  const std::vector<int> o= argsort(x);
  std::vector<T> yo(n);
  std::vector<W> cw(n);
  xs.v.resize(n), pos.resize(n);
  for(int j= 0; j < n; ++j) xs.v[j]= x[o[j]], yo[j]= y[o[j]], cw[j]= w[o[j]], pos[o[j]]= j;
  xs.build();
  ys.v= compress(yo), ys.build();
  sigma= ys.v.size();
  yid.resize(n);
  for(int j= 0; j < n; ++j) yid[j]= u32(std::make_unsigned_t<T>(yo[j]));
  while((u32(1) << lg) < sigma) ++lg;
  nw= n / 64 + 1;
  bw.assign(std::size_t(nw) * lg, 0), bc.assign(std::size_t(nw) * lg, 0), zero.assign(lg, 0);
  bit.assign(std::size_t(lg + 1) * (n + 1), W());
  make(lg, cw);
  std::vector<u32> cur(yid), nxt(n);
  std::vector<W> nx(n);
  for(int h= lg - 1; h >= 0; --h) {
   u64* B= bw.data() + std::size_t(h) * nw;
   u32* C= bc.data() + std::size_t(h) * nw;
   for(int j= 0; j * 64 < n; ++j) {
    const int m= std::min(64, n - j * 64);
    u64 v= 0;
    for(int t= 0; t < m; ++t) v|= u64(cur[j * 64 + t] >> h & 1) << t;
    B[j]= v;
   }
   u32 ones= 0;
   for(int j= 0; j < nw; ++j) C[j]= ones, ones+= __builtin_popcountll(B[j]);
   const int z= n - ones;
   zero[h]= z;
   int i0= 0, i1= z;
   for(int i= 0; i < n; ++i) {
    const u32 bt= cur[i] >> h & 1;
    const int j= bt ? i1 : i0;
    nxt[j]= cur[i], nx[j]= cw[i];
    i1+= bt, i0+= bt ^ 1;
   }
   cur.swap(nxt), cw.swap(nx);
   make(h, cw);
  }
 }
 int size() const { return n; }
 // 点 i (作るときの添字) の重みに v を足す。
 void add(int i, W v) {
  i= pos[i];
  upd(lg, i, v);
  const u32 t= yid[i];
  for(int h= lg - 1; h >= 0; --h) {
   const int r1= rank1(h, i);
   i= t >> h & 1 ? zero[h] + r1 : i - r1;
   upd(h, i, v);
  }
 }
 // x が [xl, xr)、y が [yl, yr) に入る点の重みの和。xl >= xr や yl >= yr なら 0。
 [[gnu::always_inline]] W sum(T xl, T xr, T yl, T yr) const {
  if(!(xl < xr) || !(yl < yr)) return W();
  const u32 a= ys.template rank<false>(yl), b= ys.template rank<false>(yr);
  return a < b ? sum_range_rank(xs.template rank<false>(xl), xs.template rank<false>(xr), a, b) : W();
 }
};
