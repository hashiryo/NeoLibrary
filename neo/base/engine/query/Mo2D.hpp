#pragma once
#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <vector>
// 二次元の Mo's algorithm (二次元莫队)。部分行列 [xl, xr) x [yl, yr) のクエリを集め、4 つの端点を動かす。
//   Mo2D mo;
//   mo.query(xl, xr, yl, yr);
//   mo.run(add_x, erase_x, add_y, erase_y, out);
// add_x(x, yl, yr) は行 x を列 [yl, yr) の範囲で足す。erase_x は外す。add_y(y, xl, xr) は列 y を行 [xl, xr)
// の範囲で足す。行や列を足すときは、そのときのもう一方の範囲が渡るので、その範囲を走査して状態を更新する。
// 並べ方は xl と yl のブロックで、xr と yr を交互に並べる。計算量は 4 端点の移動回数と、1 回の移動で走査する
// 範囲の長さで決まる。
struct Mo2D {
 std::vector<int> qxl_, qxr_, qyl_, qyr_;
 int nx_= 0, ny_= 0;
 void query(int xl, int xr, int yl, int yr) {
  assert(0 <= xl && xl <= xr && 0 <= yl && yl <= yr);
  qxl_.push_back(xl), qxr_.push_back(xr), qyl_.push_back(yl), qyr_.push_back(yr);
  if(nx_ < xr) nx_= xr;
  if(ny_ < yr) ny_= yr;
 }
 int size() const { return (int)qxl_.size(); }
 bool empty() const { return qxl_.empty(); }
 template <class AX, class EX, class AY, class EY, class Out>
 void run(const AX& add_x, const EX& erase_x, const AY& add_y, const EY& erase_y, const Out& out) const {
  const int q= (int)qxl_.size();
  if(q == 0) return;
  const int root= std::max(1, (int)std::sqrt((double)q));
  int bx= nx_ / root, by= ny_ / root;
  if(bx < 1) bx= 1;
  if(by < 1) by= 1;
  std::vector<int> ord(q);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b) {
   const int axb= qxl_[a] / bx, bxb= qxl_[b] / bx;
   if(axb != bxb) return axb < bxb;
   const int ayb= qyl_[a] / by, byb= qyl_[b] / by;
   if(ayb != byb) return (axb & 1) ? ayb > byb : ayb < byb;
   if(qxr_[a] != qxr_[b]) return (ayb & 1) ? qxr_[a] > qxr_[b] : qxr_[a] < qxr_[b];
   return (qxr_[a] & 1) ? qyr_[a] > qyr_[b] : qyr_[a] < qyr_[b];
  });
  int xl= 0, xr= 0, yl= 0, yr= 0;
  for(int id: ord) {
   while(xl > qxl_[id]) add_x(--xl, yl, yr);
   while(xr < qxr_[id]) add_x(xr++, yl, yr);
   while(yl > qyl_[id]) add_y(--yl, xl, xr);
   while(yr < qyr_[id]) add_y(yr++, xl, xr);
   while(xl < qxl_[id]) erase_x(xl++, yl, yr);
   while(xr > qxr_[id]) erase_x(--xr, yl, yr);
   while(yl < qyl_[id]) erase_y(yl++, xl, xr);
   while(yr > qyr_[id]) erase_y(--yr, xl, xr);
   out(id);
  }
 }
};
