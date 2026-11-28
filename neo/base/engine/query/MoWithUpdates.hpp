#pragma once
#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <vector>
// 時間軸つきの Mo's algorithm (带修莫队)。クエリ (l, r, t) を集める。t は適用済みの更新の数で、状態は更新
// 0..t-1 を適用したあとの列の区間 [l, r)。add と erase は区間の端点、apply(t, l, r) と undo(t, l, r) は
// 時刻の移動で、l と r はそのときの区間 (更新の位置が区間に入るかを渡す側で判定する)。
//   MoWithUpdates mo;
//   mo.query(l, r, t);
//   mo.run(add_left, add_right, erase_left, erase_right, apply, undo, out);
// 並べ方は l、r、t の 3 本で、ブロックの幅は n^{2/3}。計算量は O((N + Q) N^{2/3} T) 程度 (N=Q)。
struct MoWithUpdates {
 std::vector<int> ql_, qr_, qt_;
 int n_= 0;
 void query(int l, int r, int t) {
  assert(0 <= l && l <= r && t >= 0);
  ql_.push_back(l), qr_.push_back(r), qt_.push_back(t);
  if(n_ < r) n_= r;
 }
 int size() const { return (int)ql_.size(); }
 bool empty() const { return ql_.empty(); }
 template <class AL, class AR, class EL, class ER, class AP, class UN, class Out>
 void run(const AL& add_left, const AR& add_right, const EL& erase_left, const ER& erase_right, const AP& apply, const UN& undo, const Out& out) const {
  const int q= (int)ql_.size();
  if(q == 0) return;
  int bs= 1;
  if(n_ > 1) bs= (int)std::pow((double)n_, 2.0 / 3.0);
  if(bs < 1) bs= 1;
  std::vector<int> ord(q);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b) {
   const int alb= ql_[a] / bs, blb= ql_[b] / bs;
   if(alb != blb) return alb < blb;
   const int arb= qr_[a] / bs, brb= qr_[b] / bs;
   if(arb != brb) return (alb & 1) ? arb > brb : arb < brb;
   return (arb & 1) ? qt_[a] > qt_[b] : qt_[a] < qt_[b];
  });
  int l= 0, r= 0, t= 0;
  for(int id: ord) {
   const int ql= ql_[id], qr= qr_[id], qt= qt_[id];
   while(l > ql) add_left(--l);
   while(r < qr) add_right(r++);
   while(l < ql) erase_left(l++);
   while(r > qr) erase_right(--r);
   while(t < qt) apply(t++, l, r);
   while(t > qt) undo(--t, l, r);
   out(id);
  }
 }
};
