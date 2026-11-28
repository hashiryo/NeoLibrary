#pragma once
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
// 回滚する Mo's algorithm (回滚莫队)。区間 [l, r) のクエリを集め、追加だけをして、クエリごとに一時的に足した
// ぶんを undo() で巻き戻す。削除で答えを戻せないもの (最大値、最頻値など) に使う。
//   RollbackMo mo;
//   mo.query(l, r);
//   mo.run(add, undo, out);
// add(i) は i を足す。undo() は直前に足したものを 1 つ外す (LIFO)。各ブロックの終わりに、そのブロックで
// 足したものを全部 undo() して空に戻す。状態は run の前と各ブロックの始まりで空にしておく。
// 計算量は O(N sqrt(Q) T) で、T は add か undo 1 回の計算量。
struct RollbackMo {
 std::vector<int> ql_, qr_;
 int n_= 0;
 void query(int l, int r) {
  ql_.push_back(l), qr_.push_back(r);
  if(n_ < r) n_= r;
 }
 int size() const { return (int)ql_.size(); }
 bool empty() const { return ql_.empty(); }
 template <class Add, class Undo, class Out>
 void run(const Add& add, const Undo& undo, const Out& out) const {
  const int q= (int)ql_.size();
  if(q == 0) return;
  int bs= n_ / std::max(1, (int)std::sqrt((double)q));
  if(bs < 1) bs= 1;
  std::vector<int> ord(q);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b) {
   const int ab= ql_[a] / bs, bb= ql_[b] / bs;
   if(ab != bb) return ab < bb;
   return qr_[a] < qr_[b];
  });
  for(int i= 0; i < q;) {
   const int bb= ql_[ord[i]] / bs;
   const int be= std::min(n_, (bb + 1) * bs);
   int j= i;
   while(j < q && ql_[ord[j]] / bs == bb) ++j;
   int rcur= be, kept= 0;
   for(int k= i; k < j; ++k) {
    const int id= ord[k], l= ql_[id], r= qr_[id];
    if(r <= be) {
     int cnt= 0;
     for(int x= l; x < r; ++x) add(x), ++cnt;
     out(id);
     while(cnt--) undo();
    } else {
     while(rcur < r) add(rcur++), ++kept;
     int cnt= 0;
     for(int x= l; x < be; ++x) add(x), ++cnt;
     out(id);
     while(cnt--) undo();
    }
   }
   while(kept--) undo();
   i= j;
  }
 }
};
