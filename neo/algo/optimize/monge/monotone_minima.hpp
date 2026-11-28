#pragma once
#include <algorithm>
#include <vector>
// 行の最小の列が行について単調な行列の、各行の最小の列を分割統治で求める。value の評価は O(h log w) 回。
// value(r, c) は行 r 列 c の値を T で返す。有効でない (r, c) は T の最大値を返すこと。
// 行の最小の列が行について単調なのは、行列が totally monotone のとき (Monge や逆 Monge を含む)。
// procon-judge の yosupo-min-plus-convolution-convex-arbitrary の核を、値の関数を受ける形にしたもの。
// 記録は algo-notes の notes/dp-optimization.md にある。
template <class T, class F> std::vector<int> monotone_minima(int h, int w, const F& value) {
 std::vector<int> best(h);
 if(h == 0 || w == 0) return best;
 auto rec= [&](auto& rec, int h1, int h2, int w1, int w2) -> void {
  if(h1 == h2) return;
  const int hm= (h1 + h2) / 2;
  int bj= w1;
  T bv= value(hm, w1);
  for(int j= w1 + 1; j < w2; ++j) {
   const T v= value(hm, j);
   if(v < bv) bv= v, bj= j;
  }
  best[hm]= bj;
  rec(rec, h1, hm, w1, bj + 1);
  rec(rec, hm + 1, h2, bj, w2);
 };
 rec(rec, 0, h, 0, w);
 return best;
}
// 行 r の有効な列が [lo(r), hi(r)) のときの、各行の最小の列。lo と hi は広義単調増加で、どの行も
// lo(r) < hi(r)。行の最小の列が行について単調であること (Monge など)。有効でないセルを評価しないので、
// 値の関数の中で範囲を判定して T の最大値を返す monotone_minima より速い。dp の各層のように有効な列が
// 階段になる使い方で使う。procon-judge の cses-2086 の核 (値だけを返す形と対) を写したもの。
// 記録は algo-notes の notes/dp-optimization.md にある。
template <class T, class F, class Lo, class Hi> std::vector<int> monotone_minima_bounded(int h, const F& value, const Lo& lo, const Hi& hi) {
 std::vector<int> best(h);
 if(h == 0) return best;
 auto rec= [&](auto& rec, int h1, int h2, int w1, int w2) -> void {
  if(h1 == h2) return;
  const int hm= (h1 + h2) / 2;
  const int wl= std::max(w1, lo(hm)), wr= std::min(w2, hi(hm));
  int bj= wl;
  T bv= value(hm, wl);
  for(int j= wl + 1; j < wr; ++j) {
   const T v= value(hm, j);
   if(v < bv) bv= v, bj= j;
  }
  best[hm]= bj;
  rec(rec, h1, hm, w1, bj + 1);
  rec(rec, hm + 1, h2, bj, w2);
 };
 rec(rec, 0, h, lo(0), hi(h - 1));
 return best;
}
// 上の関数の、各行の最小の値だけを返す形。呼ぶ側が添字から値を引き直さないので、dp の各層のように
// 最小の値がそのまま次の層に要るときはこちらが速い (procon-judge の cses-2086 で x64-gcc 56.4 対 59.2 ms)。
template <class T, class F, class Lo, class Hi> std::vector<T> monotone_minima_bounded_values(int h, const F& value, const Lo& lo, const Hi& hi) {
 std::vector<T> best(h);
 if(h == 0) return best;
 auto rec= [&](auto& rec, int h1, int h2, int w1, int w2) -> void {
  if(h1 == h2) return;
  const int hm= (h1 + h2) / 2;
  const int wl= std::max(w1, lo(hm)), wr= std::min(w2, hi(hm));
  int bj= wl;
  T bv= value(hm, wl);
  for(int j= wl + 1; j < wr; ++j) {
   const T v= value(hm, j);
   if(v < bv) bv= v, bj= j;
  }
  best[hm]= bv;
  rec(rec, h1, hm, w1, bj + 1);
  rec(rec, hm + 1, h2, bj, w2);
 };
 rec(rec, 0, h, lo(0), hi(h - 1));
 return best;
}
