#pragma once
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
