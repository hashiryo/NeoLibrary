#pragma once
#include <numeric>
#include <vector>
// 行の最小の列が行について単調な行列の、各行の最小の列を SMAWK で求める。value の評価は O(h + w) 回。
// value(r, c) は行 r 列 c の値を T で返す。有効でない (r, c) は T の最大値を返すこと。
// 列を減らしながら奇数番目の行で再帰し、偶数番目の行は両隣の答えの間だけを見る。
// - 列を減らす段は、同じ (行, 山の頂点) の値を覚えて 2 回目の評価を省く。
// - 偶数番目の行の走査は、最良の値を覚えて 1 候補 1 評価にする。
// procon-judge の yosupo-min-plus-convolution-convex-arbitrary の核を、値の関数を受ける形にしたもの。
// 記録は algo-notes の notes/dp-optimization.md にある。
template <class T, class F> std::vector<int> smawk(int h, int w, const F& value) {
 std::vector<int> best(h);
 if(h == 0 || w == 0) return best;
 std::vector<int> rows(h), cols(w);
 std::iota(rows.begin(), rows.end(), 0), std::iota(cols.begin(), cols.end(), 0);
 auto rec= [&](auto& rec, const std::vector<int>& rows, const std::vector<int>& cols) -> void {
  if(rows.empty()) return;
  // 列を減らす。各行の最小になりえない列を落とす。
  std::vector<int> stk;
  stk.reserve(cols.size());
  int crow= -1, ctop= -1;
  T cval= 0;
  for(int col : cols) {
   while(!stk.empty()) {
    const int row= rows[(int)stk.size() - 1], top= stk.back();
    T topv;
    if(row == crow && top == ctop) topv= cval;
    else topv= cval= value(row, top), crow= row, ctop= top;
    if(value(row, col) <= topv) stk.pop_back();
    else break;
   }
   if((int)stk.size() < (int)rows.size()) stk.push_back(col);
  }
  // 奇数番目の行を解く。
  std::vector<int> odd;
  odd.reserve(rows.size() / 2);
  for(size_t i= 1; i < rows.size(); i+= 2) odd.push_back(rows[i]);
  rec(rec, odd, stk);
  // 偶数番目の行は、両隣の答えの列の間だけを見る。
  size_t pos= 0;
  for(size_t i= 0; i < rows.size(); i+= 2) {
   const int hm= rows[i];
   size_t lp= 0;
   if(i > 0) {
    while(stk[pos] != best[rows[i - 1]]) ++pos;
    lp= pos;
   } else pos= 0;
   size_t hp= stk.size() - 1;
   if(i + 1 < rows.size()) {
    size_t q= pos;
    while(stk[q] != best[rows[i + 1]]) ++q;
    hp= q;
   }
   int bj= stk[lp];
   T bv= value(hm, bj);
   for(size_t t= lp + 1; t <= hp; ++t) {
    const int j= stk[t];
    const T v= value(hm, j);
    if(v < bv) bv= v, bj= j;
   }
   best[hm]= bj;
  }
 };
 rec(rec, rows, cols);
 return best;
}
