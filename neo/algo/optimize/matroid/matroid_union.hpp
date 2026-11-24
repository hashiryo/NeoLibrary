#pragma once
#include <vector>
#include "neo/algo/optimize/matroid/matroid_intersection.hpp"
// マトロイドの和。台集合 [0, n) が同じ k 個のマトロイド M_0, ..., M_{k-1} について、互いに素な独立集合 I_0, ..., I_{k-1}
// の和集合が最大になる分け方を求める。要素 x を入れた集合の番号 owner[x] (どれにも入れない要素は -1) を返す。
// - matroid_union(ms) は、和集合の要素数が最大になる分け方を返す。
// - matroid_union(ms, order) は、order の順に要素を入れ、入らない要素は捨てる。和もマトロイドなので、order を重みの大きい
//   順に並べると、和集合の重みの和が最大になる (辺素な 2 本の全域森の重みの和の最小なら、重みの小さい順)。
// ms の型は同じで、matroid_intersection に渡す型 (size、build、can_add と、can_exchange か circuit) か、さらに add(x)
// (I + x が独立のとき、中の構造を作り直さずに x を足し、can_add だけを正しく保つ) を持つ型。add を持つ型は、そのまま入る
// 要素ごとの build を省ける。ms の各マトロイドは、中の構造を書き換えられる。
// Edmonds のマトロイドの分割の算法。要素を 1 つずつ入れ、どの M_j にもそのまま入らなければ、交換の列を幅優先で探す。要素 u
// を I_j に入れると、I_j + u の閉路の要素 y が押し出されるので、y をほかの I_j' に入れ直す。どれかにそのまま入る要素に着いた
// ら、列に沿って入れ替える。列が見つからない要素は捨てる。
// 速くする工夫が 2 つある。見つからなかった探索でたどった要素の集合 S は、S の要素をどの I_j に入れても S の要素しか押し出
// さず、そのまま入る要素も無い。入れ替えで I_j と S の共通部分は変わらないので、この性質はずっと保たれ、S は以後のどの探索
// の最短の列にも現れない。そこで以後の探索では S の要素を飛ばす (外すと 10 倍から 100 倍遅くなった)。また order を渡さない
// ときは、入れ替えで各 I_j の張る空間が変わらず、そのまま入る要素で広がるだけなので、一度どの I_j にもそのまま入らなかった
// 要素は以後もそのまま入らない。そこでそのまま入る要素を先に全部入れてから、残りを探索で入れる。
// procon-judge の self-matroid-forest-union と self-matroid-xor-union で書き比べた形 (中身は
// problems/_shared/matroid/union_algos.hpp)。記録は algo-notes の notes/matroid.md にある。
namespace matroid_internal {
template <class M>
concept HasAdd= requires(M& m, int x) { m.add(x); };
// matroid_union が使う操作の包み。circuit が無い型には、I の要素を全部 can_exchange で調べて閉路を作る。
template <class M> class UnionAdapter {
 M& m_;
 const std::vector<int>* I_= nullptr;
public:
 explicit UnionAdapter(M& m): m_(m) {}
 void build(const std::vector<int>& I) { I_= &I, m_.build(I); }
 [[gnu::always_inline]] bool can_add(int x) { return m_.can_add(x); }
 void add(int x) {
  if constexpr(HasAdd<M>) m_.add(x);
 }
 template <class F> void circuit(int x, F&& f) {
  if constexpr(HasCircuit<M>) m_.circuit(x, f);
  else
   for(int y: *I_)
    if(m_.can_exchange(y, x)) f(y);
 }
};
template <class M> std::vector<int> union_impl(std::vector<M>& ms, const std::vector<int>* order) {
 const int k= ms.size(), n= k ? ms[0].size() : 0;
 std::vector<UnionAdapter<M>> as;
 as.reserve(k);
 for(auto& m: ms) as.emplace_back(m);
 std::vector<int> owner(n, -1), pos(n, -1), par(n), q;
 std::vector<std::vector<int>> I(k);
 // st[j]: 0 は I[j] で build 済み、1 は build のあと add だけした (can_add は正しい)、2 は作り直しが要る
 std::vector<char> st(k, 2), frozen(n);
 std::vector<unsigned> vis(n, 0);
 unsigned stamp= 0;
 auto move= [&](int x, int j) {
  if(owner[x] >= 0) {
   auto& v= I[owner[x]];
   const int last= v.back();
   v[pos[x]]= last, pos[last]= pos[x], v.pop_back(), st[owner[x]]= 2;
  }
  owner[x]= j;
  if(j >= 0) pos[x]= I[j].size(), I[j].push_back(x), st[j]= 2;
 };
 auto direct= [&](int s) {
  for(int j= 0; j < k; ++j) {
   if(st[j] == 2) as[j].build(I[j]), st[j]= 0;
   if(as[j].can_add(s)) {
    if constexpr(HasAdd<M>) {
     owner[s]= j, pos[s]= I[j].size(), I[j].push_back(s), as[j].add(s);
     if(st[j] == 0) st[j]= 1;
    } else move(s, j);
    return true;
   }
  }
  return false;
 };
 auto search= [&](int s) {
  for(int j= 0; j < k; ++j)
   if(st[j] != 0) as[j].build(I[j]), st[j]= 0;
  ++stamp;
  vis[s]= stamp, par[s]= -1;
  q.assign(1, s);
  int got= -1, got_j= -1;
  for(size_t h= 0; h < q.size() && got < 0; ++h) {
   const int u= q[h];
   auto visit= [&](int y) {
    if(got >= 0 || vis[y] == stamp || frozen[y]) return;
    vis[y]= stamp, par[y]= u, q.push_back(y);
    for(int j= 0; j < k; ++j)
     if(j != owner[y] && as[j].can_add(y)) {
      got= y, got_j= j;
      return;
     }
   };
   for(int j= 0; j < k && got < 0; ++j)
    if(j != owner[u]) as[j].circuit(u, visit);
  }
  if(got < 0) {
   for(int y: q) frozen[y]= 1;
   return;
  }
  for(int x= got, j= got_j; x >= 0;) {
   const int old= owner[x], p= par[x];
   move(x, j), j= old, x= p;
  }
 };
 if(order) {
  for(int s: *order)
   if(!direct(s)) search(s);
 } else {
  std::vector<int> rest;
  for(int s= 0; s < n; ++s)
   if(!direct(s)) rest.push_back(s);
  for(int s: rest) search(s);
 }
 return owner;
}
}  // namespace matroid_internal
template <class M> std::vector<int> matroid_union(std::vector<M>& ms) { return matroid_internal::union_impl(ms, nullptr); }
template <class M> std::vector<int> matroid_union(std::vector<M>& ms, const std::vector<int>& order) { return matroid_internal::union_impl(ms, &order); }
