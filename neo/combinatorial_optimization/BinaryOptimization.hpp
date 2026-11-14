#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <limits>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>
#include "neo/combinatorial_optimization/MaxFlow.hpp"
// 燃やす埋める。0 か 1 の値を取る変数 x_0, ..., x_{n-1} について、定数、1 変数の項、2 変数の項の和の最小を求める。項は
// add で値の表として足し、solve() で最小を、argmin() で最小を与える割り当てを取り出す。
// 2 変数の項は劣モジュラ (c00 + c11 <= c01 + c10) でなくてよい。そうでない項は、変数の 0 と 1 を入れ替える反転を探して
// 劣モジュラに直す (劣モジュラな項は両端を同じ向きに、優モジュラな項は逆向きに反転する、という制約を 2 色に塗る)。格子の
// 市松模様や、行と列のような二部の形は、問題文のまま書ける。どう反転しても直らないときは solve() が nullopt を返す。
// 中は最小カット (MaxFlow) で、変数ごとに頂点 1 つ、1 次でない 2 変数の項ごとに両向きの辺 1 本を張る。2 変数の項は、
// 2 つの向きの辺と 1 次の部分に分け、1 次の部分の増分ができるだけ 0 に近くなる分け方を選ぶ。禁止のための大きな値は
// 辺の容量に入り、定数の項に漏れない。
// 費用は整数の型。和は中で 128 bit で数え、禁止したい組には大きな費用 (1e18 や 2^62) を置いてよい。最小カットに入りえ
// ない大きな容量は、自明なカット (全部 0 か全部 1) の値 + 1 に縮めて 64 bit で流す。自明なカットの値が 2^61 以上のとき
// は、128 bit の容量で流す。ただし、GNU の拡張が無く __int128 が整数の型と見なされない環境 (-std=c++20 など) では、
// 64 bit で 2^61 まで流し、最小カットがそれを超えると正しい値にならない。最大化は費用の符号を反転して渡す。
template <class Cost> class BinaryOptimization {
 static_assert(std::is_integral_v<Cost>, "費用は整数の型");
 using i64= long long;
 using i128= __int128;
 static constexpr bool WIDE= std::is_integral_v<i128> && std::numeric_limits<i128>::is_specialized;  // MaxFlow<__int128> が組めるか
 using W= std::conditional_t<WIDE, i128, i64>;
public:
 explicit BinaryOptimization(int n= 0): n_(n), a_(n, {0, 0}) {}
 // 変数を 1 つ足し、番号を返す。
 int add_variable() { return a_.push_back({0, 0}), n_++; }
 int size() const { return n_; }
 // 定数の項。
 [[gnu::always_inline]] void add(Cost c) { base_+= c; }
 // x_i が 0 なら c0、1 なら c1。
 [[gnu::always_inline]] void add(int i, Cost c0, Cost c1) {
  assert(0 <= i && i < n_);
  a_[i][0]+= c0, a_[i][1]+= c1;
 }
 // (x_i, x_j) が (0, 0) なら c00、(0, 1) なら c01、(1, 0) なら c10、(1, 1) なら c11。
 [[gnu::always_inline]] void add(int i, int j, Cost c00, Cost c01, Cost c10, Cost c11) {
  assert(0 <= i && i < n_ && 0 <= j && j < n_);
  if(i == j) return add(i, c00, c11);
  const i128 w= i128(c00) + c11 - c01 - c10;
  if(w == 0) {  // 1 変数の項 2 つに分かれる
   a_[i][0]+= c00, a_[i][1]+= c10, a_[j][1]+= i128(c01) - c00;
   return;
  }
  anti_|= w > 0;
  ts_.push_back({i, j, {c00, c01, c10, c11}});
 }
 // 最小。どう反転しても劣モジュラにならない 2 変数の項の組があれば nullopt。
 std::optional<Cost> solve() {
  mf_.reset(), mfw_.reset();
  const int n= n_;
  flip_.assign(n, 0);
  if(anti_ && !decide_flips()) return std::nullopt;
  // y_v = x_v xor flip_v とした 1 変数の表 u_[v][y]。2 変数の項の 1 次の部分もここに足す。y_v = 1 を t の側とする。
  u_.resize(n);
  for(int v= 0; v < n; ++v) u_[v]= {a_[v][flip_[v]], a_[v][flip_[v] ^ 1]};
  for(const auto& e: ts_) {
   const auto [t00, p, q, l, r]= split(e);
   u_[e.i][0]+= t00, u_[e.i][1]+= t00 + p, u_[e.j][1]+= q;
  }
  // 最小カットの上界として、全部 0、全部 1、1 変数の項だけで最小になる割り当て (強制した変数が多いときに効く) のカットの
  // 値を見る。
  i128 cst= base_, ts= 0, tt= 0, tg= 0;
  for(int v= 0; v < n; ++v) {
   const i128 m= std::min(u_[v][0], u_[v][1]);
   cst+= m, u_[v][0]-= m, u_[v][1]-= m, ts+= u_[v][1], tt+= u_[v][0];
  }
  for(const auto& e: ts_)
   if(const bool yi= u_[e.i][1] == 0, yj= u_[e.j][1] == 0; yi != yj) {
    const auto [t00, p, q, l, r]= split(e);
    tg+= yi ? l : r;
   }
  const i128 ub= std::min({ts, tt, tg}), lim= i128(1) << 61;  // 最小カットは ub 以下
  nsol_= n;
  i128 f;
  if(ub < lim) f= run(mf_.emplace(n + 2), ub + 1, false);
  else if constexpr(WIDE) f= run(mfw_.emplace(n + 2), ub + 1, false);
  else f= run(mf_.emplace(n + 2), lim, true);
  const i128 ans= cst + f;
  if(ans > std::numeric_limits<Cost>::max()) return std::numeric_limits<Cost>::max();
  if(ans < std::numeric_limits<Cost>::min()) return std::numeric_limits<Cost>::min();
  return Cost(ans);
 }
 // 直前の solve の最小を与える割り当て (x_i は 0 か 1)。solve が nullopt を返したあとは空。
 std::vector<int> argmin() {
  if(!mf_ && !mfw_) return {};
  const std::vector<bool> cut= mf_ ? mf_->min_cut(nsol_) : mfw_->min_cut(nsol_);
  std::vector<int> x(nsol_);
  for(int v= 0; v < nsol_; ++v) x[v]= int(!cut[v]) ^ flip_[v];
  return x;
 }
private:
 struct Term {
  int i, j;
  Cost c[4];  // c[(x_i << 1) | x_j]
 };
 int n_, nsol_= 0;
 i128 base_= 0;
 std::vector<std::array<i128, 2>> a_, u_;
 std::vector<Term> ts_;
 bool anti_= false;
 std::vector<unsigned char> flip_;
 std::optional<MaxFlow<i64>> mf_;
 std::optional<MaxFlow<W>> mfw_;
 // 反転したあとの表 t を、t00 + p y_i + q y_j + l [y_i = 1, y_j = 0] + r [y_i = 0, y_j = 1] に分ける (l, r >= 0)。
 // p はできるだけ 0 に近くとる。{t00, p, q, l, r} を返す。
 std::array<i128, 5> split(const Term& e) const {
  const int fi= flip_[e.i] << 1, fj= flip_[e.j];
  const i128 t00= e.c[fi | fj], t01= e.c[fi | (fj ^ 1)], t10= e.c[(fi ^ 2) | fj], t11= e.c[(fi ^ 2) | (fj ^ 1)];
  const i128 p= std::clamp(i128(0), t11 - t01, t10 - t00), q= t11 - t00 - p;
  return {t00, p, q, t10 - t00 - p, t01 - t00 - q};
 }
 template <class C> i128 run(MaxFlow<C>& mf, i128 cmax, bool limited) {
  auto cl= [cmax](i128 c) { return C(c < cmax ? c : cmax); };
  const int n= nsol_, s= n, t= n + 1;
  for(int v= 0; v < n; ++v)
   if(u_[v][1] > 0) mf.add_edge(s, v, cl(u_[v][1]));  // y_v = 1 で u_[v][1]
   else if(u_[v][0] > 0) mf.add_edge(v, t, cl(u_[v][0]));
  for(const auto& e: ts_) {
   const auto [t00, p, q, l, r]= split(e);
   assert(l >= 0 && r >= 0 && l + r > 0);
   mf.add_edge(e.j, e.i, cl(l), cl(r));  // j から i へ l (y_i = 1 かつ y_j = 0)、i から j へ r
  }
  return limited ? mf.flow(s, t, C(cmax)) : mf.flow(s, t);
 }
 // 劣モジュラな項は両端の反転を同じに、優モジュラな項は違うものにする。重み付きの union-find で 2 色に塗る。
 bool decide_flips() {
  const int n= n_;
  std::vector<int> up(n, -1);
  std::vector<unsigned char> pr(n, 0);
  auto find= [&](int v, unsigned char& p) {
   for(p= 0; up[v] >= 0; v= up[v]) p^= pr[v];
   return v;
  };
  for(const auto& e: ts_) {
   const unsigned char want= i128(e.c[0]) + e.c[3] - e.c[1] - e.c[2] > 0;
   unsigned char pi, pj;
   int ri= find(e.i, pi), rj= find(e.j, pj);
   if(ri == rj) {
    if((pi ^ pj) != want) return false;
    continue;
   }
   if(up[ri] > up[rj]) std::swap(ri, rj);
   up[ri]+= up[rj], up[rj]= ri, pr[rj]= pi ^ pj ^ want;
  }
  for(int v= 0; v < n; ++v) find(v, flip_[v]);
  return true;
 }
};
