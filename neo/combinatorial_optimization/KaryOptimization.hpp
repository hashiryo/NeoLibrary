#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <limits>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>
#include "neo/combinatorial_optimization/MaxFlow.hpp"
// k 値の燃やす埋める。変数 x_i が 0 から k_i - 1 の値を取るとき、定数、1 変数の項、2 変数の項の和の最小を求める。項は
// add で値の表 (か、値を返す関数) として足し、solve() で最小を、argmin() で最小を与える割り当てを取り出す。
// 1 変数の項は何でもよい。2 変数の項 f(a, b) は Monge (f(a, b) + f(a + 1, b + 1) <= f(a, b + 1) + f(a + 1, b)) である
// 必要がある。逆向きの不等号の項 (逆 Monge) は、変数の値の並びを逆にする反転を探して Monge に直す (Monge の項は両端を同じ
// 向きに、逆 Monge の項は逆向きに反転する、という制約を 2 色に塗る)。どちらでもない項があるか、どう反転しても直らない
// ときは solve() が nullopt を返す。禁止は大きな費用で書けるが、表は有限の値のまま Monge であること。たとえば
// x_i - x_j >= 2 を禁じるなら、INF [a - b >= 2] は Monge でないので、INF max(0, a - b - 1) と書く。
// 中は最小カット (MaxFlow) で、変数 i に頂点 k_i - 1 個 ([x_i >= a] が 1 なら t の側) を鎖に並べ、1 変数の表の値を鎖の
// 辺の容量に、並びを保つ大きな容量を逆向きに載せる (Ishikawa の形)。2 変数の項は、表の 2 階の差分が 0 でない升目ごとに
// 両向きの辺 1 本で、2 つの向きへの分け方は 1 次の部分の増分ができるだけ 0 に近くなるように選ぶ。禁止のための大きな値は
// 辺の容量に入り、定数の項に漏れない。2 変数の項を全組に足す必要はなく、足した項の分だけの辺になる。
// 費用は整数の型。和は中で 128 bit で数える。最小カットに入りえない大きな容量は、自明なカット (全部 0 か全部 k_i - 1) の
// 値 + 1 に縮めて 64 bit で流す。自明なカットの値が 2^61 以上のときは、128 bit の容量で流す。ただし、GNU の拡張が無く
// __int128 が整数の型と見なされない環境 (-std=c++20 など) では、64 bit で 2^61 まで流し、最小カットがそれを超えると正しい
// 値にならない。最大化は費用の符号を反転して渡す。
template <class Cost> class KaryOptimization {
 static_assert(std::is_integral_v<Cost>, "費用は整数の型");
 using i64= long long;
 using i128= __int128;
 static constexpr bool WIDE= std::is_integral_v<i128> && std::numeric_limits<i128>::is_specialized;  // MaxFlow<__int128> が組めるか
 using W= std::conditional_t<WIDE, i128, i64>;
public:
 // 変数 i は 0 から k[i] - 1 の値を取る (k[i] >= 1)。
 explicit KaryOptimization(const std::vector<int>& k): k_(k), off_(k.size() + 1, 0) {
  for(std::size_t i= 0; i < k.size(); ++i) assert(k[i] >= 1), off_[i + 1]= off_[i] + k[i];
  th_.assign(off_.back(), 0);
 }
 explicit KaryOptimization(int n, int k): KaryOptimization(std::vector<int>(n, k)) {}
 KaryOptimization(): KaryOptimization(std::vector<int>{}) {}
 // 0 から k - 1 の値を取る変数を 1 つ足し、番号を返す。
 int add_variable(int k) {
  assert(k >= 1);
  k_.push_back(k), off_.push_back(off_.back() + k), th_.resize(off_.back(), 0);
  return int(k_.size()) - 1;
 }
 int size() const { return k_.size(); }
 // 定数の項。
 void add(Cost c) { base_+= c; }
 // x_i が a なら c[a] (c の大きさは k_i)。
 void add(int i, const std::vector<Cost>& c) {
  assert(0 <= i && i < size() && int(c.size()) == k_[i]);
  for(int a= 0; a < k_[i]; ++a) th_[off_[i] + a]+= c[a];
 }
 // x_i が a なら f(a)。
 template <class F>
  requires std::is_invocable_v<F&, int>
 void add(int i, F&& f) {
  assert(0 <= i && i < size());
  for(int a= 0; a < k_[i]; ++a) th_[off_[i] + a]+= Cost(f(a));
 }
 // (x_i, x_j) が (a, b) なら c[a k_j + b] (c の大きさは k_i k_j)。
 void add(int i, int j, const std::vector<Cost>& c) {
  assert(0 <= i && i < size() && 0 <= j && j < size() && c.size() == std::size_t(k_[i]) * k_[j]);
  add_pair(i, j, c.data());
 }
 // (x_i, x_j) が (a, b) なら f(a, b)。
 template <class F>
  requires std::is_invocable_v<F&, int, int>
 void add(int i, int j, F&& f) {
  assert(0 <= i && i < size() && 0 <= j && j < size());
  buf_.resize(std::size_t(k_[i]) * k_[j]);
  for(int a= 0; a < k_[i]; ++a)
   for(int b= 0; b < k_[j]; ++b) buf_[std::size_t(a) * k_[j] + b]= Cost(f(a, b));
  add_pair(i, j, buf_.data());
 }
 // 最小。Monge でも逆 Monge でもない 2 変数の項があるか、どう反転しても Monge にならない項の組があれば nullopt。
 std::optional<Cost> solve() {
  mf_.reset(), mfw_.reset();
  if(bad_) return std::nullopt;
  const int n= size();
  rev_.assign(n, 0);
  if(anti_ && !decide_reversals()) return std::nullopt;
  nid_.assign(n + 1, 0);
  int kmax= 1;
  for(int i= 0; i < n; ++i) nid_[i + 1]= nid_[i] + k_[i] - 1, kmax= std::max(kmax, k_[i]);
  // x'_i (x_i を反転したもの) で見た 1 変数の表 u_。2 変数の項の 1 次の部分もここに足す。
  u_.resize(off_[n]);
  for(int i= 0; i < n; ++i)
   for(int a= 0; a < k_[i]; ++a) u_[off_[i] + a]= th_[off_[i] + (rev_[i] ? k_[i] - 1 - a : a)];
  std::vector<i128> ls(kmax), rs(kmax);
  for(const auto& p: ps_) {
   const int ki= k_[p.i], kj= k_[p.j];
   std::fill(ls.begin(), ls.begin() + ki, 0), std::fill(rs.begin(), rs.begin() + kj, 0);
   for(int a= 1; a < ki; ++a)
    for(int b= 1; b < kj; ++b) {
     const auto [l, r]= cell(p, a, b);
     ls[a]+= l, rs[b]+= r;
    }
   const i128 c00= at(p, 0, 0);
   i128 acc= 0;
   for(int a= 0; a < ki; ++a) acc+= ls[a], u_[off_[p.i] + a]+= at(p, a, 0) - acc;
   acc= 0;
   for(int b= 0; b < kj; ++b) acc+= rs[b], u_[off_[p.j] + b]+= at(p, 0, b) - c00 - acc;
  }
  // 最小カットの上界として、全部 0、全部 k_i - 1、1 変数の項だけで最小になる割り当て (強制した変数が多いときに効く) の
  // カットの値を見る。
  i128 cst= base_, ts= 0, tt= 0, tg= 0;
  std::vector<int> g(n);
  for(int i= 0; i < n; ++i) {
   i128* u= u_.data() + off_[i];
   const int a= std::min_element(u, u + k_[i]) - u;
   const i128 m= u[a];
   cst+= m, g[i]= a;
   for(int b= 0; b < k_[i]; ++b) u[b]-= m;
   ts+= u[k_[i] - 1], tt+= u[0];
  }
  for(const auto& p: ps_)
   for(int a= 1; a < k_[p.i]; ++a)
    for(int b= 1; b < k_[p.j]; ++b)
     if(const bool ai= g[p.i] >= a, bj= g[p.j] >= b; ai != bj) {
      const auto [l, r]= cell(p, a, b);
      tg+= ai ? l : r;
     }
  const i128 ub= std::min({ts, tt, tg}), lim= i128(1) << 61;  // 最小カットは ub 以下
  const int nv= nid_[n];
  nsol_= n;
  i128 f;
  if(ub < lim) f= run(mf_.emplace(nv + 2), ub + 1, false);
  else if constexpr(WIDE) f= run(mfw_.emplace(nv + 2), ub + 1, false);
  else f= run(mf_.emplace(nv + 2), lim, true);
  const i128 ans= cst + f;
  if(ans > std::numeric_limits<Cost>::max()) return std::numeric_limits<Cost>::max();
  if(ans < std::numeric_limits<Cost>::min()) return std::numeric_limits<Cost>::min();
  return Cost(ans);
 }
 // 直前の solve の最小を与える割り当て (0 <= x_i < k_i)。solve が nullopt を返したあとは空。
 std::vector<int> argmin() {
  if(!mf_ && !mfw_) return {};
  const int s= nid_[nsol_];
  const std::vector<bool> cut= mf_ ? mf_->min_cut(s) : mfw_->min_cut(s);
  std::vector<int> x(nsol_);
  for(int i= 0; i < nsol_; ++i) {
   int c= 0;
   for(int v= nid_[i]; v < nid_[i + 1]; ++v) c+= !cut[v];
   x[i]= rev_[i] ? k_[i] - 1 - c : c;
  }
  return x;
 }
private:
 struct Pair {
  int i, j;
  bool anti;
  std::size_t pos;  // tab_ の中の表の先頭
 };
 std::vector<int> k_, off_, nid_;
 int nsol_= 0;
 i128 base_= 0;
 std::vector<i128> th_, u_;
 std::vector<Pair> ps_;
 std::vector<Cost> tab_, buf_;
 bool anti_= false, bad_= false;
 std::vector<unsigned char> rev_;
 std::optional<MaxFlow<i64>> mf_;
 std::optional<MaxFlow<W>> mfw_;
 void add_pair(int i, int j, const Cost* c) {
  const int ki= k_[i], kj= k_[j];
  if(i == j) {
   assert(ki == kj);
   for(int a= 0; a < ki; ++a) th_[off_[i] + a]+= c[std::size_t(a) * kj + a];
   return;
  }
  bool neg= false, pos= false;
  for(int a= 1; a < ki; ++a)
   for(int b= 1; b < kj; ++b) {
    const i128 w= i128(c[std::size_t(a) * kj + b]) - c[std::size_t(a - 1) * kj + b] - c[std::size_t(a) * kj + b - 1] + c[std::size_t(a - 1) * kj + b - 1];
    neg|= w < 0, pos|= w > 0;
   }
  if(neg && pos) {
   bad_= true;
   return;
  }
  if(!neg && !pos) {  // f(a, b) = f(a, 0) + f(0, b) - f(0, 0) に分かれる
   for(int a= 0; a < ki; ++a) th_[off_[i] + a]+= c[std::size_t(a) * kj];
   for(int b= 0; b < kj; ++b) th_[off_[j] + b]+= i128(c[b]) - c[0];
   return;
  }
  anti_|= pos;
  ps_.push_back({i, j, pos, tab_.size()});
  tab_.insert(tab_.end(), c, c + std::size_t(ki) * kj);
 }
 // 反転したあとの表の値 f'(a, b)。
 i128 at(const Pair& p, int a, int b) const {
  const int ki= k_[p.i], kj= k_[p.j];
  if(rev_[p.i]) a= ki - 1 - a;
  if(rev_[p.j]) b= kj - 1 - b;
  return tab_[p.pos + std::size_t(a) * kj + b];
 }
 // 升目 (a, b) の 2 階の差分を、[x'_i >= a, x'_j < b] に掛かる l と [x'_i < a, x'_j >= b] に掛かる r に分ける (どちらも 0 以上)。
 // l は、その升目の 2 行 2 列で見た x'_i の側の増分 f'(a, b - 1) - f'(a - 1, b - 1) にできるだけ近くとる。
 std::pair<i128, i128> cell(const Pair& p, int a, int b) const {
  const i128 t00= at(p, a - 1, b - 1), t01= at(p, a - 1, b), t10= at(p, a, b - 1), t11= at(p, a, b);
  const i128 lam= t01 + t10 - t00 - t11, l= std::clamp(t10 - t00, i128(0), lam);
  assert(lam >= 0);
  return {l, lam - l};
 }
 template <class C> i128 run(MaxFlow<C>& mf, i128 cmax, bool limited) {
  auto cl= [cmax](i128 c) { return C(c < cmax ? c : cmax); };
  const int n= nsol_, s= nid_[n], t= s + 1;
  for(int i= 0; i < n; ++i) {
   // 頂点 nid_[i] + a - 1 は [x'_i >= a] で、x'_i = a のときに頂点 a + 1 から a への辺 (頂点 k_i は s、頂点 0 は t) が切れる。
   // 逆向きの容量 cmax は、頂点 a が s の側で a + 1 が t の側になる (並びが崩れる) カットを除く。
   const int k= k_[i], b= nid_[i] - 1;
   const i128* u= u_.data() + off_[i];
   if(k == 1) continue;
   for(int a= 0; a < k; ++a) {
    const bool inner= 0 < a && a < k - 1;
    if(u[a] > 0 || inner) mf.add_edge(a == k - 1 ? s : b + a + 1, a == 0 ? t : b + a, cl(u[a]), inner ? C(cmax) : C(0));
   }
  }
  for(const auto& p: ps_) {
   const int bi= nid_[p.i] - 1, bj= nid_[p.j] - 1;
   for(int a= 1; a < k_[p.i]; ++a)
    for(int b= 1; b < k_[p.j]; ++b) {
     const auto [l, r]= cell(p, a, b);
     if(l > 0 || r > 0) mf.add_edge(bj + b, bi + a, cl(l), cl(r));
    }
  }
  return limited ? mf.flow(s, t, C(cmax)) : mf.flow(s, t);
 }
 // Monge の項は両端の反転を同じに、逆 Monge の項は違うものにする。重み付きの union-find で 2 色に塗る。
 bool decide_reversals() {
  const int n= size();
  std::vector<int> up(n, -1);
  std::vector<unsigned char> pr(n, 0);
  auto find= [&](int v, unsigned char& p) {
   for(p= 0; up[v] >= 0; v= up[v]) p^= pr[v];
   return v;
  };
  for(const auto& e: ps_) {
   const unsigned char want= e.anti;
   unsigned char pi, pj;
   int ri= find(e.i, pi), rj= find(e.j, pj);
   if(ri == rj) {
    if((pi ^ pj) != want) return false;
    continue;
   }
   if(up[ri] > up[rj]) std::swap(ri, rj);
   up[ri]+= up[rj], up[rj]= ri, pr[rj]= pi ^ pj ^ want;
  }
  for(int v= 0; v < n; ++v) find(v, rev_[v]);
  return true;
 }
};
