#pragma once
#include <cassert>
#include <limits>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>
#include "neo/combinatorial_optimization/MaxFlow.hpp"
// 流量に下限と上限のある網の流れ。add_edge(from, to, lower, upper) で辺を足し、add_supply(v, x) で頂点の湧き出しを足す。
// solve() は、どの辺も下限と上限の間にあり、どの頂点でも出る量から入る量を引いたものが湧き出しに等しい流れ (b-flow) が
// あるかを返す。max_flow(s, t) と min_flow(s, t) は、s と t だけはその釣り合いを外し、s から t へ F (0 以上) だけ流す
// 流れのうち、F の最大と最小を返す。流れが無ければ nullopt。解いたあとは get_edge(i) と edges() で辺ごとの流量を読める。
// 負の下限は、to から from へ流れてよい量を表す (無向辺は add_edge(u, v, -c, c))。
// 中身は MaxFlow を使う決まった形。辺ごとに、下限と上限の間で 0 にいちばん近い量を先に流しておき、残りの幅を容量にした
// 辺を MaxFlow に張る。先に流した分で頂点に出る過不足は、超頂点 S から過剰な頂点へ、不足の頂点から超頂点 T への辺で
// 埋め、S から T へ全部流れるかを見る。s と t のあるときは、t から s へ容量の無限大の辺 (補助の辺) を足して同じことを
// する。最大は、そのまま s から t へ流せるだけ流す (補助の辺の逆向きの弧は s から t へ直接つながるので、最大流は
// 必ずそれを使い切り、増えた量が F の最大になる)。最小は、補助の辺の容量を今の流量 f0 に下げて止め、t から s へ f0 まで
// 戻す。下限も湧き出しも無ければ (先に流す量が全部 0)、超頂点を使わずに MaxFlow で s から t へ流すだけになる。
// 辺は MaxFlow に直接溜め、最初の問い合わせは、その網に超頂点と補助の辺を足して流す。2 回目からの問い合わせと、問い合わせ
// のあとの辺や頂点の追加では、溜めた辺の元の容量から MaxFlow を組み直す (前の流れは捨てる)。
// Cap は符号付きの整数の型。辺ごとの upper - lower、下限の絶対値の和と湧き出しの絶対値の和、s から t への流量の答えが
// Cap に収まること。
template <class Cap> class MaxFlowLowerBound {
 static_assert(std::is_integral_v<Cap> && std::is_signed_v<Cap>, "容量は符号付きの整数の型");
public:
 struct Edge {
  int from, to;
  Cap lower, upper, flow;  // flow は from から to へ正、逆向きに流れていれば負
 };
 explicit MaxFlowLowerBound(int n= 0): n_(n), ex_(n), mf_(n) {}
 int add_vertex() { return reset(), ok_= false, ex_.push_back(0), mf_.add_vertex(), n_++; }
 // from から to へ、流量が lower 以上 upper 以下の辺。番号を返す。
 int add_edge(int from, int to, Cap lower, Cap upper) {
  assert(0 <= from && from < n_ && 0 <= to && to < n_ && lower <= upper);
  reset(), ok_= false;
  const Cap b= lower > 0 ? lower : upper < 0 ? upper : Cap(0);  // 先に流しておく量
  base_.push_back(b), ex_[from]-= b, ex_[to]+= b;
  return mf_.add_edge(from, to, upper - b, b - lower);
 }
 // 頂点 v の湧き出しに x を足す (負なら吸い込み)。
 void add_supply(int v, Cap x) { ok_= false, ex_[v]+= x; }
 // 下限と上限と湧き出しを満たす流れ (b-flow) があれば true。
 bool solve() { return ok_= saturate(-1, -1); }
 // 下限と上限と湧き出しを満たし、s から t へ F (0 以上) だけ流す流れのうち、F の最大。無ければ nullopt。
 std::optional<Cap> max_flow(int s, int t) {
  assert(0 <= s && s < n_ && 0 <= t && t < n_);
  if(s == t) return solve() ? std::optional<Cap>(Cap(0)) : std::nullopt;
  if(!(ok_= saturate(s, t))) return std::nullopt;
  return used_= true, mf_.flow(s, t);
 }
 // F の最小。無ければ nullopt。
 std::optional<Cap> min_flow(int s, int t) {
  assert(0 <= s && s < n_ && 0 <= t && t < n_);
  if(s == t) return solve() ? std::optional<Cap>(Cap(0)) : std::nullopt;
  if(!(ok_= saturate(s, t))) return std::nullopt;
  if(aux_ < 0) return Cap(0);  // 先に流した流れのままで F = 0
  const Cap f0= mf_.get_edge(aux_).flow;
  mf_.change_cap(aux_, f0, 0, t, s);  // 流量は新しい容量に収まるので、流れは動かない
  return f0 - mf_.flow(t, s, f0);
 }
 // 直前の solve、max_flow、min_flow の流れ。流れが無かったときと、そのあとに辺や頂点や湧き出しを足したときは、flow は
 // 下限になる。
 Edge get_edge(int i) {
  const typename MaxFlow<Cap>::Edge e= mf_.get_edge(i);
  const Cap b= base_[i];
  return {e.from, e.to, Cap(b - e.rev_cap), Cap(b + e.cap), Cap(ok_ ? b + e.flow : b - e.rev_cap)};
 }
 std::vector<Edge> edges() {
  std::vector<Edge> ret(base_.size());
  for(int i= 0; i < (int)ret.size(); ++i) ret[i]= get_edge(i);
  return ret;
 }
 int size() const { return n_; }
private:
 int n_;
 std::vector<Cap> ex_, base_;  // ex_: 先に流した分と湧き出しで頂点に残る過剰 (負なら不足)
 MaxFlow<Cap> mf_;
 bool used_= false, ok_= false;  // used_: mf_ に超頂点や流れがある
 int aux_= -1;                   // 補助の辺の番号 (足していなければ -1)
 // mf_ を、溜めた辺の元の容量だけの網に戻す。
 void reset() {
  if(!used_) return;
  MaxFlow<Cap> g(n_);
  for(int i= 0; i < (int)base_.size(); ++i) {
   const typename MaxFlow<Cap>::Edge e= mf_.get_edge(i);
   g.add_edge(e.from, e.to, e.cap, e.rev_cap);
  }
  mf_= std::move(g), used_= false;
 }
 // 過不足を超頂点 S と T への辺で埋めて S から T へ流し、全部流れれば true。s >= 0 なら補助の辺 t → s も足す。過不足が
 // 無ければ、何も足さずに true を返す (aux_ は -1)。
 bool saturate(int s, int t) {
  reset(), aux_= -1;
  Cap pos= 0, neg= 0;
  for(int v= 0; v < n_; ++v) (ex_[v] > 0 ? pos : neg)+= ex_[v];
  if(pos + neg != 0) return false;
  if(pos == 0) return true;
  used_= true;
  const int S= mf_.add_vertex(), T= mf_.add_vertex();
  for(int v= 0; v < n_; ++v)
   if(ex_[v] > 0) mf_.add_edge(S, v, ex_[v]);
   else if(ex_[v] < 0) mf_.add_edge(v, T, -ex_[v]);
  if(s >= 0) aux_= mf_.add_edge(t, s, std::numeric_limits<Cap>::max());
  return mf_.flow(S, T) == pos;
 }
};
