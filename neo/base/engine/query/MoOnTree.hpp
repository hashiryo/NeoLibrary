#pragma once
#include <algorithm>
#include <bit>
#include <utility>
#include <vector>
#include "neo/base/engine/query/Mo.hpp"
#include "neo/base/engine/query/MoWithUpdates.hpp"
// 木のパスの Mo's algorithm (树上莫队)。木を括弧列 (Euler tour で各頂点を 2 回並べた長さ 2N の列) に落とし、
// パス u-v をその区間に対応させる。区間の端点を 1 つ動かすと、その位置の頂点が「入る」と「出る」を交互に
// 繰り返すので、add と erase ではなく toggle(v, entering) を渡す。entering が真なら v が区間に入り、偽なら出る。
//   MoOnTree mo(n);
//   mo.add_edge(u, v);
//   mo.query(u, v);
//   mo.build(root);             // 辺とクエリを足したあとに 1 回
//   mo.run(toggle, out);
// パスに LCA が入らない形の区間では、答えるときだけ LCA を足して外す。LCA は二分探索 (前計算 O(N log N))。
// 計算量は O((N + Q) sqrt(N) T log N) 程度で、T は toggle 1 回の計算量。木が連結であることを前提にする。
struct MoOnTree {
 int n_= 0;
 std::vector<std::pair<int, int>> edges_;
 std::vector<int> qu_, qv_;
 std::vector<int> tour_, in_, out_, dep_, par_, up_;
 int lg_= 0;
 std::vector<int> ql_, qr_, extra_;
 std::vector<char> on_;
 MoOnTree()= default;
 explicit MoOnTree(int n): n_(n) {}
 int size() const { return n_; }
 void add_edge(int u, int v) {
  edges_.emplace_back(u, v);
  if(n_ <= u) n_= u + 1;
  if(n_ <= v) n_= v + 1;
 }
 void query(int u, int v) {
  qu_.push_back(u), qv_.push_back(v);
  if(n_ <= u) n_= u + 1;
  if(n_ <= v) n_= v + 1;
 }
 int lca(int u, int v) const {
  if(dep_[u] < dep_[v]) std::swap(u, v);
  int d= dep_[u] - dep_[v];
  for(int k= 0; d; ++k, d>>= 1)
   if(d & 1) u= up_[k * n_ + u];
  if(u == v) return u;
  for(int k= lg_ - 1; k >= 0; --k)
   if(up_[k * n_ + u] != up_[k * n_ + v]) u= up_[k * n_ + u], v= up_[k * n_ + v];
  return up_[u];
 }
 // 括弧列と LCA を作り、クエリを区間 [ql, qr) と、足す LCA (extra、不要なら -1) に直す。
 void build(int root= 0) {
  if(n_ == 0) return;
  std::vector<int> deg(n_, 0);
  for(auto [u, v]: edges_) ++deg[u], ++deg[v];
  std::vector<int> start(n_ + 1, 0);
  for(int i= 0; i < n_; ++i) start[i + 1]= start[i] + deg[i];
  std::vector<int> adj(start[n_]), pos(start);
  for(auto [u, v]: edges_) adj[pos[u]++]= v, adj[pos[v]++]= u;
  in_.resize(n_), out_.resize(n_), dep_.resize(n_), par_.resize(n_);
  tour_.clear(), tour_.reserve(2 * n_);
  std::vector<int> st, it(n_);
  par_[root]= root, dep_[root]= 0, in_[root]= 0, tour_.push_back(root), st.push_back(root);
  for(int i= 0; i < n_; ++i) it[i]= start[i];
  while(!st.empty()) {
   const int v= st.back();
   if(it[v] < start[v + 1]) {
    const int u= adj[it[v]++];
    if(u == par_[v]) continue;
    par_[u]= v, dep_[u]= dep_[v] + 1, in_[u]= (int)tour_.size(), tour_.push_back(u);
    st.push_back(u);
   } else {
    out_[v]= (int)tour_.size(), tour_.push_back(v);
    st.pop_back();
   }
  }
  lg_= (int)std::bit_width((unsigned)std::max(n_, 1));
  up_.assign((std::size_t)lg_ * n_, root);
  for(int v= 0; v < n_; ++v) up_[v]= par_[v];
  for(int k= 1; k < lg_; ++k)
   for(int v= 0; v < n_; ++v) up_[(std::size_t)k * n_ + v]= up_[(std::size_t)(k - 1) * n_ + up_[(std::size_t)(k - 1) * n_ + v]];
  const int q= (int)qu_.size();
  ql_.resize(q), qr_.resize(q), extra_.resize(q);
  for(int i= 0; i < q; ++i) {
   int u= qu_[i], v= qv_[i];
   if(in_[u] > in_[v]) std::swap(u, v);
   const int w= lca(u, v);
   if(w == u) ql_[i]= in_[u], qr_[i]= in_[v] + 1, extra_[i]= -1;
   else ql_[i]= out_[u], qr_[i]= in_[v] + 1, extra_[i]= w;
  }
 }
 bool is_on(int v) const { return on_[v]; }
 template <class Toggle, class Out>
 void run(const Toggle& toggle, const Out& out) {
  Mo mo;
  for(int i= 0; i < (int)ql_.size(); ++i) mo.query(ql_[i], qr_[i]);
  on_.assign(n_, 0);
  auto tg= [&](int p) {
   const int v= tour_[p];
   toggle(v, on_[v] == 0);
   on_[v]^= 1;
  };
  auto ov= [&](int id) {
   const int w= extra_[id];
   if(w >= 0 && !on_[w]) {
    toggle(w, true), on_[w]= 1, out(id), toggle(w, false), on_[w]= 0;
   } else out(id);
  };
  mo.run(tg, tg, ov);
 }
 template <class Toggle, class Out>
 void run_hilbert(const Toggle& toggle, const Out& out) {
  Mo mo;
  for(int i= 0; i < (int)ql_.size(); ++i) mo.query(ql_[i], qr_[i]);
  on_.assign(n_, 0);
  auto tg= [&](int p) {
   const int v= tour_[p];
   toggle(v, on_[v] == 0);
   on_[v]^= 1;
  };
  auto ov= [&](int id) {
   const int w= extra_[id];
   if(w >= 0 && !on_[w]) {
    toggle(w, true), on_[w]= 1, out(id), toggle(w, false), on_[w]= 0;
   } else out(id);
  };
  mo.run_hilbert(tg, tg, ov);
 }
};

// 木のパスに時間軸を足した Mo's algorithm (带修改树上莫队)。クエリ (u, v, t) を集める。t は適用済みの更新の数。
// apply(t, on) と undo(t, on) の on は、そのときどの頂点が区間に入っているか (on[v] が真なら入っている)。
// 更新の頂点が on なら答えも直す。
struct MoOnTreeWithUpdates {
 MoOnTree tree_;
 std::vector<int> qt_;
 MoOnTreeWithUpdates()= default;
 explicit MoOnTreeWithUpdates(int n): tree_(n) {}
 int size() const { return tree_.n_; }
 void add_edge(int u, int v) { tree_.add_edge(u, v); }
 void query(int u, int v, int t) {
  tree_.query(u, v);
  qt_.push_back(t);
 }
 void build(int root= 0) { tree_.build(root); }
 bool is_on(int v) const { return tree_.on_[v]; }
 template <class Toggle, class Apply, class Undo, class Out>
 void run(const Toggle& toggle, const Apply& apply, const Undo& undo, const Out& out) {
  MoWithUpdates mo;
  for(int i= 0; i < (int)tree_.ql_.size(); ++i) mo.query(tree_.ql_[i], tree_.qr_[i], qt_[i]);
  tree_.on_.assign(tree_.n_, 0);
  auto tg= [&](int p) {
   const int v= tree_.tour_[p];
   toggle(v, tree_.on_[v] == 0);
   tree_.on_[v]^= 1;
  };
  auto ap= [&](int t, int, int) { apply(t, tree_.on_); };
  auto un= [&](int t, int, int) { undo(t, tree_.on_); };
  auto ov= [&](int id) {
   const int w= tree_.extra_[id];
   if(w >= 0 && !tree_.on_[w]) {
    toggle(w, true), tree_.on_[w]= 1, out(id), toggle(w, false), tree_.on_[w]= 0;
   } else out(id);
  };
  mo.run(tg, tg, tg, tg, ap, un, ov);
 }
};
