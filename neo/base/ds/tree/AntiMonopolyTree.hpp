#pragma once
// オフラインの部分木加算と部分木和。頂点 1 個 1 節のクラスタ木 (top tree の一種) で、親・部分木の大きさ・
// 親辺の重み・部分木和・遅延を持つ。2/3 の大きさ平衡 (anti-monopoly) を保ち、splay を使わない。
// 辺の重みは削除時刻の符号反転で、閉路では重みの大きい辺を外す。部分木は辺 (v, p) を切った成分で、
// イベント列の再生で切る・足す・聞く・戻す。Library Checker の dynamic_tree_subtree_add_subtree_sum の
// 首位 (anon123 の 401546) の核を写した。記録は algo-notes の notes/dynamic-tree.md。
//   FoldableAntiMonopolyTree<T>: link / erase / cut / add_vertex / add_all / fold を持つ。
//   OfflineSubtree<T>: link / cut / add_all / fold を積んで solve。fold の答えは fold の順に返る。
#include <algorithm>
#include <cassert>
#include <climits>
#include <cstdint>
#include <concepts>
#include <unordered_map>
#include <utility>
#include <vector>

template <class T>
concept AmTreeValue = std::default_initializable<T> && std::constructible_from<T, int> && requires(T a, T b) {
 a + b, a - b, a* b, a+= b, a-= b;
};

template <class T>
 requires AmTreeValue<T>
class FoldableAntiMonopolyTree {
 int n;
 std::vector<int> parent_, size_, weight_;
 std::vector<T> sum_, lazy_;
 [[gnu::always_inline]] void promote(int v) {
  int p = parent_[v];
  int len = size_[v];
  T x = lazy_[p];
  size_[p] -= len;
  sum_[p] -= sum_[v] + x * T(len);
  parent_[v] = parent_[p];
  lazy_[v] += x;
  sum_[v] += x * T(len);
  if (weight_[v] < weight_[p]) {
   x = lazy_[v];
   lazy_[p] -= x;
   sum_[p] -= x * T(size_[p]);
   size_[v] += size_[p];
   sum_[v] += sum_[p] + x * T(size_[p]);
   std::swap(weight_[v], weight_[p]);
   parent_[p] = v;
  }
 }
 [[gnu::always_inline]] void upward_maintain(int v) {
  while (~parent_[v]) {
   int p = parent_[v];
   if (3LL * size_[v] <= 2LL * size_[p]) {
    v = p;
    continue;
   }
   promote(v);
  }
 }
 [[gnu::always_inline]] int root(int v) const {
  while (~parent_[v]) v = parent_[v];
  return v;
 }
 [[gnu::always_inline]] T potential(int v) const {
  T x = T();
  for (; ~v; v = parent_[v]) x += lazy_[v];
  return x;
 }

 public:
 explicit FoldableAntiMonopolyTree(const std::vector<T>& a)
     : n(int(a.size())), parent_(n, -1), size_(n, 1), weight_(n, INT_MAX), sum_(a), lazy_(n, T()) {}
 int size() const { return n; }
 // 経路 u-v の最大の重みと、その節点。非連結なら {INT_MAX, -1}。
 std::pair<int, int> path_max(int u, int v) {
  upward_maintain(u), upward_maintain(v);
  int max_w = INT_MIN, t = -1;
  while (u != v) {
   if (size_[u] > size_[v]) std::swap(u, v);
   if (!~parent_[u]) return {INT_MAX, -1};
   if (weight_[u] > max_w) max_w = weight_[u], t = u;
   u = parent_[u];
  }
  return {max_w, t};
 }
 // v の親辺を切る。
 void cut(int v) {
  int len = size_[v];
  T x = T();
  for (int p = parent_[v]; ~p; p = parent_[p]) x += lazy_[p];
  int p = parent_[v];
  if (~p) {
   int dsize = -len;
   T dsum = T() - (sum_[v] + lazy_[p] * T(len));
   for (;;) {
    size_[p] += dsize;
    sum_[p] += dsum;
    int q = parent_[p];
    if (!~q) break;
    dsum += lazy_[q] * T(dsize);
    p = q;
   }
  }
  lazy_[v] += x;
  sum_[v] += x * T(len);
  parent_[v] = -1;
  weight_[v] = INT_MAX;
 }
 // 辺 (u, v) を重み w で入れる。閉路なら経路の最大重みの辺を外す。連結したら真。
 bool link(int u, int v, int w) {
  if (u == v) return false;
  auto [max_w, t] = path_max(u, v);
  bool merged = max_w == INT_MAX;
  if (!merged) {
   if (w >= max_w) return false;
   cut(t);
  }
  T pu = potential(u), pv = potential(v);
  int du = 0, dv = 0;
  T su = T(), sv = T();
  for (;;) {
   if (w >= weight_[u]) {
    int p = parent_[u];
    su += lazy_[p] * T(du);
    size_[p] += du;
    sum_[p] += su;
    pu -= lazy_[u];
    u = p;
   } else if (w >= weight_[v]) {
    int p = parent_[v];
    sv += lazy_[p] * T(dv);
    size_[p] += dv;
    sum_[p] += sv;
    pv -= lazy_[v];
    v = p;
   } else {
    if (size_[u] > size_[v]) {
     std::swap(u, v), std::swap(du, dv), std::swap(su, sv), std::swap(pu, pv);
    }
    int len = size_[u];
    T value = sum_[u];
    T old_lazy = lazy_[u];
    int old_parent = parent_[u];
    int old_weight = weight_[u];
    if (~old_parent) {
     T next_sum = su + lazy_[old_parent] * T(du);
     du -= len;
     su = next_sum - value - lazy_[old_parent] * T(len);
     size_[old_parent] += du;
     sum_[old_parent] += su;
    }
    T delta = pu - old_lazy - pv;
    lazy_[u] += delta;
    sum_[u] += delta * T(len);
    T value_v = sum_[u] + lazy_[v] * T(len);
    dv += len;
    sv += value_v;
    size_[v] += len;
    sum_[v] += value_v;
    parent_[u] = v;
    weight_[u] = w;
    w = old_weight;
    if (!~old_parent) {
     for (; ~parent_[v]; v = parent_[v]) {
      int p = parent_[v];
      sv += lazy_[p] * T(dv);
      size_[p] += dv;
      sum_[p] += sv;
     }
     return merged;
    }
    pu -= old_lazy;
    u = old_parent;
   }
  }
 }
 // 辺 (u, v) の重み w の辺を外す。経路の最大重みの辺がそれより大きいときだけ切る。
 bool erase(int u, int v, int w) {
  auto [max_w, t] = path_max(u, v);
  if (!~t || w >= max_w) return false;
  cut(t);
  return true;
 }
 void add_vertex(int v, const T& x) {
  upward_maintain(v);
  for (; ~v; v = parent_[v]) sum_[v] += x;
 }
 // v の成分全体に x を足す。
 void add_all(int v, const T& x) {
  upward_maintain(v);
  int r = root(v);
  lazy_[r] += x;
  sum_[r] += x * T(size_[r]);
 }
 // v の成分の和。
 T fold(int v) {
  upward_maintain(v);
  return sum_[root(v)];
 }
};

// イベント列。kind 0 = link、1 = cut、2 = 成分全体への加算、3 = 成分の和。
template <class T>
class OfflineSubtree {
 struct Event {
  int kind, u, v, w;
  T x;
 };
 int n;
 std::unordered_map<unsigned long long, int> active;
 std::vector<Event> ev;
 static unsigned long long key(int u, int v) {
  if (u > v) std::swap(u, v);
  return (unsigned long long)(unsigned)(u + 1) << 32 | (unsigned)(v + 1);
 }

 public:
 explicit OfflineSubtree(int n_): n(n_) { ev.reserve(1024); }
 void link(int u, int v) {
  active[key(u, v)] = int(ev.size());
  ev.push_back({0, u, v, INT_MIN, T()});
 }
 void cut(int u, int v) {
  auto it = active.find(key(u, v));
  assert(it != active.end());
  ev[it->second].w = -int(ev.size());
  active.erase(it);
  ev.push_back({1, u, v, 1, T()});
 }
 void add_all(int v, T x) { ev.push_back({2, v, 0, 0, std::move(x)}); }
 void fold(int v) { ev.push_back({3, v, 0, 0, T()}); }
 // a は頂点ごとの初期値。積んだ操作を再生し、fold の答えを順に返す。
 std::vector<T> solve(const std::vector<T>& a) {
  FoldableAntiMonopolyTree<T> amt(a);
  std::vector<T> ans;
  for (int i = 0; i < int(ev.size()); ++i) {
   const Event& e = ev[i];
   if (e.kind == 0) amt.link(e.u, e.v, e.w);
   else if (e.kind == 1) amt.erase(e.u, e.v, -(i + 1));
   else if (e.kind == 2) amt.add_all(e.u, e.x);
   else ans.push_back(amt.fold(e.u));
  }
  return ans;
 }
};
