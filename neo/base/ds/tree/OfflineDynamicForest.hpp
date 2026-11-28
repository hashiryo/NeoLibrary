#pragma once
// オフラインの動的連結性。全クエリを先に積み、辺の削除時刻を重みにした最大全域森を 1 本持って再生する。
// 木は親ポインタと部分木の大きさ・部分木和だけで、splay を使わない。2/3 の大きさ平衡 (anti-monopoly) を保ち、
// 期限の短い辺ほど上に来る。削除は owner[時刻] で切る節を引くだけで、代替辺を探さない。G は可換群
// (using T、static T ti()、static T op(T, T)、static T inv(T)、using commute)。
// 使い方は link / cut / add / query を積んで solve。query は成分の和で、答えは query の順に返る。
// Library Checker の dynamic_graph_vertex_add_component_sum の首位 (nandhagk の 402761) の核を写した。
// 記録は algo-notes の notes/dynamic-tree.md。
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <utility>
#include <vector>
#include "neo/base/core/algebra/concepts.hpp"

namespace offline_dynamic_forest_internal {
template <class G>
concept CommutativeGroup = Commutative<G> && requires(typename G::T a) { G::inv(a); };

// 辺の寿命の表。2 の冪の大きさの開番地法で引く (std::unordered_map より軽い)。
class DeathTable {
 std::size_t mask = 3, shift = 62;
 std::vector<unsigned long long> keys{4, 0};
 std::vector<unsigned> vals{4, 0};
 static constexpr unsigned long long TOMB = ~0ull;
 std::size_t slot(unsigned long long k) const { return std::size_t((k * 0x9E3779B97F4A7C15ull) >> shift); }

 public:
 explicit DeathTable(std::size_t links) {
  std::size_t cap = 4;
  while (cap < links * 2) cap <<= 1, --shift;
  mask = cap - 1;
  keys.assign(cap, 0);
  vals.assign(cap, 0);
 }
 void insert(unsigned long long k, unsigned v) {
  std::size_t h = slot(k), tomb = mask + 1;
  while (keys[h]) {
   if (keys[h] == TOMB && tomb == mask + 1) tomb = h;
   h = (h + 1) & mask;
  }
  if (tomb <= mask) h = tomb;
  keys[h] = k, vals[h] = v;
 }
 unsigned erase(unsigned long long k) {
  std::size_t h = slot(k);
  while (keys[h] != k) {
   assert(keys[h] != 0);
   h = (h + 1) & mask;
  }
  keys[h] = TOMB;
  return vals[h];
 }
};

template <class G>
class ExpiryForest {
 using T = typename G::T;
 struct Node {
  unsigned parent = 0, deadline = 0, size = 1;
  T sum;
 };
 std::vector<Node> nd;
 std::vector<unsigned> owner;
 unsigned horizon;
 void unindex(unsigned x) {
  const unsigned d = nd[x].deadline;
  if (d && d < horizon) {
   assert(owner[d] == x);
   owner[d] = 0;
  }
 }
 void index(unsigned x) {
  const unsigned d = nd[x].deadline;
  if (d && d < horizon) {
   assert(!owner[d]);
   owner[d] = x;
  }
 }
 [[gnu::always_inline]] void rotate(unsigned x) {
  Node& a = nd[x];
  const unsigned y = a.parent;
  Node& b = nd[y];
  b.size -= a.size;
  b.sum = G::op(b.sum, G::inv(a.sum));
  if (a.deadline <= b.deadline) a.parent = b.parent;
  else {
   a.size += b.size;
   a.sum = G::op(a.sum, b.sum);
   unindex(x), unindex(y);
   std::swap(a.deadline, b.deadline);
   a.parent = b.parent;
   b.parent = x;
   index(x), index(y);
  }
 }
 [[gnu::always_inline]] void maintain(unsigned x) {
  while (nd[x].parent) {
   const unsigned y = nd[x].parent;
   if ((unsigned long long)nd[x].size * 3 <= (unsigned long long)nd[y].size * 2) x = y;
   else rotate(x);
  }
 }
 [[gnu::always_inline]] void reroot(unsigned x, unsigned watch) {
  while (nd[x].parent && nd[x].parent != watch) rotate(x);
 }

 public:
 ExpiryForest(const std::vector<T>& values, unsigned horizon_): nd(values.size() + 1), owner(horizon_), horizon(horizon_) {
  nd[0].size = 0;
  nd[0].sum = G::ti();
  for (unsigned i = 0; i < values.size(); ++i) nd[i + 1].sum = values[i];
 }
 void link(unsigned u, unsigned v, unsigned deadline) {
  ++u, ++v;
  assert(u != v && deadline && deadline <= horizon);
  maintain(u), maintain(v);
  reroot(v, 0);
  reroot(u, v);
  if (nd[u].parent == v) {
   if (deadline > nd[u].deadline) {
    unindex(u);
    nd[u].deadline = deadline;
    index(u);
   }
  } else {
   assert(!nd[u].parent);
   if (nd[u].size > nd[v].size) std::swap(u, v);
   nd[u].parent = v;
   nd[u].deadline = deadline;
   index(u);
   nd[v].size += nd[u].size;
   nd[v].sum = G::op(nd[v].sum, nd[u].sum);
  }
 }
 // 期限 deadline の辺を切る。
 void expire(unsigned deadline) {
  assert(deadline < horizon);
  const unsigned x = owner[deadline];
  if (!x) return;
  const unsigned count = nd[x].size;
  const T sum = nd[x].sum;
  for (unsigned p = nd[x].parent; p; p = nd[p].parent) {
   nd[p].size -= count;
   nd[p].sum = G::op(nd[p].sum, G::inv(sum));
  }
  nd[x].parent = 0;
  nd[x].deadline = 0;
  owner[deadline] = 0;
 }
 void add(unsigned v, const T& delta) {
  ++v;
  maintain(v);
  for (; v; v = nd[v].parent) nd[v].sum = G::op(nd[v].sum, delta);
 }
 T fold(unsigned v) {
  ++v;
  maintain(v);
  while (nd[v].parent) v = nd[v].parent;
  return nd[v].sum;
 }
};
}  // namespace offline_dynamic_forest_internal

template <class G>
class OfflineDynamicForest {
 static_assert(offline_dynamic_forest_internal::CommutativeGroup<G>, "OfflineDynamicForest requires a commutative group");
 using T = typename G::T;
 struct Event {
  int kind, u, v;
  T x;
 };
 int n;
 std::vector<T> init;
 std::vector<Event> ev;
 static unsigned long long key(int u, int v) {
  if (u > v) std::swap(u, v);
  return (unsigned long long)(unsigned)(u + 1) << 32 | (unsigned)(v + 1);
 }

 public:
 OfflineDynamicForest(int n_, const std::vector<T>& a): n(n_), init(a) {}
 // 辺 (u, v) を張る。単純グラフで、同じ辺が生きている間は呼ばないこと。
 void link(int u, int v) { ev.push_back({0, u, v, T()}); }
 // 辺 (u, v) を消す。生きている辺を指すこと。
 void cut(int u, int v) { ev.push_back({1, u, v, T()}); }
 // 頂点 v の値を x だけ動かす。
 void add(int v, T x) { ev.push_back({2, v, 0, std::move(x)}); }
 // v の連結成分の和を記録する。
 void query(int v) { ev.push_back({3, v, 0, T()}); }
 // 積んだ操作を先頭から再生し、query の答えを順に返す。
 std::vector<T> solve() {
  const unsigned q = unsigned(ev.size());
  std::vector<unsigned> death(q, q);
  std::size_t links = 0;
  for (const Event& e : ev) links += (e.kind == 0);
  offline_dynamic_forest_internal::DeathTable open(links);
  for (unsigned i = 0; i < q; ++i) {
   const Event& e = ev[i];
   if (e.kind > 1) continue;
   const unsigned long long k = key(e.u, e.v);
   if (e.kind == 0) open.insert(k, i);
   else death[open.erase(k)] = i;
  }
  offline_dynamic_forest_internal::ExpiryForest<G> forest(init, q);
  std::vector<T> ans;
  for (unsigned i = 0; i < q; ++i) {
   const Event& e = ev[i];
   if (e.kind == 0) forest.link(e.u, e.v, death[i]);
   else if (e.kind == 1) forest.expire(i);
   else if (e.kind == 2) forest.add(e.u, e.x);
   else ans.push_back(forest.fold(e.u));
  }
  return ans;
 }
};
