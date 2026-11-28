#pragma once
#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>
#include "neo/algo/optimize/line/LiChaoTree.hpp"
#include "neo/base/ds/range/SegmentTree.hpp"
// 永続 Li Chao 木。LiChaoTree と同じく、x の範囲 [lo, hi) と、関数の型 F (既定は直線 y = a x + b) を渡す。
// add / add_line / add_segment が新しい版を返し、query はその版の最小 (MAX なら最大) を引く。版の写しは根 1 個。
// 追加は、LiChaoTree と同じ descend と cover を、通る節を写しながら行う。全体に足す関数は O(log)、線分は O(log^2) 節を写す。
// 節は型ごとに 1 本の置き場に取り、根は -1 を空とする。pool_empty() が真を返したら、生きている版を渡して gc を呼ぶ。
template <class T= long long, bool MAX= false, class F= LiChaoLine<T, MAX>> class PersistentLiChaoTree {
 struct Node {
  F f;
  int c[2];
  [[no_unique_address]] li_chao_tree_internal::MarkOf<F> has;
 };
 using Pool= std::vector<Node, segment_tree_internal::Alloc<Node>>;
 static inline Pool pool_= [] {
  Pool p;
  return p;
 }();
 static inline std::size_t gc_at_= std::size_t(1) << 22;
 long long lo_= 0, hi_= 0;
 int depth_= 0, root_= -1;
 static bool wins(const T& u, const T& v) { return li_chao_tree_internal::wins<T, MAX>(u, v); }
 static int new_node(Node x) {
  const int i= pool_.size();
  pool_.push_back(std::move(x));
  return i;
 }
 static int make(const F& f, bool has) {
  return new_node(Node{has ? f : li_chao_tree_internal::vacant(f), {-1, -1}, {has}});
 }
 // 節 q は写したもの。q の下へ f を下ろす。道の子は写す。
 void descend_into(int q, long long l, int d, F f) const {
  for(;;) {
   if(!pool_[q].has.get()) {
    pool_[q].f= f, pool_[q].has.set();
    return;
   }
   const long long xl= std::min(lo_ + l, hi_ - 1), xr= std::min(lo_ + l + (1LL << d) - 1, hi_ - 1);
   const bool bl= wins(f(xl), pool_[q].f(xl)), br= wins(f(xr), pool_[q].f(xr));
   if(bl == br) {
    if(bl) pool_[q].f= f;
    return;
   }
   const long long m= l + (1LL << (d - 1));
   const bool bm= wins(f(lo_ + m), pool_[q].f(lo_ + m));
   if(bm) std::swap(pool_[q].f, f);
   const int dir= bl == bm;
   if(dir) l= m;
   --d;
   const int c= pool_[q].c[dir];
   if(c < 0) {
    const int nk= make(f, true);
    pool_[q].c[dir]= nk;
    return;
   }
   const int nc= new_node(pool_[c]);
   pool_[q].c[dir]= nc;
   q= nc;
  }
 }
 // 節 q は写したもの。区間 [sl, sr) へ f を足す。
 void cover_into(int q, long long l, int d, long long sl, long long sr, const F& f) const {
  const long long r= l + (1LL << d);
  if(sl <= l && r <= sr) return descend_into(q, l, d, f);
  const long long m= l + (1LL << (d - 1));
  for(int dir= 0; dir < 2; ++dir) {
   const long long cl= dir ? m : l, cr= dir ? r : m;
   if(cr <= sl || sr <= cl) continue;
   const int c= pool_[q].c[dir];
   int nc= c;
   if(c < 0) nc= make(f, false);
   else nc= new_node(pool_[c]);
   pool_[q].c[dir]= nc;
   cover_into(nc, cl, d - 1, sl, sr, f);
  }
 }
 static void gc_roots(const std::vector<int*>& roots) {
  const std::size_t before= pool_.size();
  std::vector<char> alive(before, 0);
  std::vector<int> st;
  for(int* r: roots) st.push_back(*r);
  while(!st.empty()) {
   const int c= st.back();
   st.pop_back();
   if(c < 0 || alive[c]) continue;
   alive[c]= 1;
   st.push_back(pool_[c].c[0]);
   st.push_back(pool_[c].c[1]);
  }
  std::vector<int> map(before, -1);
  int k= 0;
  for(std::size_t i= 0; i < before; ++i)
   if(alive[i]) map[i]= k++;
  for(std::size_t i= 0; i < before; ++i)
   if(alive[i]) {
    Node x= pool_[i];
    x.c[0]= x.c[0] < 0 ? -1 : map[x.c[0]];
    x.c[1]= x.c[1] < 0 ? -1 : map[x.c[1]];
    pool_[map[i]]= x;
   }
  pool_.resize(k);
  for(int* r: roots) *r= *r < 0 ? -1 : map[*r];
  const std::size_t live= k;
  gc_at_= std::max<std::size_t>(std::size_t(1) << 22, (2 * live > before ? 4 : 2) * live);
 }
public:
 PersistentLiChaoTree()= default;
 PersistentLiChaoTree(long long lo, long long hi): lo_(lo), hi_(hi) {
  while((1LL << depth_) < hi - lo) ++depth_;
 }
 // 関数 f を全体に足した新しい版。
 PersistentLiChaoTree add(const F& f) const {
  PersistentLiChaoTree res= *this;
  if(res.root_ < 0) res.root_= make(f, true);
  else {
   res.root_= new_node(pool_[root_]);
   res.descend_into(res.root_, 0, depth_, f);
  }
  return res;
 }
 // 関数 f を x in [l, r) に足した新しい版。
 PersistentLiChaoTree add(const F& f, long long l, long long r) const {
  PersistentLiChaoTree res= *this;
  long long sl= std::max(l - lo_, 0LL), sr= std::min(r - lo_, 1LL << depth_);
  if(sl >= sr) return res;
  if(res.root_ < 0) res.root_= make(f, false);
  else res.root_= new_node(pool_[root_]);
  res.cover_into(res.root_, 0, depth_, sl, sr, f);
  return res;
 }
 // 直線 y = a x + b を足した新しい版 (F が LiChaoLine のとき)。
 PersistentLiChaoTree add_line(T a, T b) const { return add(F{a, b}); }
 // 線分 y = a x + b (x in [l, r)) を足した新しい版 (F が LiChaoLine のとき)。
 PersistentLiChaoTree add_segment(long long l, long long r, T a, T b) const { return add(F{a, b}, l, r); }
 // x での最小 (MAX なら最大)。
 [[gnu::always_inline]] T query(long long x) const {
  T res= li_chao_tree_internal::worst<T, MAX>;
  const long long off= x - lo_;
  int k= root_, d= depth_;
  while(k >= 0) {
   if(pool_[k].has.get())
    if(const T v= pool_[k].f(x); wins(v, res)) res= v;
   if(d == 0) break;
   --d, k= pool_[k].c[off >> d & 1];
  }
  return res;
 }
 static void reserve(std::size_t n) {
  if(pool_.capacity() < n) pool_.reserve(std::max(n, 2 * pool_.capacity()));
 }
 static bool pool_empty() { return pool_.size() >= gc_at_; }
 static void gc(std::vector<PersistentLiChaoTree>& vs) {
  std::vector<int*> roots;
  roots.reserve(vs.size());
  for(auto& v: vs) roots.push_back(&v.root_);
  gc_roots(roots);
 }
 template <class... Ps> static void gc(PersistentLiChaoTree& v, Ps&... vs) {
  std::vector<int*> roots{&v.root_, &vs.root_...};
  gc_roots(roots);
 }
};
