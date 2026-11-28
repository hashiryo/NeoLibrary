#pragma once
#include <algorithm>
#include <cstddef>
#include <vector>
#include "neo/base/ds/range/SegmentTree.hpp"
// 永続配列。1 つの節に M 個の子を持たせ、添字 i を M 進法で下ります (i = 0 は節の値、i > 0 は
// 子 (i - 1) % M を (i - 1) / M で下る)。set は新しい版を返し、根から葉への道の節を写す。版の写しは根 1 個。
// 値の無い節は T() を返す (T は既定構築できること)。節は型ごとに 1 本の置き場に取り、0 番を空の節にする。
// M を大きくすると木が浅くなるが、節の写しが M 個の子の分だけ重くなる。union-find の親の配列には 64 を使う。
namespace persistent_array_internal {
template <class T, int M> struct Node {
 T v;
 int ch[M];
};
}
template <class T, int M= 8> class PersistentArray {
 static_assert(M > 1);
 using Node= persistent_array_internal::Node<T, M>;
 using Pool= std::vector<Node, segment_tree_internal::Alloc<Node>>;
 static constexpr std::size_t MIN_NODES= (std::size_t(1) << 28) / sizeof(Node);
 static inline Pool pool_= [] {
  Pool p;
  Node z{};
  for(int i= 0; i < M; ++i) z.ch[i]= 0;
  p.push_back(z);
  return p;
 }();
 static inline std::size_t gc_at_= std::max<std::size_t>(std::size_t(1) << 18, MIN_NODES);
 int n_= 0, root_= 0;
 static int new_node(Node x) {
  const int i= pool_.size();
  pool_.push_back(x);
  return i;
 }
 void set_nc(int i, const T& x) {
  int q= root_;
  if(q == 0) q= root_= new_node(pool_[0]);
  for(;;) {
   if(i == 0) {
    pool_[q].v= x;
    return;
   }
   const int d= (i - 1) % M;
   int c= pool_[q].ch[d];
   if(c == 0) c= pool_[q].ch[d]= new_node(pool_[0]);
   q= c, i= (i - 1) / M;
  }
 }
 int set_rec(int p, int i, const T& x) const {
  const int q= new_node(pool_[p]);
  if(i == 0) {
   pool_[q].v= x;
   return q;
  }
  const int d= (i - 1) % M;
  const int c= set_rec(pool_[p].ch[d], (i - 1) / M, x);
  pool_[q].ch[d]= c;
  return q;
 }
 static void gc_roots(const std::vector<int*>& roots) {
  const std::size_t before= pool_.size();
  std::vector<char> alive(before, 0);
  std::vector<int> st;
  for(int* r: roots) st.push_back(*r);
  while(!st.empty()) {
   const int c= st.back();
   st.pop_back();
   if(c == 0 || alive[c]) continue;
   alive[c]= 1;
   for(int i= 0; i < M; ++i) st.push_back(pool_[c].ch[i]);
  }
  std::vector<int> map(before, 0);
  int k= 1;
  for(std::size_t i= 1; i < before; ++i)
   if(alive[i]) map[i]= k++;
  for(std::size_t i= 1; i < before; ++i)
   if(alive[i]) {
    Node x= pool_[i];
    for(int j= 0; j < M; ++j) x.ch[j]= map[x.ch[j]];
    pool_[map[i]]= x;
   }
  pool_.resize(k);
  for(int* r: roots) *r= map[*r];
  const std::size_t live= k;
  gc_at_= std::max<std::size_t>(MIN_NODES, (2 * live > before ? 4 : 2) * live);
 }
public:
 PersistentArray()= default;
 explicit PersistentArray(int n, const T& v= T()): n_(n) {
  for(int i= n; i--;) set_nc(i, v);
 }
 explicit PersistentArray(const std::vector<T>& a): n_(a.size()) {
  for(int i= n_; i--;) set_nc(i, a[i]);
 }
 int size() const { return n_; }
 T get(int i) const {
  int p= root_;
  while(i > 0) {
   if(p == 0) return T();
   p= pool_[p].ch[(i - 1) % M], i= (i - 1) / M;
  }
  return p == 0 ? T() : pool_[p].v;
 }
 PersistentArray set(int i, const T& x) const {
  PersistentArray res= *this;
  res.root_= set_rec(root_, i, x);
  return res;
 }
 static void reserve(std::size_t n) {
  if(pool_.capacity() < n + 1) pool_.reserve(std::max(n + 1, 2 * pool_.capacity()));
 }
 static bool pool_empty() { return pool_.size() >= gc_at_; }
 static void gc(std::vector<PersistentArray>& vs) {
  std::vector<int*> roots;
  roots.reserve(vs.size());
  for(auto& v: vs) roots.push_back(&v.root_);
  gc_roots(roots);
 }
 template <class... Ps> static void gc(PersistentArray& v, Ps&... vs) {
  std::vector<int*> roots{&v.root_, &vs.root_...};
  gc_roots(roots);
 }
};
