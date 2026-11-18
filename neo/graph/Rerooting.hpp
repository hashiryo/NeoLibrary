#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <memory>
#include <new>
#include <utility>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 全方位木 DP。木 (森でもよい) の頂点 v ごとに、v を根にしたときの v の値を O(N) で求める。関数の受け方は今の Library と同じ。
//   Rerooting<T> dp(n, edges, put_edge, op, ui, put_vertex);
//     put_edge(v, e, x) -> U: 辺 e の向こうの部分木の値 x (型 T) を、v から見た寄与にする
//     op(l, r) -> U:          寄与をまとめる。可換で結合的、単位元は ui
//     put_vertex(v, s) -> T:  v に集まった寄与 s から v の値を作る
//   dp[v]: v を根にしたときの v の値。begin() と end() は頂点の順。
//   dp(root, v): root を根にしたときの v の部分木の値 (root と v は同じ木の頂点)。
// 作り方は、隣接の CSR を組み、根を全部先に置いて BFS した順に並べ直す (兄弟の位置が隣り合う)。下からの段は位置の逆順に回って子の
// 寄与をまとめ、上からの段は位置の順に回って、親の側の寄与と兄弟の後ろからの積と前からの積で、子から見た親の側の値を作る。どちらの
// 段も配列を位置の順に読み書きする。dp(root, v) は、初めて呼んだときに内部の根で見た行きがけ順の番号を作り、root を含む v の子を
// 兄弟の中から二分探索で探す。procon-judge の yosupo-tree-path-composite-sum で書き比べた形 (記録は algo-notes の notes/tree-queries.md)。
namespace rerooting_internal {
// 1 MB 以上の確保は 2 MB 境界に寄せ、Linux では huge page を頼む。resize で値を埋めない (どの項目も読む前に書く)。
template <class T> struct Alloc {
 using value_type= T;
 Alloc()= default;
 template <class U> Alloc(const Alloc<U>&) {}
 T* allocate(std::size_t n) {
  const std::size_t want= n * sizeof(T);
  if(want < (std::size_t(1) << 20)) return static_cast<T*>(::operator new(want, std::align_val_t(64)));
  const std::size_t H= std::size_t(1) << 21, bytes= (want + H - 1) & ~(H - 1);
  T* p= static_cast<T*>(::operator new(bytes, std::align_val_t(H)));
#ifdef __linux__
  madvise(p, bytes, MADV_HUGEPAGE);
#endif
  return p;
 }
 void deallocate(T* p, std::size_t n) {
  if(n * sizeof(T) < (std::size_t(1) << 20)) ::operator delete(p, std::align_val_t(64));
  else ::operator delete(p, std::align_val_t(std::size_t(1) << 21));
 }
 template <class U> bool operator==(const Alloc<U>&) const { return true; }
 template <class U> void construct(U* p) noexcept(noexcept(U())) { ::new(static_cast<void*>(p)) U; }
 template <class U, class... A> void construct(U* p, A&&... a) { ::new(static_cast<void*>(p)) U(std::forward<A>(a)...); }
};
template <class T> using Vec= std::vector<T, Alloc<T>>;
// 値の配列。T が bool でも 1 要素 1 個の T を持つよう (std::vector<bool> は詰めて持ち、T& を返せない)、new T[n] で持つ。
template <class T> class Buf {
 std::unique_ptr<T[]> p_;
 std::size_t n_= 0;
public:
 Buf()= default;
 explicit Buf(std::size_t n): p_(n ? new T[n] : nullptr), n_(n) {}
 Buf(const Buf& o): Buf(o.n_) { std::copy(o.p_.get(), o.p_.get() + n_, p_.get()); }
 Buf(Buf&&) noexcept= default;
 Buf& operator=(const Buf& o) {
  if(this != &o) *this= Buf(o);
  return *this;
 }
 Buf& operator=(Buf&&) noexcept= default;
 T* data() { return p_.get(); }
 const T* data() const { return p_.get(); }
 T& operator[](std::size_t i) { return p_[i]; }
 const T& operator[](std::size_t i) const { return p_[i]; }
};
}  // namespace rerooting_internal
template <class T> class Rerooting {
 int n= 0;
 rerooting_internal::Buf<T> dp_, dp1_, dp2_;  // dp_ は頂点の順、dp1_ (内部の根で見た部分木の値) と dp2_ (親の側の値) は位置の順
 rerooting_internal::Vec<int> pos_, cst_;     // 頂点 → 位置、位置 i の子の位置は [cst_[i], cst_[i + 1])
 mutable rerooting_internal::Vec<int> tin_;   // dp(root, v) で使う。位置ごとに、行きがけ順の番号と部分木の大きさ
 void build_tin() const {
  if(!tin_.empty()) return;
  tin_.resize(2 * std::size_t(n));
  int *tin= tin_.data(), *sz= tin + n;
  for(int i= n; i--;) {
   int s= 1;
   for(int j= cst_[i], e= cst_[i + 1]; j < e; ++j) s+= sz[j];
   sz[i]= s;
  }
  for(int i= 0, t= 0; i < cst_[0]; ++i) tin[i]= t, t+= sz[i];  // 根は位置 [0, cst_[0])
  for(int i= 0; i < n; ++i)
   for(int j= cst_[i], e= cst_[i + 1], t= tin[i] + 1; j < e; ++j) tin[j]= t, t+= sz[j];
 }
public:
 Rerooting()= default;
 template <class U, class PE, class OP, class PV> Rerooting(int n_, const std::vector<std::array<int, 2>>& edges, const PE& put_edge, const OP& op, const U& ui, const PV& put_vertex): n(n_) {
  const int N= n_;  // メンバの n は確保のあとで書き換わりうると gcc が見て、負の長さの fill を警告するので、引数の写しを使う
  if(N <= 0) return;
  const int m= int(edges.size());
  assert(m < N);  // 森の辺は N - (木の数) 本
  dp_= rerooting_internal::Buf<T>(N), dp1_= rerooting_internal::Buf<T>(N), dp2_= rerooting_internal::Buf<T>(N);
  pos_.resize(N), cst_.resize(N + 1);
  rerooting_internal::Vec<int> tmp(3 * std::size_t(N) + 1 + 4 * std::size_t(m));
  int *vert= tmp.data(), *pe= vert + N, *st= pe + N, *adj= st + N + 1;  // adj[2j] = 隣、adj[2j + 1] = 辺
  int* cst= cst_.data();
  std::fill(st, st + N + 1, 0);
  for(const auto& e: edges) ++st[e[0] + 1], ++st[e[1] + 1];
  for(int v= 0; v < N; ++v) st[v + 1]+= st[v];
  for(int i= 0; i < m; ++i) {
   const int a= edges[i][0], b= edges[i][1];
   adj[2 * st[a]]= b, adj[2 * st[a] + 1]= i, ++st[a];
   adj[2 * st[b]]= a, adj[2 * st[b] + 1]= i, ++st[b];
  }
  // 並べたあと、v の隣は adj の [v ? st[v - 1] : 0, st[v])。兄弟が隣り合うよう、根を全部先に置いてから BFS する。木なら根は 0 で、
  // 森のときだけ先に連結成分ごとの根を探す (pos_ を見た印に使う)。
  int k= 0;
  if(m + 1 == N) {
   vert[k]= 0, pe[k++]= -1;
  } else {
   int* seen= pos_.data();
   std::fill(seen, seen + N, 0);
   std::vector<int> stk;
   for(int s= 0; s < N; ++s) {
    if(seen[s]) continue;
    vert[k]= s, pe[k++]= -1, seen[s]= 1, stk.assign(1, s);
    while(!stk.empty()) {
     const int v= stk.back();
     stk.pop_back();
     for(int j= v ? st[v - 1] : 0, e= st[v]; j < e; ++j)
      if(!seen[adj[2 * j]]) seen[adj[2 * j]]= 1, stk.push_back(adj[2 * j]);
    }
   }
  }
  for(int i= 0; i < N; ++i) {
   const int v= vert[i], ep= pe[i];
   cst[i]= k;
   for(int j= v ? st[v - 1] : 0, e= st[v]; j < e; ++j)
    if(adj[2 * j + 1] != ep) {
     assert(k < N);  // 閉路があると同じ頂点を何度も置く
     vert[k]= adj[2 * j], pe[k++]= adj[2 * j + 1];
    }
  }
  cst[N]= N;
  for(int i= 0; i < N; ++i) pos_[vert[i]]= i;
  rerooting_internal::Buf<U> g(N), sf(N);  // g[j] は位置 j の子から親への寄与、sf[j] は j より後ろの兄弟の寄与の積
  T *dp1= dp1_.data(), *dp2= dp2_.data();
  for(int i= N; i--;) {
   const int v= vert[i];
   U s= ui;
   for(int j= cst[i], e= cst[i + 1]; j < e; ++j) {
    U x= put_edge(v, pe[j], dp1[j]);
    s= op(s, x), g[j]= std::move(x);
   }
   dp1[i]= put_vertex(v, s);
  }
  for(int i= 0; i < N; ++i) {
   const int v= vert[i], l= cst[i], r= cst[i + 1];
   U pre= pe[i] < 0 ? ui : put_edge(v, pe[i], dp2[i]);
   U s= ui;
   for(int j= r; j-- > l;) sf[j]= s, s= op(g[j], s);
   for(int j= l; j < r; ++j) dp2[j]= put_vertex(v, op(pre, sf[j])), pre= op(pre, g[j]);
   dp_[v]= put_vertex(v, pre);
  }
 }
 std::size_t size() const { return n; }
 const T& operator[](int v) const { return dp_[v]; }
 const T* begin() const { return dp_.data(); }
 const T* end() const { return dp_.data() + n; }
 // root を根にしたときの v の部分木の値。
 const T& operator()(int root, int v) const {
  if(root == v) return dp_[v];
  build_tin();
  const int *tin= tin_.data(), *sz= tin + n;
  const int pr= pos_[root], pv= pos_[v], t= tin[pr];
  if(t < tin[pv] || tin[pv] + sz[pv] <= t) return dp1_[pv];  // root は v の部分木の外
  // root を含む v の子: 兄弟は行きがけ順の番号が増える順に並ぶので、番号が t 以下の最後の子。
  const int c= int(std::upper_bound(tin + cst_[pv], tin + cst_[pv + 1], t) - tin) - 1;
  return dp2_[c];
 }
};
