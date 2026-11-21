#pragma once
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <new>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 関数グラフ (頂点 x から to[x] へ辺が 1 本ずつ出る有向グラフ)。x, f(x), f(f(x)), ... は尾をたどってから閉路を回る。
//   FunctionalGraph(to): to[x] = f(x) (0 <= to[x] < n)。
//   jump(x, k): f^k(x)。
//   path(x, k, f): x, f(x), ..., f^k(x) の k + 1 個を、位置の半開区間 [l, r) と回数 c に分けて f(l, r, c) を呼ぶ。区間は道の順で、
//                  閉路を何周もする分は、閉路の区間 1 つに周回の数を付けて渡す。
//   to_seq(x)、to_vertex(i): 頂点と位置の対応。
//   dist_to_cycle(x): x から閉路に入るまでの辺の数 (閉路の上なら 0)。cycle_length(x): x が入る閉路の長さ。
// 位置は、閉路ごとに閉路の頂点を f の順に 1 つの区間へ並べ、閉路に入る木の頂点は、f(x) を親とした木の HL 分解で、鎖を上る
// (f をたどる) と位置が 1 ずつ増えるように並べる (重い子を先にした行きがけ順を逆にした順)。jump と path は O(log N)。
// 作り方は、入次数が 0 の頂点を葉から剥がして木の頂点を子が親より前に来る順に並べ (残りが閉路)、その順で部分木の大きさと
// 重い子を決め、逆順で位置を配る。今の Library の Period の、関数を int の配列で渡す場合を写したもの。
namespace functional_graph_internal {
// 1 MB 以上の確保は 2 MB 境界に寄せ、Linux では huge page を頼む。resize で値を埋めない。
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
 template <class U> void construct(U* p) noexcept { ::new(static_cast<void*>(p)) U; }
 template <class U, class... A> void construct(U* p, A&&... a) { ::new(static_cast<void*>(p)) U(std::forward<A>(a)...); }
};
template <class T> using Vec= std::vector<T, Alloc<T>>;
}  // namespace functional_graph_internal
class FunctionalGraph {
 int n= 0;
 functional_graph_internal::Vec<int> buf;
 // 位置 i ごとに、inf_[2i] と inf_[2i + 1] は、木の頂点なら (鎖の先頭の位置, 先頭の f の位置)、閉路の頂点なら (閉路の始まりの
 // 位置, 閉路の長さ)。dep_[i] は閉路までの辺の数、ent_[i] は入る閉路の頂点の位置、vert_[i] は頂点。頂点 x ごとに pos_[x]。
 int *inf_= nullptr, *dep_= nullptr, *ent_= nullptr, *vert_= nullptr, *pos_= nullptr;
 void bind() {
  int* p= buf.empty() ? nullptr : buf.data();
  inf_= p, dep_= p ? p + 2 * std::size_t(n) : p, ent_= p ? dep_ + n : p, vert_= p ? ent_ + n : p, pos_= p ? vert_ + n : p;
 }
public:
 FunctionalGraph()= default;
 FunctionalGraph(const FunctionalGraph& o): n(o.n), buf(o.buf) { bind(); }
 FunctionalGraph(FunctionalGraph&& o) noexcept: n(o.n), buf(std::move(o.buf)) { bind(), o.n= 0, o.bind(); }
 FunctionalGraph& operator=(const FunctionalGraph& o) {
  if(this != &o) n= o.n, buf= o.buf, bind();
  return *this;
 }
 FunctionalGraph& operator=(FunctionalGraph&& o) noexcept {
  if(this != &o) n= o.n, buf= std::move(o.buf), bind(), o.n= 0, o.bind();
  return *this;
 }
 explicit FunctionalGraph(const std::vector<int>& to): n(int(to.size())) {
  const int N= n;  // 確保のあとはメンバの n が書き換わりうると gcc が見るので、長さは写しで持つ
  buf.resize(6 * std::size_t(N)), bind();
  if(N <= 0) return;
  functional_graph_internal::Vec<int> tmp(5 * std::size_t(N));
  int *cnt= tmp.data(), *ord= cnt + N, *sz= ord + N, *hb= sz + N;  // hb[2v] = 子の部分木の大きさの最大、hb[2v + 1] = その子
  std::fill(cnt, cnt + N, 0);
  for(int x= 0; x < N; ++x) assert(0 <= to[x] && to[x] < N), ++cnt[to[x]];
  // 入次数が 0 の頂点を剥がす。剥がした頂点は cnt を -1 にし、残った頂点 (cnt > 0) が閉路。
  int k= 0;
  for(int v= 0; v < N; ++v)
   for(int u= v; cnt[u] == 0;) {
    ord[k++]= u, cnt[u]= -1, u= to[u];
    if(--cnt[u] != 0) break;
   }
  std::fill(sz, sz + N, 1), std::fill(hb, hb + 2 * std::size_t(N), 0);
  for(int i= 0; i < k; ++i) {
   const int v= ord[i], p= to[v], s= sz[v];
   sz[p]+= s;
   // 条件式で選ぶと gcc が分岐に戻すので、全ビットのマスクで選ぶ。
   const int c= -int(s > hb[2 * p]);
   hb[2 * p]= (s & c) | (hb[2 * p] & ~c), hb[2 * p + 1]= (v & c) | (hb[2 * p + 1] & ~c);
  }
  // 閉路を f の順に位置 [0, C) へ並べる。並べた頂点は cnt を -2 にする。
  int nc= 0;
  for(int v= 0; v < N; ++v) {
   if(cnt[v] <= 0) continue;
   const int cs= nc;
   int u= v;
   do pos_[u]= nc, vert_[nc]= u, dep_[nc]= 0, ent_[nc]= nc, cnt[u]= -2, ++nc, u= to[u];
   while(u != v);
   for(int i= cs; i < nc; ++i) inf_[2 * i]= cs, inf_[2 * i + 1]= nc - cs;
  }
  // 閉路の頂点 z ごとに、z に入る木の頂点 (sz[z] - 1 個) の区画を [C, N) に取る。sz[v] は v の位置を配ったあと、v の軽い子の区画の
  // 終わり (この位置の手前に置く) に書き換えて使い回す。
  for(int i= 0, t= nc; i < nc; ++i) {
   const int z= vert_[i];
   t+= sz[z] - 1, sz[z]= t;
  }
  for(int i= k; i--;) {
   const int v= ord[i], p= to[v], pp= pos_[p];
   int q;
   if(cnt[p] == -2 || hb[2 * p + 1] != v) {  // 軽い子 (閉路の頂点の子はどれも軽い): 鎖の先頭になる
    q= sz[p] - 1, sz[p]-= sz[v];
    inf_[2 * q]= q, inf_[2 * q + 1]= pp;
   } else {  // 重い子: 親のすぐ下に置き、親と同じ鎖
    q= pp - 1;
    inf_[2 * q]= inf_[2 * pp], inf_[2 * q + 1]= inf_[2 * pp + 1];
   }
   pos_[v]= q, vert_[q]= v, dep_[q]= dep_[pp] + 1, ent_[q]= dep_[pp] ? ent_[pp] : pp;
   sz[v]= q - hb[2 * v];
  }
 }
 int size() const { return n; }
 int to_seq(int x) const { return pos_[x]; }
 int to_vertex(int i) const { return vert_[i]; }
 int dist_to_cycle(int x) const { return dep_[pos_[x]]; }
 int cycle_length(int x) const { return inf_[2 * ent_[pos_[x]] + 1]; }
 // f^k(x) (k >= 0)。
 int jump(int x, long long k) const {
  assert(k >= 0);
  int p= pos_[x];
  const int d= dep_[p];
  if(k >= d) {
   const int q= ent_[p], cs= inf_[2 * q], len= inf_[2 * q + 1];
   return vert_[cs + int((q - cs + (k - d)) % len)];
  }
  // k < d なので閉路には届かない。鎖を上がるたびに組を 1 つ読む。
  for(;;) {
   const int top= inf_[2 * p];
   if(k <= top - p) return vert_[p + int(k)];
   k-= top - p + 1, p= inf_[2 * p + 1];
  }
 }
 // x, f(x), ..., f^k(x) (k >= 0) を、道の順に位置の半開区間 [l, r) と回数 c に分けて f(l, r, c) を呼ぶ。
 template <class F> void path(int x, long long k, const F& f) const {
  assert(k >= 0);
  int p= pos_[x];
  long long m= k + 1;  // まだ渡していない頂点の数
  while(dep_[p]) {
   const int top= inf_[2 * p];
   const long long len= top - p + 1;
   if(m <= len) {
    f(p, p + int(m), 1LL);
    return;
   }
   f(p, top + 1, 1LL), m-= len, p= inf_[2 * p + 1];
  }
  const int cs= inf_[2 * p], len= inf_[2 * p + 1], rest= cs + len - p;
  if(m <= rest) {
   f(p, p + int(m), 1LL);
   return;
  }
  f(p, cs + len, 1LL), m-= rest;
  if(m >= len) f(cs, cs + len, m / len);
  if(m % len) f(cs, cs + int(m % len), 1LL);
 }
};
