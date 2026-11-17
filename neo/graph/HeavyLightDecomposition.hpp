#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <new>
#include <utility>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 根付き木 (森でもよい) の HL 分解。頂点を、重い子を先にした行きがけ順の位置に並べ、道を O(log N) 個の位置の区間に、部分木を
// 1 つの区間に分ける。lca、k 個上の祖先、道の上の k 番目の頂点、深さ、距離も引ける。関数の名前は今の Library と同じ。
//   HeavyLightDecomposition(n, edges, root = 0): 無向の辺の列から。root を含む木は root を根にし、ほかの木の根は決めない。
//   HeavyLightDecomposition(par): 親の列から (根は par[v] < 0)。
// 作り方は、子が親より前に来る順を作り、その順で部分木の大きさと重い子を決め、逆順で位置を配る。無向の辺からは隣接を組まず、
// 頂点ごとに次数と隣の番号の xor を持って、次数 1 の頂点を葉から剥がす (剥がすときの xor が親)。親の列からは、p_i < i (根が 0)
// なら番号の大きい順をそのまま使い、そうでなければ子の CSR を組んで根から BFS する。
// 位置ごとに、鎖の先頭の位置と先頭の親の位置 (根の鎖は -1) を 8 byte の組で持ち、lca と道の分け方は鎖を 1 本上がるたびに組を
// 1 つ読む。頂点ごとに位置と部分木の終わりを並べて持つ。procon-judge の yosupo-lca などで書き比べた形 (記録は algo-notes の notes/tree-queries.md)。
namespace hld_internal {
using u64= unsigned long long;
// 1 MB 以上の確保は 2 MB 境界に寄せ、Linux では huge page を頼む。
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
 // resize で値を埋めない (すぐ上書きするので、0 で埋める分と初めて触る分を省く)。
 template <class U> void construct(U* p) noexcept { ::new(static_cast<void*>(p)) U; }
 template <class U, class... A> void construct(U* p, A&&... a) { ::new(static_cast<void*>(p)) U(std::forward<A>(a)...); }
};
template <class T> using Vec= std::vector<T, Alloc<T>>;
inline u64 chain(int head, int up) { return u64(unsigned(head)) | u64(unsigned(up)) << 32; }
inline int head_of(u64 c) { return int(unsigned(c)); }
inline int up_of(u64 c) { return int(c >> 32); }
}  // namespace hld_internal
class LCA;
template <class M> class PathProduct;
class HeavyLightDecomposition {
 friend class LCA;
 template <class M> friend class PathProduct;
 using u64= hld_internal::u64;
 int n= 0;
 // 配列はすべて buf の中に置く。同じ添字で引く配列の先頭が 2 MB 境界でそろわないよう、配列の間を 2 KB と 64 byte 空ける。
 hld_internal::Vec<u64> buf;
 u64* ch_= nullptr;   // 位置 → 鎖の先頭の位置 (下位 32 bit) と先頭の親の位置 (上位 32 bit、根の鎖は -1)
 int* pe_= nullptr;   // 頂点 v → pe_[2v] が位置、pe_[2v + 1] が部分木の終わり (半開区間)。部分木の問い合わせで 1 回だけ読む
 int* vert_= nullptr;  // 位置 → 頂点
 int* par_= nullptr;   // 頂点 → 親 (根は -1)
 int* dep_= nullptr;   // 頂点 → 深さ
 static constexpr std::size_t GAP= (2048 + 64) / 8;  // u64 の個数
 static std::size_t words(int n) { return 2 * std::size_t(n) + 3 * ((std::size_t(n) + 1) / 2 + GAP) + 2 * GAP; }
 void bind() {
  if(buf.empty()) {
   ch_= nullptr, pe_= vert_= par_= dep_= nullptr;
   return;
  }
  u64* p= buf.data();
  ch_= p, p+= std::size_t(n) + GAP;
  pe_= reinterpret_cast<int*>(p), p+= std::size_t(n) + GAP;
  int** a[3]= {&vert_, &par_, &dep_};
  for(int** q: a) *q= reinterpret_cast<int*>(p), p+= (std::size_t(n) + 1) / 2 + GAP;
 }
 void alloc(int n_) { n= n_, buf.resize(words(n)), bind(); }
 // ord は子が親より前に来る順で、末尾の nr 個が根 (ord[n - 1] の根が位置 0)。par_ は書いてあるとする。
 // work は 3n + 1 個の int の作業領域。
 void layout(const int* ord, int nr, int* work) {
  int *sz= work, *hb= work + n;  // hb[2v] = 子の部分木の大きさの最大、hb[2v + 1] = その子
  std::fill(sz, sz + n, 1), std::fill(hb, hb + 2 * std::size_t(n), 0);
  for(int i= 0; i < n - nr; ++i) {
   const int v= ord[i], p= par_[v], s= sz[v];
   sz[p]+= s;
   // 条件式で選ぶと gcc が分岐に戻すので、全ビットのマスクで選ぶ。
   const int c= -int(s > hb[2 * p]);
   hb[2 * p]= (s & c) | (hb[2 * p] & ~c), hb[2 * p + 1]= (v & c) | (hb[2 * p + 1] & ~c);
  }
  // 上から下へ位置を配る。sz[v] は v の位置を配ったあと、v の子を置く次の位置に書き換えて使い回す。
  int top= 0;
  for(int i= n; i-- > n - nr;) {
   const int v= ord[i], x= top;
   top+= sz[v];
   ch_[x]= hld_internal::chain(x, -1), pe_[2 * v]= x, pe_[2 * v + 1]= top, vert_[x]= v, dep_[v]= 0;
   sz[v]= x + 1 + hb[2 * v];
  }
  for(int i= n - nr; i--;) {
   const int v= ord[i], p= par_[v], pp= pe_[2 * p], nx= sz[p], s= sz[v];
   const int h= -int(hb[2 * p + 1] == v);  // v は p の子なので、hb[2p + 1] は p の本当の重い子
   const int x= ((pp + 1) & h) | (nx & ~h);
   sz[p]= nx + (s & ~h);
   const u64 hm= 0 - u64(h & 1), cp= ch_[pp];
   ch_[x]= (cp & hm) | (hld_internal::chain(x, pp) & ~hm);
   pe_[2 * v]= x, pe_[2 * v + 1]= x + s, vert_[x]= v, dep_[v]= dep_[p] + 1;
   sz[v]= x + 1 + hb[2 * v];
  }
 }
public:
 HeavyLightDecomposition()= default;
 HeavyLightDecomposition(const HeavyLightDecomposition& o): n(o.n), buf(o.buf) { bind(); }
 HeavyLightDecomposition(HeavyLightDecomposition&& o) noexcept: n(o.n), buf(std::move(o.buf)) { bind(), o.n= 0, o.bind(); }
 HeavyLightDecomposition& operator=(const HeavyLightDecomposition& o) {
  if(this != &o) n= o.n, buf= o.buf, bind();
  return *this;
 }
 HeavyLightDecomposition& operator=(HeavyLightDecomposition&& o) noexcept {
  if(this != &o) n= o.n, buf= std::move(o.buf), bind(), o.n= 0, o.bind();
  return *this;
 }
 // 無向の辺の列から。森でもよい (辺の数は n - 木の数)。root を含む木は root を根にして位置 0 から置く。
 HeavyLightDecomposition(int n_, const std::vector<std::array<int, 2>>& edges, int root= 0) {
  alloc(n_);
  if(n <= 0) return;
  assert(0 <= root && root < n);
  hld_internal::Vec<int> tmp(4 * std::size_t(n) + 1);
  int *ord= tmp.data(), *work= ord + n, *dx= work + n;  // dx[2v] = 次数、dx[2v + 1] = 隣の番号の xor (work の後ろ 2n を使う)
  std::fill(dx, dx + 2 * std::size_t(n), 0);
  for(const auto& e: edges) {
   const int a= e[0], b= e[1];
   ++dx[2 * a], dx[2 * a + 1]^= b, ++dx[2 * b], dx[2 * b + 1]^= a;
  }
  dx[2 * root]+= 2;  // root は剥がさない
  std::fill(par_, par_ + n, -1);
  int k= 0;
  for(int v= 0; v < n; ++v)
   for(int u= v; dx[2 * u] == 1;) {
    const int p= dx[2 * u + 1];
    ord[k++]= u, par_[u]= p, dx[2 * u]= 0;
    dx[2 * p + 1]^= u, --dx[2 * p];
    u= p;
   }
  // 剥がれずに残った頂点が根。root を最後に置いて位置 0 にする。
  const int nr= n - k;
  assert(int(edges.size()) == k);  // 閉路があると剥がしきれない
  for(int v= n; v--;)
   if(par_[v] < 0 && v != root) ord[k++]= v;
  ord[k]= root;
  layout(ord, nr, work);
 }
 // 親の列から (根は par[v] < 0)。根が複数なら番号の小さい根の木から順に位置を配る。
 explicit HeavyLightDecomposition(const std::vector<int>& par) {
  alloc(int(par.size()));
  if(n <= 0) return;
  hld_internal::Vec<int> tmp(5 * std::size_t(n) + 2);
  int *ord= tmp.data(), *work= ord + n;
  std::copy(par.begin(), par.end(), par_);
  // p_i < i (根が 0) なら、番号の大きい順がそのまま子が親より前に来る順。分岐の無い and の積み上げで確かめる。
  bool ok= par[0] < 0;
  for(int i= 1; i < n; ++i) ok&= unsigned(par[i]) < unsigned(i);
  int nr= 1;
  if(ok) {
   for(int i= 0; i < n; ++i) ord[i]= n - 1 - i;
  } else {
   // 子の CSR を組み、根から BFS した順を ord の後ろから前へ置く (根は番号の小さい順に ord[n - 1] から)。
   int *st= work + 3 * std::size_t(n) + 1, *cs= work;  // st は n + 1 個 (tmp の末尾)、cs は n 個
   std::fill(st, st + n + 1, 0);
   for(int i= 0; i < n; ++i)
    if(par[i] >= 0) ++st[par[i] + 1];
   for(int i= 0; i < n; ++i) st[i + 1]+= st[i];
   for(int i= 0; i < n; ++i)
    if(par[i] >= 0) cs[st[par[i]]++]= i;
   int t= n;
   for(int v= 0; v < n; ++v)
    if(par[v] < 0) ord[--t]= v;
   nr= n - t;
   for(int i= n - 1; i >= t; --i) {
    const int v= ord[i];
    for(int j= v ? st[v - 1] : 0, e= st[v]; j < e; ++j) ord[--t]= cs[j];
   }
   assert(t == 0);  // 閉路があると根から届かない頂点が残る
  }
  layout(ord, nr, work);
 }
 int to_seq(int v) const { return pe_[2 * v]; }
 int to_vertex(int i) const { return vert_[i]; }
 std::size_t size() const { return n; }
 int parent(int v) const { return par_[v]; }
 // v を含む鎖の先頭の頂点。
 int head(int v) const { return vert_[hld_internal::head_of(ch_[pe_[2 * v]])]; }
 int root(int v) const {
  for(int x= pe_[2 * v];;) {
   const u64 c= ch_[x];
   if(hld_internal::up_of(c) < 0) return vert_[hld_internal::head_of(c)];
   x= hld_internal::up_of(c);
  }
 }
 bool connected(int u, int v) const { return root(u) == root(v); }
 // u が v の部分木にあるか。
 bool in_subtree(int u, int v) const {
  const int x= pe_[2 * u];
  return pe_[2 * v] <= x && x < pe_[2 * v + 1];
 }
 int subtree_size(int v) const { return pe_[2 * v + 1] - pe_[2 * v]; }
 // v の部分木の位置の半開区間 [l, r)。
 std::pair<int, int> subtree(int v) const { return {pe_[2 * v], pe_[2 * v + 1]}; }
 int depth(int v) const { return dep_[v]; }
 // 同じ木の 2 頂点の lca。鎖の先頭の位置が大きい方を上げる。
 int lca(int u, int v) const {
  int x= pe_[2 * u], y= pe_[2 * v];
  u64 a= ch_[x], b= ch_[y];
  while(hld_internal::head_of(a) != hld_internal::head_of(b)) {
   if(hld_internal::head_of(a) > hld_internal::head_of(b)) x= hld_internal::up_of(a), a= ch_[x];
   else y= hld_internal::up_of(b), b= ch_[y];
  }
  return vert_[std::min(x, y)];
 }
 // v の k 個上の祖先 (0 <= k <= depth(v))。
 int la(int v, int k) const {
  assert(0 <= k && k <= dep_[v]);
  for(int x= pe_[2 * v];;) {
   const u64 c= ch_[x];
   const int h= hld_internal::head_of(c);
   if(x - h >= k) return vert_[x - k];
   k-= x - h + 1, x= hld_internal::up_of(c);
  }
 }
 // u から v への道の k 番目 (u が 0 番目) の頂点。道が k より短ければ -1。
 int jump(int u, int v, int k) const {
  if(!k) return u;
  if(u == v) return -1;
  if(k == 1) return in_subtree(v, u) ? la(v, dep_[v] - dep_[u] - 1) : par_[u];
  const int w= lca(u, v), d_uw= dep_[u] - dep_[w], d_vw= dep_[v] - dep_[w];
  return k > d_uw + d_vw ? -1 : k <= d_uw ? la(u, k) : la(v, d_uw + d_vw - k);
 }
 int dist(int u, int v) const { return dep_[u] + dep_[v] - 2 * dep_[lca(u, v)]; }
 // u から v への道を位置の閉区間に分け、道の順に f(l, r) を呼ぶ。u から lca へ上る区間は l > r (位置 l から r へ向かう)、
 // lca から v へ下りる区間は l <= r。edge が真なら lca を除く (辺の値を子の頂点に置く使い方)。
 template <class F> void path(int u, int v, const F& f, bool edge= false) const {
  int x= pe_[2 * u], y= pe_[2 * v];
  u64 a= ch_[x], b= ch_[y];
  int dl[64], dr[64], k= 0;  // v の側の区間は下から見つかるので、積んでおいて最後に逆順で呼ぶ
  while(hld_internal::head_of(a) != hld_internal::head_of(b)) {
   if(hld_internal::head_of(a) > hld_internal::head_of(b)) f(x, hld_internal::head_of(a)), x= hld_internal::up_of(a), a= ch_[x];
   else dl[k]= hld_internal::head_of(b), dr[k++]= y, y= hld_internal::up_of(b), b= ch_[y];
  }
  if(x < y) f(x + edge, y);
  else if(y + edge <= x) f(x, y + edge);
  while(k--) f(dl[k], dr[k]);
 }
};
