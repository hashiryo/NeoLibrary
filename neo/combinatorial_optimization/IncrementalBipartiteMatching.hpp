#pragma once
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <new>
#include <utility>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 増加路を 1 本ずつ探して大きくしていく二部マッチング。左の頂点 [0, L)、右の頂点 [0, R) と、その間の辺を持つ。頂点の番号は
// BipartiteMatching と同じく、左 a を a、右 b を L + b とする (mate)。マッチングは空で始め、augment を呼んだときだけ 1 本ずつ
// 増やす。左の頂点を足す順で答えが決まる問題 (値を小さい順に足し、増加路が見つかった値を採る貪欲法など) のための型。
// - add_edge(a, b) は辺を足す。erase_edge(a, b) は辺を 1 本抜き、組に使っていて同じ組の辺がほかに無ければ両端を空ける。
//   どちらもマッチングを大きくはしない。
// - augment(a) は空いた左 a から、augment() は空いた左の頂点のどれかから、増加路を幅優先で探し、見つかれば 1 本だけ流す。
// 隣接は 1 本の配列の上に、左の頂点ごとの区画として置く。構築で渡した辺は CSR と同じ並び (区画の大きさは辺の数ちょうど)
// にし、add_edge で区画が溢れたら、2 倍の大きさの区画を配列の末尾に取って写す。erase_edge は区画の中で末尾と入れ替える。
// 増加路が見つからなかった探索で訪ねた左の頂点の集合 D には印を残し、次からの探索で入らない。D の頂点の隣の右の頂点はどれも
// 埋まっていて、その相手も D にあるので、D からは空いた右の頂点へ届かない。D の外を通る増加路は D の頂点の隣を通らない
// ので、流したあとも D の頂点の隣と相手は変わらず、印は消さなくてよい (見つかった探索で訪ねた頂点の印だけを消す)。D の
// 頂点から辺を足したときと、D の頂点が使う辺を抜いたときは、印を全部消す。
// 2 MB 以上の配列は、Linux では 2 MB の区画 (huge page) に置き、先頭を配列ごとにずらす (MaxFlow と同じ allocator)。探索の
// 待ち行列と、頂点を見つけた左の頂点の表は 0 で埋めず、使うところまでしかページに触らない。
// procon-judge の luogu-P1640、cf-1728-f、cf-1139-e、loj-2338、aoj-3198、luogu-P1971 で書き比べた keep2_bfs (中身は
// problems/_shared/kuhn/inc_matching.hpp) を写し、配列の置き方を lib-neo-hp (problems/_shared/kuhn/neo_hp.hpp) で試して
// 入れたもの。記録は algo-notes の notes/network-flow.md にある。
namespace incremental_bipartite_matching_internal {
// 2 MB 以上の確保を、Linux では 2 MB 境界の mmap に置いて MADV_HUGEPAGE を頼み、先頭を確保ごとに 2 KB + 64 byte の倍数だけ
// ずらす allocator (k 回目は k mod 16 倍)。値が 10^6 個のように左の頂点が多いと、4 KB のページに初めて触るたびに OS が
// ページを割り当てる手間が大きい。同じ添字で引く頂点ごとの配列がどれも 2 MB 境界から始まると、同じ添字の要素がキャッシュの
// 同じ組を取り合うので、ずらす。返すときは 2 MB 境界へ切り下げて先頭に戻す。
inline unsigned big_alloc_count= 0;
template <class T> struct BigAlloc {
 using value_type= T;
 static constexpr std::size_t H= std::size_t(1) << 21, STEP= 2048 + 64;
 BigAlloc()= default;
 template <class U> BigAlloc(const BigAlloc<U>&) {}
 T* allocate(std::size_t n) {
#ifdef __linux__
  if(n * sizeof(T) >= H) {
   const std::size_t off= (big_alloc_count++ % 16) * STEP, len= (n * sizeof(T) + off + H - 1) & ~(H - 1);
   void* m= mmap(nullptr, len + H, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
   if(m == MAP_FAILED) throw std::bad_alloc();
   char *p= static_cast<char*>(m), *a= reinterpret_cast<char*>((reinterpret_cast<std::uintptr_t>(p) + H - 1) & ~(H - 1));
   if(a != p) munmap(p, a - p);
   munmap(a + len, p + H - a);
   madvise(a, len, MADV_HUGEPAGE);
   return reinterpret_cast<T*>(a + off);
  }
#endif
  return static_cast<T*>(::operator new(n * sizeof(T)));
 }
 void deallocate(T* p, [[maybe_unused]] std::size_t n) {
#ifdef __linux__
  if(n * sizeof(T) >= H) {
   char* a= reinterpret_cast<char*>(reinterpret_cast<std::uintptr_t>(p) & ~std::uintptr_t(H - 1));
   munmap(a, (n * sizeof(T) + (reinterpret_cast<char*>(p) - a) + H - 1) & ~(H - 1));
   return;
  }
#endif
  ::operator delete(p);
 }
 template <class U> bool operator==(const BigAlloc<U>&) const { return true; }
};
// construct を既定の初期化にして、std::vector<int, NoInit<int>>(n) が 0 で埋めないようにする。
template <class T> struct NoInit: BigAlloc<T> {
 template <class U> struct rebind {
  using other= NoInit<U>;
 };
 NoInit()= default;
 template <class U> NoInit(const NoInit<U>&) {}
 template <class U> void construct(U* p) { ::new(static_cast<void*>(p)) U; }
 template <class U, class... A> void construct(U* p, A&&... a) { ::new(static_cast<void*>(p)) U(std::forward<A>(a)...); }
};
}  // namespace incremental_bipartite_matching_internal
class IncrementalBipartiteMatching {
 template <class T> using Vec= std::vector<T, incremental_bipartite_matching_internal::BigAlloc<T>>;
 template <class T> using Raw= std::vector<T, incremental_bipartite_matching_internal::NoInit<T>>;
public:
 IncrementalBipartiteMatching(int l, int r, const std::vector<std::array<int, 2>>& edges= {}): L(l), beg_(l + 1, 0), len_(l, 0), adj_(edges.size()), ml_(l, -1), mr_(r, -1), pre_(l), q_(l), vis_(l, 0) {
  for(auto& e: edges) ++len_[e[0]];
  for(int i= 0; i < L; ++i) beg_[i + 1]= beg_[i] + len_[i], len_[i]= 0;
  for(auto& e: edges) adj_[beg_[e[0]] + len_[e[0]]++]= e[1];
 }
 // 左 a と右 b を結ぶ辺を足す。マッチングは変えない。
 void add_edge(int a, int b) {
  if(cap_.empty()) cap_.assign(len_.begin(), len_.end());
  if(len_[a] == cap_[a]) grow(a);
  adj_[beg_[a] + len_[a]++]= b;
  if(vis_[a] == stamp_) ++stamp_;  // a からは空いた右へ届かないと覚えていたが、辺が増えたので忘れる
 }
 // 左 a と右 b を結ぶ辺を 1 本抜く (無ければ何もしない)。組に使っていて、同じ組の辺がほかに無ければ、両端が空く。
 void erase_edge(int a, int b) {
  int* p= adj_.data() + beg_[a];
  const int j= std::find(p, p + len_[a], b) - p;
  if(j == len_[a]) return;
  p[j]= p[--len_[a]];
  if(ml_[a] != b || std::find(p, p + len_[a], b) != p + len_[a]) return;
  if(vis_[a] == stamp_) ++stamp_;
  ml_[a]= -1, mr_[b]= -1, --sz_;
 }
 // 空いた左 a から増加路を探し、見つかれば流して true。a が空いていないか、a から届かなければ false。
 [[gnu::always_inline]] bool augment(int a) {
  if(ml_[a] >= 0 || vis_[a] == stamp_) return false;
  vis_[a]= stamp_, pre_[a]= -1, q_[0]= a;
  return bfs(1);
 }
 // 空いた左の頂点のどれかから増加路を探し、見つかれば 1 本だけ流して true。false なら今のマッチングが最大。
 bool augment() {
  int t= 0;
  for(int a= 0; a < L; ++a)
   if(ml_[a] < 0 && vis_[a] != stamp_) vis_[a]= stamp_, pre_[a]= -1, q_[t++]= a;
  return bfs(t);
 }
 // マッチングの辺の数。
 [[gnu::always_inline]] int size() const { return sz_; }
 // 頂点 v (左 a は a、右 b は L + b) の相手。いなければ -1。
 [[gnu::always_inline]] int mate(int v) const {
  if(v < L) return ml_[v] < 0 ? -1 : L + ml_[v];
  return mr_[v - L];
 }
private:
 int L, sz_= 0;
 unsigned stamp_= 1;
 // 左 l の隣は adj_[beg_[l], beg_[l] + len_[l])。cap_ は区画の大きさ (初めて足すときに作る)。
 Vec<int> beg_, len_, adj_;
 std::vector<int> cap_;
 // 相手 (いなければ -1)、探索で l を見つけた左の頂点、探索の待ち行列 (pre_ と q_ は 0 で埋めない)。
 Vec<int> ml_, mr_;
 Raw<int> pre_, q_;
 // vis_[l] == stamp_ なら、訪ねたか、空いた右へ届かないと分かっている。
 Vec<unsigned> vis_;
 void grow(int a) {
  const int nb= adj_.size(), nc= std::max(4, 2 * cap_[a]);
  adj_.resize(nb + nc);
  std::copy(adj_.begin() + beg_[a], adj_.begin() + beg_[a] + len_[a], adj_.begin() + nb);
  beg_[a]= nb, cap_[a]= nc;
 }
 // q_[0, t) に積んだ左の頂点 (印を付けてある) から幅優先で探す。pre_[l] は l を見つけた左の頂点で、そのとき通った右の頂点は
 // 見つけたときの l の相手。
 bool bfs(int t) {
  const int *adj= adj_.data(), *beg= beg_.data(), *len= len_.data();
  int* ml= ml_.data();
  int* mr= mr_.data();
  int* pre= pre_.data();
  int* q= q_.data();
  unsigned* vis= vis_.data();
  const unsigned st= stamp_;
  for(int h= 0; h < t; ++h) {
   const int l= q[h];
   const int *p= adj + beg[l], *e= p + len[l];
   for(; p != e; ++p) {
    const int r= *p, m= mr[r];
    if(m < 0) {
     for(int x= l, y= r;;) {
      const int px= pre[x], old= ml[x];
      ml[x]= y, mr[y]= x;
      if(px < 0) break;
      x= px, y= old;
     }
     ++sz_;
     for(int i= 0; i < t; ++i) vis[q[i]]= 0;
     return true;
    }
    if(vis[m] != st) vis[m]= st, pre[m]= l, q[t++]= m;
   }
  }
  return false;
 }
};
