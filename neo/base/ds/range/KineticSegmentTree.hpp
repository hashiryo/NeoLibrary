#pragma once
#include <algorithm>
#include <bit>
#include <cassert>
#include <cstddef>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 各要素が直線 (x_i, y_i) で、区間 [l, r) の各 i へ a x_i + b を足すことと、区間の最大 (MAX = true) か最小を引くことを
// 交互にできる kinetic segment tree。節は部分木で勝っている直線とその値、値が変わりうるまでに足せる a の幅 (nxt) を持つ。
// 作用が幅に収まる間は節の値をずらすだけで済ませ、収まらない作用は子へ配って降り、2 つの子から取り直す (Segment Tree
// Beats と同じ形)。幅は作用のたびに減るので、a は 0 以上でないと使えない。区間の値は節の値が常に今の値なので、覆う節の
// 値を比べるだけで取れる。maspy の ds/segtree/beats_kinetic.hpp と同じ持ち方で、procon-judge の
// self-kinetic-range-add-range-max / -min で書き比べた beats を写したもの。x, y, nxt, a x + b が 64 bit に収まる範囲で
// 使う。記録は algo-notes の notes/kinetic.md にある。
namespace kinetic_segment_tree_internal {
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
};
template <class T> using Vec= std::vector<T, Alloc<T>>;
// 印。unsigned char の配列だと、書くたびに gcc がほかのメンバと重なりうるとみなして読み直す。
enum class Flag: bool { off, on };
}  // namespace kinetic_segment_tree_internal
template <bool MAX= true> class KineticSegmentTree {
 using i64= long long;
 using i128= __int128_t;
 using Flag= kinetic_segment_tree_internal::Flag;
 template <class T> using Vec= kinetic_segment_tree_internal::Vec<T>;
 static constexpr i64 INF= i64(4) << 60;
 struct Node {
  i64 x= 0, y= 0, nxt= INF;
  int idx= -1;  // -1 は空 (大きさを 2 の冪に揃えたときの余り)
 };
 int n= 0, sz= 1;
 Vec<Node> d;
 Vec<i64> la, lb;  // 子へ配っていない作用 y += la x + lb
 Vec<Flag> has;
 static bool better(const Node& A, const Node& B) {
  if(A.idx < 0) return false;
  if(B.idx < 0) return true;
  if constexpr(MAX) return A.y > B.y || (A.y == B.y && A.idx < B.idx);
  else return A.y < B.y || (A.y == B.y && A.idx < B.idx);
 }
 void pull(int k) {
  const Node& L= d[k << 1];
  const Node& R= d[k << 1 | 1];
  if(L.idx < 0) {
   d[k]= R;
   return;
  }
  if(R.idx < 0) {
   d[k]= L;
   return;
  }
  const bool lw= better(L, R);
  const Node& W= lw ? L : R;
  const Node& Z= lw ? R : L;
  Node M= W;
  M.nxt= std::min(L.nxt, R.nxt);
  i64 num, den;
  if constexpr(MAX) {
   if(W.x >= Z.x) {
    d[k]= M;
    return;
   }
   num= W.y - Z.y, den= Z.x - W.x;  // 敗者 Z の傾きが大きい。Z が勝つ最初の整数時刻まで
  } else {
   if(W.x <= Z.x) {
    d[k]= M;
    return;
   }
   num= Z.y - W.y, den= W.x - Z.x;  // 敗者 Z の傾きが小さい
  }
  const i128 t= (i128)num / den + 1;
  if(t < (i128)M.nxt) M.nxt= (i64)t;
  d[k]= M;
 }
 bool try_apply(int k, i64 a, i64 b) {
  Node& v= d[k];
  if(v.idx < 0) return true;
  if(v.nxt <= a) return false;
  v.y+= a * v.x + b;
  v.nxt-= a;
  if(k < sz) la[k]+= a, lb[k]+= b, has[k]= Flag::on;
  return true;
 }
 void push(int k) {
  if(has[k] == Flag::off) return;
  const i64 a= la[k], b= lb[k];
  has[k]= Flag::off, la[k]= 0, lb[k]= 0;
  const bool p= try_apply(k << 1, a, b), q= try_apply(k << 1 | 1, a, b);
  // 節が a を受け入れたとき a < 節の nxt <= 子の nxt なので、子も受け入れる。
  assert(p && q);
  (void)p, (void)q;
 }
 void add_rec(int k, int nl, int nr, int l, int r, i64 a, i64 b) {
  if(r <= nl || nr <= l) return;
  if(l <= nl && nr <= r && try_apply(k, a, b)) return;
  if(k >= sz) return;  // 葉は nxt = INF なので、ここへは来ない
  push(k);
  const int m= (nl + nr) >> 1;
  add_rec(k << 1, nl, m, l, r, a, b);
  add_rec(k << 1 | 1, m, nr, l, r, a, b);
  pull(k);
 }
 i64 prod_rec(int k, int nl, int nr, int l, int r) {
  if(r <= nl || nr <= l || d[k].idx < 0) return MAX ? -INF : INF;
  if(l <= nl && nr <= r) return d[k].y;
  push(k);  // 端で降りるときは子を今の時刻にそろえる
  const int m= (nl + nr) >> 1;
  const i64 A= prod_rec(k << 1, nl, m, l, r), B= prod_rec(k << 1 | 1, m, nr, l, r);
  if constexpr(MAX) return std::max(A, B);
  else return std::min(A, B);
 }
 // 葉 i の道の節の作用を配り、葉の値を今の時刻にする。道は path[0] が葉、path[t - 1] が根。
 int lift(int i, int* path) {
  int k= sz + i, t= 0;
  path[t++]= k;
  while(k > 1) k>>= 1, path[t++]= k;
  for(int u= t - 1; u >= 0; --u)
   if(path[u] < sz) push(path[u]);
  return t;
 }
public:
 KineticSegmentTree()= default;
 explicit KineticSegmentTree(const std::vector<i64>& x, const std::vector<i64>& y): n(int(x.size())), sz(int(std::bit_ceil(unsigned(std::max(n, 1))))) {
  d.assign(2 * sz, Node{}), la.assign(sz, 0), lb.assign(sz, 0), has.assign(sz, Flag::off);
  for(int i= 0; i < n; ++i) d[sz + i]= {x[i], y[i], INF, i};
  for(int k= sz - 1; k >= 1; --k) pull(k);
 }
 int size() const { return n; }
 // i in [l, r) に y_i += a x_i + b。a は 0 以上。
 [[gnu::always_inline]] void add(int l, int r, i64 a, i64 b) {
  if(l >= r) return;
  add_rec(1, 0, sz, l, r, a, b);
 }
 // (x_i, y_i) を (x, y) にする。y は今の時刻の値。
 [[gnu::always_inline]] void set(int i, i64 x, i64 y) {
  int path[32];
  const int t= lift(i, path);
  d[sz + i]= {x, y, INF, i};
  for(int u= 1; u < t; ++u) pull(path[u]);
 }
 // max(y_l, ..., y_{r - 1})。MAX = false なら最小。
 [[gnu::always_inline]] i64 prod(int l, int r) { return prod_rec(1, 0, sz, l, r); }
 [[gnu::always_inline]] i64 all_prod() { return d[1].idx < 0 ? (MAX ? -INF : INF) : d[1].y; }
 // y_i。
 i64 get(int i) {
  int path[32];
  lift(i, path);
  return d[sz + i].y;
 }
};
