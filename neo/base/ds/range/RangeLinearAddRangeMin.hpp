#pragma once
#include <algorithm>
#include <bit>
#include <cstddef>
#include <limits>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 整数の列 a に、区間 [l, r) の各 i へ b i + c を足すことと、区間の最小を引くことを交互にできる。点 (i, a_i) の下側の凸包を、
// 節ごとに左右の子の凸包を結ぶ橋 (共通の接線の 2 点) で持つ二分木 (Overmars と van Leeuwen の形) に、kinetic segment tree の
// 「最小を取る葉と、それが最小のままでいる幅 (melt)」を重ねる。区間に足すのは、要素 i を傾き i の直線と見て時刻を b 進める
// ことにあたり、作用が節の melt に収まる間は節の最小をずらすだけで済む。収まらなかった節は「最小が古い」とだけ記し、区間の最小
// で要るときに求め直す。求め直しは、まず古い子へ再帰して 2 つの子から取り直し (kinetic segment tree と同じ)、古い節を 32 個
// より多く見たらそれをやめて、橋をたどって 1 本の道で最小の葉へ下りる。道の上の橋の傾きから、その葉が最小のままでいる幅も
// ちょうど決まる。時刻を前後に行き来させると kinetic segment tree は再帰が木の全体に広がるが、ここでは 1 節が
// O(32 + log N) で済む。橋は、区間に足したときは端の道の節に「古い」印を付けるだけにし、たどるときに古い子から探し直す。
// 1 回の操作は償却で O(log^2 N)。procon-judge の yosupo-range-linear-add-range-min で書き比べた kinetic_bridge を写したもの。
namespace range_linear_add_range_min_internal {
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
}  // namespace range_linear_add_range_min_internal
class RangeLinearAddRangeMin {
 using i64= long long;
 using i128= __int128;
 using Flag= range_linear_add_range_min_internal::Flag;
 template <class T> using Vec= range_linear_add_range_min_internal::Vec<T>;
 static constexpr i64 INF= i64(1) << 62;
 static constexpr int BUDGET= 32;  // 子から取り直すときに見る古い節の数の上限
 struct Pt {
  i64 y, x;
 };
 struct Geo {
  Pt p, q;     // 橋。葉なら p = q でその点。節に積んだ作用まで足した値
  i64 lm, lk;  // 子へ配っていない作用 y += lm x + lk
 };
 struct Kin {
  i64 mn, ax;  // 部分木の最小と、それを取る葉の位置
  i64 fw, bw;  // 時刻を前へ fw、後ろへ bw まで進めても ax が最小のまま
 };
 int n= 0, sz= 1, lg= 0;
 Vec<Geo> g;
 Vec<Kin> kc;
 Vec<Flag> stale, ok;  // 橋が古いか、最小が新しいか
 static i128 cross(const Pt& a, const Pt& b, const Pt& c) { return i128(b.x - a.x) * (c.y - a.y) - i128(b.y - a.y) * (c.x - a.x); }
 void all_apply(int k, i64 m, i64 c) {
  Geo& x= g[k];
  x.p.y+= m * x.p.x + c, x.q.y+= m * x.q.x + c;
  Kin& h= kc[k];
  if(k >= sz) {
   h.mn+= m * h.ax + c;
   return;
  }
  x.lm+= m, x.lk+= c;
  if(ok[k] == Flag::on) {
   if(m >= 0 ? m <= h.fw : -m <= h.bw) h.mn+= m * h.ax + c, h.fw= std::min(INF, h.fw - m), h.bw= std::min(INF, h.bw + m);
   else ok[k]= Flag::off;
  }
 }
 void push(int k) {
  if(g[k].lm || g[k].lk) all_apply(2 * k, g[k].lm, g[k].lk), all_apply(2 * k + 1, g[k].lm, g[k].lk), g[k].lm= g[k].lk= 0;
 }
 // 子の最小がどちらも新しく、k の作用を配り終えているときに、2 つの子から取る。
 void pull(int k) {
  const Kin &L= kc[2 * k], &R= kc[2 * k + 1];
  const bool lw= L.mn <= R.mn;
  const Kin &w= lw ? L : R, &o= lw ? R : L;
  Kin& h= kc[k];
  h.mn= w.mn, h.ax= w.ax, h.fw= std::min(L.fw, R.fw), h.bw= std::min(L.bw, R.bw);
  if(w.ax > o.ax) h.fw= std::min(h.fw, (o.mn - w.mn) / (w.ax - o.ax));
  if(o.ax > w.ax) h.bw= std::min(h.bw, (o.mn - w.mn) / (o.ax - w.ax));
  ok[k]= Flag::on;
 }
 // 節 k の橋を探し直す。古い子から先に探し直す。
 void rebuild(int k) {
  push(k);
  int l= 2 * k, r= 2 * k + 1;
  if(stale[l] == Flag::on) rebuild(l);
  if(stale[r] == Flag::on) rebuild(r);
  i64 border= r;  // 右の子の最初の葉の位置
  while(border < sz) border<<= 1;
  border-= sz;
  i64 am= 0, ak= 0, bm= 0, bk= 0;  // 左と右の道で、下りるときに通り過ぎた節の作用の和
  Pt a= g[l].p, b= g[l].q, c= g[r].p, e= g[r].q;
  auto movel= [&](int f) {
   am+= g[l].lm, ak+= g[l].lk, l= 2 * l + f;
   a= g[l].p, b= g[l].q, a.y+= am * a.x + ak, b.y+= am * b.x + ak;
  };
  auto mover= [&](int f) {
   bm+= g[r].lm, bk+= g[r].lk, r= 2 * r + f;
   c= g[r].p, e= g[r].q, c.y+= bm * c.x + bk, e.y+= bm * e.x + bk;
  };
  while(l < sz || r < sz) {
   if(l < sz && cross(a, b, c) < 0) movel(0);
   else if(r < sz && cross(b, c, e) < 0) mover(1);
   else if(l >= sz) mover(0);
   else if(r >= sz) movel(1);
   else {
    // 直線 a b と直線 c e の交点が、右の子の始まりより左にあるか。
    const i128 c1= i128(b.x - a.x) * (e.y - c.y) - i128(b.y - a.y) * (e.x - c.x), c2= i128(b.x - a.x) * (b.y - c.y) - i128(b.y - a.y) * (b.x - c.x);
    if(c1 == 0 && c2 == 0 ? c.x < border : i128(c.x) * c1 + i128(e.x - c.x) * c2 < c1 * border) movel(1);
    else mover(0);
   }
  }
  g[k].p= a, g[k].q= c, stale[k]= Flag::off;
 }
 // 橋をたどって節 k の最小の葉へ下り、その葉が最小のままでいる幅とともに k に記す。
 void descend(int k) {
  if(stale[k] == Flag::on) rebuild(k);
  i64 m= 0, c= 0, fw= INF, bw= INF;
  int i= k;
  while(i < sz) {
   const Geo& x= g[i];
   const i64 W= x.q.x - x.p.x, D= (x.q.y - x.p.y) + m * W;
   // 左へ下りる (D >= 0) なら、時刻を後ろへ D / W まで戻してもそのまま。右なら、前へ (-D - 1) / W まで。
   if(D >= 0) bw= std::min(bw, D / W);
   else fw= std::min(fw, (-D - 1) / W);
   m+= x.lm, c+= x.lk, i= 2 * i + (D < 0);
  }
  kc[k]= {g[i].p.y + m * g[i].p.x + c, i - sz, fw, bw};
  ok[k]= Flag::on;
 }
 // 古い子へ再帰して 2 つの子から取り直す。古い節を budget 個より多く見たらやめて偽を返す (それまでに取り直した節は正しい)。
 bool fix(int k, int& budget) {
  if(ok[k] == Flag::on) return true;
  if(--budget < 0) return false;
  push(k);
  if(!fix(2 * k, budget) || !fix(2 * k + 1, budget)) return false;
  pull(k);
  return true;
 }
 // 節 k の最小を新しくする。k の祖先の作用は配り終えている。
 void ensure(int k) {
  if(ok[k] == Flag::on) return;
  int budget= BUDGET;
  if(!fix(k, budget)) descend(k);
 }
 [[gnu::always_inline]] void push_ends(int l, int r) {
  for(int i= lg; i >= 1; --i) {
   if(((l >> i) << i) != l) push(l >> i);
   if(((r >> i) << i) != r) push((r - 1) >> i);
  }
 }
 // 区間に足したあと、端の道の節 k の橋を古いとし、子の最小がどちらも新しければ取り直す。
 [[gnu::always_inline]] void touch(int k) {
  stale[k]= Flag::on;
  if(ok[2 * k] == Flag::on && ok[2 * k + 1] == Flag::on) pull(k);
  else ok[k]= Flag::off;
 }
public:
 RangeLinearAddRangeMin()= default;
 explicit RangeLinearAddRangeMin(const std::vector<i64>& a): n(int(a.size())) {
  sz= int(std::bit_ceil(unsigned(std::max(n, 1)))), lg= std::countr_zero(unsigned(sz));
  g.assign(2 * sz, {{0, 0}, {0, 0}, 0, 0}), kc.assign(2 * sz, {0, 0, INF, INF});
  stale.assign(2 * sz, Flag::off), ok.assign(2 * sz, Flag::on);
  for(int i= 0; i < sz; ++i) {
   const i64 y= i < n ? a[i] : 0;
   g[sz + i].p= g[sz + i].q= {y, i}, kc[sz + i]= {y, i, INF, INF};
  }
  for(int k= sz - 1; k >= 1; --k) rebuild(k), pull(k);
 }
 int size() const { return n; }
 // i in [l, r) に a[i] += b i + c。
 [[gnu::always_inline]] void add(int l, int r, i64 b, i64 c) {
  if(l >= r) return;
  l+= sz, r+= sz;
  push_ends(l, r);
  for(int x= l, y= r; x < y; x>>= 1, y>>= 1) {
   if(x & 1) all_apply(x++, b, c);
   if(y & 1) all_apply(--y, b, c);
  }
  for(int i= 1; i <= lg; ++i) {
   if(((l >> i) << i) != l) touch(l >> i);
   if(((r >> i) << i) != r) touch((r - 1) >> i);
  }
 }
 // min(a[l], ..., a[r - 1]) (l < r)。
 [[gnu::always_inline]] i64 min(int l, int r) {
  l+= sz, r+= sz;
  push_ends(l, r);
  i64 res= std::numeric_limits<i64>::max();
  for(; l < r; l>>= 1, r>>= 1) {
   if(l & 1) ensure(l), res= std::min(res, kc[l].mn), ++l;
   if(r & 1) --r, ensure(r), res= std::min(res, kc[r].mn);
  }
  return res;
 }
 // a[i]。
 i64 get(int i) const {
  const int p= sz + i;
  i64 y= g[p].p.y;
  for(int k= p >> 1; k >= 1; k>>= 1) y+= g[k].lm * i + g[k].lk;
  return y;
 }
};
