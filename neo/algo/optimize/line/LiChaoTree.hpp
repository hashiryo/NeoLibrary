#pragma once
#include <algorithm>
#include <concepts>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>
// 関数の集合に関数を足しながら、点 x での最小 (MAX が true なら最大) を引く Li Chao 木。2 つの型を置く。
// - LiChaoTree<T, MAX, F>(lo, hi): x の範囲 [lo, hi) を渡す。節は使う所だけ作る。
// - LiChaoTreePoints<T, MAX, F>(xs): 問いの x を全部先に渡す。座標を詰めた完全二分木を配列に置き、問いは葉から根へ上る。
//   線分が多いときはこちらが速い。
// 載せる関数の型 F は、x (long long) を受けて T を返す operator() を持ち、どの 2 つも高々 1 回しか交わらない族 (擬直線) で
// あること。既定は直線 y = a x + b (LiChaoLine<T, MAX>) で、add_line と add_segment で足せる。ほかの関数は add で足す。
// 値は T に収まること (直線なら a x + b)。関数が 1 つも載っていない点の答えは、最小なら T の最大値、最大なら T の最小値。
// 追加は、節の両端で新しい関数が両方勝てば置き換え、両方負ければ止め、片方だけ勝つなら中点で勝つほうを節に残して、
// 負けたほうを、それが勝つ端のある側の子へ下ろす。交わりが高々 1 回なので、両端で勝てば節の区間の全体で勝つ。
// 関数の載っていない節は、F が「どこでも答えにならない値を返す関数」を static の none() で出せるならそれを置き (直線は
// 傾き 0 で切片が T の最大値か最小値)、出せないなら節ごとの印で見分ける。印を持たないほうが、節が小さく、問いに判定が
// 要らないので速い。
// procon-judge の yosupo-line-add-get-min、yosupo-segment-add-get-min、yuki-2458、yuki-1297、yuki-409、atcoder-dp-z で書き
// 比べた核 (problems/_shared/lines/lichao_fn.hpp) を写したもの。記録は algo-notes の notes/dp-optimization.md にある。
template <class T= long long, bool MAX= false> struct LiChaoLine {
 T a, b;
 [[gnu::always_inline]] T operator()(long long x) const { return a * x + b; }
 static LiChaoLine none() { return {0, MAX ? std::numeric_limits<T>::lowest() : std::numeric_limits<T>::max()}; }
};
namespace li_chao_tree_internal {
template <class F> concept HasNone= requires {
 { F::none() } -> std::convertible_to<F>;
};
// 節に関数が載っているかの印。none() を持つ F では空の型になり、いつも「載っている」と答える。
struct Mark {
 bool v;
 bool get() const { return v; }
 void set() { v= true; }
};
struct NoMark {
 NoMark(bool= true) {}
 static constexpr bool get() { return true; }
 static constexpr void set() {}
};
template <class F> using MarkOf= std::conditional_t<HasNone<F>, NoMark, Mark>;
template <class F> F vacant(const F& f) {
 if constexpr(HasNone<F>) return F::none();
 else return f;
}
template <class T, bool MAX> [[gnu::always_inline]] inline bool wins(const T& u, const T& v) { return MAX ? v < u : u < v; }
template <class T, bool MAX> constexpr T worst= MAX ? std::numeric_limits<T>::lowest() : std::numeric_limits<T>::max();
}  // namespace li_chao_tree_internal
template <class T= long long, bool MAX= false, class F= LiChaoLine<T, MAX>> class LiChaoTree {
 struct Node {
  F f;
  int c[2];
  [[no_unique_address]] li_chao_tree_internal::MarkOf<F> has;
 };
 long long lo_, hi_;
 int depth_= 0, root_= -1;  // 根が覆う幅は 2^depth_
 std::vector<Node> t_;
 static bool wins(const T& u, const T& v) { return li_chao_tree_internal::wins<T, MAX>(u, v); }
 int make(const F& f, bool has) {
  t_.push_back({has ? f : li_chao_tree_internal::vacant(f), {-1, -1}, {has}});
  return (int)t_.size() - 1;
 }
 // 節 k (左端 l、幅 2^d) の下に関数を足す。
 void descend(int k, long long l, int d, F f) {
  for(;;) {
   Node& nd= t_[k];
   if(!nd.has.get()) {
    nd.f= f, nd.has.set();
    return;
   }
   const long long xl= std::min(lo_ + l, hi_ - 1), xr= std::min(lo_ + l + (1LL << d) - 1, hi_ - 1);
   const bool bl= wins(f(xl), nd.f(xl)), br= wins(f(xr), nd.f(xr));
   if(bl == br) {
    if(bl) nd.f= f;
    return;
   }
   const long long m= l + (1LL << (d - 1));
   const bool bm= wins(f(lo_ + m), nd.f(lo_ + m));
   if(bm) std::swap(nd.f, f);
   const int dir= bl == bm;  // 負けた関数が左端で勝つなら左 (0)、そうでなければ右 (1)
   if(dir) l= m;
   --d;
   const int c= nd.c[dir];
   if(c < 0) {
    const int nk= make(f, true);  // t_ が伸びると nd は無効になるので、添字で書き直す
    t_[k].c[dir]= nk;
    return;
   }
   k= c;
  }
 }
 void cover(int k, long long l, int d, long long sl, long long sr, const F& f) {
  const long long r= l + (1LL << d);
  if(sl <= l && r <= sr) return descend(k, l, d, f);
  const long long m= l + (1LL << (d - 1));
  for(int dir= 0; dir < 2; ++dir) {
   const long long cl= dir ? m : l, cr= dir ? r : m;
   if(cr <= sl || sr <= cl) continue;
   int c= t_[k].c[dir];
   if(c < 0) c= make(f, false), t_[k].c[dir]= c;
   cover(c, cl, d - 1, sl, sr, f);
  }
 }
public:
 // x の範囲は [lo, hi)。reserve_nodes は節の数の見込み (全体に足す関数 1 つにつき高々 1 つ)。
 LiChaoTree(long long lo, long long hi, int reserve_nodes= 0): lo_(lo), hi_(hi) {
  while((1LL << depth_) < hi - lo) ++depth_;
  t_.reserve(reserve_nodes);
 }
 // 関数 f を全体に足す。
 void add(const F& f) {
  if(root_ < 0) root_= make(f, true);
  else descend(root_, 0, depth_, f);
 }
 // 関数 f を x in [l, r) にだけ足す。
 void add(const F& f, long long l, long long r) {
  l= std::max(l - lo_, 0LL), r= std::min(r - lo_, 1LL << depth_);
  if(l >= r) return;
  if(root_ < 0) root_= make(f, false);
  cover(root_, 0, depth_, l, r, f);
 }
 // 直線 y = a x + b を足す (F が LiChaoLine のとき)。
 void add_line(T a, T b) { add(F{a, b}); }
 // 線分 y = a x + b (x in [l, r)) を足す (F が LiChaoLine のとき)。
 void add_segment(long long l, long long r, T a, T b) { add(F{a, b}, l, r); }
 // x での最小 (MAX なら最大)。
 [[gnu::always_inline]] T query(long long x) const {
  T res= li_chao_tree_internal::worst<T, MAX>;
  const long long off= x - lo_;
  int k= root_, d= depth_;
  while(k >= 0) {
   const Node& nd= t_[k];
   if(nd.has.get())
    if(const T v= nd.f(x); wins(v, res)) res= v;
   if(d == 0) break;
   --d, k= nd.c[off >> d & 1];
  }
  return res;
 }
};
template <class T= long long, bool MAX= false, class F= LiChaoLine<T, MAX>> class LiChaoTreePoints {
 struct Slot {
  F f;
  [[no_unique_address]] li_chao_tree_internal::MarkOf<F> has;
 };
 int n_= 1, m_= 0;            // 葉の数 (2 の冪) と、重複を除いた点の数
 std::vector<long long> xs_;  // 長さ n_。重複を除いて並べ、末尾を最後の値で埋める
 std::vector<Slot> t_;        // 節 1..2n-1。葉は n..2n-1
 static bool wins(const T& u, const T& v) { return li_chao_tree_internal::wins<T, MAX>(u, v); }
 void init(std::vector<long long>&& pts) {
  const int m= (int)pts.size();
  int n= 1;
  while(n < m) n<<= 1;
  n_= n, m_= m, xs_= std::move(pts);
  xs_.resize(n, m ? xs_.back() : 0);
  t_.assign(2 * (std::size_t)n, Slot{li_chao_tree_internal::vacant(F{}), {false}});
 }
 // 節 k (葉の添字の区間 [l, r)) の下に関数を足す。
 void descend(int k, int l, int r, F f) {
  for(;;) {
   Slot& e= t_[k];
   if(!e.has.get()) {
    e.f= f, e.has.set();
    return;
   }
   const long long xl= xs_[l], xr= xs_[r - 1];
   const bool bl= wins(f(xl), e.f(xl)), br= wins(f(xr), e.f(xr));
   if(bl == br) {
    if(bl) e.f= f;
    return;
   }
   const int mid= (l + r) >> 1;
   const bool bm= wins(f(xs_[mid]), e.f(xs_[mid]));
   if(bm) std::swap(e.f, f);
   if(bl != bm) k= 2 * k, r= mid;
   else k= 2 * k + 1, l= mid;
  }
 }
public:
 struct SortedUnique {};
 // 問いの x を全部渡す (順は問わず、重複があってよい)。
 explicit LiChaoTreePoints(std::vector<long long> xs) {
  std::sort(xs.begin(), xs.end()), xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
  init(std::move(xs));
 }
 // 並べて重複を除いた x を渡す。呼ぶ側が点の添字を持つ使い方のため。
 LiChaoTreePoints(std::vector<long long> xs, SortedUnique) { init(std::move(xs)); }
 // x 以上の最初の点の添字 (並べて重複を除いた列での位置)。
 int index(long long x) const { return (int)(std::lower_bound(xs_.begin(), xs_.begin() + m_, x) - xs_.begin()); }
 // 関数 f を全体に足す。
 void add(const F& f) { descend(1, 0, n_, f); }
 // 関数 f を、点の添字の区間 [l, r) にだけ足す。
 void add_index(const F& f, int l, int r) {
  if(l >= r) return;
  const int n= n_;
  int lk= l + n, rk= r + n, w= 1;
  for(; lk < rk; lk>>= 1, rk>>= 1, w<<= 1) {
   if(lk & 1) descend(lk, lk * w - n, lk * w - n + w, f), ++lk;
   if(rk & 1) --rk, descend(rk, rk * w - n, rk * w - n + w, f);
  }
 }
 // 関数 f を x in [l, r) にだけ足す。
 void add(const F& f, long long l, long long r) { add_index(f, index(l), index(r)); }
 // 直線 y = a x + b を足す (F が LiChaoLine のとき)。
 void add_line(T a, T b) { add(F{a, b}); }
 // 線分 y = a x + b (x in [l, r)) を足す (F が LiChaoLine のとき)。
 void add_segment(long long l, long long r, T a, T b) { add(F{a, b}, l, r); }
 // 添字 i の点での最小 (MAX なら最大)。
 [[gnu::always_inline]] T query_index(int i) const {
  const long long x= xs_[i];
  T res= li_chao_tree_internal::worst<T, MAX>;
  for(int k= i + n_; k >= 1; k>>= 1)
   if(t_[k].has.get())
    if(const T v= t_[k].f(x); wins(v, res)) res= v;
  return res;
 }
 // 先に渡した点 x での最小 (MAX なら最大)。
 [[gnu::always_inline]] T query(long long x) const { return query_index(index(x)); }
};
