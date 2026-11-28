#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>
// 区分線形な凸関数 f を持つ。折れ点と、そこで傾きが増える量を treap に載せ、最小を取る折れ点 (mn) の
// 左の部分木と右の部分木に分けて持つ。各部分木は x に足す遅延 z を持つ。ノードはインスタンスごとの
// 配列に置き、0 番を番兵にする。
// 主な操作: add_abs (a|x-x0|)、add_ramp (max(0, a(x-x0)))、add_max (max(a(x-x0), b(x-x0)))、
// add_const、add_linear (a x)、add_inf (定義域の制限)、shift (f(x) <- f(x-x0))、
// chmin_cum (f(x) <- min_{y<=x} f(y)、引数 true で min_{x<=y})、
// chmin_slide_win (f(x) <- min_{lb<=y<=ub} f(x-y))、min、argmin、f(x)、operator+= (凸関数どうしの和)。
// remove 系は subtract_max、subtract_ramp、subtract_abs。
// 今の Library の PiecewiseLinearConvex から、永続 (persistent) とクラス static のプールを外したもの。
// T は 32 bit か 64 bit の符号付き整数 (計算は 64 bit か 128 bit で行う)。
// 記録は algo-notes の notes/dp-optimization.md にある。
template <class T= long long> class PiecewiseLinearConvex {
 using D= std::conditional_t<std::is_integral_v<T>, std::conditional_t<(sizeof(T) < 8), long long, __int128_t>, long double>;
 struct Node {
  int ch[2]= {0, 0};
  T z= 0, x= 0, d= 0, a= 0;
  D s= 0;
  std::size_t sz= 0;
 };
 std::vector<Node> n_= std::vector<Node>(1);  // 0 番は番兵
 static bool lt(T a, T b) {
  if constexpr(std::is_floating_point_v<T>) return 1e-15 < b - a;
  else return a < b;
 }
 template <bool g> static bool lgt(T a, T b) {
  if constexpr(g) return lt(b, a);
  else return lt(a, b);
 }
 template <bool l> static bool lte(T a, T b) {
  if constexpr(l) return lt(a, b);
  else return !lt(b, a);
 }
 int create(T d, T x) { return n_.push_back(Node{{0, 0}, 0, x, d, d, D(x) * d, 1}), (int)n_.size() - 1; }
 void update(int t) {
  const int l= n_[t].ch[0], r= n_[t].ch[1];
  n_[t].sz= 1 + n_[l].sz + n_[r].sz, n_[t].a= n_[t].d + n_[l].a + n_[r].a, n_[t].s= D(n_[t].x) * n_[t].d + n_[l].s + n_[r].s;
 }
 void prop(int& t, T v) {
  if(!t) return;
  n_[t].z+= v, n_[t].s+= D(v) * n_[t].a, n_[t].x+= v;
 }
 void push(int t) {
  if(n_[t].z != 0) prop(n_[t].ch[0], n_[t].z), prop(n_[t].ch[1], n_[t].z), n_[t].z= 0;
 }
 template <bool r> int join_(int t, int a, int b) {
  push(a);
  if constexpr(r) b= join<0>(b, t, n_[a].ch[0]);
  else b= join<0>(n_[a].ch[1], t, b);
  if(n_[n_[a].ch[r]].sz * 4 >= n_[b].sz) return n_[a].ch[!r]= b, update(a), a;
  return n_[a].ch[!r]= n_[b].ch[r], update(a), n_[b].ch[r]= a, update(b), b;
 }
 template <bool b= 1> int join(int l, int t, int r) {
  if constexpr(!b) (void)0;
  if(n_[l].sz > n_[r].sz * 4) return join_<0>(t, l, r);
  if(n_[r].sz > n_[l].sz * 4) return join_<1>(t, r, l);
  return n_[t].ch[0]= l, n_[t].ch[1]= r, update(t), t;
 }
 std::array<int, 3> split(int t, T x) {
  if(!t) return {0, 0, 0};
  push(t);
  if(lt(n_[t].x, x)) {
   auto [a, b, c]= split(n_[t].ch[1], x);
   return {join(n_[t].ch[0], t, a), b, c};
  } else if(lt(x, n_[t].x)) {
   auto [a, b, c]= split(n_[t].ch[0], x);
   return {a, b, join(c, t, n_[t].ch[1])};
  }
  return {n_[t].ch[0], t, n_[t].ch[1]};
 }
 int unite(int l, int r) {
  if(!l) return r;
  if(!r) return l;
  push(l);
  auto [a, b, c]= split(r, n_[l].x);
  return n_[l].d+= n_[b].d, join<0>(unite(a, n_[l].ch[0]), l, unite(n_[l].ch[1], c));
 }
 int insert(int t, T x, T d) {
  if(!t) return create(d, x);
  push(t);
  if(lt(x, n_[t].x)) return join<0>(insert(n_[t].ch[0], x, d), t, n_[t].ch[1]);
  if(lt(n_[t].x, x)) return join<0>(n_[t].ch[0], t, insert(n_[t].ch[1], x, d));
  return n_[t].d+= d, update(t), t;
 }
 int erase(int t, T x, T d) {
  assert(t);
  push(t);
  if(lt(x, n_[t].x)) return join(erase(n_[t].ch[0], x, d), t, n_[t].ch[1]);
  if(lt(n_[t].x, x)) return join(n_[t].ch[0], t, erase(n_[t].ch[1], x, d));
  assert(!lt(n_[t].d, d));
  if(lt(d, n_[t].d)) return n_[t].d-= d, update(t), t;
  if(n_[t].ch[1]) {
   auto [a, s]= pop<0>(n_[t].ch[1]);
   return join(n_[t].ch[0], s, a);
  }
  return n_[t].ch[0];
 }
 template <bool r> std::pair<int, int> pop(int t) {
  if(push(t); !n_[t].ch[r]) return {n_[t].ch[!r], t};
  auto [a, s]= pop<r>(n_[t].ch[r]);
  if constexpr(r) return {join(n_[t].ch[!r], t, a), s};
  else return {join(a, t, n_[t].ch[!r]), s};
 }
 template <bool r> int cut(int t, T x) {
  if(!t) return t;
  if(push(t); lgt<r>(n_[t].x, x)) return cut<r>(n_[t].ch[!r], x);
  if(lgt<r>(x, n_[t].x)) {
   if constexpr(r) return join(n_[t].ch[0], t, cut<1>(n_[t].ch[1], x));
   else return join(cut<0>(n_[t].ch[0], x), t, n_[t].ch[1]);
  }
  return n_[t].ch[!r];
 }
 template <bool r> D calc_y(int t, T x, T ol, D ou) {
  for(; t;) {
   if(push(t); lgt<r>(n_[t].x, x)) t= n_[t].ch[!r];
   else {
    ol+= n_[n_[t].ch[!r]].a, ou+= n_[n_[t].ch[!r]].s;
    if(!lgt<r>(x, n_[t].x)) break;
    ol+= n_[t].d, ou+= D(n_[t].x) * n_[t].d, t= n_[t].ch[r];
   }
  }
  return D(x) * ol - ou;
 }
 template <bool r> std::array<int, 3> split(int t, T p, T& ol, D& ou) {
  push(t);
  T s= ol + n_[n_[t].ch[!r]].a;
  if(lt(p, s)) {
   auto [a, b, c]= split<r>(n_[t].ch[!r], p, ol, ou);
   if constexpr(r) return {a, b, join(c, t, n_[t].ch[r])};
   else return {join(n_[t].ch[r], t, a), b, c};
  }
  ol= s + n_[t].d;
  if(lt(ol, p)) {
   ou+= n_[n_[t].ch[!r]].s + D(n_[t].x) * n_[t].d;
   auto [a, b, c]= split<r>(n_[t].ch[r], p, ol, ou);
   if constexpr(r) return {join(n_[t].ch[!r], t, a), b, c};
   else return {a, b, join(c, t, n_[t].ch[!r])};
  }
  ou+= n_[n_[t].ch[!r]].s;
  return {n_[t].ch[0], t, n_[t].ch[1]};
 }
 template <bool l, bool r> std::pair<int, int> split_cum(int t, T p, T& ol, D& ou) {
  push(t);
  T s= ol + n_[n_[t].ch[!r]].a;
  if(lte<l>(p, s)) {
   auto [c, b]= split_cum<l, r>(n_[t].ch[!r], p, ol, ou);
   if constexpr(l) {
    if constexpr(r) return {join(c, t, n_[t].ch[r]), b};
    else return {join(n_[t].ch[r], t, c), b};
   } else return {c, b};
  }
  ol= s + n_[t].d;
  if(lte<!l>(ol, p)) {
   ou+= n_[n_[t].ch[!r]].s + D(n_[t].x) * n_[t].d;
   auto [a, b]= split_cum<l, r>(n_[t].ch[r], p, ol, ou);
   if constexpr(l) return {a, b};
   else {
    if constexpr(r) return {join(n_[t].ch[!r], t, a), b};
    else return {join(a, t, n_[t].ch[!r]), b};
   }
  }
  ou+= n_[n_[t].ch[!r]].s;
  return {n_[t].ch[!r ^ l], t};
 }
 int mn= 0, lr[2]= {0, 0};
 bool bf[2]= {false, false};
 T o[2]= {0, 0}, rem= 0, bx[2]= {0, 0};
 D y= 0;
 D calc_y(T x) {
  if(!mn) return 0;
  if(lt(x, n_[mn].x)) return -calc_y<0>(lr[0], x, o[0], D(n_[mn].x) * o[0]);
  if(lt(n_[mn].x, x)) return calc_y<1>(lr[1], x, o[1], D(n_[mn].x) * o[1]);
  return 0;
 }
 void slope_eval(bool neg) {
  T p= neg ? -rem : rem, ol= o[neg];
  if(p <= ol) o[neg]-= p, o[!neg]+= p, y+= D(n_[mn].x) * rem;
  else {
   D ou= D(n_[mn].x) * ol;
   auto [a, b, c]= neg ? split<1>(lr[neg], p, ol, ou) : split<0>(lr[neg], p, ol, ou);
   o[neg]= ol - p, ol-= n_[b].d, ou+= D(n_[b].x) * (o[!neg]= p - ol);
   if(neg) y-= ou, lr[!neg]= join(lr[!neg], mn, a), lr[neg]= c;
   else y+= ou, lr[!neg]= join(c, mn, lr[!neg]), lr[neg]= a;
   mn= b;
  }
  rem= 0;
 }
 template <bool l, bool neg> void slope_eval_cum() {
  T p= neg ? -rem : rem, ol= o[neg];
  if(lte<l>(p, ol)) o[neg]-= p, o[!neg]+= p, y+= D(n_[mn].x) * rem;
  else {
   D ou= D(n_[mn].x) * ol;
   auto [a, b]= split_cum<l, neg>(lr[neg], p, ol, ou);
   o[neg]= ol - p, ol-= n_[b].d, ou+= D(n_[b].x) * (o[!neg]= p - ol);
   if constexpr(l) lr[neg]= a;
   else {
    if constexpr(neg) lr[!neg]= join(lr[!neg], mn, a);
    else lr[!neg]= join(a, mn, lr[!neg]);
   }
   if constexpr(neg) y-= ou;
   else y+= ou;
   mn= b;
  }
  rem= 0;
 }
 template <bool r> void add_inf(T x0) {
  if(bf[r] && !lgt<r>(bx[r], x0)) return;
  if(assert(!bf[!r] || !lgt<r>(bx[!r], x0)), bf[r]= true, bx[r]= x0; !mn) return;
  if(lgt<r>(x0, n_[mn].x)) return lr[r]= cut<r>(lr[r], x0), void();
  D q= n_[lr[!r]].s + D(n_[mn].x) * o[!r];
  T v= o[!r] + n_[lr[!r]].a;
  lr[!r]= cut<r>(lr[!r], x0);
  if(!r) y-= q, rem+= v;
  else y+= q, rem-= v;
  if(lr[!r]) std::tie(lr[r], mn)= pop<!r>(lr[!r]), lr[!r]= 0;
  else mn= lr[r]= 0;
  o[r]= n_[mn].d, o[!r]= 0;
 }
 void prop(T x) {
  n_[mn].x+= x;
 }
public:
 // f(x) := 0
 PiecewiseLinearConvex() {}
 // (既定で reserve しない。木の DP は PLC を 1 頂点につき 1 個作って移すので、確保を
 // 先に取ると捨てるぶんが積み上がる。chain 2 万で reserve(256) が 8 倍の代金だった。)
 // f(x) := sum max(0, a(x-x0))
 explicit PiecewiseLinearConvex(const std::vector<std::pair<T, T>>& ramps) {
  int m= ramps.size();
  if(!m) return;
  std::vector<std::pair<T, T>> w(m);
  int s= 0, t= 0;
  for(auto [d, x]: ramps) {
   if(lt(d, 0)) y-= D(d) * x, rem+= d, d= -d;
   if(!lt(0, d)) continue;
   w[s++]= {d, x};
  }
  std::sort(w.begin(), w.begin() + s, [](auto a, auto b) { return a.second < b.second; });
  for(int i= 0; i < s; ++i) {
   if(t && !lt(w[t - 1].second, w[i].second) && !lt(w[i].second, w[t - 1].second)) w[t - 1].first+= w[i].first;
   else w[t++]= w[i];
  }
  if(!t) return;
  mn= create(w[0].first, w[0].second), o[1]= n_[mn].d, lr[1]= build(w.begin() + 1, w.begin() + t);
 }
 // f(x) += c
 void add_const(D c) { y+= c; }
 // f(x) += a x
 void add_linear(T a) { rem+= a; }
 // f(x) += max(a(x-x0), b(x-x0)), (a < b)
 void add_max(T a, T b, T x0) {
  assert(lt(a, b));
  if(bf[0] && x0 <= bx[0]) y-= D(b) * x0, rem+= b;
  else if(bf[1] && bx[1] <= x0) y-= D(a) * x0, rem+= a;
  else if(T c= b - a; mn) {
   if(lt(n_[mn].x, x0)) lr[1]= insert(lr[1], x0, c), y-= D(a) * x0, rem+= a;
   else if(lt(x0, n_[mn].x)) lr[0]= insert(lr[0], x0, c), y-= D(b) * x0, rem+= b;
   else n_[mn].d+= c, o[1]+= c, y-= D(a) * x0, rem+= a;
  } else mn= create(c, x0), y-= D(a) * x0, rem+= a, o[0]= 0, o[1]= c;
 }
 // f(x) += max(0, a(x-x0))
 void add_ramp(T a, T x0) {
  if(lt(0, a)) add_max(0, a, x0);
  else if(lt(a, 0)) add_max(a, 0, x0);
 }
 // f(x) += a|x-x0|
 void add_abs(T a, T x0) {
  if(assert(!lt(a, 0)); lt(0, a)) add_max(-a, a, x0);
 }
 // right=false: x < x0 で f(x) += inf、right=true: x0 < x で f(x) += inf
 void add_inf(bool right= false, T x0= 0) { return right ? add_inf<1>(x0) : add_inf<0>(x0); }
 // f(x) -= max(a(x-x0), b(x-x0)), (a < b)
 void subtract_max(T a, T b, T x0) {
  assert(lt(a, b));
  if(bf[0] && x0 <= bx[0]) y+= D(b) * x0, rem-= b;
  else if(bf[1] && bx[1] <= x0) y+= D(a) * x0, rem-= a;
  else {
   assert(mn);
   T c= b - a;
   if(lt(n_[mn].x, x0)) lr[1]= erase(lr[1], x0, c), y+= D(a) * x0, rem-= a;
   else if(lt(x0, n_[mn].x)) lr[0]= erase(lr[0], x0, c), y+= D(b) * x0, rem-= b;
   else {
    assert(!lt(n_[mn].d, c));
    if(lt(c, n_[mn].d)) n_[mn].d-= c, o[1]+= o[0] - c, y+= D(o[0] + a) * x0, rem-= o[0] + a, o[0]= 0;
    else if(lr[1]) y+= D(o[0] + a) * x0, rem-= o[0] + a, o[0]= 0, std::tie(lr[1], mn)= pop<0>(lr[1]), o[1]= n_[mn].d;
    else if(lr[0]) y-= D(o[1] - b) * x0, rem+= o[1] - b, o[1]= 0, std::tie(lr[0], mn)= pop<1>(lr[0]), o[0]= n_[mn].d;
    else mn= 0, y+= D(a + o[0]) * x0, rem-= a + o[0], o[0]= o[1]= 0;
   }
  }
 }
 // f(x) -= max(0, a(x-x0))
 void subtract_ramp(T a, T x0) {
  if(lt(0, a)) subtract_max(0, a, x0);
  else if(lt(a, 0)) subtract_max(a, 0, x0);
 }
 // f(x) -= a|x-x0|
 void subtract_abs(T a, T x0) {
  if(assert(!lt(a, 0)); lt(0, a)) subtract_max(-a, a, x0);
 }
 // f(x) <- f(x-x0)
 void shift(T x0) {
  if(bx[0]+= x0, bx[1]+= x0, y-= D(rem) * x0; mn) prop(x0), prop(lr[0], x0), prop(lr[1], x0);
 }
 // rev=false: f(x) <- min_{y<=x} f(y)、rev=true: f(x) <- min_{x<=y} f(y)
 void chmin_cum(bool rev= false) {
  if(bf[0] && bf[1] && !lt(bx[0], bx[1])) y+= D(rem) * bx[0], rem= 0;
  else if(bool r= lt(rem, 0); r || lt(0, rem)) {
   T u= (r ? -rem : rem) - o[r] - n_[lr[r]].a;
   if(!lt(u, 0)) {
    if(r ^ rev) {
     if(lt(0, u) && bf[r]) {
      D q= n_[lr[r]].s + D(n_[mn].x) * o[r] + D(u) * bx[r];
      if(r ? y-= q : y+= q; mn) lr[!r]= join(lr[0], mn, lr[1]);
      o[!r]= u, rem= 0, mn= create(u, bx[r]), lr[!rev]= 0, o[!rev]= 0;
     }
    } else {
     assert(bf[r]);
     D q= n_[lr[r]].s + D(n_[mn].x) * o[r] + D(u) * bx[r];
     (r ? y-= q : y+= q), rem= 0, mn= lr[0]= lr[1]= 0, o[0]= o[1]= 0;
    }
    bf[!rev]= false;
    return;
   }
   if((r ^ rev)) r ? slope_eval_cum<0, 1>() : slope_eval_cum<0, 0>();
   else r ? slope_eval_cum<1, 1>() : slope_eval_cum<1, 0>();
   n_[mn].d= o[rev];
  } else if(mn) {
   if(!lt(0, o[rev])) {
    if(lr[rev]) std::tie(lr[rev], mn)= rev ? pop<0>(lr[rev]) : pop<1>(lr[rev]), o[rev]= n_[mn].d;
    else mn= 0;
   } else n_[mn].d= o[rev];
  }
  bf[!rev]= false, lr[!rev]= 0, o[!rev]= 0;
 }
 // f(x) <- min_{lb<=y<=ub} f(x-y)。(lb <= ub)
 void chmin_slide_win(T lb, T ub) {
  assert(lb <= ub);
  if(bf[0] && bf[1] && !lt(bx[0], bx[1])) y+= D(rem) * bx[0], rem= 0;
  else {
   if(bool r= lt(rem, 0); r || lt(0, rem)) {
    T u= (r ? -rem : rem) - o[r] - n_[lr[r]].a;
    if(lt(0, u)) {
     T b[2]= {lb, ub};
     if(bf[r]) {
      D q= n_[lr[r]].s + D(n_[mn].x) * o[r] + D(u) * bx[r];
      if(r ? y-= q : y+= q; mn) lr[!r]= join(lr[0], mn, lr[1]), prop(lr[!r], b[!r]);
      lr[r]= 0, rem= 0, o[!r]= u, o[r]= 0, mn= create(u, bx[r] + b[!r]);
     } else {
      y-= D(rem) * b[!r];
      if(mn) prop(b[!r]), prop(lr[0], b[!r]), prop(lr[1], b[!r]);
     }
     bx[0]+= lb, bx[1]+= ub;
     return;
    }
    slope_eval(r);
   }
   if(mn) {
    if(!lt(0, o[0])) prop(ub);
    else if(!lt(0, o[1])) prop(lb);
    else lr[1]= join<0>(0, create(o[1], n_[mn].x), lr[1]), prop(lb), n_[mn].d= o[0], o[1]= 0;
    prop(lr[0], lb), prop(lr[1], ub);
   }
  }
  bx[0]+= lb, bx[1]+= ub;
 }
 // f(x)。定義域の外は nullopt
 std::optional<D> operator()(T x) {
  if(bf[0] && x < bx[0]) return std::nullopt;
  if(bf[1] && bx[1] < x) return std::nullopt;
  return calc_y(x) + D(rem) * x + y;
 }
 // 最小値。下に有界でなければ nullopt
 std::optional<D> min() {
  bool r= lt(rem, 0);
  if(!r && !lt(0, rem)) return y;
  T u= (r ? -rem : rem) - o[r] - n_[lr[r]].a;
  if(lt(0, u)) {
   if(!bf[r]) return std::nullopt;
   D q= n_[lr[r]].s + D(n_[mn].x) * o[r] + D(u) * bx[r];
   return r ? y - q : y + q;
  }
  return slope_eval(r), y;
 }
 // 最小を取る x の区間 [l, r]。解が無ければ {1, 0}
 std::array<T, 2> argmin() {
  if(bool r= lt(rem, 0); r || lt(0, rem)) {
   if(lt(o[r] + n_[lr[r]].a, (r ? -rem : rem))) {
    if(!bf[r]) return {1, 0};
    return {bx[r], bx[r]};
   }
   slope_eval(r);
  }
  std::array<T, 2> ret= {bx[0], bx[1]};
  int t= mn;
  if(!t) return ret;
  bool r= lt(0, o[0]);
  if(r && lt(0, o[1])) ret[0]= ret[1]= n_[t].x;
  else if(ret[!r]= n_[t].x, t= lr[r]; t) {
   for(; push(t), n_[t].ch[!r];) t= n_[t].ch[!r];
   ret[r]= n_[t].x;
  } else if(!bf[r]) return {1, 0};
  return ret;
 }
 std::size_t size() const { return n_[lr[0]].sz + n_[lr[1]].sz + !!mn; }
 // f += g。g は値で受けるので、右辺値 (`f += dfs(...)` の形) ではノードを移し、this を
 // コピーしない。木の DP で累積に += を続けても O(N^2) にならない。
 PiecewiseLinearConvex& operator+=(PiecewiseLinearConvex g) { return absorb(std::move(g)), *this; }
 PiecewiseLinearConvex operator+(PiecewiseLinearConvex g) const {
  PiecewiseLinearConvex ret= *this;
  ret.absorb(std::move(g));
  return ret;
 }
private:
 // g を this に併合する。足し算は可換だが、unite の引数の向きは mylib と同じ (this が
 // 左) に保つ。向きを変えると木の形が変わり、平らな最小のときに argmin が返す区間が
 // 変わる。ノードの移し直しは常に少ないほうの側だけにする (chain のような形の木の DP で、
 // 深さのぶんコピーが積み上がって O(N^2) になるのを防ぐ)。
 void absorb(PiecewiseLinearConvex&& g) {
  // g の定義域を this に、this の定義域を g に。mylib と同じく、this の定義域は
  // 動かす前のものを g に使う。
  const bool tbf0= bf[0], tbf1= bf[1];
  const T tbx0= bx[0], tbx1= bx[1];
  if(g.bf[0]) add_inf(false, g.bx[0]);
  if(g.bf[1]) add_inf(true, g.bx[1]);
  if(tbf0) g.add_inf(false, tbx0);
  if(tbf1) g.add_inf(true, tbx1);
  y+= g.y, rem+= g.rem;
  if(!g.mn) return;
  // this の定義域 (bf/bx) は、bx が shift の累積も持つので、土台を移すときに控えて戻す。
  const bool fbf0= bf[0], fbf1= bf[1];
  const T fbx0= bx[0], fbx1= bx[1];
  if(g.n_.size() > n_.size()) {
   // g の配列を土台にし、this のノードをその後ろへ移す。添字を付け直してから
   // 土台を this に移し、unite は (this, g) の順で呼ぶ。
   const bool t_has= mn != 0;
   const int t_l0= lr[0], t_mn= mn, t_l1= lr[1], g_l0= g.lr[0], g_mn= g.mn, g_l1= g.lr[1];
   const T to0= o[0], go0= g.o[0];
   const D ysum= y;
   const T remsum= rem;
   const int base= (int)g.n_.size(), off= base - 1;
   for(std::size_t i= 1; i < n_.size(); ++i) g.n_.push_back(n_[i]);
   for(int i= base; i < (int)g.n_.size(); ++i)
    for(int c= 0; c < 2; ++c) g.n_[i].ch[c]= g.n_[i].ch[c] ? g.n_[i].ch[c] + off : 0;
   const int tl0= t_l0 ? t_l0 + off : 0, tmn= t_mn ? t_mn + off : 0, tl1= t_l1 ? t_l1 + off : 0;
   *this= std::move(g);
   y= ysum, rem= remsum;
   bf[0]= fbf0, bf[1]= fbf1, bx[0]= fbx0, bx[1]= fbx1;
   if(!t_has) return;
   y+= n_[tl0].s + D(n_[tmn].x) * to0 + n_[g_l0].s + D(n_[g_mn].x) * go0;
   rem-= to0 + n_[tl0].a + go0 + n_[g_l0].a;
   int t= unite(join(tl0, tmn, tl1), join(g_l0, g_mn, g_l1));
   std::tie(lr[1], mn)= pop<0>(t), lr[0]= 0, o[0]= 0, o[1]= n_[mn].d;
   return;
  }
  // this の配列を土台にし、g のノードをその後ろへ移す。reserve を毎回ちょうどの
  // 大きさで呼ぶと、1 ノードずつ併合する形で配列を作り直し続けて O(N^2) になる。
  const int base= (int)n_.size(), off= base - 1;
  for(std::size_t i= 1; i < g.n_.size(); ++i) n_.push_back(g.n_[i]);
  for(int i= base; i < (int)n_.size(); ++i)
   for(int c= 0; c < 2; ++c) n_[i].ch[c]= n_[i].ch[c] ? n_[i].ch[c] + off : 0;
  const int gmn= g.mn + off, gl0= g.lr[0] ? g.lr[0] + off : 0, gl1= g.lr[1] ? g.lr[1] + off : 0;
  if(!mn) {
   mn= gmn, lr[0]= gl0, lr[1]= gl1, o[0]= g.o[0], o[1]= g.o[1];
   return;
  }
  y+= n_[lr[0]].s + D(n_[mn].x) * o[0] + n_[gl0].s + D(n_[gmn].x) * g.o[0];
  rem-= o[0] + n_[lr[0]].a + g.o[0] + n_[gl0].a;
  int t= unite(join(lr[0], mn, lr[1]), join(gl0, gmn, gl1));
  std::tie(lr[1], mn)= pop<0>(t), lr[0]= 0, o[0]= 0, o[1]= n_[mn].d;
 }
 template <class Iter> int build(Iter bg, Iter ed) {
  if(bg == ed) return 0;
  auto md= bg + (ed - bg) / 2;
  int t= create(md->first, md->second);
  return n_[t].ch[0]= build(bg, md), n_[t].ch[1]= build(md + 1, ed), update(t), t;
 }
};
