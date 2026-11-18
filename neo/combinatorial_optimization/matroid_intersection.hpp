#pragma once
#include <algorithm>
#include <array>
#include <concepts>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>
// マトロイドの交差。台集合 [0, n) が同じ 2 つのマトロイド M1、M2 の共通独立集合を求める。
// - matroid_intersection(m1, m2) は、要素数が最大の共通独立集合 (要素の番号の列) を返す。
// - weighted_matroid_intersection(m1, m2, w, k) は、要素数 j ごとの重みの和の最大 best[j] (j = 0, 1, ...) と、要素数が
//   最大 (k >= 0 なら k で打ち切る) で重みの和が best.back() の集合を返す。
// マトロイドの型は次を持つ。例は neo/combinatorial_optimization/matroids.hpp にある。
//   int size() const                    台集合の大きさ n
//   void build(const std::vector<int>& I)  独立集合 I に合わせて中の構造を作り直す
//   bool can_add(int x)                 x が I に無いとき、I + x が独立か
// さらに、次の 2 つのうち書きやすいほう (両方でもよい) を持つ。足りないほうはここで補う。
//   bool can_exchange(int y, int x)     y が I にあり x が I に無いとき、I - y + x が独立か
//   void circuit(int x, F f)            I + x が独立でないとき、I + x のただ 1 つの閉路に入る I の要素 y ごとに f(y) を呼ぶ
// can_exchange が無い型には、x の閉路を build ごとに覚えて、y が入っているかを答える。circuit が無い型には、I の要素を
// 全部 can_exchange で調べて閉路を作る (重み付きで増加路 1 本ごとに O(n |I|) 回の判定になる)。
//
// 重みなしは、2 引数の判定で交換グラフの層を BFS で作り、層に沿った最短の増加路を DFS で次々に探す (Cunningham の方法)。
// M2 への問いは、どれも x を固定して y を変える形にしている。can_exchange を持たない型は M2 の側に回す。
// 重み付きは、閉路の列挙で交換グラフを作り、要素の長さ (I に無い x は -w(x)、I にある y は w(y)) の和と辺の数の組で最小
// の路を Bellman-Ford (待ち行列) で求めて、1 本ずつ増やす。
// procon-judge の self-matroid-graphic-partition、self-matroid-xor-partition、self-matroid-weighted-graphic-partition で
// 書き比べた形 (中身は problems/_shared/matroid/flex.hpp)。記録は algo-notes の notes/matroid.md にある。
namespace matroid_internal {
template <class M>
concept HasExchange= requires(M& m, int y, int x) {
 { m.can_exchange(y, x) } -> std::convertible_to<bool>;
};
template <class M>
concept HasCircuit= requires(M& m, int x) { m.circuit(x, [](int) {}); };
template <class M> class Adapter {
 static_assert(HasExchange<M> || HasCircuit<M>, "a matroid needs can_exchange or circuit");
 M& m_;
 const std::vector<int>* I_= nullptr;
 // circuit だけの型の覚え書き。cache_[x] は stamp_[x] == ver_ のとき今の I での x の閉路。mark_[y] == tag_ は、y が owner_ の
 // 閉路に入っていること。
 std::vector<std::vector<int>> cache_;
 std::vector<unsigned> stamp_, mark_;
 unsigned ver_= 1, tag_= 0;
 int owner_= -1;
public:
 explicit Adapter(M& m): m_(m) {
  if constexpr(!HasExchange<M>) cache_.resize(m.size()), stamp_.assign(m.size(), 0), mark_.assign(m.size(), 0);
 }
 int size() const { return m_.size(); }
 void build(const std::vector<int>& I) {
  I_= &I, m_.build(I);
  if constexpr(!HasExchange<M>) ++ver_, owner_= -1;
 }
 [[gnu::always_inline]] bool can_add(int x) { return m_.can_add(x); }
 [[gnu::always_inline]] bool can_exchange(int y, int x) {
  if constexpr(HasExchange<M>) return m_.can_exchange(y, x);
  else {
   if(m_.can_add(x)) return true;
   if(owner_ != x) {
    owner_= x, ++tag_;
    if(stamp_[x] != ver_) stamp_[x]= ver_, cache_[x].clear(), m_.circuit(x, [&](int y) { cache_[x].push_back(y); });
    for(int z: cache_[x]) mark_[z]= tag_;
   }
   return mark_[y] == tag_;
  }
 }
 template <class F> void circuit(int x, F&& f) {
  if constexpr(HasCircuit<M>) m_.circuit(x, f);
  else
   for(int y: *I_)
    if(m_.can_exchange(y, x)) f(y);
 }
};
// 重みなし。b への問いは x を固定する形だけ。a への問いも、層を作るところと終点に届かない要素を除くところは x を固定する
// 順に回す (増加路の DFS で I の要素 u から次の層の x を探すところだけ x が変わる)。
template <class A, class B> std::vector<int> intersect(A& a, B& b) {
 const int n= a.size();
 std::vector<char> in(n), dead(n), good(n);
 std::vector<int> I, rest_out, rest_in, it, stk;
 auto rebuild= [&] {
  I.erase(std::remove_if(I.begin(), I.end(), [&](int y) { return !in[y]; }), I.end());
  for(int v: stk)
   if(in[v]) I.push_back(v);
  a.build(I), b.build(I);
 };
 rebuild();
 std::vector<std::vector<int>> L;
 while(true) {
  // 層を作る。偶数の層は I に無い要素、奇数の層は I の要素。L[0] は M1 に足せる要素で、M2 に足せる要素が現れた層で止める。
  L.assign(1, {});
  rest_out.clear();
  for(int x= 0; x < n; ++x)
   if(!in[x]) (a.can_add(x) ? L[0] : rest_out).push_back(x);
  rest_in= I;
  int top= -1;
  for(int d= 0; !L[d].empty(); d+= 2) {
   if(std::any_of(L[d].begin(), L[d].end(), [&](int x) { return b.can_add(x); })) {
    top= d;
    break;
   }
   L.emplace_back(), L.emplace_back();
   for(int x: L[d])
    for(size_t k= 0; k < rest_in.size();)
     if(const int y= rest_in[k]; b.can_exchange(y, x)) L[d + 1].push_back(y), rest_in[k]= rest_in.back(), rest_in.pop_back();
     else ++k;
   for(size_t k= 0; k < rest_out.size();) {
    const int x= rest_out[k];
    if(std::any_of(L[d + 1].begin(), L[d + 1].end(), [&](int y) { return a.can_exchange(y, x); })) L[d + 2].push_back(x), rest_out[k]= rest_out.back(), rest_out.pop_back();
    else ++k;
   }
  }
  if(top < 0) break;
  // 今の I で終点まで届かない要素を、後ろの層から除く。
  {
   auto& T= L[top];
   T.erase(std::remove_if(T.begin(), T.end(), [&](int x) { return !b.can_add(x); }), T.end());
   for(int d= top - 1; d >= 0; --d) {
    auto& C= L[d];
    const auto& nx= L[d + 1];
    if(d % 2 == 0) C.erase(std::remove_if(C.begin(), C.end(), [&](int u) { return std::none_of(nx.begin(), nx.end(), [&](int v) { return b.can_exchange(v, u); }); }), C.end());
    else {
     for(int y: C) good[y]= 0;
     for(int x: nx)
      for(int y: C)
       if(!good[y] && a.can_exchange(y, x)) good[y]= 1;
     C.erase(std::remove_if(C.begin(), C.end(), [&](int y) { return !good[y]; }), C.end());
    }
   }
  }
  // 層に沿った最短の増加路を、今の I での判定で次々に探す。行き止まりの要素はこの段階の残りで使わない。
  std::fill(dead.begin(), dead.end(), 0);
  it.assign(top + 1, 0);
  bool any= false;
  for(int s: L[0]) {
   if(dead[s] || in[s] || !a.can_add(s)) continue;
   stk.assign(1, s), it[0]= 0;
   while(!stk.empty()) {
    const int d= (int)stk.size() - 1, u= stk.back();
    if(d == top) {
     if(b.can_add(u)) {
      for(int v: stk) in[v]^= 1, dead[v]= 1;
      rebuild(), any= true;
      break;
     }
     dead[u]= 1, stk.pop_back();
     continue;
    }
    const auto& nx= L[d + 1];
    bool moved= false;
    while(it[d] < (int)nx.size()) {
     const int v= nx[it[d]++];
     if(dead[v]) continue;
     if(d % 2 == 0 ? b.can_exchange(v, u) : a.can_exchange(u, v)) {
      stk.push_back(v), it[d + 1]= 0, moved= true;
      break;
     }
    }
    if(!moved) dead[u]= 1, stk.pop_back();
   }
  }
  stk.clear();
  if(!any) break;
 }
 return I;
}
}  // namespace matroid_internal
template <class M1, class M2> std::vector<int> matroid_intersection(M1&& m1, M2&& m2) {
 using T1= std::remove_cvref_t<M1>;
 using T2= std::remove_cvref_t<M2>;
 matroid_internal::Adapter<T1> a(m1);
 matroid_internal::Adapter<T2> b(m2);
 if constexpr(!matroid_internal::HasExchange<T1> && matroid_internal::HasExchange<T2>) return matroid_internal::intersect(b, a);
 else return matroid_internal::intersect(a, b);
}
template <class M1, class M2, class T> std::pair<std::vector<T>, std::vector<int>> weighted_matroid_intersection(M1&& m1, M2&& m2, const std::vector<T>& w, int k= -1) {
 matroid_internal::Adapter<std::remove_cvref_t<M1>> a(m1);
 matroid_internal::Adapter<std::remove_cvref_t<M2>> b(m2);
 const int n= a.size();
 const T INF= std::numeric_limits<T>::max() / 4;
 std::vector<char> in(n), src(n), snk(n), inq(n);
 std::vector<int> I, hop(n), par(n), start(n + 1), to, path, q(n + 1);
 std::vector<T> d(n), best(1, T());
 std::vector<std::array<int, 2>> es;
 a.build(I), b.build(I);
 while(k < 0 || (int)I.size() < k) {
  // 交換グラフ。y ∈ C1(x) なら y → x、y ∈ C2(x) なら x → y。
  es.clear();
  for(int x= 0; x < n; ++x) {
   src[x]= snk[x]= 0;
   if(in[x]) continue;
   if(a.can_add(x)) src[x]= 1;
   else a.circuit(x, [&](int y) { es.push_back({y, x}); });
   if(b.can_add(x)) snk[x]= 1;
   else b.circuit(x, [&](int y) { es.push_back({x, y}); });
  }
  std::fill(start.begin(), start.end(), 0);
  for(auto& [u, v]: es) ++start[u + 1];
  for(int u= 0; u < n; ++u) start[u + 1]+= start[u];
  to.resize(es.size());
  for(auto& [u, v]: es) to[start[u]++]= v;
  for(int u= n; u > 0; --u) start[u]= start[u - 1];
  start[0]= 0;
  // 長さの和と辺の数の組で最小の路。待ち行列は長さ n + 1 の輪 (同じ要素は同時に 1 つしか入らない)。始点へ入る辺は使わない。
  std::fill(d.begin(), d.end(), INF);
  int qh= 0, qt= 0;
  for(int x= 0; x < n; ++x)
   if(src[x]) d[x]= -w[x], hop[x]= 0, par[x]= -1, inq[x]= 1, q[qt++]= x;
  while(qh != qt) {
   const int u= q[qh];
   qh= qh == n ? 0 : qh + 1, inq[u]= 0;
   for(int j= start[u]; j < start[u + 1]; ++j) {
    const int v= to[j];
    if(src[v]) continue;
    const T nd= d[u] + (in[v] ? w[v] : -w[v]);
    if(nd < d[v] || (nd == d[v] && hop[u] + 1 < hop[v])) {
     d[v]= nd, hop[v]= hop[u] + 1, par[v]= u;
     if(!inq[v]) inq[v]= 1, q[qt]= v, qt= qt == n ? 0 : qt + 1;
    }
   }
  }
  int t= -1;
  for(int x= 0; x < n; ++x)
   if(snk[x] && d[x] < INF && (t < 0 || d[x] < d[t] || (d[x] == d[t] && hop[x] < hop[t]))) t= x;
  if(t < 0) break;
  path.clear();
  for(int v= t; v >= 0; v= par[v]) path.push_back(v);
  for(int v: path) in[v]^= 1;
  I.erase(std::remove_if(I.begin(), I.end(), [&](int y) { return !in[y]; }), I.end());
  for(int v: path)
   if(in[v]) I.push_back(v);
  a.build(I), b.build(I);
  best.push_back(best.back() - d[t]);
 }
 return {best, I};
}
