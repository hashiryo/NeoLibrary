#pragma once
#include <algorithm>
#include <array>
#include <vector>
// マトロイドの交差 (matroid_intersection.hpp) に渡すマトロイドの型。どれも次を持つ。
//   int size() const                       台集合の大きさ
//   void build(const std::vector<int>& I)  独立集合 I に合わせて中の構造を作り直す
//   bool can_add(int x)                    x が I に無いとき、I + x が独立か
// さらに、2 引数の交換の判定 can_exchange(y, x) (I - y + x が独立か) と、閉路の列挙 circuit(x, f) (I + x の閉路に入る
// I の要素 y ごとに f(y)) のうち、型にとって書きやすいほうを持つ。
// 型は、独立と分かっている要素を足すだけの add(x) (x が I に無く、I + x が独立のとき、中の構造を作り直さずに x を足し、
// can_add だけを正しく保つ。次の build まで can_exchange と circuit は使えない) も持てる。持つ型は、マトロイドの和
// (matroid_union.hpp) が、そのまま入る要素ごとの build を省く。GraphicMatroid、PartitionMatroid、XorMatroid が持つ。
// - GraphicMatroid: 辺の集合で、森になるもの。両方を持つ。
// - PartitionMatroid: 要素ごとの組で、組 c からは cap[c] 個までのもの。両方を持つ。
// - XorMatroid: 64 bit のベクトルの集合で、xor で線形独立なもの。両方を持つ。
// - TransversalMatroid: 二部グラフの左の頂点の集合で、マッチングで覆えるもの。circuit だけを持つ。
// - CographicMatroid: 辺の集合で、取り除いても各連結成分が連結のままのもの (グラフ的マトロイドの双対)。can_exchange
//   だけを持つ。

// グラフ的マトロイド。build で森を根付き木にし、深さ、親への辺、オイラー順の入りと出を持つ。can_exchange は、y の下側の
// 頂点の部分木に x の両端のちょうど一方が入るかで O(1)。circuit は x の両端から深いほうを親へたどる。
class GraphicMatroid {
 int n_;
 std::vector<std::array<int, 2>> e_;
 std::vector<int> start_, adj_, comp_, dep_, pe_, pv_, tin_, tout_, it_, stk_;
 // add のあとの連結。comp_ の値 (成分の根) の上の Union-Find で、added_ のときだけ使う。
 std::vector<int> cu_;
 bool added_= false;
 int cfind(int v) {
  while(cu_[v] != v) v= cu_[v]= cu_[cu_[v]];
  return v;
 }
public:
 GraphicMatroid(int n, const std::vector<std::array<int, 2>>& edges): n_(n), e_(edges), start_(n + 1), comp_(n), dep_(n), pe_(n), pv_(n), tin_(n), tout_(n), it_(n), cu_(n) {}
 int size() const { return e_.size(); }
 void build(const std::vector<int>& I) {
  added_= false;
  std::fill(start_.begin(), start_.end(), 0);
  for(int y: I) ++start_[e_[y][0] + 1], ++start_[e_[y][1] + 1];
  for(int v= 0; v < n_; ++v) start_[v + 1]+= start_[v];
  adj_.resize(2 * I.size());
  for(int v= 0; v < n_; ++v) it_[v]= start_[v];
  for(int y: I) adj_[it_[e_[y][0]]++]= y, adj_[it_[e_[y][1]]++]= y;
  std::fill(comp_.begin(), comp_.end(), -1);
  for(int r= 0, t= 0; r < n_; ++r) {
   if(comp_[r] >= 0) continue;
   comp_[r]= r, dep_[r]= 0, pe_[r]= -1, pv_[r]= -1, tin_[r]= t++, it_[r]= start_[r];
   stk_.assign(1, r);
   while(!stk_.empty()) {
    const int v= stk_.back();
    if(it_[v] == start_[v + 1]) {
     tout_[v]= t, stk_.pop_back();
     continue;
    }
    const int y= adj_[it_[v]++], u= e_[y][0] ^ e_[y][1] ^ v;
    if(comp_[u] >= 0) continue;
    comp_[u]= r, dep_[u]= dep_[v] + 1, pe_[u]= y, pv_[u]= v, tin_[u]= t++, it_[u]= start_[u];
    stk_.push_back(u);
   }
  }
 }
 [[gnu::always_inline]] bool can_add(int x) {
  if(!added_) return comp_[e_[x][0]] != comp_[e_[x][1]];
  return cfind(comp_[e_[x][0]]) != cfind(comp_[e_[x][1]]);
 }
 void add(int x) {
  if(!added_) {
   added_= true;
   for(int v= 0; v < n_; ++v) cu_[v]= v;
  }
  cu_[cfind(comp_[e_[x][0]])]= cfind(comp_[e_[x][1]]);
 }
 [[gnu::always_inline]] bool can_exchange(int y, int x) {
  if(can_add(x)) return true;
  const int c= dep_[e_[y][0]] > dep_[e_[y][1]] ? e_[y][0] : e_[y][1];
  auto inside= [&](int v) { return tin_[c] <= tin_[v] && tin_[v] < tout_[c]; };
  return inside(e_[x][0]) != inside(e_[x][1]);
 }
 template <class F> void circuit(int x, F&& f) const {
  int u= e_[x][0], v= e_[x][1];
  while(u != v) {
   if(dep_[u] < dep_[v]) std::swap(u, v);
   f(pe_[u]), u= pv_[u];
  }
 }
};
// 分割マトロイド。要素 x は組 part[x] に属し (-1 なら制限なし)、組 c からは cap[c] 個まで選べる。閉路は、x と同じ組で
// I に入っている要素の全部。
class PartitionMatroid {
 std::vector<int> part_, cap_, cnt_, start_, mem_;
public:
 PartitionMatroid(const std::vector<int>& part, const std::vector<int>& cap): part_(part), cap_(cap), cnt_(cap.size()), start_(cap.size() + 1) {}
 int size() const { return part_.size(); }
 void build(const std::vector<int>& I) {
  std::fill(cnt_.begin(), cnt_.end(), 0);
  for(int y: I)
   if(part_[y] >= 0) ++cnt_[part_[y]];
  for(size_t c= 0; c < cnt_.size(); ++c) start_[c + 1]= start_[c] + cnt_[c];
  mem_.resize(start_.back());
  for(int y: I)
   if(part_[y] >= 0) mem_[start_[part_[y] + 1] - cnt_[part_[y]]--]= y;
  for(size_t c= 0; c < cnt_.size(); ++c) cnt_[c]= start_[c + 1] - start_[c];
 }
 [[gnu::always_inline]] bool can_add(int x) const { return part_[x] < 0 || cnt_[part_[x]] < cap_[part_[x]]; }
 void add(int x) {
  if(part_[x] >= 0) ++cnt_[part_[x]];
 }
 [[gnu::always_inline]] bool can_exchange(int y, int x) const { return can_add(x) || part_[y] == part_[x]; }
 template <class F> void circuit(int x, F&& f) const {
  for(int k= start_[part_[x]]; k < start_[part_[x] + 1]; ++k) f(mem_[k]);
 }
};
// GF(2) の線形マトロイド (64 bit のベクトル)。build で I の基底を掃き出し、各行が I のどの要素の xor かを 64 bit の印で
// 持つ (|I| <= 64)。x を基底で消して 0 になれば従属で、残った印が閉路。消した結果は覚えておき、基底が変わるまで使い回す。
// add は行を足すだけなので、従属だった要素の結果はそのまま使え、独立だった要素の結果だけが古くなる (add の回数の印を
// 別に付けて見分ける)。build は I が今の並びと同じなら何もしない (add した順の I を渡されたときがこれ)。
class XorMatroid {
 std::vector<unsigned long long> v_, memo_;
 std::vector<int> idx_, elem_;
 std::vector<unsigned> stamp_, astamp_;
 unsigned long long row_[64]= {}, mask_[64]= {};
 unsigned now_= 1, adds_= 0;
 // x を消したあとの印。独立なら ~0
 unsigned long long reduce(int x) {
  if(stamp_[x] == now_ && (memo_[x] != ~0ULL || astamp_[x] == adds_)) return memo_[x];
  stamp_[x]= now_;
  unsigned long long a= v_[x], m= 0;
  while(a) {
   const int b= 63 - __builtin_clzll(a);
   if(!row_[b]) return astamp_[x]= adds_, memo_[x]= ~0ULL;
   a^= row_[b], m^= mask_[b];
  }
  return memo_[x]= m;
 }
 void insert(int x) {
  const int k= elem_.size();
  elem_.push_back(x), idx_[x]= k;
  unsigned long long a= v_[x], m= 1ULL << k;
  while(a) {
   const int b= 63 - __builtin_clzll(a);
   if(!row_[b]) {
    row_[b]= a, mask_[b]= m;
    break;
   }
   a^= row_[b], m^= mask_[b];
  }
 }
public:
 explicit XorMatroid(const std::vector<unsigned long long>& v): v_(v), memo_(v.size()), idx_(v.size(), -1), stamp_(v.size(), 0), astamp_(v.size(), 0) {}
 int size() const { return v_.size(); }
 void build(const std::vector<int>& I) {
  if(I == elem_) return;
  for(int y: elem_) idx_[y]= -1;
  elem_.clear();
  std::fill(row_, row_ + 64, 0);
  for(int y: I) insert(y);
  ++now_;
 }
 bool can_add(int x) { return reduce(x) == ~0ULL; }
 void add(int x) { insert(x), ++adds_; }
 bool can_exchange(int y, int x) {
  const unsigned long long m= reduce(x);
  return m == ~0ULL || (m >> idx_[y] & 1);
 }
 template <class F> void circuit(int x, F&& f) {
  for(unsigned long long m= reduce(x); m; m&= m - 1) f(elem_[__builtin_ctzll(m)]);
 }
};
// 横断マトロイド。左の頂点 l 個、右の頂点 r 個の二部グラフで、左の頂点の集合のうちマッチングで覆えるものを独立とする
// (要素は左の頂点)。build は前のマッチングのうち I に残る組を残し、I の空いた頂点から増加路を探して覆い直す。そのあと、
// 空いた右の頂点から交互路を逆にたどって、空いた右の頂点へ届く右の頂点 (good) に印を付ける。can_add(x) は x の隣に
// good があるか。circuit は x から交互路を幅優先でたどり、届いた右の頂点の相手 (I の要素) を列挙する。
class TransversalMatroid {
 int l_, r_;
 std::vector<int> ls_, la_, rs_, ra_, ml_, mr_, q_, par_;
 std::vector<char> in_, good_;
 std::vector<unsigned> vis_;
 unsigned tag_= 0;
 bool augment(int s) {
  ++tag_;
  q_.assign(1, s);
  for(size_t h= 0; h < q_.size(); ++h) {
   const int u= q_[h];
   for(int k= ls_[u]; k < ls_[u + 1]; ++k) {
    const int r= la_[k];
    if(vis_[r] == tag_) continue;
    vis_[r]= tag_, par_[r]= u;
    if(mr_[r] < 0) {
     for(int c= r; c >= 0;) {
      const int a= par_[c], nc= ml_[a];
      ml_[a]= c, mr_[c]= a, c= nc;
     }
     return true;
    }
    q_.push_back(mr_[r]);
   }
  }
  return false;
 }
public:
 TransversalMatroid(int l, int r, const std::vector<std::array<int, 2>>& edges): l_(l), r_(r), ls_(l + 1), rs_(r + 1), ml_(l, -1), mr_(r, -1), par_(r), in_(l), good_(r), vis_(r, 0) {
  for(auto& [a, b]: edges) ++ls_[a + 1], ++rs_[b + 1];
  for(int i= 0; i < l; ++i) ls_[i + 1]+= ls_[i];
  for(int i= 0; i < r; ++i) rs_[i + 1]+= rs_[i];
  la_.resize(edges.size()), ra_.resize(edges.size());
  std::vector<int> pl(ls_.begin(), ls_.end() - 1), pr(rs_.begin(), rs_.end() - 1);
  for(auto& [a, b]: edges) la_[pl[a]++]= b, ra_[pr[b]++]= a;
 }
 int size() const { return l_; }
 void build(const std::vector<int>& I) {
  std::fill(in_.begin(), in_.end(), 0);
  for(int y: I) in_[y]= 1;
  for(int a= 0; a < l_; ++a)
   if(!in_[a] && ml_[a] >= 0) mr_[ml_[a]]= -1, ml_[a]= -1;
  for(int y: I)
   if(ml_[y] < 0) augment(y);
  std::fill(good_.begin(), good_.end(), 0);
  q_.clear();
  for(int r= 0; r < r_; ++r)
   if(mr_[r] < 0) good_[r]= 1, q_.push_back(r);
  for(size_t h= 0; h < q_.size(); ++h)
   for(int k= rs_[q_[h]]; k < rs_[q_[h] + 1]; ++k)
    if(const int a= ra_[k]; ml_[a] >= 0 && !good_[ml_[a]]) good_[ml_[a]]= 1, q_.push_back(ml_[a]);
 }
 bool can_add(int x) const {
  for(int k= ls_[x]; k < ls_[x + 1]; ++k)
   if(good_[la_[k]]) return true;
  return false;
 }
 template <class F> void circuit(int x, F&& f) {
  ++tag_;
  q_.assign(1, x);
  for(size_t h= 0; h < q_.size(); ++h)
   for(int k= ls_[q_[h]]; k < ls_[q_[h] + 1]; ++k)
    if(const int r= la_[k]; vis_[r] != tag_) vis_[r]= tag_, f(mr_[r]), q_.push_back(mr_[r]);
 }
};
// グラフ的マトロイドの双対。辺の集合 I のうち、G から I を取り除いても連結成分の数が変わらないものを独立とする。build は
// G - I の DFS の森を作り、lowlink で橋を求め、オイラー順の入りと出を持つ。can_add(x) は x が G - I の橋でないこと。
// can_exchange(y, x) は、x が橋なら、x の下側の頂点の部分木に y の両端のちょうど一方が入るか。
class CographicMatroid {
 int n_;
 std::vector<std::array<int, 2>> e_;
 std::vector<int> start_, adj_, tin_, tout_, low_, pe_, it_, stk_, child_;
 std::vector<char> in_, bridge_;
public:
 CographicMatroid(int n, const std::vector<std::array<int, 2>>& edges): n_(n), e_(edges), start_(n + 1), tin_(n), tout_(n), low_(n), pe_(n), it_(n), child_(edges.size()), in_(edges.size()), bridge_(edges.size()) {}
 int size() const { return e_.size(); }
 void build(const std::vector<int>& I) {
  std::fill(in_.begin(), in_.end(), 0);
  for(int y: I) in_[y]= 1;
  const int m= e_.size();
  std::fill(start_.begin(), start_.end(), 0);
  for(int x= 0; x < m; ++x)
   if(!in_[x]) ++start_[e_[x][0] + 1], ++start_[e_[x][1] + 1];
  for(int v= 0; v < n_; ++v) start_[v + 1]+= start_[v];
  adj_.resize(start_[n_]);
  for(int v= 0; v < n_; ++v) it_[v]= start_[v];
  for(int x= 0; x < m; ++x)
   if(!in_[x]) adj_[it_[e_[x][0]]++]= x, adj_[it_[e_[x][1]]++]= x;
  std::fill(bridge_.begin(), bridge_.end(), 0);
  std::fill(tin_.begin(), tin_.end(), -1);
  for(int r= 0, t= 0; r < n_; ++r) {
   if(tin_[r] >= 0) continue;
   tin_[r]= low_[r]= t++, pe_[r]= -1, it_[r]= start_[r];
   stk_.assign(1, r);
   while(!stk_.empty()) {
    const int v= stk_.back();
    if(it_[v] == start_[v + 1]) {
     tout_[v]= t, stk_.pop_back();
     if(pe_[v] >= 0) {
      const int p= e_[pe_[v]][0] ^ e_[pe_[v]][1] ^ v;
      low_[p]= std::min(low_[p], low_[v]);
      if(low_[v] > tin_[p]) bridge_[pe_[v]]= 1, child_[pe_[v]]= v;
     }
     continue;
    }
    const int x= adj_[it_[v]++];
    if(x == pe_[v]) continue;
    const int u= e_[x][0] ^ e_[x][1] ^ v;
    if(tin_[u] >= 0) low_[v]= std::min(low_[v], tin_[u]);
    else tin_[u]= low_[u]= t++, pe_[u]= x, it_[u]= start_[u], stk_.push_back(u);
   }
  }
 }
 [[gnu::always_inline]] bool can_add(int x) const { return !bridge_[x]; }
 [[gnu::always_inline]] bool can_exchange(int y, int x) const {
  if(!bridge_[x]) return true;
  const int c= child_[x];
  auto inside= [&](int v) { return tin_[c] <= tin_[v] && tin_[v] < tout_[c]; };
  return inside(e_[y][0]) != inside(e_[y][1]);
 }
};
