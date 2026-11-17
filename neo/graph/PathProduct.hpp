#pragma once
#include <cassert>
#include <vector>
#include "neo/graph/HeavyLightDecomposition.hpp"
// 木の頂点に値を置き、1 点の代入と、u から v への道の積 (道の順の op、可換でなくてよい) を O(log N) で答える。構造体 M は
// SegmentTree と同じモノイド (using T、static T ti()、static T op(T a, T b))。木は HeavyLightDecomposition で渡す (中身は写すので、
// 作ったあとに HLD を捨ててよい)。
// HLD の鎖ごとに二分木を組む (中国の競プロで全局平衡二叉树と呼ばれる形で、Library Checker の vertex_set_path_composite の首位
// chaihf の 270199 と同じ組み方)。軽い子を先に配る行きがけ順の位置 vp を作ると、鎖の隣どうしの vp の差は、上の頂点と軽い子の
// 部分木の大きさの和になる。鎖の頂点の優先度を「vp と鎖の次の頂点の vp の xor の最上位 bit」として、優先度の大きい方を根にした
// 位置の順の二分木を組むと、ぶら下がる部分木の大きい頂点ほど根に近くなる。軽い辺で上の鎖へ移ると、入った頂点の重みは下の鎖の
// 全体の大きさ以上なので、鎖の二分木を上がる段数を道に沿って足すと O(log N) に収まる。
// 節は自分の値 v と、部分木の前向きの積 f (位置の小さい順) と後ろ向きの積 b、左の部分木と自分の積 tf と tb を持つ。鎖の頂点 k から
// 鎖の先頭までの積は、k から二分木の根へ上がり、位置が k 以下の祖先の tf か tb を集めて求める。1 点の代入は、二分木の根まで
// 上がりながら節ごとに op を 4 回する。セグメント木 2 本で書く形と比べると、軽い辺を多く通る木では速く、一直線の木で代入が
// 多いと遅い (記録は algo-notes の notes/tree-queries.md)。
template <class M> class PathProduct {
 using T= typename M::T;
 using u64= hld_internal::u64;
 struct Node {
  T v, f, b, tf, tb;
  int l, r, p;
 };
 int n= 0;
 std::vector<Node> nd;           // 位置 + 1 → 節。nd[0] は番兵で、f と b は ti()
 std::vector<unsigned char> pr;  // 位置 + 1 → 優先度
 std::vector<u64> ch;            // 位置 → 鎖の先頭の位置と先頭の親の位置 (HLD と同じ)
 std::vector<int> pos;           // 頂点 → 位置
 [[gnu::always_inline]] void maintain(int u) {
  Node& z= nd[u];
  z.tf= M::op(nd[z.l].f, z.v), z.f= M::op(z.tf, nd[z.r].f);
  z.tb= M::op(z.v, nd[z.l].b), z.b= M::op(nd[z.r].b, z.tb);
 }
 // 鎖の先頭から節 k までの、後ろ向きの積 (k から先頭へ) と前向きの積。
 T bprefix(int k) const {
  T s= M::ti();
  for(int u= k; u; u= nd[u].p)
   if(u <= k) s= M::op(s, nd[u].tb);
  return s;
 }
 T fprefix(int k) const {
  T s= M::ti();
  for(int u= k; u; u= nd[u].p)
   if(u <= k) s= M::op(nd[u].tf, s);
  return s;
 }
 // 同じ鎖の節 [l, r] の、後ろ向きの積 (r から l へ) と前向きの積。l と r から二分木を上がり、出会った節で合わせる。
 // 優先度は祖先ほど大きく、等しいのは右の子だけなので、等しいときは r の側を上げる。
 T brange(int l, int r) const {
  T al= M::ti(), ar= M::ti();
  int u= l, v= r;
  while(u != v) {
   if(pr[u] < pr[v]) {
    if(l <= u) al= M::op(M::op(nd[nd[u].r].b, nd[u].v), al);
    u= nd[u].p;
   } else {
    if(v <= r) ar= M::op(ar, nd[v].tb);
    v= nd[v].p;
   }
  }
  return M::op(M::op(ar, nd[u].v), al);
 }
 T frange(int l, int r) const {
  T al= M::ti(), ar= M::ti();
  int u= l, v= r;
  while(u != v) {
   if(pr[u] < pr[v]) {
    if(l <= u) al= M::op(al, M::op(nd[u].v, nd[nd[u].r].f));
    u= nd[u].p;
   } else {
    if(v <= r) ar= M::op(nd[v].tf, ar);
    v= nd[v].p;
   }
  }
  return M::op(M::op(al, nd[u].v), ar);
 }
 template <class F> void build(const HeavyLightDecomposition& h, const F& init) {
  n= int(h.size());
  nd.assign(n + 1, Node{M::ti(), M::ti(), M::ti(), M::ti(), M::ti(), 0, 0, 0});
  pr.assign(n + 1, 0), ch.resize(n), pos.resize(n);
  std::vector<int> sz(n), vp(n), nxt(n);
  for(int i= 0; i < n; ++i) {
   const int v= h.to_vertex(i), hv= h.head(v), hp= h.to_seq(hv), up= h.parent(hv);
   pos[v]= i, sz[i]= h.subtree_size(v);
   ch[i]= hld_internal::chain(hp, up < 0 ? -1 : h.to_seq(up));
  }
  // 軽い子を先に配る位置 vp (1 始まり) を、HLD の位置の小さい順 (親が子より前) に決める。HLD の位置では重い子は親のすぐ後ろ。
  for(int i= 0; i < n; ++i) {
   const int hp= hld_internal::head_of(ch[i]), up= hld_internal::up_of(ch[i]);
   if(hp != i) vp[i]= vp[i - 1] + sz[i - 1] - sz[i];  // 重い子は親の部分木の最後
   else if(up < 0) vp[i]= i + 1;                       // 根は HLD と同じ位置
   else vp[i]= nxt[up], nxt[up]+= sz[i];               // 軽い子は親の次から順に
   nxt[i]= vp[i] + 1;
  }
  // 優先度と、鎖ごとの二分木。節の番号は位置 + 1。鎖の頂点を上から順に、右の背骨を上がりながら差し込む。
  for(int i= 0; i < n; ++i) {
   const int x= i + 1, heavy= sz[i] > 1 ? sz[i + 1] : 0;
   const unsigned a= unsigned(vp[i]), b= unsigned(vp[i] + sz[i] - heavy);
   pr[x]= (unsigned char)(31 - __builtin_clz(a ^ b));
   nd[x].v= init(h.to_vertex(i)), nd[x].l= nd[x].r= nd[x].p= 0;
   if(hld_internal::head_of(ch[i]) != i) {
    int q= x - 1, l= 0;
    while(q && pr[q] < pr[x]) l= q, q= nd[q].p;
    nd[q].r= x, nd[x].p= q, nd[x].l= l, nd[l].p= x;
   }
  }
  nd[0].l= nd[0].r= nd[0].p= 0;
  // 帰りがけ順に作る。右の子の無い節 x から、位置が x 以下の祖先へ上がりながら作ると、どの節も左右を作り終えてから作られる。
  for(int x= 1; x <= n; ++x)
   if(nd[x].r == 0)
    for(int u= x; u && u <= x; u= nd[u].p) maintain(u);
 }
public:
 PathProduct()= default;
 // a は頂点ごとの値。
 PathProduct(const HeavyLightDecomposition& h, const std::vector<T>& a) {
  assert(a.size() == h.size());
  build(h, [&](int v) { return a[v]; });
 }
 // 値をすべて ti() にする。
 explicit PathProduct(const HeavyLightDecomposition& h) {
  build(h, [](int) { return M::ti(); });
 }
 int size() const { return n; }
 T get(int v) const { return nd[pos[v] + 1].v; }
 void set(int v, T x) {
  int u= pos[v] + 1;
  nd[u].v= x;
  for(; u; u= nd[u].p) maintain(u);
 }
 // op(a[u], ..., a[v]) (u から v への道の順)。
 T prod(int u, int v) const {
  int x= pos[u], y= pos[v];
  u64 a= ch[x], b= ch[y];
  T up= M::ti(), dn= M::ti();
  while(hld_internal::head_of(a) != hld_internal::head_of(b)) {
   if(hld_internal::head_of(a) > hld_internal::head_of(b)) up= M::op(up, bprefix(x + 1)), x= hld_internal::up_of(a), a= ch[x];
   else dn= M::op(fprefix(y + 1), dn), y= hld_internal::up_of(b), b= ch[y];
  }
  if(x >= y) up= M::op(up, brange(y + 1, x + 1));
  else dn= M::op(frange(x + 1, y + 1), dn);
  return M::op(up, dn);
 }
};
