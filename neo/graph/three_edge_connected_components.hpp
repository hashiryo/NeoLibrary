#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 無向グラフの三辺連結成分分解。comp[v] は頂点 v の成分の番号 (0 以上 count 未満)。自己ループと多重辺があってもよい。
// DFS の前に次数 1 以下の頂点を待ち行列で剥がす (剥がした頂点の辺は橋なので、どれも 1 頂点の成分)。頂点ごとに次数と残っている隣接の
// 番号の XOR を 1 つの配列に並べて持ち、隣接を組まずに辺の列だけで、分岐なしに剥がす。剥がしたあとの 2-core の頂点に 0 から番号を
// 振り直して CSR を組み (端が剥がれた辺と自己ループは、辺の番号で回す 8 つの番兵に書いて捨てる)、Tsin の吸収の形の 1 回の DFS で
// 分ける。2-core の次数 2 の頂点 (鎖の頂点) はどれも 1 頂点の成分で、鎖を 1 本の辺に縮めても (直列の縮約) ほかの頂点の成分は変わら
// ないので、DFS には次数 3 以上の頂点 (核) だけを載せる。核の頂点から鎖の頂点に当たったら、記録 {CSR の始まり, 終わり, 最初の 2 つの
// 隣接} の XOR で鎖の向こう端の核の頂点まで歩き、鎖の両端の頂点に両端の核の頂点の XOR を印として残して、向こう端を隣接とみなす。
// 反対の端から当たったときは印から引く。核を含まない閉路には触らない。DFS の頂点ごとの状態は 32 byte の構造体 1 つにまとめ、親への
// 木の辺は 1 本だけ飛ばす。作業領域は 2 MB 以上なら、Linux では 2 MB 境界の mmap に置いて MADV_HUGEPAGE を頼む。
// procon-judge の yosupo-three-edge-connected-components で書き比べた tsin_chain_hp を写したもの。
struct ThreeEdgeConnectedComponents {
 int count;              // 成分の数
 std::vector<int> comp;  // 頂点ごとの成分の番号 (0 以上 count 未満)
};
namespace three_edge_connected_components_internal {
using u32= unsigned;
constexpr u32 NIL= ~u32(0), CM= NIL - 1;  // CM: 歩いた鎖の端の頂点の in
// u32 の n 個分の作業領域。0 では埋めない。
class Buffer {
 static constexpr std::size_t H= std::size_t(1) << 21, P= 4096;
 u32* p= nullptr;
#ifdef __linux__
 std::size_t mapped= 0;  // mmap で取った大きさ
#endif
public:
 explicit Buffer(std::size_t n) {
#ifdef __linux__
  if(n * sizeof(u32) >= H) {
   const std::size_t bytes= (n * sizeof(u32) + P - 1) & ~(P - 1);
   void* m= mmap(nullptr, bytes + H, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
   if(m != MAP_FAILED) {
    char *q= static_cast<char*>(m), *a= reinterpret_cast<char*>((reinterpret_cast<std::uintptr_t>(q) + H - 1) & ~(H - 1));
    if(a != q) munmap(q, a - q);
    munmap(a + bytes, q + H - a);
    madvise(a, bytes, MADV_HUGEPAGE);
    p= reinterpret_cast<u32*>(a), mapped= bytes;
    return;
   }
  }
#endif
  p= new u32[n];
 }
 ~Buffer() {
#ifdef __linux__
  if(mapped) {
   munmap(p, mapped);
   return;
  }
#endif
  delete[] p;
 }
 Buffer(const Buffer&)= delete;
 Buffer& operator=(const Buffer&)= delete;
 u32* get() const { return p; }
};
// v に n 個分の領域を取り、Linux ではその領域のページに huge page を頼む (領域の内側にある 2 MB の区切りだけが載る)。
template <class T> inline void reserve_huge(std::vector<T>& v, std::size_t n) {
 v.reserve(n);
#ifdef __linux__
 const std::uintptr_t b= reinterpret_cast<std::uintptr_t>(v.data()) & ~std::uintptr_t(4095), e= (reinterpret_cast<std::uintptr_t>(v.data() + n) + 4095) & ~std::uintptr_t(4095);
 if(b < e) madvise(reinterpret_cast<void*>(b), e - b, MADV_HUGEPAGE);
#endif
}
}
inline ThreeEdgeConnectedComponents three_edge_connected_components(int n, const std::vector<std::array<int, 2>>& edges) {
 using namespace three_edge_connected_components_internal;
 struct R {
  u32 s, e, n0, n1;
 };
 struct F {
  u32 v, i, e, par;
 };
 // DFS の頂点ごとの状態。in は行きがけ順 (NIL なら未訪問、鎖の端なら CM で、low に鎖の両端の核の頂点の XOR を持つ)、out は帰りの
 // 行きがけ順、deg は部分木から出る辺の数の見積もり (符号なしで、負なら大きい値)、path は下へ伸びるまだ吸収していない頂点の列、
 // par は吸収した祖先。
 struct S {
  u32 in, out, low, deg, path, par, pad0, pad1;
 };
 const u32 N= n;
 const std::size_t m= edges.size();
 auto W= [](std::size_t x) { return (x + 15) & ~std::size_t(15); };
 Buffer buf(W(2 * std::size_t(N + 1)) + 4 * W(N + 10) + W(2 * m + 2) + W(4 * std::size_t(N) + 4) + W(8 * std::size_t(N) + 8) + W(4 * std::size_t(N) + 4));
 u32* cur= buf.get();
 auto take= [&](std::size_t x) {
  u32* r= cur;
  cur+= W(x);
  return r;
 };
 u32 *dg= take(2 * std::size_t(N + 1)), *q= take(N + 1), *off= take(N + 10), *adj= take(2 * m + 2), *id= take(N + 1), *ord= take(N + 1);
 R* rec= reinterpret_cast<R*>(take(4 * std::size_t(N) + 4));
 S* st= reinterpret_cast<S*>(take(8 * std::size_t(N) + 8));
 F* fs= reinterpret_cast<F*>(take(4 * std::size_t(N) + 4));
 // dg[2 v] が次数、dg[2 v + 1] が残っている隣接の XOR。自己ループは足す値を 0 にして捨てる。16 本先の辺の端を先読みする。
 std::memset(dg, 0, 2 * std::size_t(N) * sizeof(u32));
 for(std::size_t i= 0; i < m; ++i) {
  const auto& f= edges[std::min(i + 16, m - 1)];
  __builtin_prefetch(dg + 2 * f[0]), __builtin_prefetch(dg + 2 * f[1]);
  const u32 a= edges[i][0], b= edges[i][1], d= a != b;
  dg[2 * a]+= d, dg[2 * b]+= d, dg[2 * a + 1]^= b & -d, dg[2 * b + 1]^= a & -d;
 }
 // 次数 1 の頂点を剥がす。積むかどうかは分岐にせず、毎回 q の末尾に書いて、次数がちょうど 1 になったときだけ末尾を進める。出した頂点の
 // 次数が 0 なら XOR も 0 で、隣接は頂点 0 になり、0 を足し引きするだけで済む。16 個先と 8 個先の頂点を先読みする。
 u32 qt= 0;
 for(u32 v= 0; v < N; ++v) q[qt]= v, qt+= dg[2 * v] == 1;
 for(u32 qh= 0; qh < qt; ++qh) {
  __builtin_prefetch(dg + 2 * q[std::min(qh + 16, qt - 1)]);
  __builtin_prefetch(dg + 2 * dg[2 * q[std::min(qh + 8, qt - 1)] + 1]);
  const u32 v= q[qh], live= dg[2 * v], u= dg[2 * v + 1];
  dg[2 * v]= 0, dg[2 * u + 1]^= v & -live;
  const u32 du= dg[2 * u] - live;
  dg[2 * u]= du, q[qt]= u, qt+= live & (du == 1);
 }
 // 2-core の頂点に番号を振り、CSR を組む。q は使い終わったので、2-core の番号から頂点への表にする。off[x] は x までの次数の累積で、
 // 辺を置くたびに減らして始まりにする。捨てる辺は番兵 n2 + 1 から n2 + 8 に回して書き、番兵の位置は減らさない。
 u32* orig= q;
 u32 s= 0, n2= 0;
 for(u32 v= 0; v < N; ++v) {
  const u32 d= dg[2 * v], c= -u32(d >= 2);
  id[v]= n2 | ~c, orig[n2]= v, s+= d, off[n2]= s, n2-= c;
 }
 for(u32 j= 0; j <= 8; ++j) off[n2 + j]= s;
 for(std::size_t i= 0; i < m; ++i) {
  const auto& f= edges[std::min(i + 16, m - 1)];
  __builtin_prefetch(id + f[0]), __builtin_prefetch(id + f[1]);
  const auto& g= edges[std::min(i + 8, m - 1)];
  __builtin_prefetch(off + std::min(id[g[0]], n2)), __builtin_prefetch(off + std::min(id[g[1]], n2));
  const u32 ia= id[edges[i][0]], ib= id[edges[i][1]], sk= -u32((ia == NIL) | (ib == NIL) | (ia == ib)), dm= n2 + 1 + (u32(i) & 7);
  const u32 a= (ia & ~sk) | (dm & sk), b= (ib & ~sk) | (dm & sk), l= 1 + sk;
  const u32 pa= off[a] - l;
  adj[pa]= b, off[a]= pa;
  const u32 pb= off[b] - l;
  adj[pb]= a, off[b]= pb;
 }
 for(u32 x= 0; x < n2; ++x) {
  const u32 b= off[x], e= off[x + 1];
  rec[x]= {b, e, adj[b], adj[b + 1]};  // 2-core の頂点の次数は 2 以上
  st[x].in= NIL;
 }
 // Tsin の DFS。いまの頂点の隣接の位置 i が s0 なら n0、s0 + 1 なら n1、それ以外は adj[i]。子から帰ったら s0 を NIL - 1 にして adj
 // から読む。子 w から帰ったとき、木の辺 (v, w) が 2 本の辺の切断に入る (w の道が空で deg[w] <= 1) なら w の道はそこで切れ、そうでなければ
 // w の道と v の道のうち low の小さい方を v の道にし、もう片方を v に吸収する。子孫 w から v への後ろ向きの辺を見たら、v の道のうち w の
 // 祖先にあたる頂点を v に吸収する。
 u32 cnt= 0;
 for(u32 r= 0; r < n2; ++r) {
  if(st[r].in != NIL) continue;
  const R c= rec[r];
  if(c.e - c.s == 2) continue;  // 根は核の頂点だけ
  u32 v= r, i= c.s, e= c.e, s0= c.s, n0= c.n0, n1= c.n1, par= NIL;
  F* sp= fs;
  st[v]= {cnt, 0, NIL, 0, NIL, v, 0, 0}, ord[cnt++]= v;
  for(;;) {
   if(i != e) {
    u32 w= i == s0 ? n0 : i == s0 + 1 ? n1 : adj[i];
    ++i;
    u32 iw= st[w].in;
    if(iw == CM) w= st[w].low ^ v, iw= st[w].in;
    else if(iw == NIL) {
     R rw= rec[w];
     if(rw.e - rw.s == 2) {
      u32 prev= v, now= w;
      do {
       const u32 nx= rw.n0 ^ rw.n1 ^ prev;
       prev= now, now= nx, rw= rec[now];
      } while(rw.e - rw.s == 2);
      st[w].in= CM, st[w].low= v ^ now, st[prev].in= CM, st[prev].low= v ^ now;
      w= now, iw= st[w].in;
     }
    }
    if(w == v) continue;  // 1 つの核の頂点に戻る鎖 (縮めると自己ループ)
    if(w == par) {
     par= NIL;
     continue;
    }
    if(iw != NIL) {
     S& sv= st[v];
     if(iw < sv.in) ++sv.deg, sv.low= std::min(sv.low, iw);
     else {
      --sv.deg;
      u32 u= sv.path;
      while(u != NIL && st[u].in <= iw && iw < st[u].out) st[u].par= v, sv.deg+= st[u].deg, u= st[u].path;
      sv.path= u;
     }
     continue;
    }
    *sp++= {v, i, e, par};
    const R cw= rec[w];
    par= v, v= w, i= s0= cw.s, e= cw.e, n0= cw.n0, n1= cw.n1;
    st[v]= {cnt, 0, NIL, 0, NIL, v, 0, 0}, ord[cnt++]= v;
   } else {
    st[v].out= cnt;
    if(sp == fs) break;
    u32 w= v;
    --sp, v= sp->v, i= sp->i, e= sp->e, par= sp->par, s0= NIL - 1;
    S& sv= st[v];
    if(st[w].path == NIL && st[w].deg <= 1) {
     sv.deg+= st[w].deg, sv.low= std::min(sv.low, st[w].low);
     continue;
    }
    if(st[w].deg == 0) w= st[w].path;
    if(st[w].low < sv.low) sv.low= st[w].low, std::swap(w, sv.path);
    while(w != NIL) st[w].par= v, sv.deg+= st[w].deg, w= st[w].path;
   }
  }
 }
 // 成分の番号。核の頂点は行きがけ順に、吸収した祖先の番号を写す (吸収する祖先は先に訪れている)。DFS に載らない 2-core の頂点と剥がした
 // 頂点は 1 頂点の成分。番号は off に置く (off は使い終わった)。
 u32* cn= off;
 u32 k= 0;
 for(u32 t= 0; t < cnt; ++t) {
  const u32 x= ord[t], y= st[x].par;
  cn[x]= y == x ? k++ : cn[y];
 }
 for(u32 x= 0; x < n2; ++x)
  if(st[x].in >= CM) cn[x]= k++;
 ThreeEdgeConnectedComponents res;
 reserve_huge(res.comp, N);
 res.comp.assign(N, 0);
 for(u32 v= 0; v < N; ++v) {
  const u32 x= id[v];
  res.comp[v]= int(x == NIL ? k++ : cn[x]);
 }
 res.count= int(k);
 return res;
}
