#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 無向グラフの二辺連結成分分解。comp[v] は頂点 v の成分の番号 (0 以上 count 未満)。自己ループと多重辺があってもよい。
// DFS の前に次数 1 以下の頂点を待ち行列で剥がす。剥がした頂点の辺は橋なので、剥がした頂点はどれも 1 頂点の成分。頂点ごとに次数と
// 残っている隣接の番号の XOR を 1 つの配列に並べて持ち、次数が 1 になった頂点の残りの隣接は XOR の値そのものなので、隣接を組まずに
// 辺の列だけで剥がせる。積むかどうかは分岐にせず、待ち行列の先の頂点を先読みする。剥がしたあとの 2-core の頂点に 0 から番号を
// 振り直し、2-core の CSR を組む (端が剥がれた辺と自己ループは、辺の番号で回す 8 つの番兵に書いて捨てる)。DFS は頂点ごとの記録
// {CSR の始まり, 終わり, 最初の隣接, 2 番目の隣接} で回す path-based (Gabow) で、親への木の辺は 1 本だけ飛ばし (平行な辺は後ろ向きの
// 辺)、次数 2 の頂点では親でない方の隣接 n0 ^ n1 ^ 親 だけを見る。成分の番号は、2-core の成分を DFS で閉じた順に振り、そのあとに
// 剥がした頂点を頂点の番号の順に振る。作業領域は 2 MB 以上なら、Linux では 2 MB 境界の mmap に置いて MADV_HUGEPAGE を頼む。
// procon-judge の yosupo-two-edge-connected-components で書き比べた gabow_crec_aos_d8pf_hp を写したもの。
struct TwoEdgeConnectedComponents {
 int count;              // 成分の数
 std::vector<int> comp;  // 頂点ごとの成分の番号 (0 以上 count 未満)
};
namespace two_edge_connected_components_internal {
using u32= unsigned;
constexpr u32 NIL= ~u32(0);
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
}
inline TwoEdgeConnectedComponents two_edge_connected_components(int n, const std::vector<std::array<int, 2>>& edges) {
 using namespace two_edge_connected_components_internal;
 struct R {
  u32 s, e, n0, n1;
 };
 struct F {
  u32 v, i, e, par, pv;
 };
 const u32 N= n;
 const std::size_t m= edges.size();
 auto W= [](std::size_t x) { return (x + 15) & ~std::size_t(15); };
 Buffer buf(W(2 * std::size_t(N + 1)) + 5 * W(N + 1) + W(N + 10) + W(2 * m + 2) + W(5 * std::size_t(N)) + W(4 * std::size_t(N)));
 u32* cur= buf.get();
 auto take= [&](std::size_t x) {
  u32* r= cur;
  cur+= W(x);
  return r;
 };
 // dg[2 v] が次数、dg[2 v + 1] が残っている隣接の XOR。id は頂点から 2-core の番号 (剥がした頂点は NIL)、q は剥がしの待ち行列で、
 // 使い終わったら 2-core の番号から頂点への表にする。off は番兵 8 つの分だけ多めに取る。鍵は DFS のスタックに載っている間はその位置
 // (1 から)、成分に入ったら NIL - (成分の番号)。
 u32 *dg= take(2 * std::size_t(N + 1)), *id= take(N + 1), *q= take(N + 1), *off= take(N + 10), *adj= take(2 * m + 2);
 u32 *st= take(N + 1), *bs= take(N + 1), *key= take(N + 1);
 F* fs= reinterpret_cast<F*>(take(5 * std::size_t(N)));
 R* rec= reinterpret_cast<R*>(take(4 * std::size_t(N)));
 auto D= [&](u32 v) -> u32& { return dg[2 * v]; };
 auto X= [&](u32 v) -> u32& { return dg[2 * v + 1]; };
 std::memset(dg, 0, 2 * std::size_t(N) * sizeof(u32));
 // 次数と XOR。自己ループは足す値を 0 にして分岐なしで捨てる。16 本先の辺の端を先読みする。
 for(std::size_t i= 0; i < m; ++i) {
  const auto& f= edges[std::min(i + 16, m - 1)];
  __builtin_prefetch(&D(f[0])), __builtin_prefetch(&D(f[1]));
  const u32 a= edges[i][0], b= edges[i][1], d= a != b;
  D(a)+= d, D(b)+= d, X(a)^= b & -d, X(b)^= a & -d;
 }
 // 次数 1 の頂点を剥がす。積むかどうかは分岐にせず、毎回 q の末尾に書いて、次数がちょうど 1 になったときだけ末尾を進める。
 // 出した頂点の次数が 0 なら XOR も 0 で、隣接は頂点 0 になり、0 を足し引きするだけで済む。16 個先と 8 個先の頂点を先読みする。
 u32 qt= 0;
 for(u32 v= 0; v < N; ++v) q[qt]= v, qt+= D(v) == 1;
 for(u32 qh= 0; qh < qt; ++qh) {
  const u32 a= q[std::min(qh + 16, qt - 1)];
  __builtin_prefetch(&D(a));
  const u32 ub= X(q[std::min(qh + 8, qt - 1)]);
  __builtin_prefetch(&D(ub));
  const u32 v= q[qh], live= D(v), u= X(v);
  D(v)= 0, X(u)^= v & -live;
  const u32 du= D(u) - live;
  D(u)= du, q[qt]= u, qt+= live & (du == 1);
 }
 // 2-core の頂点に番号を振り、CSR を組む。off[x] は x までの次数の累積で、辺を置くたびに減らして始まりにする。端が剥がれた辺と
 // 自己ループは、番兵 n2 + 1 から n2 + 8 に辺の番号で回して書き、番兵の位置は減らさない (同じ番地の読み書きが辺ごとに直列に
 // つながらない)。番兵の書く位置は 2-core の辺の終わりで、off[n2] もそこを指す。16 本先の辺の端の番号と 8 本先の辺の端の off を
 // 先読みする。
 u32* orig= q;
 u32 s= 0, n2= 0;
 for(u32 v= 0; v < N; ++v) {
  const u32 d= D(v), c= -u32(d >= 2);
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
 std::memset(key, 0, (n2 + 1) * sizeof(u32));
 for(u32 x= 0; x < n2; ++x) {
  const u32 b= off[x], e= off[x + 1];
  rec[x]= {b, e, adj[b], adj[b + 1]};  // 2-core の頂点の次数は 2 以上
 }
 // path-based。bs は成分の根の候補の位置のスタック。いまの頂点の隣接の位置 i が s0 なら n0、s0 + 1 なら n1、それ以外は adj[i]。
 // 子から帰ったら s0 を NIL - 1 にして adj から読む。pv はいまの頂点のスタックでの位置。無向グラフの DFS には横向きの辺が無く、
 // 成分に入った頂点への辺は子孫への辺なので、成分に入った頂点の鍵をどの位置より大きい値にしておけば、何も捨てない。
 u32 top= 0, bt= 0, nc= 0;
 for(u32 r= 0; r < n2; ++r) {
  if(key[r]) continue;
  const R c= rec[r];
  u32 v= r, i= c.s, e= c.e, s0= c.s, n0= c.n0, n1= c.n1, par= NIL;
  st[top++]= v, key[v]= top, bs[bt++]= top;
  u32 pv= top;
  F* sp= fs;
  for(;;) {
   if(i != e) {
    const u32 w= i == s0 ? n0 : i == s0 + 1 ? n1 : adj[i];
    ++i;
    if(w == par) {
     par= NIL;
     continue;
    }
    const u32 kw= key[w];
    if(!kw) {
     *sp++= {v, i, e, par, pv};
     const R cw= rec[w];
     par= v, v= w, i= s0= cw.s, e= cw.e, n0= cw.n0, n1= cw.n1;
     const u32 m2= -u32(e - i == 2);
     n0^= (n1 ^ par) & m2, e-= m2 & 1, par|= m2;
     st[top++]= v, key[v]= pv= top, bs[bt++]= top;
    } else
     while(kw < bs[bt - 1]) --bt;
   } else {
    if(pv == bs[bt - 1]) {
     --bt;
     const u32 k= NIL - nc++;
     while(top >= pv) key[st[--top]]= k;
    }
    if(sp == fs) break;
    --sp, v= sp->v, i= sp->i, e= sp->e, par= sp->par, pv= sp->pv, s0= NIL - 1;
   }
  }
 }
 // 剥がした頂点は、2-core の成分のあとに頂点の番号の順に番号を振る。key[n2] は 0 なので、剥がした頂点 (id = NIL) で引いてもよい。
 TwoEdgeConnectedComponents res{0, std::vector<int>(N)};
 u32 k= nc;
 for(u32 v= 0; v < N; ++v) {
  const u32 x= id[v], p= -u32(x == NIL);
  res.comp[v]= int(((NIL - key[std::min(x, n2)]) & ~p) | (k & p)), k-= p;
 }
 res.count= int(k);
 return res;
}
