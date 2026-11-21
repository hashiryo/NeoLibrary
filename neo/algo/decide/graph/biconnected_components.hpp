#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 無向グラフの二重連結成分 (関節点で分けた成分) 分解。成分 k の頂点は vertex[start[k]] から vertex[start[k + 1] - 1]。関節点は
// 複数の成分に入り、辺を持たない頂点はそれだけで 1 つの成分。自己ループは無く、多重辺はあってもよい。
// DFS の前に次数 1 の頂点を待ち行列で剥がす。剥がした頂点 v と、そのとき残っていた唯一の隣接 u は、それだけで 1 つの成分 {v, u}
// (辺 vu は橋)。頂点ごとに次数と残っている隣接の番号の XOR を 1 つの配列に並べて持ち、u は XOR の値そのものなので、隣接を組まずに
// 辺の列だけで剥がせる。積むかどうか、成分を置くかどうかは分岐にせず、待ち行列の先の頂点を先読みする。剥がしたあとの 2-core の
// 頂点に 0 から番号を振り直し、2-core の CSR を組む (端が剥がれた辺は、辺の番号で回す 8 つの番兵に書いて捨てる)。DFS は頂点ごとの
// 記録 {CSR の始まり, 終わり, 最初の隣接, 2 番目の隣接} で回す Hopcroft–Tarjan で、頂点のスタックから成分を切り出す。下りた直後は
// 親への辺を番号の比較で飛ばし (low[v] >= ord[親] の判定は変わらないので、多重辺もすべて飛ばしてよい)、次数 2 の頂点では親で
// ない方の隣接 n0 ^ n1 ^ 親 だけを見る。low はフレームに持つ。作業領域は 2 MB 以上なら、Linux では 2 MB 境界の mmap に置いて
// MADV_HUGEPAGE を頼む。procon-judge の yosupo-biconnected-components で書き比べた tarjan_crec_aos_d8pf_hp を写したもの。
struct BiconnectedComponents {
 int count;                // 成分の数
 std::vector<int> start;   // 長さ count + 1
 std::vector<int> vertex;  // 成分 k の頂点は vertex[start[k]] から vertex[start[k + 1] - 1]
};
namespace biconnected_components_internal {
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
inline BiconnectedComponents biconnected_components(int n, const std::vector<std::array<int, 2>>& edges) {
 using namespace biconnected_components_internal;
 struct R {
  u32 s, e, n0, n1;
 };
 struct F {
  u32 v, i, e, lv, ov;
 };
 const u32 N= n;
 const std::size_t m= edges.size();
 auto W= [](std::size_t x) { return (x + 15) & ~std::size_t(15); };
 Buffer buf(W(2 * std::size_t(N + 1)) + 4 * W(N + 1) + W(N + 10) + W(2 * m + 2) + 2 * W(N + m + 2) + W(5 * std::size_t(N)) + W(4 * std::size_t(N)));
 u32* cur= buf.get();
 auto take= [&](std::size_t x) {
  u32* r= cur;
  cur+= W(x);
  return r;
 };
 // dg[2 v] が次数、dg[2 v + 1] が残っている隣接の XOR。id は頂点から 2-core の番号 (剥がした頂点は NIL)、q は剥がしの待ち行列で、
 // 使い終わったら 2-core の番号から頂点への表にする。off は番兵 8 つの分だけ多めに取る。成分 k の頂点は bv[bo[k]] から
 // bv[bo[k + 1] - 1] で、bv の nv0 から後ろ (DFS が置いた分) は 2-core の番号。
 u32 *dg= take(2 * std::size_t(N + 1)), *id= take(N + 1), *q= take(N + 1), *off= take(N + 10), *adj= take(2 * m + 2);
 u32 *ord= take(N + 1), *st= take(N + 1), *bv= take(N + m + 2), *bo= take(N + m + 2);
 F* fs= reinterpret_cast<F*>(take(5 * std::size_t(N)));
 R* rec= reinterpret_cast<R*>(take(4 * std::size_t(N)));
 auto D= [&](u32 v) -> u32& { return dg[2 * v]; };
 auto X= [&](u32 v) -> u32& { return dg[2 * v + 1]; };
 std::memset(dg, 0, 2 * std::size_t(N) * sizeof(u32));
 // 次数と XOR。16 本先の辺の端を先読みする。
 for(std::size_t i= 0; i < m; ++i) {
  const auto& f= edges[std::min(i + 16, m - 1)];
  __builtin_prefetch(&D(f[0])), __builtin_prefetch(&D(f[1]));
  const u32 a= edges[i][0], b= edges[i][1];
  ++D(a), ++D(b), X(a)^= b, X(b)^= a;
 }
 // 辺を持たない頂点は 1 頂点の成分。成分を置くかどうかも分岐にせず、毎回書いて、置くときだけ末尾を進める。
 u32 nv= 0, nb= 0, qt= 0;
 bo[0]= 0;
 for(u32 v= 0; v < N; ++v) {
  const u32 d= D(v);
  bv[nv]= v, nv+= d == 0, bo[nb + 1]= nv, nb+= d == 0, q[qt]= v, qt+= d == 1;
 }
 // 次数 1 の頂点を剥がす。積むかどうかは分岐にせず、毎回 q の末尾に書いて、次数がちょうど 1 になったときだけ末尾を進める。
 // 出した頂点の次数が 0 (残りの隣接がもう剥がれた) なら XOR も 0 で、隣接は頂点 0 になり、0 を足し引きするだけで済む。成分
 // {v, u} を置くかどうかも次数の値 (0 か 1) で選ぶ。16 個先と 8 個先の頂点を先読みする。
 for(u32 qh= 0; qh < qt; ++qh) {
  const u32 a= q[std::min(qh + 16, qt - 1)];
  __builtin_prefetch(&D(a));
  const u32 ub= X(q[std::min(qh + 8, qt - 1)]);
  __builtin_prefetch(&D(ub));
  const u32 v= q[qh], live= D(v), u= X(v);
  bv[nv]= v, bv[nv + 1]= u, nv+= 2 * live, bo[nb + 1]= nv, nb+= live;
  D(v)= 0, X(u)^= v & -live;
  const u32 du= D(u) - live;
  D(u)= du, q[qt]= u, qt+= live & (du == 1);
 }
 // 2-core の頂点に番号を振り、CSR を組む。off[x] は x までの次数の累積で、辺を置くたびに減らして始まりにする。端が剥がれた辺は、
 // 番兵 n2 + 1 から n2 + 8 に辺の番号で回して書き、番兵の位置は減らさない (同じ番地の読み書きが辺ごとに直列につながらない)。
 // 番兵の書く位置は 2-core の辺の終わりで、off[n2] もそこを指す。16 本先の辺の端の番号と 8 本先の辺の端の off を先読みする。
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
  const u32 ia= id[edges[i][0]], ib= id[edges[i][1]], sk= -u32((ia == NIL) | (ib == NIL)), dm= n2 + 1 + (u32(i) & 7);
  const u32 a= (ia & ~sk) | (dm & sk), b= (ib & ~sk) | (dm & sk), l= 1 + sk;
  const u32 pa= off[a] - l;
  adj[pa]= b, off[a]= pa;
  const u32 pb= off[b] - l;
  adj[pb]= a, off[b]= pb;
 }
 const u32 nv0= nv;
 std::memset(ord, 0, n2 * sizeof(u32));
 for(u32 x= 0; x < n2; ++x) {
  const u32 b= off[x], e= off[x + 1];
  rec[x]= {b, e, adj[b], adj[b + 1]};  // 2-core の頂点の次数は 2 以上
 }
 // Hopcroft–Tarjan。いまの頂点の隣接の位置 i が s0 なら n0、s0 + 1 なら n1、それ以外は adj[i]。子から帰ったら s0 を NIL - 1 に
 // して adj から読み、親への辺も後ろ向きの辺として扱う (判定は変わらない)。lv はいまの頂点の lowlink、ov は訪問順。子 c から帰って
 // low[c] >= ord[いまの頂点] なら、スタックの c から上といまの頂点が 1 つの成分。
 u32 idx= 0, top= 0;
 for(u32 r= 0; r < n2; ++r) {
  if(ord[r]) continue;
  const R c= rec[r];
  u32 v= r, i= c.s, e= c.e, s0= c.s, n0= c.n0, n1= c.n1, par= NIL, lv= ++idx, ov= lv;
  F* sp= fs;
  ord[v]= ov, st[top++]= v;
  for(;;) {
   if(i != e) {
    const u32 w= i == s0 ? n0 : i == s0 + 1 ? n1 : adj[i];
    ++i;
    if(w == par) continue;
    const u32 ow= ord[w];
    if(!ow) {
     *sp++= {v, i, e, lv, ov};
     const R cw= rec[w];
     par= v, v= w, i= s0= cw.s, e= cw.e, n0= cw.n0, n1= cw.n1;
     const u32 m2= -u32(e - i == 2);
     n0^= (n1 ^ par) & m2, e-= m2 & 1, par|= m2;
     lv= ov= ord[v]= ++idx, st[top++]= v;
    } else if(ow < lv) lv= ow;
   } else {
    if(sp == fs) {
     --top;  // 根。子の成分はどれも子から帰るときに切り出してあるので、根だけが残っている
     break;
    }
    const u32 lc= lv, ch= v;
    --sp, v= sp->v, i= sp->i, e= sp->e, lv= sp->lv, ov= sp->ov, s0= NIL - 1, par= NIL;
    if(lc >= ov) {
     u32 x;
     do x= st[--top], bv[nv++]= x;
     while(x != ch);
     bv[nv++]= v, bo[++nb]= nv;
    } else if(lc < lv) lv= lc;
   }
  }
 }
 BiconnectedComponents res{int(nb), std::vector<int>(bo, bo + nb + 1), std::vector<int>(nv)};
 for(u32 j= 0; j < nv0; ++j) res.vertex[j]= int(bv[j]);
 for(u32 j= nv0; j < nv; ++j) res.vertex[j]= int(orig[bv[j]]);
 return res;
}
