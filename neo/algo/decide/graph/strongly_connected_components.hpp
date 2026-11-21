#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 有向グラフの強連結成分分解。comp[v] は頂点 v の成分の番号で、成分を縮めた DAG のトポロジカル順に振る (どの辺 (a, b) でも
// comp[a] <= comp[b])。DFS の前に入次数 0 の頂点を待ち行列で剥がす。剥がした頂点はどれも 1 頂点の成分で、剥がした順が
// トポロジカル順の先頭に来る。剥がされない頂点から剥がした頂点への辺は無いので、残りは剥がした頂点に触れずに、path-based
// (Gabow) の反復の DFS で分ける。頂点ごとに、鍵、最初の隣接、残りの隣接の連結リストの先頭を、別々の 4 byte の配列に置く。
// DFS が下りる 1 歩は最初の隣接の読み込み 1 段で済み、配列が小さいので L2 にも載りやすい。隣接は辺の列を 1 回なめて組み、
// 同じ走査で入次数も数える。剥がしでは待ち行列の先の頂点を先読みする。作業領域は 2 MB 以上なら、Linux では 2 MB 境界の mmap
// に置いて MADV_HUGEPAGE を頼む。procon-judge の yosupo-scc で書き比べた gabow_split_trimpf_fused_hp を写したもの。
struct StronglyConnectedComponents {
 int count;              // 成分の数
 std::vector<int> comp;  // 頂点ごとの成分の番号 (0 以上 count 未満)
};
namespace scc_internal {
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
inline StronglyConnectedComponents strongly_connected_components(int n, const std::vector<std::array<int, 2>>& edges) {
 using namespace scc_internal;
 struct F {
  u32 v, e;  // 頂点と、残りの隣接の連結リストのうちまだ見ていない先頭
 };
 const u32 N= n;
 const std::size_t m= edges.size();
 auto W= [](std::size_t x) { return (x + 15) & ~std::size_t(15); };
 Buffer buf(5 * W(N) + W(2 * m + 2) + W(2 * std::size_t(N)) + W(N + 33) + W(N + 1));
 u32* cur= buf.get();
 auto take= [&](std::size_t x) {
  u32* r= cur;
  cur+= W(x);
  return r;
 };
 // ky: 鍵、fst: 最初の隣接、hd: 残りの隣接の連結リストの先頭、nx[2 e] と nx[2 e + 1]: 辺 e の行き先と次の辺。
 // 辺 m は剥がしの番兵で、行き先が N (deg[N] は 0 から減らすので 0 に戻らず、積まれない)、次が NIL。
 u32 *ky= take(N), *fst= take(N), *hd= take(N), *nx= take(2 * m + 2), *st= take(N), *bs= take(N);
 F* fs= reinterpret_cast<F*>(take(2 * std::size_t(N)));
 u32 *q= take(N + 33), *deg= take(N + 1);
 std::memset(ky, 0, N * sizeof(u32)), std::memset(fst, 0xff, N * sizeof(u32)), std::memset(hd, 0xff, N * sizeof(u32));
 std::memset(deg, 0, (N + 1) * sizeof(u32));
 nx[2 * m]= N, nx[2 * m + 1]= NIL;
 // その頂点の最初の辺かどうかは、分岐にせずマスクで選ぶ。最初の辺でも nx の辺の場所には書く (使われない)。
 for(std::size_t i= 0; i < m; ++i) {
  const u32 a= edges[i][0], b= edges[i][1], f= fst[a], h= hd[a], k= -u32(f == NIL);
  nx[2 * i]= b, nx[2 * i + 1]= h;
  fst[a]= (b & k) | (f & ~k), hd[a]= (h & k) | (u32(i) & ~k);
  ++deg[b];
 }
 // 入次数 0 の頂点を剥がす。剥がした頂点の鍵は NIL - (剥がした順)。積むかどうかは分岐にせず、毎回 q の末尾に書いて、
 // 入次数が 0 になったときだけ末尾を進める。16 個先の頂点の記録と、8 個先の頂点の行き先の入次数と残りの隣接を先読みする。
 u32 qt= 0;
 for(u32 v= 0; v < N; ++v) q[qt]= v, qt+= deg[v] == 0;
 for(u32 qh= 0; qh < qt; ++qh) {
  const u32 a= q[std::min(qh + 16, qt - 1)];
  __builtin_prefetch(fst + a), __builtin_prefetch(hd + a), __builtin_prefetch(ky + a);
  const u32 b= q[std::min(qh + 8, qt - 1)], hb= hd[b];
  __builtin_prefetch(deg + std::min(fst[b], N)), __builtin_prefetch(nx + 2 * (hb == NIL ? 0 : hb));
  const u32 v= q[qh];
  ky[v]= NIL - qh;
  u32 w= std::min(fst[v], N);
  q[qt]= w, qt+= --deg[w] == 0;
  for(u32 e= hd[v]; e != NIL; e= nx[2 * e + 1]) w= nx[2 * e], q[qt]= w, qt+= --deg[w] == 0;
 }
 const u32 P= qt;
 // path-based。鍵は st に載っている間はその位置 (1 から)、成分に入ったら成分の番号 (2N から下がる。どの位置より大きい)。
 // bs は根の候補の位置のスタック。nf は、いまの頂点の最初の隣接のうち、まだ見ていないもの。
 u32 top= 0, bt= 0, cc= 2 * N;
 for(u32 r= 0; r < N; ++r) {
  if(ky[r]) continue;
  u32 v= r, nf= fst[r], ce= hd[r];
  F* sp= fs;
  st[top++]= v, ky[v]= top, bs[bt++]= top;
  for(;;) {
   u32 w;
   if(nf != NIL) w= nf, nf= NIL;
   else if(ce != NIL) w= nx[2 * ce], ce= nx[2 * ce + 1];
   else {
    const u32 p= ky[v];
    if(p == bs[bt - 1]) {
     --bt;
     while(top >= p) ky[st[--top]]= cc;
     --cc;
    }
    if(sp == fs) break;
    --sp, v= sp->v, ce= sp->e;
    continue;
   }
   const u32 kw= ky[w], fw= fst[w];
   if(!kw) *sp++= {v, ce}, v= w, nf= fw, ce= hd[w], st[top++]= v, ky[v]= top, bs[bt++]= top;
   else
    while(kw < bs[bt - 1]) --bt;
  }
 }
 // 剥がした頂点は剥がした順、残りは剥がした数のあとに DFS で見つけた成分を並べる (DFS は成分をトポロジカル順の逆に、
 // 大きい番号から振るので、番号の小さい順がトポロジカル順になる)。
 StronglyConnectedComponents res{int(P + (2 * N - cc)), std::vector<int>(N)};
 const u32 base= cc + 1;
 for(u32 v= 0; v < N; ++v) {
  const u32 x= ky[v];
  res.comp[v]= int(x > 2 * N ? NIL - x : x - base + P);
 }
 return res;
}
