#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 2-SAT。n 個の変数と節 a ∨ b の列を受け、すべての節を満たす値を求める。リテラルは、変数 i (0 から n - 1) の肯定が i、否定が
// ~i (= -1 - i)。満たせるなら各変数の値 (0 か 1) の列を、満たせなければ空の列を返す。
// 含意グラフ (節 a ∨ b から辺 ¬a → b と ¬b → a) の強連結成分で解く。辺の列は作らず、リテラル x の出現の並び occ[x] (x を含む
// 節のもう一方のリテラル) の CSR を 1 つだけ組む。含意グラフでは、x から出る辺の行き先が occ[¬x]、x の入次数が |occ[x]| なので、
// 節を 1 回なめて出現を数えれば、CSR の区切りと剥がしの入次数の両方が手に入る (数えた配列をそのまま入次数に使う)。リテラル
// ごとの配列は鍵の 4 byte だけで、10^6 リテラルでも作業領域を小さく保つ。入次数 0 のリテラルを待ち行列で剥がしてから、残りを
// path-based (Gabow) の反復の DFS で分ける。成分の番号はトポロジカル順で、x が真になるのは x の成分が ¬x の成分より後ろの
// とき。作業領域は 2 MB 以上なら、Linux では 2 MB 境界の mmap に置いて MADV_HUGEPAGE を頼む。procon-judge の yosupo-two-sat で
// 書き比べた occ_gabow_trim_hp を写したもの。
namespace two_sat_internal {
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
// リテラルを頂点にする。肯定の i は 2 i、否定の ~i は 2 i + 1。
inline u32 vertex(int l) { return l >= 0 ? 2 * u32(l) : 2 * u32(~l) + 1; }
}
inline std::vector<int> two_sat(int n, const std::vector<std::array<int, 2>>& clauses) {
 using namespace two_sat_internal;
 struct F {
  u32 v, i, e;
 };
 const u32 N= 2 * u32(n);
 const std::size_t m= clauses.size();
 auto W= [](std::size_t x) { return (x + 15) & ~std::size_t(15); };
 Buffer buf(4 * W(N + 1) + W(N) + W(N + 33) + W(2 * m + 1) + W(3 * std::size_t(N)));
 u32* cur= buf.get();
 auto take= [&](std::size_t x) {
  u32* r= cur;
  cur+= W(x);
  return r;
 };
 // deg: 出現の数 = 入次数 (deg[N] は剥がしの番兵)。off と adj: 出現の CSR (adj は範囲の外を 1 語読んでよいように多めに取る)。
 u32 *deg= take(N + 1), *off= take(N + 1), *adj= take(2 * m + 1), *key= take(N), *q= take(N + 33), *st= take(N + 1), *bs= take(N + 1);
 F* fs= reinterpret_cast<F*>(take(3 * std::size_t(N)));
 std::memset(deg, 0, (N + 1) * sizeof(u32));
 for(auto& c: clauses) ++deg[vertex(c[0])], ++deg[vertex(c[1])];
 // off[x] は x までの出現の数の累積。出現を後ろから置いて、off[x] を x の始まりにする。
 if(N) {
  off[0]= deg[0];
  for(u32 x= 1; x < N; ++x) off[x]= off[x - 1] + deg[x];
 }
 off[N]= 2 * m;
 for(std::size_t j= m; j--;) {
  const u32 a= vertex(clauses[j][0]), b= vertex(clauses[j][1]);
  adj[--off[a]]= b, adj[--off[b]]= a;
 }
 std::memset(key, 0, N * sizeof(u32));
 // 入次数 0 の頂点を剥がす。剥がした頂点の鍵は NIL - (剥がした順)。x を剥がすと、x の出辺の行き先 occ[x ^ 1] の入次数が減る。
 // 積むかどうかは分岐にせず、毎回 q の末尾に書いて、入次数が 0 になったときだけ末尾を進める。
 u32 qt= 0;
 for(u32 x= 0; x < N; ++x) q[qt]= x, qt+= deg[x] == 0;
 for(u32 qh= 0; qh < qt; ++qh) {
  const u32 x= q[qh];
  key[x]= NIL - qh;
  for(u32 i= off[x ^ 1], e= off[(x ^ 1) + 1]; i < e; ++i) {
   const u32 z= adj[i];
   q[qt]= z, qt+= --deg[z] == 0;
  }
 }
 const u32 P= qt;
 // path-based。鍵は st に載っている間はその位置 (1 から)、成分に入ったら成分の番号 (2N から下がる。どの位置より大きい)。
 // 頂点 v の出辺の範囲は off[v ^ 1] から off[(v ^ 1) + 1]。
 u32 top= 0, bt= 0, cc= 2 * N;
 for(u32 r= 0; r < N; ++r) {
  if(key[r]) continue;
  u32 v= r, i= off[r ^ 1], e= off[(r ^ 1) + 1];
  F* sp= fs;
  st[top++]= v, key[v]= top, bs[bt++]= top;
  for(;;) {
   if(i != e) {
    const u32 w= adj[i++], kw= key[w];
    if(!kw) *sp++= {v, i, e}, v= w, i= off[w ^ 1], e= off[(w ^ 1) + 1], st[top++]= v, key[v]= top, bs[bt++]= top;
    else
     while(kw < bs[bt - 1]) --bt;
   } else {
    const u32 p= key[v];
    if(p == bs[bt - 1]) {
     --bt;
     while(top >= p) key[st[--top]]= cc;
     --cc;
    }
    if(sp == fs) break;
    --sp, v= sp->v, i= sp->i, e= sp->e;
   }
  }
 }
 // 成分の番号 (トポロジカル順)。剥がした頂点は剥がした順、残りは剥がした数のあとに DFS で見つけた成分を並べる。
 const u32 base= cc + 1;
 auto comp= [&](u32 v) {
  const u32 x= key[v];
  return x > 2 * N ? NIL - x : x - base + P;
 };
 std::vector<int> val(n);
 for(int i= 0; i < n; ++i) {
  const u32 a= comp(2 * i), b= comp(2 * i + 1);
  if(a == b) return {};
  val[i]= a > b;
 }
 return val;
}
