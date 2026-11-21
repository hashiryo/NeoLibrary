#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 有向グラフの支配木。s を根とし、返す列の v 番目は頂点 v の直接支配点 (v = s なら s、s から届かなければ -1)。自己ループと多重辺が
// あってもよい。SEMI-NCA で求める。s から反復の DFS で届く頂点に行きがけ順の番号を振り、以後は番号の上で計算する。DFS の隣接は、
// 頂点ごとに最初の出辺の行き先と残りの出辺の連結リストの先頭を別々の配列に置き、辺の列を 1 回なめて組む (DFS が下りる 1 歩が
// 読み込み 1 段で済む)。届く頂点どうしの辺だけで逆向きの CSR を番号の上で組み直す (1 回目の走査で引いた端の番号の組を残し、2 回目の
// 走査は引き直さない)。準支配点 (sdom) は番号の大きい順に、逆向きの辺の始点の eval の最小で求める。森の各点には森の親と道の上の
// sdom の最小の値を並べて持ち、eval は道を 1 回登りながら各点の親を祖父に付け替える経路半減で、最小の値を返す。直接支配点は、番号の
// 小さい順に、木の親から sdom より大きいあいだ直接支配点をたどって決める。作業領域は 2 MB 以上なら、Linux では 2 MB 境界の mmap に
// 置いて MADV_HUGEPAGE を頼む。procon-judge の yosupo-dominator-tree で書き比べた snca_halve_hp を写したもの。
namespace dominator_tree_internal {
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
// v に n 個分の領域を取り、Linux ではその領域のページに huge page を頼む (領域の内側にある 2 MB の区切りだけが載る)。
template <class T> inline void reserve_huge(std::vector<T>& v, std::size_t n) {
 v.reserve(n);
#ifdef __linux__
 const std::uintptr_t b= reinterpret_cast<std::uintptr_t>(v.data()) & ~std::uintptr_t(4095), e= (reinterpret_cast<std::uintptr_t>(v.data() + n) + 4095) & ~std::uintptr_t(4095);
 if(b < e) madvise(reinterpret_cast<void*>(b), e - b, MADV_HUGEPAGE);
#endif
}
struct AL {
 u32 anc, lab;  // 森の親と、道の上の sdom の最小の値
};
// v から森の根の手前までの sdom の最小。各点の値は「その点から親の手前までの最小」のままにして、親を祖父に付け替える。
[[gnu::always_inline]] inline u32 eval(u32 v, AL* al) {
 if(al[v].anc == NIL) return v;  // 根 (まだ処理していない頂点と s) の sdom は番号そのもの
 u32 x= v, res= NIL;
 for(;;) {
  const u32 a= al[x].anc;
  if(al[a].anc == NIL) return std::min(res, al[x].lab);
  if(al[a].lab < al[x].lab) al[x].lab= al[a].lab;
  al[x].anc= al[a].anc, res= std::min(res, al[x].lab), x= al[x].anc;
  if(al[x].anc == NIL) return res;
 }
}
}
inline std::vector<int> dominator_tree(int n, int s, const std::vector<std::array<int, 2>>& edges) {
 using namespace dominator_tree_internal;
 const u32 N= n;
 const std::size_t m= edges.size();
 auto W= [](std::size_t x) { return (x + 15) & ~std::size_t(15); };
 Buffer buf(9 * W(N + 2) + 2 * W(2 * m + 2) + W(m + 1) + W(3 * std::size_t(N)) + W(2 * std::size_t(N) + 2));
 u32* cur= buf.get();
 auto take= [&](std::size_t x) {
  u32* r= cur;
  cur+= W(x);
  return r;
 };
 // num[v] は頂点 v の番号 (届かなければ NIL)、vert[i] は番号 i の頂点、par[i] は木の親の番号、idom[i] は直接支配点の番号。
 u32 *ro= take(N + 2), *ra= take(m + 1), *num= take(N + 1), *vert= take(N + 1), *idom= take(N + 1), *par= take(N + 1), *fs= take(3 * std::size_t(N));
 u32 *fst= take(N + 1), *hd= take(N + 1), *nx= take(2 * m + 2);
 // 最初の行き先 fst と、残りの出辺の連結リストの先頭 hd。その頂点の最初の辺かどうかはマスクで選ぶ。
 std::memset(fst, 0xff, N * sizeof(u32)), std::memset(hd, 0xff, N * sizeof(u32)), std::memset(num, 0xff, N * sizeof(u32));
 for(std::size_t j= 0; j < m; ++j) {
  const u32 a= edges[j][0], b= edges[j][1], f= fst[a], h= hd[a], k= -u32(f == NIL);
  nx[2 * j]= b, nx[2 * j + 1]= h;
  fst[a]= (b & k) | (f & ~k), hd[a]= (h & k) | (u32(j) & ~k);
 }
 // 反復の DFS。フレームは {頂点, 残りの連結リストのうちまだ見ていない先頭, 頂点の番号}。nf はいまの頂点の最初の行き先のうち
 // まだ見ていないもの。
 u32 R= 0;
 {
  u32* sp= fs;
  num[s]= R, par[R]= 0, vert[R++]= s;
  u32 v= s, iv= 0, nf= fst[s], ce= hd[s];
  for(;;) {
   u32 w;
   if(nf != NIL) w= nf, nf= NIL;
   else if(ce != NIL) w= nx[2 * ce], ce= nx[2 * ce + 1];
   else {
    if(sp == fs) break;
    sp-= 3, v= sp[0], ce= sp[1], iv= sp[2];
    continue;
   }
   const u32 nw= num[w], fw= fst[w], hw= hd[w];
   if(nw != NIL) continue;
   sp[0]= v, sp[1]= ce, sp[2]= iv, sp+= 3;
   num[w]= R, par[R]= iv, vert[R]= w, iv= R++, v= w, nf= fw, ce= hw;
  }
 }
 // 逆向きの CSR を番号の上で組む。届かない頂点が端の辺は番兵 R に数えて末尾に捨てる。
 std::memset(ro, 0, (R + 1) * sizeof(u32));
 u32* tp= nx;  // 隣接の連結リストは使い終わった
 for(std::size_t j= 0; j < m; ++j) {
  const u32 a= num[edges[j][0]], b= num[edges[j][1]], ok= (a != NIL) & (b != NIL);
  tp[2 * j]= a, tp[2 * j + 1]= ok ? b : R;
  ++ro[ok ? b : R];
 }
 for(u32 k= 1; k <= R; ++k) ro[k]+= ro[k - 1];
 for(std::size_t j= m; j--;) ra[--ro[tp[2 * j + 1]]]= tp[2 * j];
 // sdom。sdom は値を fst に置く (fst は使い終わった)。
 AL* al= reinterpret_cast<AL*>(take(2 * std::size_t(R) + 2));
 u32* sdom= fst;
 for(u32 k= 0; k < R; ++k) al[k]= {NIL, k};
 for(u32 w= R; --w >= 1;) {
  u32 sw= w;
  for(u32 j= ro[w]; j < ro[w + 1]; ++j) sw= std::min(sw, eval(ra[j], al));
  sdom[w]= sw, al[w]= {par[w], sw};
 }
 idom[0]= 0;
 for(u32 w= 1; w < R; ++w) {
  u32 d= par[w];
  while(d > sdom[w]) d= idom[d];
  idom[w]= d;
 }
 std::vector<int> res;
 reserve_huge(res, N);
 res.assign(N, -1);
 for(u32 k= 0; k < R; ++k) res[vert[k]]= int(vert[idom[k]]);
 return res;
}
