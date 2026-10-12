#pragma once
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 両端に文字を足し引きできる列の回文木 (Wang、Yang、Zhu の double-ended eertree、arXiv 2210.02292 の表面を記録する方法)。回文の頂点
// ごとに、長さ、suffix link、親 (両端を除いた回文)、末尾の文字、cnt (最長の回文接尾辞として現れる位置の数)、suffix link の木での子の数、
// 文字ごとの遷移と direct link (自分の真の回文接尾辞のうち、自分の中で直前が c である最長のもの。Rubinchik と Shur) を持つ。列の各位置
// には、そこから始まる表面とそこで終わる表面 (ほかの回文の出現の接頭辞にも接尾辞にもならない回文の出現) を持つ。足すときは、最長の回文
// 接尾辞 (接頭辞) から direct link を 1 回引いて伸ばせる回文を見つけ、無ければ頂点を作り、direct link は suffix link の先の表を写して
// 1 か所を直す。除くときは、その回文が一意 (子が無く cnt が 1) なら頂点を消す。どの操作にも繰り返しが無く、O(K) で済む。表と列は
// 1 つの 0 で埋まった領域 (2 MB 以上なら、Linux では 2 MB 境界の mmap に置いて huge page を頼む) に置き、EVEN を 0 番にして、表面の
// 初めの値 (EVEN) を埋める手間を省く。procon-judge の yosupo-palindromes-in-deque で書き比べた dq_tk2 を写したもの。
namespace palindrome_deque_internal {
// bytes byte の 0 で埋まった領域。2 MB 以上なら、Linux では 2 MB 境界の mmap に置いて huge page を頼む。
class Zeroed {
 void* p= nullptr;
#ifdef __linux__
 std::size_t mapped= 0;
#endif
public:
 explicit Zeroed(std::size_t bytes) {
#ifdef __linux__
  constexpr std::size_t H= std::size_t(1) << 21;
  if(bytes >= H) {
   const std::size_t b= (bytes + 4095) & ~std::size_t(4095);
   void* m= mmap(nullptr, b + H, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
   if(m != MAP_FAILED) {
    char *q= static_cast<char*>(m), *a= reinterpret_cast<char*>((reinterpret_cast<std::uintptr_t>(q) + H - 1) & ~(H - 1));
    if(a != q) munmap(q, a - q);
    munmap(a + b, q + H - a);
    madvise(a, b, MADV_HUGEPAGE);
    p= a, mapped= b;
    return;
   }
  }
#endif
  p= std::calloc(bytes ? bytes : 1, 1);
 }
 ~Zeroed() {
#ifdef __linux__
  if(mapped) {
   munmap(p, mapped);
   return;
  }
#endif
  std::free(p);
 }
 Zeroed(const Zeroed&)= delete;
 Zeroed& operator=(const Zeroed&)= delete;
 void* get() const { return p; }
};
struct Node {
 int len, link, par, cnt, linkcnt, ch;
};
}
// 文字は 0 以上 K 未満の整数。n は足す操作 (push_front と push_back) の回数の上限。
template <int K= 26> class PalindromeDeque {
 static_assert(1 <= K && K <= 256);
 using Node= palindrome_deque_internal::Node;
 static constexpr int EVEN= 0, ODD= 1;  // 子の 0 は「無し」(子は長さ 1 以上なので EVEN にならない)
 static std::size_t bytes(int n) {
  const std::size_t cap= std::size_t(n) + 2, P= 2 * std::size_t(n) + 3;
  return cap * K * 2 * sizeof(int) + cap * sizeof(Node) + 2 * P * sizeof(int) + P;
 }
 palindrome_deque_internal::Zeroed buf;
 int *nxt, *dlt, *pre, *suf;  // 遷移、direct link、各位置から始まる表面、各位置で終わる表面
 Node* nd;
 unsigned char* data;
 int start, end, m= 2, live= 0;
 int* nx(int v) const { return nxt + std::size_t(v) * K; }
 int* dl(int v) const { return dlt + std::size_t(v) * K; }
 // u から c で伸ばした回文の頂点を返す (無ければ作る)。ext(w) は、回文 w の出現の外側の文字が c か。pc(L) は、新しい頂点の中で長さ L
 // の回文の接尾辞 (接頭辞) の外側に隣る文字。
 template <class Ext, class Pc> int grow(int u, int c, Ext ext, Pc pc) {
  int& e= nx(u)[c];
  if(e) return e;
  const int x= m++;
  int lx= EVEN;
  if(u != ODD) {
   const int lu= nd[u].link;
   const int w= lu == ODD || ext(lu) ? lu : dl(lu)[c];
   lx= nx(w)[c];
  }
  nd[x]= {nd[u].len + 2, lx, u, 0, 0, c};
  ++nd[lx].linkcnt;
  std::memcpy(dl(x), dl(lx), sizeof(int) * K);
  dl(x)[pc(nd[lx].len)]= lx;
  ++live;
  return e= x;
 }
public:
 explicit PalindromeDeque(int n): buf(bytes(n)) {
  const std::size_t cap= std::size_t(n) + 2, P= 2 * std::size_t(n) + 3;
  nxt= static_cast<int*>(buf.get()), dlt= nxt + cap * K;
  nd= reinterpret_cast<Node*>(dlt + cap * K);
  pre= reinterpret_cast<int*>(nd + cap), suf= pre + P;
  data= reinterpret_cast<unsigned char*>(suf + P);
  nd[ODD]= {-1, ODD, ODD, 0, 0, 0}, nd[EVEN]= {0, ODD, ODD, 0, 0, 0};
  for(int c= 0; c < K; ++c) dl(ODD)[c]= dl(EVEN)[c]= ODD;
  start= end= n + 1;
 }
 PalindromeDeque(const PalindromeDeque&)= delete;
 PalindromeDeque& operator=(const PalindromeDeque&)= delete;
 void push_back(int c) {
  const int e= end;
  const int v= start == end ? EVEN : suf[e - 1];
  auto ext= [&](int w) {
   const int p= e - 1 - nd[w].len;
   return p >= start && data[p] == c;
  };
  const int u= ext(v) ? v : dl(v)[c];
  const int x= grow(u, c, ext, [&](int L) { return L == 0 ? c : int(data[e - L]); });
  data[e]= (unsigned char)c, end= e + 1;
  pre[e]= EVEN, pre[e + 1 - nd[x].len]= x, suf[e]= x;
  const int lx= nd[x].link, lq= e - nd[x].len + nd[lx].len;
  if(nd[lx].len >= 1 && suf[lq] == lx) suf[lq]= EVEN;
  ++nd[x].cnt;
 }
 void push_front(int c) {
  const int b= start - 1;
  const int v= start == end ? EVEN : pre[start];
  auto ext= [&](int w) {
   const int p= start + nd[w].len;
   return p < end && data[p] == c;
  };
  const int u= ext(v) ? v : dl(v)[c];
  const int x= grow(u, c, ext, [&](int L) { return L == 0 ? c : int(data[b + L]); });
  data[b]= (unsigned char)c, start= b;
  pre[b]= x, suf[b]= EVEN, suf[b + nd[x].len - 1]= x;
  const int lx= nd[x].link, lq= b + nd[x].len - nd[lx].len;
  if(nd[lx].len >= 1 && pre[lq] == lx) pre[lq]= EVEN;
  ++nd[x].cnt;
 }
 // 空でない列で呼ぶ。
 void pop_back() {
  const int e= end - 1, v= suf[e];
  if(nd[v].linkcnt == 0 && nd[v].cnt == 1) nx(nd[v].par)[nd[v].ch]= 0, --nd[nd[v].link].linkcnt, --live;
  --nd[v].cnt;
  const int lv= nd[v].link, st= e - nd[v].len + 1;
  if(nd[lv].len >= 1 && nd[suf[st + nd[lv].len - 1]].len < nd[lv].len) suf[st + nd[lv].len - 1]= lv, pre[st]= lv;
  else pre[st]= EVEN;
  end= e;
 }
 // 空でない列で呼ぶ。
 void pop_front() {
  const int b= start, v= pre[b];
  if(nd[v].linkcnt == 0 && nd[v].cnt == 1) nx(nd[v].par)[nd[v].ch]= 0, --nd[nd[v].link].linkcnt, --live;
  --nd[v].cnt;
  const int lv= nd[v].link, en= b + nd[v].len - 1;
  if(nd[lv].len >= 1 && nd[pre[en - nd[lv].len + 1]].len < nd[lv].len) pre[en - nd[lv].len + 1]= lv, suf[en]= lv;
  else suf[en]= EVEN;
  start= b + 1;
 }
 int size() const { return end - start; }
 // 列に現れる異なる (空でない) 回文の数。
 int distinct() const { return start == end ? 0 : live; }
 // 最長の回文接頭辞の長さ (空の列なら 0)。
 int longest_prefix() const { return start == end ? 0 : nd[pre[start]].len; }
 // 最長の回文接尾辞の長さ (空の列なら 0)。
 int longest_suffix() const { return start == end ? 0 : nd[suf[end - 1]].len; }
};
