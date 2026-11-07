#pragma once
#include <cstddef>
#include <cstdint>
#include <memory>
#include <new>
#include <string>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 文字列の多重集合を持つトライ。頂点ごとに、その頂点を通る要素の数と、その頂点で終わる要素の数を持つ。頂点の大きさは子の数で変える
// (adaptive radix tree の Node4 と Node256 の考え方)。子が 4 つまでの頂点は 32 byte で、子の文字 4 つを 1 語に詰めて、語の中の byte の
// 一致で探す。5 つ目の子が来たら K 個の子の欄を持つ頂点に替えて、親の欄を書き換える。根は初めから K 個の欄。深い頂点はほとんど子が
// 1 つなので、K 個の欄を全部の頂点に持つ形より表がずっと小さい。頂点の列は 2 MB 以上なら Linux で huge page を頼む。
// procon-judge の self-string-prefix-count で書き比べた trie_art を写したもの。
namespace trie_internal {
template <class T> struct HugeAlloc {
 using value_type= T;
 static constexpr std::size_t H= std::size_t(1) << 21, P= 4096;
 HugeAlloc()= default;
 template <class U> HugeAlloc(const HugeAlloc<U>&) {}
 T* allocate(std::size_t n) {
#ifdef __linux__
  if(n * sizeof(T) >= H) {
   const std::size_t bytes= (n * sizeof(T) + P - 1) & ~(P - 1);
   void* m= mmap(nullptr, bytes + H, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
   if(m == MAP_FAILED) throw std::bad_alloc();
   char *p= static_cast<char*>(m), *a= reinterpret_cast<char*>((reinterpret_cast<std::uintptr_t>(p) + H - 1) & ~(H - 1));
   if(a != p) munmap(p, a - p);
   munmap(a + bytes, p + H - a);
   madvise(a, bytes, MADV_HUGEPAGE);
   return reinterpret_cast<T*>(a);
  }
#endif
  return std::allocator<T>().allocate(n);
 }
 void deallocate(T* p, std::size_t n) {
#ifdef __linux__
  if(n * sizeof(T) >= H) {
   munmap(p, (n * sizeof(T) + P - 1) & ~(P - 1));
   return;
  }
#endif
  std::allocator<T>().deallocate(p, n);
 }
 template <class U> bool operator==(const HugeAlloc<U>&) const { return true; }
 template <class U> bool operator!=(const HugeAlloc<U>&) const { return false; }
};
}
// 文字は Base から Base + K - 1 まで (K <= 255)。
template <int K= 26, char Base= 'a'> class Trie {
 static_assert(1 <= K && K <= 255);
 using u32= std::uint32_t;
 struct Small {
  u32 keys;  // 使っていない byte は 0xff
  int cnt, end, ch[4];
  u32 num;
 };
 struct Big {
  int ch[K], cnt, end;
 };
 template <class T> using Vec= std::vector<T, trie_internal::HugeAlloc<T>>;
 // 頂点の参照は、正なら小さい頂点、負なら大きい頂点 (-1 - 番号)、0 は子なし。根は大きい頂点の 0 番。
 Vec<Small> sm;
 Vec<Big> bg;
 int& cnt(int r) { return r > 0 ? sm[r].cnt : bg[-r - 1].cnt; }
 int& end(int r) { return r > 0 ? sm[r].end : bg[-r - 1].end; }
 int child(int r, int c) const {
  if(r < 0) return bg[-r - 1].ch[c];
  const Small& v= sm[r];
  const u32 x= v.keys ^ (0x01010101u * u32(c)), z= (x - 0x01010101u) & ~x & 0x80808080u;
  return z ? v.ch[__builtin_ctz(z) >> 3] : 0;
 }
 // 頂点 r に文字 c の子を足して、その参照を返す。r が 4 つの子を持つ小さい頂点なら大きい頂点に替え、r と、r を指す親の欄 (親の参照 pr
 // と欄の番号 pk) を書き換える。
 int add_child(int& r, int c, int pr, int pk) {
  const int w= sm.size();
  sm.push_back(Small{0xffffffffu, 0, 0, {0, 0, 0, 0}, 0});
  if(r < 0) return bg[-r - 1].ch[c]= w;
  Small& v= sm[r];
  if(v.num < 4) {
   v.keys= (v.keys & ~(0xffu << (8 * v.num))) | (u32(c) << (8 * v.num));
   v.ch[v.num++]= w;
   return w;
  }
  Big b{};
  for(int k= 0; k < 4; ++k) b.ch[(v.keys >> (8 * k)) & 0xff]= v.ch[k];
  b.cnt= v.cnt, b.end= v.end, b.ch[c]= w;
  bg.push_back(b);
  r= -int(bg.size());
  (pr < 0 ? bg[-pr - 1].ch[pk] : sm[pr].ch[pk])= r;
  return w;
 }
 void add(const std::string& s, int d) {
  int r= -1, pr= 0, pk= 0;
  bg[0].cnt+= d;
  for(char ch: s) {
   const int c= (unsigned char)(ch - Base);
   int w= child(r, c);
   if(!w) w= add_child(r, c, pr, pk);
   pr= r;
   if(r < 0) pk= c;
   else {
    const Small& v= sm[r];
    pk= 0;
    while(v.ch[pk] != w) ++pk;
   }
   r= w;
   cnt(r)+= d;
  }
  end(r)+= d;
 }
 // 文字列 s の頂点 (無ければ 0)
 int walk(const std::string& s) const {
  int r= -1;
  for(char ch: s)
   if(!(r= child(r, (unsigned char)(ch - Base)))) return 0;
  return r;
 }
public:
 Trie(): sm(1), bg(1) {}
 // s を 1 つ足す。
 void insert(const std::string& s) { add(s, 1); }
 // s を 1 つ除く。s は今の多重集合にあること。
 void erase(const std::string& s) { add(s, -1); }
 // 要素の数
 int size() const { return bg[0].cnt; }
 // s と等しい要素の数
 int count(const std::string& s) const {
  const int r= walk(s);
  return r ? (r > 0 ? sm[r].end : bg[-r - 1].end) : 0;
 }
 // p で始まる要素の数。p が空なら要素の数。
 int count_prefix(const std::string& p) const {
  const int r= walk(p);
  return r ? (r > 0 ? sm[r].cnt : bg[-r - 1].cnt) : 0;
 }
};
