#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <new>
#include <string>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 接尾辞オートマトン。文字を 1 つずつ末尾に足して作り、状態ごとに suffix link、最長の長さ、遷移を持つ。遷移の持ち方は遷移の数で変える
// (adaptive radix tree の Node4 と Node256 の考え方)。状態は 32 byte で、4 本までの遷移を、文字 4 つを詰めた 1 語と行き先 4 つで持ち、
// 文字は語の中の byte の一致で探す。5 本目の遷移が来たら、K 個の欄の表を別の列に取ってそこへ移す。状態の遷移はほとんど 1〜2 本なので、
// K 個の欄を全部の状態に持つ形より表がずっと小さい。異なる部分文字列の数は、文字を足すたびに len[cur] - len[link[cur]] だけ増える。
// 列は 2 MB 以上なら Linux で huge page を頼む。procon-judge の self-string-online-distinct-substrings で書き比べた sam_art を写したもの。
namespace suffix_automaton_internal {
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
template <int K= 26, char Base= 'a'> class SuffixAutomaton {
 static_assert(1 <= K && K <= 255);
 using u32= std::uint32_t;
 struct St {
  int link, len;
  u32 keys;  // 使っていない byte は 0xff
  int to[4];
  int big;  // K 個の欄の表の番号 (無ければ -1)
 };
 template <class T> using Vec= std::vector<T, suffix_automaton_internal::HugeAlloc<T>>;
 Vec<St> st;
 Vec<std::array<int, K>> bigs;
 int lst= 0;
 long long total= 0;
 // keys の中で byte c の位置 (無ければ -1)
 static int find(u32 keys, int c) {
  const u32 x= keys ^ (0x01010101u * u32(c)), z= (x - 0x01010101u) & ~x & 0x80808080u;
  return z ? __builtin_ctz(z) >> 3 : -1;
 }
 // 状態 v の文字 c の遷移の先 (無ければ 0。根は遷移の先にならない)
 int get(int v, int c) const {
  const St& s= st[v];
  if(s.big >= 0) return bigs[s.big][c];
  const int k= find(s.keys, c);
  return k >= 0 ? s.to[k] : 0;
 }
 void set(int v, int c, int t) {
  St& s= st[v];
  if(s.big >= 0) {
   bigs[s.big][c]= t;
   return;
  }
  int k= find(s.keys, c);
  if(k < 0) k= find(s.keys, 0xff);
  if(k >= 0) {
   s.keys= (s.keys & ~(0xffu << (8 * k))) | (u32(c) << (8 * k)), s.to[k]= t;
   return;
  }
  std::array<int, K> b{};
  for(int j= 0; j < 4; ++j) b[(s.keys >> (8 * j)) & 0xff]= s.to[j];
  b[c]= t;
  s.big= bigs.size();
  bigs.push_back(b);
 }
public:
 SuffixAutomaton() { st.push_back(St{-1, 0, 0xffffffffu, {0, 0, 0, 0}, -1}); }
 explicit SuffixAutomaton(const std::string& s): SuffixAutomaton() {
  st.reserve(2 * s.size() + 1);
  for(char c: s) push(c);
 }
 // 文字 ch を末尾に足す。
 void push(char ch) {
  const int c= (unsigned char)(ch - Base), cur= st.size();
  st.push_back(St{0, st[lst].len + 1, 0xffffffffu, {0, 0, 0, 0}, -1});
  int p= lst;
  while(p != -1 && !get(p, c)) set(p, c, cur), p= st[p].link;
  if(p != -1) {
   const int q= get(p, c);
   if(st[p].len + 1 == st[q].len) {
    st[cur].link= q;
   } else {
    const int cl= st.size();
    St x= st[q];
    x.len= st[p].len + 1;
    if(x.big >= 0) {
     const std::array<int, K> b= bigs[x.big];
     x.big= bigs.size();
     bigs.push_back(b);
    }
    st.push_back(x);
    while(p != -1 && get(p, c) == q) set(p, c, cl), p= st[p].link;
    st[q].link= st[cur].link= cl;
   }
  }
  lst= cur;
  total+= st[cur].len - st[st[cur].link].len;
 }
 // 状態の数。状態 0 は空の文字列 (根)。
 int size() const { return st.size(); }
 // 今の文字列全体を受理する状態
 int last() const { return lst; }
 // 状態 v の文字列の最長の長さ
 int len(int v) const { return st[v].len; }
 // 状態 v の suffix link (根は -1)
 int link(int v) const { return st[v].link; }
 // 状態 v から文字 ch で移る先 (無ければ -1)
 int next(int v, char ch) const {
  const int t= get(v, (unsigned char)(ch - Base));
  return t ? t : -1;
 }
 // 今の文字列の、異なる空でない部分文字列の数
 long long distinct_substrings() const { return total; }
};
