#pragma once
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 文字列の列のトライと suffix link を求める。2 つの作りを入力の形で選ぶ。(a) 遷移の表を幅 k (使われている文字の種類) で全部持ち、幅優先
// の順に、取り出した頂点の子の欄を待ち行列に積んでから、空いた欄を suffix link の行で埋める (子の有無で分岐しないよう、子の無い欄でも
// 待ち行列の末尾に書いて末尾は子のときだけ進める)。(b) 1 つ目の子を頂点に直接持ち、枝分かれする頂点にだけ幅 k の行を割り当てる。行の
// 無い頂点は 0 の行を指して子を分岐なしで引き、行のある頂点の行は取り出したときに遷移で埋め、suffix link をたどる途中で行のある頂点に
// 着いたらその欄で止まる。(a) は suffix link の先が散らばる入力 (でたらめな文字列を多く並べたもの) で速く、(b) はほぼ 1 本道のトライ
// (周期的な文字列やフィボナッチ文字列) と文字の種類が多い入力で速く、表も小さい。文字の種類が 9 以上か、文字列が 2 本以下か、先頭から
// 最大 65536 個の長さ 12 の部分文字列の種類が少なければ (b) を使う。作業の配列は huge page を頼んだ 1 つの領域に置き、答えの列は領域を
// 取ってから huge page を頼む。procon-judge の yosupo-aho-corasick で書き比べた ac_mix5 を写したもの。
namespace aho_corasick_internal {
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
// v に n 個分の領域を取り (中身はまだ書かない)、Linux ではその領域のページに huge page を頼む。vector の領域は 2 MB 境界から始まるとは
// 限らないので、huge page に載るのは領域の内側にある 2 MB の区切りだけになる。
template <class T> inline void reserve_huge(std::vector<T>& v, std::size_t n) {
 v.reserve(n);
#ifdef __linux__
 const std::uintptr_t b= reinterpret_cast<std::uintptr_t>(v.data()) & ~std::uintptr_t(4095), e= (reinterpret_cast<std::uintptr_t>(v.data() + n) + 4095) & ~std::uintptr_t(4095);
 if(b < e) madvise(reinterpret_cast<void*>(b), e - b, MADV_HUGEPAGE);
#endif
}
struct Coded {
 std::array<unsigned char, 256> code;
 int k;
 std::size_t total;
};
// 現れる byte に小さい順に 0 から番号を振り、種類の数と文字の数の和を 1 回の走査で求める。
inline Coded remap(const std::vector<std::string>& ss) {
 Coded r;
 std::array<bool, 256> used{};
 std::size_t total= 0;
 for(const auto& s: ss) {
  total+= s.size();
  for(const char c: s) used[(unsigned char)c]= true;
 }
 int k= 0;
 for(int c= 0; c < 256; ++c) r.code[c]= used[c] ? (unsigned char)k++ : 0;
 r.k= k, r.total= total;
 return r;
}
// 長さ 12 の部分文字列の種類が少ないか (k <= 8 のとき)。1 文字 3 bit に詰めた 36 bit の値を 2^16 bit の表に散らし、立った bit が 1024
// 以下なら少ないとみる。長さ 12 以上の文字列が無ければ false。
inline bool repetitive(const std::vector<std::string>& ss, const Coded& cd) {
 std::vector<std::uint64_t> bits(1 << 10, 0);
 constexpr std::uint64_t MASK= (std::uint64_t(1) << 36) - 1;
 int seen= 0;
 for(const auto& s: ss) {
  if(s.size() < 12) continue;
  std::uint64_t x= 0;
  for(std::size_t i= 0; i < s.size() && seen < 65536; ++i) {
   x= (x << 3 | cd.code[(unsigned char)s[i]]) & MASK;
   if(i >= 11) {
    const std::uint64_t h= (x * 0x9E3779B97F4A7C15ull) >> 48;
    bits[h >> 6]|= std::uint64_t(1) << (h & 63), ++seen;
   }
  }
  if(seen >= 65536) break;
 }
 if(!seen) return false;
 int distinct= 0;
 for(const std::uint64_t b: bits) distinct+= std::popcount(b);
 return distinct <= 1024;
}
// (a) 遷移の表を全部持つ作り。
inline void full(const std::vector<std::string>& ss, const Coded& cd, std::vector<int>& par, std::vector<int>& link, std::vector<int>& word) {
 const auto& code= cd.code;
 const int k= cd.k;
 const std::size_t cap= cd.total + 1;
 Zeroed buf((cap * k + cap) * sizeof(int));  // 遷移の表 (cap 行、幅 k) と幅優先の待ち行列 (cap 個)
 int *ch= static_cast<int*>(buf.get()), *q= ch + cap * k;
 // push_back は末尾の位置を読み書きする待ちが要素ごとに連なるので、先に 0 で埋めて添字で書く
 reserve_huge(par, cap), reserve_huge(word, ss.size());
 par.resize(cap), word.resize(ss.size());
 int n= 1;
 for(std::size_t i= 0; i < ss.size(); ++i) {
  int v= 0;
  for(const char x: ss[i]) {
   int& w= ch[std::size_t(v) * k + code[(unsigned char)x]];
   if(!w) w= n++, par[w]= v;
   v= w;
  }
  word[i]= v;
 }
 par.resize(n);
 reserve_huge(link, n);
 link.assign(n, 0);
 int head= 0, tail= 0;
 for(int c= 0; c < k; ++c) {  // 根の子の suffix link は 0 のまま
  const int w= ch[c];
  q[tail]= w, tail+= w != 0;
 }
 // 頂点 v を取り出したとき、v の suffix link l は v より浅く、l の行は遷移で埋め終えている。v の子 w の suffix link は行 l の同じ文字の
 // 欄で、v の行の空いた欄には行 l の欄を写す。子の無い欄 (w = 0) でも link[0] と待ち行列の末尾に書く。link[0] は最後に 0 に戻す。
 while(head < tail) {
  const int v= q[head++];
  int* rv= ch + std::size_t(v) * k;
  const int* rl= ch + std::size_t(link[v]) * k;
  for(int c= 0; c < k; ++c) {
   const int w= rv[c], d= rl[c];
   link[w]= d;
   q[tail]= w, tail+= w != 0;
   rv[c]= w | (d & -int(w == 0));
  }
 }
 link[0]= 0;
}
// (b) 1 つ目の子を頂点に直接持つ作り。ow は 1 つ目の子、oc はその文字 (255 なら子が無い)、rb は行の先頭 (0 なら行が無く、0 の行を指す)、
// lk は作業用の suffix link。
struct Node {
 int ow, rb, lk;
 unsigned char oc, pad[3];
};
inline void hybrid(const std::vector<std::string>& ss, const Coded& cd, std::vector<int>& par, std::vector<int>& link, std::vector<int>& word) {
 const auto& code= cd.code;
 const int k= cd.k;
 const std::size_t cap= cd.total + 1;
 const std::size_t rows_cap= (ss.size() + 3) * std::size_t(k);  // 0 の行、根の行、枝分かれする頂点 (文字列の本数未満) の行
 Zeroed buf(cap * sizeof(Node) + (rows_cap + cap) * sizeof(int));
 Node* nd= static_cast<Node*>(buf.get());
 int *rows= reinterpret_cast<int*>(nd + cap), *q= rows + rows_cap;
 reserve_huge(par, cap), reserve_huge(word, ss.size());
 par.resize(cap), word.resize(ss.size());
 constexpr unsigned char NONE= 255;
 nd[0].oc= NONE, nd[0].rb= k;
 int n= 1, nrb= 2 * k;  // 次に割り当てる行の先頭
 for(std::size_t i= 0; i < ss.size(); ++i) {
  int v= 0;
  for(const char x: ss[i]) {
   const int c= code[(unsigned char)x];
   const int a= nd[v].ow, b= rows[nd[v].rb + c];
   int w= nd[v].oc == c ? a : b;  // 行が無ければ 0 の行を読むので、分岐なしで引ける
   if(!w) {
    w= n++;
    nd[w].oc= NONE, par[w]= v;
    if(nd[v].oc == NONE && !nd[v].rb) {
     nd[v].oc= (unsigned char)c, nd[v].ow= w;
    } else {
     if(!nd[v].rb) nd[v].rb= nrb, nrb+= k, rows[nd[v].rb + nd[v].oc]= nd[v].ow;
     rows[nd[v].rb + c]= w;
    }
   }
   v= w;
  }
  word[i]= v;
 }
 par.resize(n);
 // delta(u, c) は、u に行があればその欄、無ければ 1 つ目の子の文字が c ならその子、そうでなければ delta(lk[u], c)。行を持つ頂点 v は、
 // 取り出したときに行の空いた欄を delta(v, c) で埋めるので、v より浅い頂点の行はすべて遷移の行になっている。
 auto delta= [&](int u, int c) {
  while(true) {
   const Node& e= nd[u];
   if(e.rb) return rows[e.rb + c];
   if(e.oc == c) return e.ow;
   u= e.lk;
  }
 };
 int head= 0, tail= 0;
 for(int c= 0; c < k; ++c) {  // 根の子の suffix link は 0 のまま
  const int w= rows[k + c];
  q[tail]= w, tail+= w != 0;
 }
 while(head < tail) {
  const int v= q[head++], lv= nd[v].lk;
  if(const int rb= nd[v].rb) {
   int* r= rows + rb;
   // 子の無い欄 (w = 0) でも nd[0].lk と待ち行列の末尾に書く。根の lk は読まない
   if(const int rl= nd[lv].rb) {
    for(int c= 0; c < k; ++c) {
     const int w= r[c], d= rows[rl + c];
     nd[w].lk= d;
     q[tail]= w, tail+= w != 0;
     r[c]= w | (d & -int(w == 0));
    }
   } else {
    for(int c= 0; c < k; ++c) {
     const int w= r[c], d= delta(lv, c);
     nd[w].lk= d;
     q[tail]= w, tail+= w != 0;
     r[c]= w | (d & -int(w == 0));
    }
   }
  } else if(nd[v].oc != NONE) {
   const int w= nd[v].ow;
   nd[w].lk= delta(lv, nd[v].oc), q[tail++]= w;
  }
 }
 reserve_huge(link, n);
 link.resize(n);
 link[0]= 0;
 for(int v= 1; v < n; ++v) link[v]= nd[v].lk;
}
}
// 頂点 0 が根で、頂点の番号は文字列を 0 番から順にトライへ入れたときにできる順。par[v] は親、link[v] は suffix link (v の文字列の真の
// 接尾辞のうちトライにある最長のもの)、word[i] は i 番目の文字列の頂点。par[0] と link[0] は 0。
struct AhoCorasick {
 std::vector<int> par, link, word;
};
// 文字列の列 ss のトライと suffix link。文字は unsigned char として扱う。O(L k) (L は文字の数の和、k は文字の種類)。
inline AhoCorasick aho_corasick(const std::vector<std::string>& ss) {
 using namespace aho_corasick_internal;
 AhoCorasick r;
 const Coded cd= remap(ss);
 if(cd.k > 8 || ss.size() <= 2 || repetitive(ss, cd)) hybrid(ss, cd, r.par, r.link, r.word);
 else full(ss, cd, r.par, r.link, r.word);
 return r;
}
