#pragma once
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <numeric>
#include <string>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 回文木 (Rubinchik と Shur の eertree、2015) を左から作る。位置 i の文字 x について、直前の最長の回文接尾辞から suffix link をたどり、
// s[i - len - 1] = x となる回文 u を見つける。u の x の遷移が無ければ頂点を作り、その suffix link は u の suffix link から同じように
// たどった回文の x の遷移 (長さ 1 なら長さ 0 の根)。たどる回数は全体で 2n 以下。遷移は、文字の種類 k が 32 以下なら頂点ごとの幅 k の
// 表で持つ (10^6 文字で 256 種類なら 1 GB になるので、それより多ければ、1 つ目の子を頂点に直接持ち、2 つ目以降の子を開番地法のハッシュ
// 表に置く)。頂点の列と表は 0 で埋まった 1 つの領域 (2 MB 以上なら、Linux では 2 MB 境界の mmap に置いて huge page を頼む) に置き、
// 答えの列は領域を取ってから huge page を頼む。procon-judge の yosupo-eertree で書き比べた eer_lean4 を写し、文字の種類が多いときの
// 持ち方を足したもの。
namespace eertree_internal {
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
struct Out {
 std::vector<int> len, par, link, last;
};
// 幅 k の表で持つ形。
struct DenseNode {
 int len, link;
};
template <class T, class Code> void dense(const T* t, int n, int k, Code code, Out& out) {
 const int cap= n + 2;
 Zeroed buf(std::size_t(cap) * sizeof(DenseNode) + std::size_t(cap) * k * sizeof(int));  // 頂点の列と遷移の表 (cap 行、幅 k)
 DenseNode* nd= static_cast<DenseNode*>(buf.get());
 int* nx= reinterpret_cast<int*>(nd + cap);
 nd[0]= {-1, 0}, nd[1]= {0, 0};
 int cnt= 2, last= 1;
 for(int i= 0; i < n; ++i) {
  const T x= t[i];
  const int c= code(i);
  int u= last;
  while(true) {
   const int p= i - nd[u].len - 1;
   if(p >= 0 && t[p] == x) break;
   u= nd[u].link;
  }
  int& e= nx[std::size_t(u) * k + c];
  if(!e) {
   const int w= cnt++, L= nd[u].len + 2;
   int l= 1;
   if(L > 1) {
    int v= nd[u].link;
    while(true) {
     const int p= i - nd[v].len - 1;
     if(p >= 0 && t[p] == x) break;
     v= nd[v].link;
    }
    l= nx[std::size_t(v) * k + c];
   }
   nd[w]= {L, l};
   out.len.push_back(L), out.par.push_back(u), out.link.push_back(l);
   e= w;
  }
  last= e;
  out.last[i]= last;
 }
}
// 1 つ目の子を頂点に直接持ち、2 つ目以降の子を (頂点, 文字) を鍵にした開番地法のハッシュ表に置く形。oc の下位 31 bit は 1 つ目の子の
// 文字の番号 (NONE なら子が無い)、最上位 bit はハッシュ表に子があるか。
struct SparseNode {
 int len, link, ow;
 unsigned oc;
};
template <class T, class Code> void sparse(const T* t, int n, Code code, Out& out) {
 constexpr unsigned NONE= 0x7fffffffu, MULTI= 0x80000000u;
 const int cap= n + 2;
 int bits= 1;
 while((std::size_t(1) << bits) < 2 * std::size_t(cap)) ++bits;
 const std::size_t hcap= std::size_t(1) << bits;
 Zeroed buf(std::size_t(cap) * sizeof(SparseNode) + hcap * (sizeof(std::uint64_t) + sizeof(int)));
 SparseNode* nd= static_cast<SparseNode*>(buf.get());
 std::uint64_t* key= reinterpret_cast<std::uint64_t*>(nd + cap);  // (頂点 + 1) << 32 | 文字。0 は空き
 int* val= reinterpret_cast<int*>(key + hcap);
 auto slot= [&](int u, unsigned c) {
  const std::uint64_t kk= std::uint64_t(u + 1) << 32 | c;
  std::size_t h= std::size_t((kk * 0x9E3779B97F4A7C15ull) >> (64 - bits));
  while(key[h] && key[h] != kk) h= (h + 1) & (hcap - 1);
  return h;
 };
 auto child= [&](int u, unsigned c) {
  const SparseNode& e= nd[u];
  if((e.oc & NONE) == c) return e.ow;
  if(!(e.oc & MULTI)) return 0;
  const std::size_t h= slot(u, c);
  return key[h] ? val[h] : 0;
 };
 nd[0]= {-1, 0, 0, NONE}, nd[1]= {0, 0, 0, NONE};
 int cnt= 2, last= 1;
 for(int i= 0; i < n; ++i) {
  const T x= t[i];
  const unsigned c= unsigned(code(i));
  int u= last;
  while(true) {
   const int p= i - nd[u].len - 1;
   if(p >= 0 && t[p] == x) break;
   u= nd[u].link;
  }
  int e= child(u, c);
  if(!e) {
   const int w= cnt++, L= nd[u].len + 2;
   int l= 1;
   if(L > 1) {
    int v= nd[u].link;
    while(true) {
     const int p= i - nd[v].len - 1;
     if(p >= 0 && t[p] == x) break;
     v= nd[v].link;
    }
    l= child(v, c);
   }
   nd[w]= {L, l, 0, NONE};
   out.len.push_back(L), out.par.push_back(u), out.link.push_back(l);
   SparseNode& f= nd[u];
   if((f.oc & NONE) == NONE) {
    f.ow= w, f.oc= (f.oc & MULTI) | c;
   } else {
    const std::size_t h= slot(u, c);
    key[h]= std::uint64_t(u + 1) << 32 | c, val[h]= w, f.oc|= MULTI;
   }
   e= w;
  }
  last= e;
  out.last[i]= last;
 }
}
template <class T, class Code> void run(const T* t, int n, int k, Code code, Out& out) {
 const std::size_t cap= std::size_t(n) + 2;
 reserve_huge(out.len, cap), reserve_huge(out.par, cap), reserve_huge(out.link, cap), reserve_huge(out.last, n);
 out.len= {-1, 0}, out.par= {0, 0}, out.link= {0, 0};
 out.last.resize(n);
 if(k <= 32) dense(t, n, std::max(k, 1), code, out);
 else sparse(t, n, code, out);
}
}
// 回文木。頂点 0 は長さ -1 の根、頂点 1 は長さ 0 の根で、頂点 2 以降が異なる回文 (初めて現れて終わる位置の順)。len[v] は長さ、par[v] は
// 両端の文字を除いた回文の頂点 (長さ 1 の回文なら 0)、link[v] は suffix link (自分より短い最長の回文接尾辞の頂点、長さ 1 なら 1)。根の
// par と link は 0。last[i] は s[0, i] の最長の回文接尾辞の頂点。
struct Eertree {
 std::vector<int> len, par, link, last;
};
// 文字列 s の回文木。文字は unsigned char として比べる。O(n)。
inline Eertree eertree(const std::string& s) {
 const int n= s.size();
 const unsigned char* t= reinterpret_cast<const unsigned char*>(s.data());
 std::array<unsigned char, 256> code{};
 std::array<bool, 256> used{};
 for(int i= 0; i < n; ++i) used[t[i]]= true;
 int k= 0;
 for(int c= 0; c < 256; ++c)
  if(used[c]) code[c]= (unsigned char)k++;
 eertree_internal::Out out;
 eertree_internal::run(t, n, k, [&](int i) { return int(code[t[i]]); }, out);
 return {std::move(out.len), std::move(out.par), std::move(out.link), std::move(out.last)};
}
// 列 s の回文木。値は == と < で比べる。O(n log n) (番号の振り直しの整列)。
template <class T> Eertree eertree(const std::vector<T>& s) {
 const int n= s.size();
 std::vector<int> ord(n);
 std::iota(ord.begin(), ord.end(), 0);
 std::sort(ord.begin(), ord.end(), [&](int x, int y) { return s[x] < s[y]; });
 std::vector<int> code(n);
 int k= 0;
 for(int i= 0; i < n; ++i) {
  if(i > 0 && s[ord[i - 1]] < s[ord[i]]) ++k;
  code[ord[i]]= k;
 }
 if(n > 0) ++k;
 eertree_internal::Out out;
 eertree_internal::run(s.data(), n, k, [&](int i) { return int(code[i]); }, out);
 return {std::move(out.len), std::move(out.par), std::move(out.link), std::move(out.last)};
}
