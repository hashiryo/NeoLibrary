#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// s の接頭辞と t の部分文字列の最長共通部分列の長さを、Tiskin の seaweed の組み替えで求める。s の行を 1 つずつ足しながら、t の各列の
// 下端に出る糸の出どころ h[j] (上端の列 x から入った糸なら x、左端から入った糸なら -1) を求める。行 i では左端から -1 の糸を入れ、列を
// 左から順に見て、s[i] = t[j] か h[j] が今の横の糸より小さければ 2 本を入れ替える。a 行のあとで、s[0, a) と t[b, c) の答えは
// D(b, c) - b (D(b, c) = #{j < c : h[j] < b}。j < b なら h[j] <= j < b なので、その分はちょうど b 個)。h[j] < b となる j の集合は、
// h[j] = -1 の位置に、x < b の糸 x が下端に出る位置 pos[x] を足したもの。各行で b0 = 64 k ごとに集合を 64 bit の語の列と各語より前の
// 位置の数 cum で控えておくと、D(b, c) = cum[c / 64] + (語 c / 64 のうち c より下の popcount) + #{x in [b0, b) : pos[x] < c}
// (b0 は b を 64 の倍数に切り下げた値) で、最後の項は pos の 64 個を AVX2 の比較 4 回で数える。控えの語は、h の 64 個を 1 回読んで
// 控えごとに比べてまとめる。表は 1 つの 0 で埋まった領域 (2 MB 以上なら、Linux では 2 MB 境界の mmap に置いて huge page を頼む) に
// 置く。procon-judge の yosupo-prefix-substring-lcs で書き比べた psl_ckpt を写し、行の幅と控えの数を t の長さで決めるようにしたもの。
namespace prefix_substring_lcs_internal {
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
struct Entry {
 std::uint64_t word;
 std::uint32_t cum, pad;
};
}
// s[0, a) と t[b, c) の最長共通部分列の長さを O(1) で答える。|t| < 32767。構築は O(|s| |t| + |s| (|t| / 64)^2) で、表は
// (|s| + 1) ((|t| / 64 + 1)^2 * 16 + (|t| / 64 + 2) * 128) バイト (|s| = |t| = 1000 で 6 MB ほど)。
class PrefixSubstringLCS {
 using Entry= prefix_substring_lcs_internal::Entry;
 static constexpr std::int16_t BIG= 32767;
 int m, n, nk, pw;  // nk は控えの数 (各控えの語の数も同じ)、pw は pos の行の幅
 std::size_t ck_row;
 prefix_substring_lcs_internal::Zeroed buf;
 std::int16_t* pos;
 Entry* ck;
 template <class T> void build(const T* s, const T* t) {
  std::vector<std::int16_t> h(std::size_t(nk) * 64, BIG);  // n 以降は BIG で埋め、どの b0 よりも小さくならないようにする
  for(int j= 0; j < n; ++j) h[j]= std::int16_t(j);
  for(int a= 0; a <= m; ++a) {
   if(a > 0) {
    const T x= s[a - 1];
    int v= -1;
    for(int j= 0; j < n; ++j) {
     const int hj= h[j];
     const bool sw= (x == t[j]) | (hj < v);
     h[j]= std::int16_t(sw ? v : hj);
     v= sw ? hj : v;
    }
   }
   std::int16_t* pr= pos + std::size_t(a) * pw;
   for(int x= 0; x < pw; ++x) pr[x]= BIG;  // x >= n と右端に出た糸は BIG、pw - 1 番は h[j] = -1 の捨て場
   for(int j= 0; j < n; ++j) {
    const int hj= h[j];
    pr[hj >= 0 ? hj : pw - 1]= std::int16_t(j);
   }
   // 控え k の語 e の bit j は h[64 e + j] < 64 k。
   Entry* ca= ck + std::size_t(a) * ck_row;
   for(int e= 0; e < nk; ++e) {
    const __m256i* hp= reinterpret_cast<const __m256i*>(h.data() + 64 * e);
    const __m256i h0= _mm256_loadu_si256(hp), h1= _mm256_loadu_si256(hp + 1), h2= _mm256_loadu_si256(hp + 2), h3= _mm256_loadu_si256(hp + 3);
    for(int k= 0; k < nk; ++k) {
     const __m256i vb= _mm256_set1_epi16(std::int16_t(64 * k));
     const __m256i lo= _mm256_permute4x64_epi64(_mm256_packs_epi16(_mm256_cmpgt_epi16(vb, h0), _mm256_cmpgt_epi16(vb, h1)), 0xd8);
     const __m256i hi= _mm256_permute4x64_epi64(_mm256_packs_epi16(_mm256_cmpgt_epi16(vb, h2), _mm256_cmpgt_epi16(vb, h3)), 0xd8);
     ca[std::size_t(k) * nk + e].word= std::uint64_t(std::uint32_t(_mm256_movemask_epi8(lo))) | std::uint64_t(std::uint32_t(_mm256_movemask_epi8(hi))) << 32;
    }
   }
   for(int k= 0; k < nk; ++k)
    for(int e= 0, acc= 0; e < nk; ++e) ca[std::size_t(k) * nk + e].cum= std::uint32_t(acc), acc+= std::popcount(ca[std::size_t(k) * nk + e].word);
  }
 }
 static std::size_t bytes(int m, int nk, int pw) { return std::size_t(m + 1) * (std::size_t(pw) * sizeof(std::int16_t) + std::size_t(nk) * nk * sizeof(Entry)); }
public:
 template <class T> PrefixSubstringLCS(const std::vector<T>& s, const std::vector<T>& t): m(s.size()), n(t.size()), nk((n >> 6) + 1), pw((nk + 1) * 64), ck_row(std::size_t(nk) * nk), buf(bytes(m, nk, pw)) {
  pos= static_cast<std::int16_t*>(buf.get());
  ck= reinterpret_cast<Entry*>(pos + std::size_t(m + 1) * pw);
  build(s.data(), t.data());
 }
 PrefixSubstringLCS(const std::string& s, const std::string& t): m(s.size()), n(t.size()), nk((n >> 6) + 1), pw((nk + 1) * 64), ck_row(std::size_t(nk) * nk), buf(bytes(m, nk, pw)) {
  pos= static_cast<std::int16_t*>(buf.get());
  ck= reinterpret_cast<Entry*>(pos + std::size_t(m + 1) * pw);
  build(s.data(), t.data());
 }
 PrefixSubstringLCS(const PrefixSubstringLCS&)= delete;
 PrefixSubstringLCS& operator=(const PrefixSubstringLCS&)= delete;
 // s[0, a) と t[b, c) の最長共通部分列の長さ。0 <= a <= |s|、0 <= b <= c <= |t|。
 int query(int a, int b, int c) const {
  const int k0= b >> 6, b0= k0 << 6;
  const Entry e= ck[std::size_t(a) * ck_row + std::size_t(k0) * nk + (c >> 6)];
  const int d= int(e.cum) + std::popcount(e.word & ((std::uint64_t(1) << (c & 63)) - 1));
  const std::int16_t* pp= pos + std::size_t(a) * pw + b0;
  const __m256i vc= _mm256_set1_epi16(std::int16_t(c)), vl= _mm256_set1_epi16(std::int16_t(b - b0)), l16= _mm256_set1_epi16(16);
  __m256i ln= _mm256_setr_epi16(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
  int cnt= 0;
  for(int i= 0; i < 4; ++i) {
   const __m256i p= _mm256_loadu_si256(reinterpret_cast<const __m256i*>(pp + 16 * i));
   cnt+= std::popcount(unsigned(_mm256_movemask_epi8(_mm256_and_si256(_mm256_cmpgt_epi16(vc, p), _mm256_cmpgt_epi16(vl, ln)))));
   ln= _mm256_add_epi16(ln, l16);
  }
  return d + (cnt >> 1) - b;
 }
};
