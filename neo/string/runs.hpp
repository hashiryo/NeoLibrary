#pragma once
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstring>
#include <iterator>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
#include "neo/data_structure/FastSet.hpp"
#include "neo/data_structure/RangeMin.hpp"
#include "neo/string/suffix_array.hpp"
// 列のすべての run を Runs 定理 (Bannai、I、Inenaga、Nakashima、Takeda、Tsuruta 2017) で求める。どの run も、2 つの順序 (文字の順と
// その逆) のどちらかで、最長の Lyndon 接頭辞として現れる長さ t の根を持つ。位置 i の最長の Lyndon 接頭辞の終わりは、i より後ろで
// 順位が i より小さい (逆の順では大きい) 最初の位置なので、接尾辞配列を順位の小さい順に読み、順位が今より小さい位置の集合と大きい
// 位置の集合 (FastSet) から、i より後ろで最初の位置を引く。逆の順では列の終わりを最も大きいものとして扱い、終わりまで届く根は見ない
// (原始的でない根になりうるうえ、その run は文字の順の側で見つかる)。根 [i, j) は前後の 8 byte の排他的論理和で一致の長さを数え、
// 伸ばした長さが周期の 2 倍以上なら run とする。同じ周期の前に見つけた run の中に収まる根は飛ばすので、周期ごとには始まりの順に
// 見つかり、周期で数え分けるだけで並ぶ。8 byte とも一致したときのループで比べた長さが全体で 32n を超えたら、LCE を LCP 配列と
// RangeMin で引く形に切り替えるので、最悪でも O(n)。接尾辞配列、2 つの Lyndon 配列、写した列、run の列は 1 つの作業領域に置き、
// 1 MB 以上なら、Linux では 2 MB 境界に寄せて huge page を頼む。procon-judge の yosupo-runenumerate で書き比べた lyndon_fs_hp2 を
// 写したもの。
namespace runs_internal {
using u64= unsigned long long;
template <class C> constexpr int W= 8 / sizeof(C);                    // 8 byte に入る要素の数
template <class C> constexpr int SH= 3 + std::countr_zero(sizeof(C));  // 一致した bit の数から要素の数への shift
template <class C> inline u64 load(const C* p) {
 u64 x;
 std::memcpy(&x, p, 8);
 return x;
}
// t[i + k] = t[j + k] が続く長さ。e 個までは一致していると分かっているときに、lim 個を上限として求める。
template <class C> inline int forward(const C* t, int i, int j, int e, int lim) {
 for(; e < lim; e+= W<C>)
  if(const u64 x= load(t + i + e) ^ load(t + j + e)) return std::min(e + (std::countr_zero(x) >> SH<C>), lim);
 return lim;
}
// t[i - 1 - k] = t[j - 1 - k] が続く長さ。b 個までは一致していると分かっているときに、lim 個を上限とし、stop 個以上と分かったところで
// 止める (止めたときの値は下限)。
template <class C> inline int backward(const C* t, int i, int j, int b, int lim, int stop) {
 for(; b < lim && b < stop; b+= W<C>)
  if(const u64 y= load(t + i - W<C> - b) ^ load(t + j - W<C> - b)) return std::min(b + (std::countl_zero(y) >> SH<C>), lim);
 return std::min(b, lim);
}
// bytes byte の作業領域。0 では埋めない。hp なら 2 MB 境界に寄せて huge page を頼む。
class Buffer {
 static constexpr std::size_t H= std::size_t(1) << 21;
 std::unique_ptr<char[]> raw;
 char* p= nullptr;
public:
 Buffer(std::size_t bytes, bool hp) {
  const std::size_t a= hp ? H : 64, b= (bytes + a - 1) & ~(a - 1);
  raw.reset(new char[b + a]);
  p= reinterpret_cast<char*>((reinterpret_cast<std::uintptr_t>(raw.get()) + a - 1) & ~std::uintptr_t(a - 1));
#ifdef __linux__
  if(hp) madvise(p, b, MADV_HUGEPAGE);
#endif
 }
 char* get() const { return p; }
};
// 列の任意の 2 つの接尾辞の LCP を O(1) で引く。
class Lce {
 std::vector<int> isa;
 RangeMin rm;
public:
 template <class C> Lce(const C* s, int n, const int* sa): isa(n) {
  const std::vector<int> lcp= suffix_array_internal::lcp_phi(s, n, sa);
  for(int i= 0; i < n; ++i) isa[sa[i]]= i;
  rm= RangeMin(std::vector<unsigned>(lcp.begin(), lcp.end()));
 }
 int operator()(int i, int j) const {  // i != j、どちらも n 未満
  const auto [a, b]= std::minmax(isa[i], isa[j]);
  return int(rm.min(a, b));
 }
};
// 比べた長さが予算を超えたあとの道。初めて呼ばれたときに LCE の表を作り、根 [i, j) の run を {始まり, 終わり} で返す (run でなければ
// 始まりを -1 にする)。根を伸ばすループの中へ展開させないよう noinline にする。
template <class C> struct Slow {
 const C* s;
 int n, K;
 const int* sa;
 std::unique_ptr<Lce> fwd, bwd;
 [[gnu::noinline]] std::pair<int, int> operator()(int i, int j, int p) {
  if(!fwd) {
   const std::vector<C> r(std::make_reverse_iterator(s + n), std::make_reverse_iterator(s));
   std::vector<int> rsa(n);
   suffix_array_internal::sa_is(r.data(), n, K, rsa.data());
   fwd= std::make_unique<Lce>(s, n, sa), bwd= std::make_unique<Lce>(r.data(), n, rsa.data());
  }
  const int e= j < n ? (*fwd)(i, j) : 0, b= i > 0 ? (*bwd)(n - i, n - j) : 0;  // 反転した列で、n - i と n - j から始まる接尾辞
  return e + b >= p ? std::pair<int, int>(i - b, j + e) : std::pair<int, int>(-1, 0);
 }
};
// 値が 0 以上 K 未満の列 s[0, n) の run を (周期, 始まり, 終わり) の昇順で返す。
template <class C> std::vector<std::array<int, 3>> solve(const C* s, int n, int K) {
 if(n < 2) return {};
 auto up8= [](std::size_t x) { return (x + 7) & ~std::size_t(7); };
 const std::size_t half= std::size_t(n) / 2 + 1, words= up8(std::size_t(n) * 4);
 const std::size_t o_n0= words, o_n1= o_n0 + words, o_last= o_n1 + words, o_t= o_last + up8(half * 4),
                   o_found= o_t + up8((std::size_t(n) + 2 * W<C>) * sizeof(C)), bytes= o_found + std::size_t(n) * 12;
 Buffer buf(bytes, bytes >= (std::size_t(1) << 20));
 char* w= buf.get();
 int* sa= reinterpret_cast<int*>(w);
 int *nxt0= reinterpret_cast<int*>(w + o_n0), *nxt1= reinterpret_cast<int*>(w + o_n1);  // 2 つの順序の Lyndon 根の終わり
 int* last= reinterpret_cast<int*>(w + o_last);  // 周期ごとの、最後に見つけた run の終わり
 C* tb= reinterpret_cast<C*>(w + o_t);           // 前後に 8 byte の余白を付けた s の写し
 auto* found= reinterpret_cast<std::array<int, 3>*>(w + o_found);  // run の数は n 未満
 std::memset(tb, 0, 8), std::memcpy(tb + W<C>, s, std::size_t(n) * sizeof(C)), std::memset(tb + W<C> + n, 0, 8);
 const C* t= tb + W<C>;
 suffix_array_internal::sa_is(t, n, K, sa);
 {
  FastSet smaller(n), larger(n, [](int) { return true; });
  for(int k= 0; k < n; ++k) {
   const int i= sa[k];
   larger.erase(i);
   const int a= smaller.next(i + 1), b= larger.next(i + 1);
   nxt0[i]= a < 0 ? n : a, nxt1[i]= b < 0 ? n : b;
   smaller.insert(i);
  }
 }
 std::memset(last, 0, half * 4);
 long long budget= 32LL * n;
 Slow<C> slow{t, n, K, sa, nullptr, nullptr};
 int nf= 0;
 auto extend= [&](int i, int j, int p) __attribute__((always_inline)) {
  if(__builtin_expect(budget < 0, 0)) return slow(i, j, p);
  int e= std::countr_zero(load(t + i) ^ load(t + j)) >> SH<C>, b= std::countl_zero(load(t + i - W<C>) ^ load(t + j - W<C>)) >> SH<C>;
  if(__builtin_expect(e < W<C> && b < W<C>, 1)) {  // ほとんどの根は初めの 8 byte で決まる
   e= std::min(e, n - j), b= std::min(b, i);
   return e + b >= p ? std::pair<int, int>(i - b, j + e) : std::pair<int, int>(-1, 0);
  }
  e= e < W<C> ? std::min(e, n - j) : forward(t, i, j, W<C>, n - j);
  b= b < W<C> ? std::min(b, i) : backward(t, i, j, W<C>, i, p - e);
  budget-= e + b;
  if(e + b < p) return std::pair<int, int>(-1, 0);
  if(b >= W<C>) {
   const int b0= b;
   b= backward(t, i, j, b, i, i);
   budget-= b - b0;
  }
  return std::pair<int, int>(i - b, j + e);
 };
 auto check= [&](int i, int j) __attribute__((always_inline)) {
  const int p= j - i;
  if(2 * p > n || i + p <= last[p]) return;  // 長さが 2p に届かない根と、同じ周期の前に見つけた run の中に収まる根
  const auto [l, r]= extend(i, j, p);
  if(l < 0) return;
  found[nf++]= {p, l, r};
  last[p]= r;
 };
 for(int i= 0; i < n; ++i) {
  check(i, nxt0[i]);
  if(const int j= nxt1[i]; j < n) check(i, j);
 }
 int* cnt= last;
 std::memset(cnt, 0, half * 4);
 for(int k= 0; k < nf; ++k) ++cnt[found[k][0] - 1];
 for(int p= 0, acc= 0; p < int(half); ++p) {
  const int c= cnt[p];
  cnt[p]= acc, acc+= c;
 }
 std::vector<std::array<int, 3>> res(nf);
 for(int k= 0; k < nf; ++k) res[cnt[found[k][0] - 1]++]= found[k];
 return res;
}
}
// 文字列 s のすべての run を (最小の周期 t, 始まり l, 終わり r) の昇順で返す。文字は unsigned char として比べる。
inline std::vector<std::array<int, 3>> runs(const std::string& s) { return runs_internal::solve(reinterpret_cast<const unsigned char*>(s.data()), int(s.size()), 256); }
// 比べられる値の列 s のすべての run を返す。値を座標圧縮してから解く。
template <class T> std::vector<std::array<int, 3>> runs(const std::vector<T>& s) {
 const int n= s.size();
 std::vector<int> idx(n), t(n);
 std::iota(idx.begin(), idx.end(), 0);
 std::sort(idx.begin(), idx.end(), [&](int a, int b) { return s[a] < s[b]; });
 int K= 0;
 for(int i= 0; i < n; ++i) t[idx[i]]= K+= i && s[idx[i - 1]] < s[idx[i]];
 if(K < 256) {
  const std::vector<unsigned char> c(t.begin(), t.end());
  return runs_internal::solve(c.data(), n, K + 1);
 }
 return runs_internal::solve(t.data(), n, K + 1);
}
