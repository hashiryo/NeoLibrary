#pragma once
#include <algorithm>
#include <cstddef>
#include <new>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 0 以上 n 未満の整数の集合。64 分木の bitset で、段 h の i 番目の bit は、段 h - 1 の i 番目の語が 0 でないことを表す。
// 段は 1 本の配列に根から順に並べる。前後の要素を探すときは、葉の語を読む前に 1 つ上の段の bit を見て、葉の語が 0 なら
// 葉を読まずに上の段から探す (n = 10^7 で葉の段は 1.25 MB あり、L2 に収まらないことが多いが、1 つ上の段は 19.5 KB で
// L1 に載る)。入れるときと除くときは、上の段が変わらなくなったところで止める。procon-judge の yosupo-predecessor-problem
// で書き比べた fs_td を写したもの。1 MB 以上の配列は 2 MB 境界に置いて huge page を頼む (Linux のみ)。値の範囲が 2^30 の
// yosupo-set-xor-min (配列 128 MB) で 2 割速くなり、n = 10^7 (1.25 MB) では 0〜3% の差だった。
namespace fast_set_internal {
template <class T> struct Alloc {
 using value_type= T;
 Alloc()= default;
 template <class U> Alloc(const Alloc<U>&) {}
 T* allocate(std::size_t n) {
  const std::size_t want= n * sizeof(T);
  if(want < (std::size_t(1) << 20)) return static_cast<T*>(::operator new(want, std::align_val_t(64)));
  const std::size_t H= std::size_t(1) << 21, bytes= (want + H - 1) & ~(H - 1);
  T* p= static_cast<T*>(::operator new(bytes, std::align_val_t(H)));
#ifdef __linux__
  madvise(p, bytes, MADV_HUGEPAGE);
#endif
  return p;
 }
 void deallocate(T* p, std::size_t n) {
  if(n * sizeof(T) < (std::size_t(1) << 20)) ::operator delete(p, std::align_val_t(64));
  else ::operator delete(p, std::align_val_t(std::size_t(1) << 21));
 }
 template <class U> bool operator==(const Alloc<U>&) const { return true; }
};
}  // namespace fast_set_internal
class FastSet {
 using u64= unsigned long long;
 int n, lg, nw[8];
 int off[8];  // 段 h (0 が葉) の先頭
 std::vector<u64, fast_set_internal::Alloc<u64>> a;
 void layout(int m) {
  n= m, lg= 0;
  for(;; m= (m + 63) / 64) {
   nw[lg++]= (m + 63) / 64;
   if(m <= 64) break;
  }
  int tot= 0;
  for(int h= lg; h--;) off[h]= tot, tot+= nw[h];
  a.assign(tot, 0);
 }
 // 葉の語 w が 0 でないか。
 bool live(int w) const { return lg == 1 || (a[off[1] + w / 64] >> (w % 64) & 1); }
 void build_upper() {
  for(int h= 1; h < lg; ++h)
   for(int i= 0; i < nw[h - 1]; ++i) a[off[h] + i / 64]|= u64(a[off[h - 1] + i] != 0) << (i % 64);
 }
 // 0 でない w の立っている bit の位置 j のうち、j xor t が最小のもの (t は 6 bit)。上の bit から順に、t と同じ側に立っている
 // bit があればそちらへ絞る。
 static int pick(u64 w, unsigned t) {
  constexpr u64 K[6]= {0x5555555555555555, 0x3333333333333333, 0x0F0F0F0F0F0F0F0F, 0x00FF00FF00FF00FF, 0x0000FFFF0000FFFF, 0x00000000FFFFFFFF};
  for(int k= 5; k >= 0; --k) {
   const u64 c= w & (K[k] ^ (0 - u64(t >> k & 1)));
   w= c ? c : w;
  }
  return __builtin_ctzll(w);
 }
 // 根から、段ごとに子の番号 j と t の 6 bit の xor が最小の子へ下りる。
 int xor_pick(unsigned t) const {
  if(a.empty() || !a[off[lg - 1]]) return -1;
  int i= 0;
  for(int h= lg - 1; h >= 0; --h) i= i * 64 + pick(a[off[h] + i], t >> (6 * h) & 63);
  return i;
 }
public:
 // 空の集合。
 explicit FastSet(int n= 0) { layout(n); }
 // f(i) が真の i からなる集合。葉の語ごとに 64 個を調べてから書く。
 template <class F> FastSet(int n, F f) {
  layout(n);
  for(int w= 0; w < nw[0]; ++w) {
   u64 x= 0;
   for(int b= 0, e= std::min(64, n - w * 64); b < e; ++b) x|= u64(bool(f(w * 64 + b))) << b;
   a[off[0] + w]= x;
  }
  build_upper();
 }
 // bits[i / 64] の下から i % 64 番目の bit が 1 の i からなる集合。n 以上の位置の bit は 0 にしておく。
 FastSet(int n, const std::vector<u64>& bits) {
  layout(n);
  std::copy(bits.begin(), bits.begin() + std::min<std::size_t>(bits.size(), nw[0]), a.begin() + off[0]);
  build_upper();
 }
 bool contains(int i) const { return live(i / 64) && (a[off[0] + i / 64] >> (i % 64) & 1); }
 void insert(int i) {
  for(int h= 0; h < lg; ++h, i/= 64) {
   u64& w= a[off[h] + i / 64];
   const u64 was= w;
   w|= u64(1) << (i % 64);
   if(was) return;
  }
 }
 void erase(int i) {
  for(int h= 0; h < lg; ++h, i/= 64) {
   u64& w= a[off[h] + i / 64];
   w&= ~(u64(1) << (i % 64));
   if(w) return;
  }
 }
 // i 以上で最小の要素。無ければ -1。
 int next(int i) const {
  if(i < 0) i= 0;
  if(i >= n) return -1;
  if(live(i / 64))
   if(const u64 d= a[off[0] + i / 64] >> (i % 64)) return i + __builtin_ctzll(d);
  i= i / 64 + 1;
  for(int h= 1; h < lg; ++h) {
   if(i / 64 == nw[h]) break;
   const u64 d= a[off[h] + i / 64] >> (i % 64);
   if(!d) {
    i= i / 64 + 1;
    continue;
   }
   i+= __builtin_ctzll(d);
   for(int g= h - 1; g >= 0; --g) i= i * 64 + __builtin_ctzll(a[off[g] + i]);
   return i;
  }
  return -1;
 }
 // y xor x が最小になる要素 y。空なら -1。
 int xor_min(int x) const { return xor_pick(unsigned(x)); }
 // y xor x が最大になる要素 y。空なら -1。
 int xor_max(int x) const { return xor_pick(~unsigned(x)); }
 // i 以下で最大の要素。無ければ -1。
 int prev(int i) const {
  if(i >= n) i= n - 1;
  if(i < 0) return -1;
  if(live(i / 64))
   if(const u64 d= a[off[0] + i / 64] << (63 - i % 64)) return i - __builtin_clzll(d);
  i= i / 64 - 1;
  for(int h= 1; h < lg; ++h) {
   if(i == -1) break;
   const u64 d= a[off[h] + i / 64] << (63 - i % 64);
   if(!d) {
    i= i / 64 - 1;
    continue;
   }
   i-= __builtin_clzll(d);
   for(int g= h - 1; g >= 0; --g) i= i * 64 + 63 - __builtin_clzll(a[off[g] + i]);
   return i;
  }
  return -1;
 }
};
