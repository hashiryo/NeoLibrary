#pragma once
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <new>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
#include "neo/sort/compress.hpp"
// 静的な整数の列の wavelet matrix。区間の k 番目、値が x 未満の個数、[lo, hi) の個数、ちょうど x の個数、
// x 未満の最大と x 以上の最小を引く。値は内部で順位に詰める (neo/sort/compress.hpp の compress)。Compress = false なら、
// 値 (0 以上 2^31 未満) をそのまま順位として使い、詰める手間と、クエリのたびに値から順位を引く手間を省く。
// 1 段で 2 bit ずつ分ける 4 分の形で、段の数は 2 分の半分になる。段ごとに 64 個ずつのブロックを、記号 (0 から 3) の
// 上の bit の列と下の bit の列と、ブロックより前の記号ごとの数の 32 byte に詰める。k 番目は、記号ごとの数と k を
// 比べたマスクで記号を選ぶ (分岐なし)。値から順位への変換は、値の範囲を最大 65536 個の区間に分けた表で狭めてから
// cmov の二分探索で探す。2 MB 以上の配列は、Linux では 2 MB 境界の mmap に置いて MADV_HUGEPAGE を頼む。
// procon-judge の yosupo-range-kth-smallest などで書き比べた wm4 を写したもの。
namespace wavelet_matrix_internal {
using u32= unsigned;
using u64= unsigned long long;
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
template <class T> using Vec= std::vector<T, HugeAlloc<T>>;
// c が真なら全部の bit が 1 のマスク。GCC は同じ条件の選択をまとめて分岐に戻すので、空の asm で出どころを隠す。
[[gnu::always_inline]] inline int mask(bool c) {
 int m= -int(c);
 asm("" : "+r"(m));
 return m;
}
[[gnu::always_inline]] inline int sel(int m, int a, int b) { return b ^ ((a ^ b) & m); }
// 昇順の v[0, n) のうち x 未満 (LE なら x 以下) のものの数。clang にも cmov を選ばせる。
template <bool LE, class T> [[gnu::always_inline]] inline int lower_cm(const T* v, int n, T x) {
 if(n == 0) return 0;
 const T* b= v;
 while(n > 1) {
  const int half= n / 2;
  const bool go= LE ? !(x < b[half]) : b[half] < x;
#ifdef __clang__
  b= __builtin_unpredictable(go) ? b + half : b;
#else
  b= go ? b + half : b;
#endif
  n-= half;
 }
 return int(b - v) + (LE ? !(x < *b) : *b < x);
}
// 昇順の列 v (重複してよい) のうち x 未満 (LE なら x 以下) のものの数を引く表。値の範囲 [v の最小, v の最大] を最大
// 65536 個の区間に分け、区間を表で決めてから、その中を cmov の二分探索で探す。範囲の外の x は表を引かない。
template <class T> struct Lookup {
 using U= std::make_unsigned_t<T>;
 std::vector<T> v;
 U lo= 0;
 int sh= 0;
 std::vector<int> t;  // t[k] は v のうち ((v - lo) >> sh) < k のものの数
 // v を入れてから呼ぶ。
 void build() {
  if(v.empty()) return;
  lo= U(v[0]);
  const U d= U(U(v.back()) - lo);
  while((d >> sh) >= (U(1) << 16)) ++sh;
  const int K= int(d >> sh) + 1, m= v.size();
  t.resize(K + 1);
  for(int i= 0, k= 0; k <= K; ++k) {
   while(i < m && int(U(U(v[i]) - lo) >> sh) < k) ++i;
   t[k]= i;
  }
 }
 template <bool LE> [[gnu::always_inline]] int rank(T x) const {
  const int m= v.size();
  if(m == 0 || (LE ? x < v[0] : !(v[0] < x))) return 0;
  if(LE ? !(x < v[m - 1]) : v[m - 1] < x) return m;
  const int k= int(U(U(x) - lo) >> sh), a= t[k];
  return a + lower_cm<LE>(v.data() + a, t[k + 1] - a, x);
 }
};
// 値と順位の対応。Compress なら値を compress で順位に詰め、異なる値の昇順と値から順位を引く表を持つ。Compress でなければ、
// 値 (0 以上 2^31 未満) をそのまま順位とし、sigma は値の最大 + 1。
template <class T, bool Compress> struct Ranks {
 Lookup<T> val;
 u32 sigma= 0;
 // a の値を順位に書き換える。
 void build(std::vector<T>& a) {
  if constexpr(Compress) {
   val.v= compress(a), val.build();
   sigma= val.v.size();
  } else
   for(const T& x: a) sigma= std::max(sigma, u32(x) + 1);
 }
 // 値が x 未満 (LE なら x 以下) の順位の数。
 template <bool LE> [[gnu::always_inline]] u32 rank(T x) const {
  if constexpr(Compress) return val.template rank<LE>(x);
  else {
   if constexpr(std::is_signed_v<T>)
    if(x < T(0)) return 0;
   return u64(x) >= sigma ? sigma : u32(x) + LE;
  }
 }
 [[gnu::always_inline]] T value(u32 r) const {
  if constexpr(Compress) return val.v[r];
  else return T(r);
 }
};
// 順位に書き換えた a を、段ごとの処理で読む型の列にする。4 byte の整数ならそのまま、ほかは u32 に移す。
template <class T> inline auto rank_list(std::vector<T>&& a) {
 if constexpr(sizeof(T) == 4) return std::move(a);
 else {
  std::vector<u32> c(a.size());
  for(std::size_t i= 0; i < a.size(); ++i) c[i]= u32(std::make_unsigned_t<T>(a[i]));
  return c;
 }
}
}
template <class T= int, bool Compress= true> class WaveletMatrix {
 using u32= unsigned;
 using u64= unsigned long long;
 struct alignas(32) Blk {
  u64 h, l;
  u32 c[4];
 };
 int n= 0, levels= 1, nb= 0;
 u32 sigma= 0;
 wavelet_matrix_internal::Ranks<T, Compress> val;
 wavelet_matrix_internal::Vec<Blk> b;  // 段 d (下から d 番目の 2 bit) のブロック j は b[d * nb + j]
 std::vector<std::array<int, 4>> z;    // z[d][s] は段 d の記号 s より小さい記号の数
 [[gnu::always_inline]] static int rank_eq(const Blk& k, int i, int s) {
  const u64 m= (u64(1) << (i & 63)) - 1, hs= -u64(s >> 1), ls= -u64(s & 1);
  return k.c[s] + __builtin_popcountll(~(k.h ^ hs) & ~(k.l ^ ls) & m);
 }
 // 記号 s 未満の [0, i) の数 (s = 0 なら 0)。
 [[gnu::always_inline]] static int rank_lt(const Blk& k, int i, int s) {
  const u64 m= (u64(1) << (i & 63)) - 1, m2= -u64(s >= 2), m3= -u64(s == 3), m1= -u64(s == 1);
  const u64 a= m2 | (m1 & ~k.l), b3= m3 & ~k.l;
  const u32 k1= -u32(s >= 1), k2= -u32(s >= 2), k3= -u32(s >= 3);
  return int((k.c[0] & k1) + (k.c[1] & k2) + (k.c[2] & k3)) + __builtin_popcountll(((~k.h & a) | (k.h & b3)) & m);
 }
 [[gnu::always_inline]] u32 kth_rank(int l, int r, int k) const {
  using wavelet_matrix_internal::mask, wavelet_matrix_internal::sel;
  u32 res= 0;
  for(int d= levels - 1; d >= 0; --d) {
   const Blk *B= b.data() + std::size_t(d) * nb, &bl= B[l >> 6], &br= B[r >> 6];
   const u64 ml= (u64(1) << (l & 63)) - 1, mr= (u64(1) << (r & 63)) - 1;
   const int l0= bl.c[0] + __builtin_popcountll(~bl.h & ~bl.l & ml), l1= bl.c[1] + __builtin_popcountll(~bl.h & bl.l & ml);
   const int l2= bl.c[2] + __builtin_popcountll(bl.h & ~bl.l & ml), l3= l - l0 - l1 - l2;
   const int r0= br.c[0] + __builtin_popcountll(~br.h & ~br.l & mr), r1= br.c[1] + __builtin_popcountll(~br.h & br.l & mr);
   const int r2= br.c[2] + __builtin_popcountll(br.h & ~br.l & mr), r3= r - r0 - r1 - r2;
   const int a0= r0 - l0, a01= a0 + (r1 - l1), a012= a01 + (r2 - l2);
   const int g1= mask(k >= a0), g2= mask(k >= a01), g3= mask(k >= a012);
   const std::array<int, 4>& Z= z[d];
   k-= sel(g3, a012, sel(g2, a01, a0 & g1));
   l= sel(g3, Z[3] + l3, sel(g2, Z[2] + l2, sel(g1, Z[1] + l1, l0)));
   r= sel(g3, Z[3] + r3, sel(g2, Z[2] + r2, sel(g1, Z[1] + r1, r0)));
   res= res << 2 | u32((g1 & 1) + (g2 & 1) + (g3 & 1));
  }
  return res;
 }
 // [l, r) の順位 x 未満の数。x は 0 以上 sigma 以下。
 [[gnu::always_inline]] int count_rank(int l, int r, u32 x) const {
  if(x >= sigma) return r - l;
  int c= 0;
  for(int d= levels - 1; d >= 0; --d) {
   const int s= x >> (2 * d) & 3, zs= z[d][s];
   const Blk *B= b.data() + std::size_t(d) * nb, &bl= B[l >> 6], &br= B[r >> 6];
   c+= rank_lt(br, r, s) - rank_lt(bl, l, s);
   l= zs + rank_eq(bl, l, s), r= zs + rank_eq(br, r, s);
  }
  return c;
 }
 // cur は順位の列。
 template <class R> void build(std::vector<R> cur) {
  int lg= 1;
  while((u32(1) << lg) < sigma) ++lg;
  levels= (lg + 1) / 2, nb= n / 64 + 1;
  b.assign(std::size_t(nb) * levels, Blk{});
  z.assign(levels, {});
  std::vector<R> nxt(n);
  for(int d= levels - 1; d >= 0; --d) {
   const int sh= 2 * d;
   Blk* B= b.data() + std::size_t(d) * nb;
   int cnt[4]= {};
   for(int j= 0; j < nb; ++j) {
    for(int s= 0; s < 4; ++s) B[j].c[s]= cnt[s];
    const int m= std::max(0, std::min(64, n - j * 64));
    u64 hi= 0, lo= 0;
    for(int t= 0; t < m; ++t) {
     const u32 sy= u32(cur[j * 64 + t]) >> sh & 3;
     hi|= u64(sy >> 1) << t, lo|= u64(sy & 1) << t, ++cnt[sy];
    }
    B[j].h= hi, B[j].l= lo;
   }
   z[d]= {0, cnt[0], cnt[0] + cnt[1], cnt[0] + cnt[1] + cnt[2]};
   int pos[4]= {z[d][0], z[d][1], z[d][2], z[d][3]};
   for(int i= 0; i < n; ++i) nxt[pos[u32(cur[i]) >> sh & 3]++]= cur[i];
   cur.swap(nxt);
  }
 }
public:
 WaveletMatrix(): WaveletMatrix(std::vector<T>()) {}
 // Compress = false なら、a の値はどれも 0 以上 2^31 未満。
 explicit WaveletMatrix(std::vector<T> a): n(a.size()) {
  val.build(a);
  sigma= val.sigma;
  build(wavelet_matrix_internal::rank_list(std::move(a)));
 }
 int size() const { return n; }
 // [l, r) の小さい方から k 番目 (0 始まり)。0 <= k < r - l。
 [[gnu::always_inline]] T kth_smallest(int l, int r, int k) const { return val.value(kth_rank(l, r, k)); }
 // [l, r) の大きい方から k 番目 (0 始まり)。0 <= k < r - l。
 [[gnu::always_inline]] T kth_largest(int l, int r, int k) const { return kth_smallest(l, r, r - l - 1 - k); }
 // [l, r) のうち x 未満の個数。
 [[gnu::always_inline]] int count_lt(int l, int r, T x) const { return count_rank(l, r, val.template rank<false>(x)); }
 // [l, r) のうち x 以下の個数。
 [[gnu::always_inline]] int count_le(int l, int r, T x) const { return count_rank(l, r, val.template rank<true>(x)); }
 // [l, r) のうち値が [lo, hi) に入る個数。lo >= hi なら 0。
 [[gnu::always_inline]] int count(int l, int r, T lo, T hi) const { return lo < hi ? count_lt(l, r, hi) - count_lt(l, r, lo) : 0; }
 // [l, r) のうちちょうど x の個数。
 [[gnu::always_inline]] int freq(int l, int r, T x) const {
  const u32 id= val.template rank<false>(x);
  if(id == sigma || !(val.value(id) == x)) return 0;
  for(int d= levels - 1; d >= 0; --d) {
   const int s= id >> (2 * d) & 3, zs= z[d][s];
   const Blk* B= b.data() + std::size_t(d) * nb;
   l= zs + rank_eq(B[l >> 6], l, s), r= zs + rank_eq(B[r >> 6], r, s);
  }
  return r - l;
 }
 // [l, r) のうち x 未満の最大。無ければ nullopt。
 [[gnu::always_inline]] std::optional<T> max_lt(int l, int r, T x) const {
  const int c= count_lt(l, r, x);
  return c ? std::optional<T>(kth_smallest(l, r, c - 1)) : std::nullopt;
 }
 // [l, r) のうち x 以下の最大。無ければ nullopt。
 [[gnu::always_inline]] std::optional<T> max_le(int l, int r, T x) const {
  const int c= count_le(l, r, x);
  return c ? std::optional<T>(kth_smallest(l, r, c - 1)) : std::nullopt;
 }
 // [l, r) のうち x 以上の最小。無ければ nullopt。
 [[gnu::always_inline]] std::optional<T> min_ge(int l, int r, T x) const {
  const int c= count_lt(l, r, x);
  return c < r - l ? std::optional<T>(kth_smallest(l, r, c)) : std::nullopt;
 }
 // [l, r) のうち x より大きい最小。無ければ nullopt。
 [[gnu::always_inline]] std::optional<T> min_gt(int l, int r, T x) const {
  const int c= count_le(l, r, x);
  return c < r - l ? std::optional<T>(kth_smallest(l, r, c)) : std::nullopt;
 }
};
