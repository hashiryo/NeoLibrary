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
// 1 段で 3 bit ずつ分ける 8 分の形で、段の数は 2 分の 1/3 になる。段ごとに 64 個ずつのブロックを、記号 (0 から 7) の
// 3 本の bit の列と、ブロックより前の記号 t 未満の数 (t = 0..8) の 64 byte に詰める。k 番目は、両端で記号 1 から 7 未満の数を
// 7 個ずつ求め、区間の中の数と k を比べた 7 つの比較の和で記号を選ぶ (分岐なし)。[lo, hi) の個数は、2 つの境の記号が同じ段を
// 1 回で下り、違う段から下は 2 本を同じループで下りる。値から順位への変換は、値の範囲を最大 65536 個の区間に分けた表で狭めてから
// cmov の二分探索で探す。2 MB 以上の配列は、Linux では 2 MB 境界の mmap に置いて MADV_HUGEPAGE を頼む。
// procon-judge の yosupo-range-kth-smallest などで書き比べた wm8kp (k 番目) と wm8_rng ([lo, hi) の個数) を写したもの。
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
 struct alignas(64) Blk {
  u64 b0, b1, b2;  // 記号の bit 0、1、2
  u32 c[9];        // c[t] はブロックより前の記号 t 未満の数
 };
 // 段の記号と、その 3 本の bit を全部の bit に広げたマスク。
 struct Sym {
  int s;
  u64 a0, a1, a2;
  [[gnu::always_inline]] explicit Sym(int s_): s(s_), a0(-u64(s_ & 1)), a1(-u64(s_ >> 1 & 1)), a2(-u64(s_ >> 2)) {}
  // 記号 s 未満の印と、記号 s の印。
  [[gnu::always_inline]] u64 lt(const Blk& k) const { return (~k.b2 & a2) | (~(k.b2 ^ a2) & ((~k.b1 & a1) | (~(k.b1 ^ a1) & ~k.b0 & a0))); }
  [[gnu::always_inline]] u64 eq(const Blk& k) const { return ~(k.b2 ^ a2) & ~(k.b1 ^ a1) & ~(k.b0 ^ a0); }
 };
 int n= 0, levels= 1, nb= 0;
 u32 sigma= 0;
 wavelet_matrix_internal::Ranks<T, Compress> val;
 wavelet_matrix_internal::Vec<Blk> b;  // 段 d (下から d 番目の 3 bit) のブロック j は b[d * nb + j]
 std::vector<std::array<int, 8>> z;    // z[d][s] は段 d の記号 s より小さい記号の数
 [[gnu::always_inline]] static u64 low(int i) { return (u64(1) << (i & 63)) - 1; }
 // [0, i) の記号 s 未満の数と、記号 s の数。k はブロック i / 64。
 [[gnu::always_inline]] static int rank_lt(const Blk& k, int i, const Sym& y) { return int(k.c[y.s]) + __builtin_popcountll(y.lt(k) & low(i)); }
 [[gnu::always_inline]] static int rank_eq(const Blk& k, int i, const Sym& y) { return int(k.c[y.s + 1] - k.c[y.s]) + __builtin_popcountll(y.eq(k) & low(i)); }
 [[gnu::always_inline]] const Blk* level(int d) const { return b.data() + std::size_t(d) * nb; }
 // ブロックより前の記号 t 未満の数 c と、ブロックの中の印 x の popcount の和。空の asm に通して、clang が popcount をベクタの表引きに、
 // 足し算を 4 個ずつのベクタにまとめないようにする (まとめると x64-clang で k 番目が 2 割遅かった)。
 [[gnu::always_inline]] static int pc(u32 c, u64 x) {
  int v= int(c) + __builtin_popcountll(x);
  asm("" : "+r"(v));
  return v;
 }
 // ブロック k の [0, i) のうち記号 t 未満の数を L[t] に書く (t = 0..8)。
 [[gnu::always_inline]] static void rank_all(const Blk& k, int i, int* L) {
  const u64 m= low(i), n2= ~k.b2, n1= ~k.b1, a10= k.b1 & k.b0, o10= ~(k.b1 | k.b0);
  L[0]= 0;
  L[1]= pc(k.c[1], n2 & o10 & m);
  L[2]= pc(k.c[2], n2 & n1 & m);
  L[3]= pc(k.c[3], n2 & ~a10 & m);
  L[4]= pc(k.c[4], n2 & m);
  L[5]= pc(k.c[5], (n2 | o10) & m);
  L[6]= pc(k.c[6], (n2 | n1) & m);
  L[7]= pc(k.c[7], ~(k.b2 & a10) & m);
  L[8]= i;
 }
 // 区間の中の記号 t 未満の数 D[t] は、配列に書かずに比べる (gcc は配列に書いた 8 個の 32 bit を 1 回のベクタで読み直し、書き込みから
 // の受け渡しが効かずに遅かった)。D[t] も空の asm に通して、clang が比較をベクタにまとめないようにする。選んだ記号の数だけを、両端の
 // 配列から引く。
 [[gnu::always_inline]] u32 kth_rank(int l, int r, int k) const {
  u32 res= 0;
  for(int d= levels - 1; d >= 0; --d) {
   const Blk* B= level(d);
   int Ll[9], Lr[9];
   rank_all(B[l >> 6], l, Ll), rank_all(B[r >> 6], r, Lr);
   int d1= Lr[1] - Ll[1], d2= Lr[2] - Ll[2], d3= Lr[3] - Ll[3], d4= Lr[4] - Ll[4], d5= Lr[5] - Ll[5], d6= Lr[6] - Ll[6], d7= Lr[7] - Ll[7];
   asm("" : "+r"(d1), "+r"(d2), "+r"(d3), "+r"(d4), "+r"(d5), "+r"(d6), "+r"(d7));
   const int sy= int(d1 <= k) + int(d2 <= k) + int(d3 <= k) + int(d4 <= k) + int(d5 <= k) + int(d6 <= k) + int(d7 <= k);
   const std::array<int, 8>& Z= z[d];
   k-= Lr[sy] - Ll[sy];
   l= Z[sy] + Ll[sy + 1] - Ll[sy];
   r= Z[sy] + Lr[sy + 1] - Lr[sy];
   res= res << 3 | u32(sy);
  }
  return res;
 }
 // [l, r) の順位 x 未満の数。x は 0 以上 sigma 以下。
 [[gnu::always_inline]] int count_rank(int l, int r, u32 x) const {
  if(x >= sigma) return r - l;
  int c= 0;
  for(int d= levels - 1; d >= 0; --d) {
   const Sym y(x >> (3 * d) & 7);
   const Blk *B= level(d), &bl= B[l >> 6], &br= B[r >> 6];
   const int zs= z[d][y.s];
   c+= rank_lt(br, r, y) - rank_lt(bl, l, y);
   l= zs + rank_eq(bl, l, y), r= zs + rank_eq(br, r, y);
  }
  return c;
 }
 // [l, r) の順位が [a, b) に入る数 (a < b)。2 つの境の記号が同じ段では、両方の記号未満の数が打ち消し合うので数えずに下りる。
 // 記号が違う段で記号 a から b 未満の数を数え、そこから下は b の側 (記号 b 未満を足す) と a の側 (記号 a 未満を引く) を同じ
 // ループで下りる。
 [[gnu::always_inline]] int count_range_rank(int l, int r, u32 a, u32 b) const {
  if(b >= sigma) return r - l - count_rank(l, r, a);
  for(int d= levels - 1; d >= 0; --d) {
   const Sym ya(a >> (3 * d) & 7), yb(b >> (3 * d) & 7);
   const Blk *B= level(d), &bl= B[l >> 6], &br= B[r >> 6];
   if(ya.s == yb.s) {
    const int zs= z[d][yb.s];
    l= zs + rank_eq(bl, l, yb), r= zs + rank_eq(br, r, yb);
    continue;
   }
   int c= rank_lt(br, r, yb) - rank_lt(bl, l, yb) - (rank_lt(br, r, ya) - rank_lt(bl, l, ya));
   int lb= z[d][yb.s] + rank_eq(bl, l, yb), rb= z[d][yb.s] + rank_eq(br, r, yb);
   int la= z[d][ya.s] + rank_eq(bl, l, ya), ra= z[d][ya.s] + rank_eq(br, r, ya);
   for(int e= d - 1; e >= 0; --e) {
    const Sym ua(a >> (3 * e) & 7), ub(b >> (3 * e) & 7);
    const Blk *E= level(e), &kl= E[lb >> 6], &kr= E[rb >> 6], &hl= E[la >> 6], &hr= E[ra >> 6];
    c+= rank_lt(kr, rb, ub) - rank_lt(kl, lb, ub) - (rank_lt(hr, ra, ua) - rank_lt(hl, la, ua));
    lb= z[e][ub.s] + rank_eq(kl, lb, ub), rb= z[e][ub.s] + rank_eq(kr, rb, ub);
    la= z[e][ua.s] + rank_eq(hl, la, ua), ra= z[e][ua.s] + rank_eq(hr, ra, ua);
   }
   return c;
  }
  return 0;
 }
 // cur は順位の列。
 template <class R> void build(std::vector<R> cur) {
  int lg= 1;
  while((u32(1) << lg) < sigma) ++lg;
  levels= (lg + 2) / 3, nb= n / 64 + 1;
  b.assign(std::size_t(nb) * levels, Blk{});
  z.assign(levels, {});
  std::vector<R> nxt(n);
  for(int d= levels - 1; d >= 0; --d) {
   const int sh= 3 * d;
   Blk* B= b.data() + std::size_t(d) * nb;
   int cnt[8]= {};
   for(int j= 0; j < nb; ++j) {
    for(int t= 0, s= 0; t <= 8; ++t) B[j].c[t]= s, s+= t < 8 ? cnt[t] : 0;
    const int m= std::max(0, std::min(64, n - j * 64));
    u64 x0= 0, x1= 0, x2= 0;
    for(int t= 0; t < m; ++t) {
     const u32 sy= u32(cur[j * 64 + t]) >> sh & 7;
     x0|= u64(sy & 1) << t, x1|= u64(sy >> 1 & 1) << t, x2|= u64(sy >> 2) << t, ++cnt[sy];
    }
    B[j].b0= x0, B[j].b1= x1, B[j].b2= x2;
   }
   int pos[8];
   for(int s= 0, acc= 0; s < 8; ++s) z[d][s]= pos[s]= acc, acc+= cnt[s];
   for(int i= 0; i < n; ++i) nxt[pos[u32(cur[i]) >> sh & 7]++]= cur[i];
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
 [[gnu::always_inline]] int count(int l, int r, T lo, T hi) const {
  if(!(lo < hi)) return 0;
  const u32 a= val.template rank<false>(lo), b= val.template rank<false>(hi);
  return a < b ? count_range_rank(l, r, a, b) : 0;
 }
 // [l, r) のうちちょうど x の個数。
 [[gnu::always_inline]] int freq(int l, int r, T x) const {
  const u32 id= val.template rank<false>(x);
  if(id == sigma || !(val.value(id) == x)) return 0;
  for(int d= levels - 1; d >= 0; --d) {
   const Sym y(id >> (3 * d) & 7);
   const Blk* B= level(d);
   const int zs= z[d][y.s];
   l= zs + rank_eq(B[l >> 6], l, y), r= zs + rank_eq(B[r >> 6], r, y);
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
