#pragma once
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <new>
#include <random>
#include <string>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 法 2^61 - 1 の多項式のローリングハッシュ。基数 B はプロセスごとに乱択し、すべてのインスタンスで同じものを使う (別の文字列の部分文字列
// どうしも比べられる)。h[i] を s[0, i) のハッシュ、pw[i] を B^i として持ち、s[l, r) のハッシュは h[r] - h[l] B^(r - l)。2 つの部分文字列が
// 等しいかは、h[a + len] - h[b + len] = (h[a] - h[b]) B^len で見る (掛け算 1 回)。表は 2 MB 以上なら Linux では 2 MB 境界の mmap に置いて
// huge page を頼み、h は 512 個ごとに 1 つずらして置く (2 MB 境界から置いた h と pw の同じ添字の項や、2 の冪おきの項が、キャッシュの
// 同じ組に集まらないように)。procon-judge の self-string-substring-equal で書き比べた rh61_skew を写したもの。
namespace rolling_hash_internal {
using u64= std::uint64_t;
constexpr u64 M= (u64(1) << 61) - 1;
inline u64 mul(u64 a, u64 b) {
 const unsigned __int128 t= (unsigned __int128)a * b;
 const u64 r= (u64(t) & M) + u64(t >> 61);
 return r >= M ? r - M : r;
}
inline u64 add(u64 a, u64 b) {
 const u64 r= a + b;
 return r >= M ? r - M : r;
}
inline u64 sub(u64 a, u64 b) { return a >= b ? a - b : a + M - b; }
inline u64 base() {
 static const u64 B= [] {
  std::random_device rd;
  u64 z= (u64(rd()) << 32 | rd()) ^ u64(std::chrono::steady_clock::now().time_since_epoch().count());
  z= (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
  z= (z ^ (z >> 27)) * 0x94D049BB133111EBull;
  z^= z >> 31;
  return z % (M - 512) + 256;
 }();
 return B;
}
inline u64 power(u64 k) {
 u64 r= 1, b= base();
 for(; k; k>>= 1, b= mul(b, b))
  if(k & 1) r= mul(r, b);
 return r;
}
// 2 MB 以上を Linux では 2 MB 境界の mmap に置いて MADV_HUGEPAGE を頼む allocator。
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
class RollingHash {
 using u64= std::uint64_t;
 int n= 0;
 std::vector<u64, rolling_hash_internal::HugeAlloc<u64>> h, pw;
 static std::size_t at(int i) { return std::size_t(i) + (std::size_t(i) >> 9); }
 template <class F> void build(F&& val) {
  using namespace rolling_hash_internal;
  const u64 B= base();
  h.resize(at(n) + 1), pw.resize(n + 1);
  h[0]= 0, pw[0]= 1;
  u64 x= 0;
  for(int i= 0; i < n; ++i) x= add(mul(x, B), val(i)), h[at(i + 1)]= x, pw[i + 1]= mul(pw[i], B);
 }
public:
 RollingHash()= default;
 explicit RollingHash(const std::string& s): n(s.size()) {
  build([&](int i) { return u64((unsigned char)s[i]); });
 }
 // 値は 2^61 - 1 で割った余りで比べる。
 template <class T> explicit RollingHash(const std::vector<T>& s): n(s.size()) {
  build([&](int i) { return u64(s[i]) % rolling_hash_internal::M; });
 }
 int size() const { return n; }
 // s[l, r) のハッシュ。0 <= l <= r <= n。
 [[gnu::always_inline]] u64 get(int l, int r) const {
  using namespace rolling_hash_internal;
  return sub(h[at(r)], mul(h[at(l)], pw[r - l]));
 }
 // s[a, a + len) と s[b, b + len) が等しいか。a + len <= n、b + len <= n。
 [[gnu::always_inline]] bool equal(int a, int b, int len) const {
  using namespace rolling_hash_internal;
  return sub(h[at(a + len)], h[at(b + len)]) == mul(sub(h[at(a)], h[at(b)]), pw[len]);
 }
 // s[a, n) と s[b, n) の最長共通接頭辞の長さ。1、2、4、... と長さを倍にして範囲を絞ってから二分探索する (LCP が L なら比べるのは
 // 2 log L 回ほど)。
 int lcp(int a, int b) const {
  using namespace rolling_hash_internal;
  const int lim= n - std::max(a, b);
  if(a == b) return lim;
  const u64 d= sub(h[at(a)], h[at(b)]);
  auto eq= [&](int m) { return sub(h[at(a + m)], h[at(b + m)]) == mul(d, pw[m]); };
  int ok= 0, ng= lim + 1;
  for(int m= 1; m <= lim; m*= 2) {
   if(!eq(m)) {
    ng= m;
    break;
   }
   ok= m;
  }
  while(ng - ok > 1) {
   const int m= (ok + ng) / 2;
   (eq(m) ? ok : ng)= m;
  }
  return ok;
 }
 // 文字列 t のハッシュ (get と同じ基数)。
 static u64 hash(const std::string& t) {
  using namespace rolling_hash_internal;
  const u64 B= base();
  u64 x= 0;
  for(unsigned char c: t) x= add(mul(x, B), c);
  return x;
 }
 template <class T> static u64 hash(const std::vector<T>& t) {
  using namespace rolling_hash_internal;
  const u64 B= base();
  u64 x= 0;
  for(const T& v: t) x= add(mul(x, B), u64(v) % M);
  return x;
 }
 // ハッシュ x の文字列の後ろに、ハッシュ y で長さ len_y の文字列をつないだもののハッシュ。O(log len_y)。
 static u64 concat(u64 x, u64 y, std::size_t len_y) {
  using namespace rolling_hash_internal;
  return add(mul(x, power(len_y)), y);
 }
};
