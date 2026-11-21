#pragma once
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <numeric>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 列の空でない相異なる部分列の数。空も含めた a[0, i) の相異なる部分列の数を dp[i] とすると、dp[i + 1] = 2 dp[i] - dp[j] (j は i より
// 前で a[i] と同じ値が最後に出た位置、無ければ引かない) で、答えは dp[n] - 1。値の幅 (最大 - 最小 + 1) が要素の数の 2 倍か 2^16 の
// 大きい方以下なら、値から直接引く表に、その値が最後に出たときの dp の値 + 1 を持って 1 回の走査で済ませる。幅が広ければ (値 - 最小,
// 添字) の組を 11 bit ずつの基数ソートで並べて前の位置を求め、全要素で同じになる桁の段 (値の OR と AND で分かる) は飛ばす。組の列は
// 0 で埋まった 1 つの領域 (2 MB 以上なら、Linux では 2 MB 境界の mmap に置いて huge page を頼む) に置き、前の位置と dp の値は並べ終えた
// 側でない組の列の領域に置く。幅が 2^32 以上 (64 bit の値) なら std::sort で番号を振り直す。procon-judge の
// yosupo-number-of-subsequences で書き比べた nos_adapt2 を写したもの (幅の狭い列の表は nos_adapt と同じ vector)。
namespace distinct_subsequences_internal {
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
// v(i) = a[i] - 最小 (0 以上 range 未満、range <= 2^32) の列から答えを求める。
template <class V> unsigned solve(V v, int n, std::uint64_t range, unsigned mod) {
 using u64= std::uint64_t;
 if(range <= std::max<u64>(2 * u64(n), 1 << 16)) {
  std::vector<unsigned> last(range, 0);
  unsigned d= 1;
  for(int i= 0; i < n; ++i) {
   unsigned& r= last[v(i)];
   const unsigned before= d;
   d*= 2;
   d-= d >= mod ? mod : 0;
   if(r) {
    d+= mod - (r - 1);
    d-= d >= mod ? mod : 0;
   }
   r= before + 1;
  }
  return d == 0 ? mod - 1 : d - 1;
 }
 unsigned all_or= 0, all_and= ~0u;
 for(int i= 0; i < n; ++i) {
  const unsigned x= v(i);
  all_or|= x, all_and&= x;
 }
 const unsigned diff= all_or ^ all_and;  // 立っている桁だけが並べる意味を持つ
 Zeroed buf(std::size_t(n) * 2 * sizeof(u64));
 u64 *p= static_cast<u64*>(buf.get()), *q= p + n;  // 上位 32 bit が値 - 最小、下位 32 bit が添字
 std::array<std::array<int, 2048>, 3> cnt{};
 for(int i= 0; i < n; ++i) {
  const unsigned x= v(i);
  p[i]= u64(x) << 32 | unsigned(i);
  ++cnt[0][x & 2047], ++cnt[1][(x >> 11) & 2047], ++cnt[2][x >> 22];
 }
 for(int j= 0; j < 3; ++j) {
  const int sh= 11 * j;
  if(!((diff >> sh) & 2047)) continue;
  auto& c= cnt[j];
  for(int b= 0, acc= 0; b < 2048; ++b) {
   const int t= c[b];
   c[b]= acc, acc+= t;
  }
  for(int i= 0; i < n; ++i) q[c[(p[i] >> 32 >> sh) & 2047]++]= p[i];
  std::swap(p, q);
 }
 // 並べ終えた p でない側 q の領域に、前の位置 prev (n 個) と dp の値 D (D[i] は a[0, i) の空も含めた相異なる部分列の数、i < n) を置く
 int* prev= reinterpret_cast<int*>(q);
 unsigned* D= reinterpret_cast<unsigned*>(q) + n;
 prev[unsigned(p[0])]= -1;
 for(int k= 1; k < n; ++k) prev[unsigned(p[k])]= (p[k] >> 32) == (p[k - 1] >> 32) ? int(unsigned(p[k - 1])) : -1;
 unsigned d= 1;
 D[0]= 1;
 for(int i= 0; i < n; ++i) {
  unsigned x= d * 2;
  x-= x >= mod ? mod : 0;
  if(prev[i] >= 0) {
   x+= mod - D[prev[i]];
   x-= x >= mod ? mod : 0;
  }
  d= x;
  if(i + 1 < n) D[i + 1]= d;
 }
 return d == 0 ? mod - 1 : d - 1;
}
}
// 列 a の空でない相異なる部分列の数を mod で割った余り。T は整数の型で、mod は 1 以上 2^31 未満。O(n) (値の幅が 2^32 以上なら
// O(n log n))。
template <class T> unsigned distinct_subsequences(const std::vector<T>& a, unsigned mod= 998244353) {
 static_assert(std::is_integral_v<T>);
 using U= std::make_unsigned_t<T>;
 const int n= a.size();
 if(n == 0) return 0;
 if(mod == 1) return 0;
 T mn= a[0], mx= a[0];
 for(const T x: a) mn= std::min(mn, x), mx= std::max(mx, x);
 const std::uint64_t range= std::uint64_t(U(mx) - U(mn)) + 1;
 if(range != 0 && range <= (std::uint64_t(1) << 32)) {
  const T* p= a.data();
  const U lo= U(mn);
  return distinct_subsequences_internal::solve([p, lo](int i) { return unsigned(U(p[i]) - lo); }, n, range, mod);
 }
 std::vector<int> ord(n);
 std::iota(ord.begin(), ord.end(), 0);
 std::sort(ord.begin(), ord.end(), [&](int x, int y) { return a[x] < a[y]; });
 std::vector<std::uint32_t> v(n);
 std::uint32_t id= 0;
 for(int k= 0; k < n; ++k) v[ord[k]]= k > 0 && a[ord[k]] != a[ord[k - 1]] ? ++id : id;
 return distinct_subsequences_internal::solve([&v](int i) { return v[i]; }, n, std::uint64_t(id) + 1, mod);
}
// 文字列 s の空でない相異なる部分列の数を mod で割った余り。文字は unsigned char として比べる。
inline unsigned distinct_subsequences(const std::string& s, unsigned mod= 998244353) {
 if(s.empty() || mod == 1) return 0;
 const unsigned char* p= reinterpret_cast<const unsigned char*>(s.data());
 return distinct_subsequences_internal::solve([p](int i) { return unsigned(p[i]); }, int(s.size()), 256, mod);
}
