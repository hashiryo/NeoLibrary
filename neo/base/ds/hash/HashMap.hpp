#pragma once
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <type_traits>
#include <utility>
#ifdef __linux__
#include <sys/mman.h>
#endif
// u64 のキーから V への連想配列。開番地法の線形探索で、{キー, 値} の組を 2 の冪の個数だけ並べ、埋まりが半分を超えたら
// 2 倍にする。キーは bit を反転して持ち、0 を空きの印にするので、calloc や mmap の 0 がそのまま空きになり、表を埋める
// 手間が要らない。反転すると 0 になるキー 2^64 - 1 だけは表の外に置く。位置は、表を作るときに steady_clock から取った
// 乱数の種を足してから splitmix64 で混ぜた値の上位 bit で決める (最後の xorshift は下位 33 bit しか変えないので省く)。
// 乱数の奇数を掛けるだけのハッシュは連続する整数のキーでまれに大きく遅くなるので使わない。除くときは、後ろの組のうち
// 動かせるものを詰め直す (墓標を置かない)。2 MB 以上の表は、Linux では 2 MB 境界の mmap に置いて MADV_HUGEPAGE を頼む。
// procon-judge の yosupo-associative-array で書き比べた lp_hp を写したもの。
namespace hash_map_internal {
using u64= unsigned long long;
constexpr std::size_t H= std::size_t(1) << 21;
// 0 で埋まった領域を取る。mmap で取ったら mapped を true にする。
inline void* alloc_zero(std::size_t bytes, bool& mapped) {
 mapped= false;
#ifdef __linux__
 if(bytes >= H) {
  void* m= mmap(nullptr, bytes + H, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
  if(m != MAP_FAILED) {
   char *p= static_cast<char*>(m), *a= reinterpret_cast<char*>((reinterpret_cast<std::uintptr_t>(p) + H - 1) & ~(H - 1));
   if(a != p) munmap(p, a - p);
   munmap(a + bytes, p + H - a);
   madvise(a, bytes, MADV_HUGEPAGE);
   mapped= true;
   return a;
  }
 }
#endif
 return std::calloc(bytes, 1);
}
inline void free_zero(void* p, [[maybe_unused]] std::size_t bytes, [[maybe_unused]] bool mapped) {
#ifdef __linux__
 if(mapped) {
  munmap(p, bytes);
  return;
 }
#endif
 std::free(p);
}
}
template <class V> class HashMap {
 static_assert(std::is_trivially_copyable_v<V>, "V must be trivially copyable (a zero-filled V must equal V())");
 using u64= unsigned long long;
 struct Slot {
  u64 k;  // キーの bit 反転。0 は空き
  V v;
 };
 Slot* t;
 int shift;
 std::size_t mask, n;
 u64 seed;
 bool mapped, has_top;
 V top;  // キー 2^64 - 1 の値
 std::size_t home(u64 key) const {
  u64 x= key + seed;
  x= (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
  return ((x ^ (x >> 27)) * 0x94d049bb133111ebull) >> shift;
 }
 void rebuild(std::size_t cap) {
  Slot* old= t;
  const std::size_t oc= mask + 1;
  const bool old_mapped= mapped;
  t= static_cast<Slot*>(hash_map_internal::alloc_zero(cap * sizeof(Slot), mapped));
  mask= cap - 1, shift= 64 - __builtin_ctzll(cap);
  for(std::size_t j= 0; j < oc; ++j)
   if(old[j].k) {
    std::size_t i= home(~old[j].k);
    while(t[i].k) i= (i + 1) & mask;
    t[i]= old[j];
   }
  hash_map_internal::free_zero(old, oc * sizeof(Slot), old_mapped);
 }
 // キーの組の位置。無ければ、入れるならここに置く空きの位置。
 std::size_t find(u64 c, u64 key) const {
  std::size_t i= home(key);
  while(t[i].k && t[i].k != c) i= (i + 1) & mask;
  return i;
 }
 void release() {
  if(t) hash_map_internal::free_zero(t, (mask + 1) * sizeof(Slot), mapped), t= nullptr;
 }
 void reset() { t= static_cast<Slot*>(hash_map_internal::alloc_zero(16 * sizeof(Slot), mapped)), shift= 60, mask= 15, n= 0, has_top= false, top= V(); }
public:
 HashMap(): seed(std::chrono::steady_clock::now().time_since_epoch().count()) { reset(); }
 ~HashMap() { release(); }
 HashMap(const HashMap&)= delete;
 HashMap& operator=(const HashMap&)= delete;
 HashMap(HashMap&& o) noexcept: t(o.t), shift(o.shift), mask(o.mask), n(o.n), seed(o.seed), mapped(o.mapped), has_top(o.has_top), top(o.top) { o.t= nullptr, o.reset(); }
 HashMap& operator=(HashMap&& o) noexcept {
  if(this != &o) HashMap(std::move(o)).swap(*this);
  return *this;
 }
 void swap(HashMap& o) noexcept {
  std::swap(t, o.t), std::swap(shift, o.shift), std::swap(mask, o.mask), std::swap(n, o.n);
  std::swap(seed, o.seed), std::swap(mapped, o.mapped), std::swap(has_top, o.has_top), std::swap(top, o.top);
 }
 int size() const { return int(n + has_top); }
 bool empty() const { return size() == 0; }
 bool contains(u64 key) const {
  const u64 c= ~key;
  if(c == 0) [[unlikely]]
   return has_top;
  return t[find(c, key)].k != 0;
 }
 // key の値。無ければ d。
 V get(u64 key, V d= V()) const {
  const u64 c= ~key;
  if(c == 0) [[unlikely]]
   return has_top ? top : d;
  const Slot& s= t[find(c, key)];
  return s.k ? s.v : d;
 }
 // key の値への参照。無ければ V() で入れる。
 V& operator[](u64 key) {
  const u64 c= ~key;
  if(c == 0) [[unlikely]] {
   if(!has_top) has_top= true, top= V();
   return top;
  }
  std::size_t i= find(c, key);
  if(!t[i].k) {
   if(n + 1 > (mask + 1) / 2) rebuild((mask + 1) * 2), i= find(c, key);
   t[i].k= c, t[i].v= V(), ++n;
  }
  return t[i].v;
 }
 void set(u64 key, V v) { (*this)[key]= v; }
 // key を除く。あったら true。
 bool erase(u64 key) {
  const u64 c= ~key;
  if(c == 0) [[unlikely]] {
   const bool had= has_top;
   has_top= false;
   return had;
  }
  std::size_t i= find(c, key);
  if(!t[i].k) return false;
  // i の後ろの組のうち、ホームから i を通って今の位置に来たものを i へ詰め、空きを後ろへ送る。
  for(std::size_t j= (i + 1) & mask; t[j].k; j= (j + 1) & mask)
   if(((j - home(~t[j].k)) & mask) >= ((j - i) & mask)) t[i]= t[j], i= j;
  t[i].k= 0, --n;
  return true;
 }
 void clear() { release(), reset(); }
 // 入っている (キー, 値) を f(key, value) で 1 つずつ渡す。順番は決まらない。
 template <class F> void for_each(F&& f) const {
  for(std::size_t j= 0; j <= mask; ++j)
   if(t[j].k) f(~t[j].k, t[j].v);
  if(has_top) f(~0ull, top);
 }
};
