#pragma once
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <new>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 末尾に足し、先頭から消し、列の全体の積を返すキュー (スライド窓の集約、SWAG)。構造体 M は SegmentTree と同じモノイド
// (using T、static T ti()、static T op(T a, T b)) で、可換でなくてよい。
// 1 本の配列の [hd, tl) に列を持ち、境 md の前の [hd, md) には各位置から md までの積を、後ろの [md, tl) には足した値を
// そのまま置き、後ろの積を ba に持つ。先頭を消すときに前が空なら、後ろの値をその場で後ろから積に書き換えて境を tl に
// 移す (2 本のスタックで持つ形と同じで、ならして O(1)。値を別の配列へ写さない)。
namespace foldable_queue_internal {
// 1 MB 以上の確保は 2 MB 境界に寄せ、Linux では huge page を頼む。
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
  const std::size_t want= n * sizeof(T);
  if(want < (std::size_t(1) << 20)) ::operator delete(p, std::align_val_t(64));
  else ::operator delete(p, std::align_val_t(std::size_t(1) << 21));
 }
 template <class U> bool operator==(const Alloc<U>&) const { return true; }
};
}  // namespace foldable_queue_internal
template <class M> class FoldableQueue {
 using T= typename M::T;
 std::vector<T, foldable_queue_internal::Alloc<T>> a;
 int hd= 0, md= 0, tl= 0;
 T ba= M::ti();
 // 末尾に空きが無いとき。生きている要素が半分以下なら先頭へ詰め、そうでなければ倍の配列へ移す。
 void grow() {
  const int k= tl - hd, c= a.size();
  if(c && 2 * k <= c) std::move(a.begin() + hd, a.begin() + tl, a.begin());
  else {
   std::vector<T, foldable_queue_internal::Alloc<T>> b(std::max(16, 2 * c));
   std::move(a.begin() + hd, a.begin() + tl, b.begin());
   a.swap(b);
  }
  md-= hd, tl-= hd, hd= 0;
 }
public:
 FoldableQueue()= default;
 int size() const { return tl - hd; }
 bool empty() const { return hd == tl; }
 void clear() { hd= md= tl= 0, ba= M::ti(); }
 void push(const T& x) {
  if(tl == int(a.size())) grow();
  a[tl++]= x, ba= M::op(ba, x);
 }
 void pop() {
  assert(hd < tl);
  if(hd == md) {
   // a[md] は書き換えてすぐ消すので数えない。
   for(int i= tl - 1; --i > md;) a[i]= M::op(a[i], a[i + 1]);
   md= tl, ba= M::ti();
  }
  ++hd;
 }
 // op(先頭, ..., 末尾)。空なら ti()。
 T all_prod() const { return hd == md ? ba : M::op(a[hd], ba); }
};
