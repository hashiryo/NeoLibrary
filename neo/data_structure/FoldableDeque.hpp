#pragma once
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <new>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 両端に足し、両端から消し、列の全体の積を返す両端キュー (両端のスライド窓の集約)。構造体 M は SegmentTree と同じモノイド
// (using T、static T ti()、static T op(T a, T b)) で、可換でなくてよい。
// 1 本の配列の [hd, tl) に列を持ち、境 md の前の [hd, md) には各位置から md までの積を、後ろの [md, tl) には md から各位置
// までの積を、値と並べて置く (2 本のスタックを向かい合わせに置いた形)。片方が空のときにそちらから消すと、境を列の真ん中へ
// 移して両方の積をその場で数え直す (2 本のスタックで半分を移す形と同じで、ならして O(1)。値を別の配列へ写さない)。
namespace foldable_deque_internal {
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
}  // namespace foldable_deque_internal
template <class M> class FoldableDeque {
 using T= typename M::T;
 struct Node {
  T x, s;
 };
 std::vector<Node, foldable_deque_internal::Alloc<Node>> a;
 int hd= 0, md= 0, tl= 0;
 // 境を m に移し、前と後ろの積を数え直す。
 void rebuild(int m) {
  md= m;
  if(hd < md) {
   a[md - 1].s= a[md - 1].x;
   for(int i= md - 1; --i >= hd;) a[i].s= M::op(a[i].x, a[i + 1].s);
  }
  if(md < tl) {
   a[md].s= a[md].x;
   for(int i= md + 1; i < tl; ++i) a[i].s= M::op(a[i - 1].s, a[i].x);
  }
 }
 // 端に空きが無いとき。生きている要素が半分以下なら同じ配列の真ん中へ、そうでなければ倍の配列の真ん中へ移す。
 void recenter() {
  const int k= tl - hd, c= a.size(), nc= c && 2 * k <= c ? c : std::max(16, 2 * c), nh= (nc - k) / 2;
  if(nc == c) {
   if(nh < hd) std::move(a.begin() + hd, a.begin() + tl, a.begin() + nh);
   else std::move_backward(a.begin() + hd, a.begin() + tl, a.begin() + nh + k);
  } else {
   std::vector<Node, foldable_deque_internal::Alloc<Node>> b(nc);
   std::move(a.begin() + hd, a.begin() + tl, b.begin() + nh);
   a.swap(b);
  }
  md+= nh - hd, tl+= nh - hd, hd= nh;
 }
public:
 FoldableDeque()= default;
 int size() const { return tl - hd; }
 bool empty() const { return hd == tl; }
 void clear() { hd= md= tl= int(a.size()) / 2; }
 void push_front(const T& x) {
  if(hd == 0) recenter();
  --hd;
  a[hd].x= x, a[hd].s= hd + 1 == md ? x : M::op(x, a[hd + 1].s);
 }
 void push_back(const T& x) {
  if(tl == int(a.size())) recenter();
  a[tl].x= x, a[tl].s= tl == md ? x : M::op(a[tl - 1].s, x);
  ++tl;
 }
 void pop_front() {
  assert(hd < tl);
  if(hd == md) rebuild(hd + (tl - hd + 1) / 2);
  ++hd;
 }
 void pop_back() {
  assert(hd < tl);
  if(md == tl) rebuild(hd + (tl - hd) / 2);
  --tl;
 }
 const T& front() const { return a[hd].x; }
 const T& back() const { return a[tl - 1].x; }
 const T& operator[](int i) const { return a[hd + i].x; }
 // op(先頭, ..., 末尾)。空なら ti()。
 T all_prod() const {
  if(hd == md) return md == tl ? M::ti() : a[tl - 1].s;
  return md == tl ? a[hd].s : M::op(a[hd].s, a[tl - 1].s);
 }
};
