#pragma once
#include <cstddef>
#include <vector>
// 傾きが単調な順に直線 y = a x + b を足しながら、点 x での最小 (MAX が true なら最大) を引く Convex Hull Trick。
// - 足す順の約束: 最小なら傾き a が増えない順、最大なら減らない順 (同じ傾きが続いてよい)。逆の順のときは、x の符号を
//   反転して (傾きを -a にして) 使う。
// - query(x) は任意の x を二分探索で引く。query_monotone(x) は、前の問いより小さくない x を、包絡の先頭から進めて引く。
// 包絡に出る直線を配列に持ち、足すたびに末尾の要らなくなった直線を除く。要るかどうかは交点の位置の比較で、切片の差と
// 傾きの差の積を __int128 で比べる。T は 64 bit までの符号付き整数で、a x + b が T に収まること。直線を 1 本も足さないうちの
// 問いは未定義。最大は、傾きと切片の符号を反転して最小で持つ。
// 傾きの順を問わないときや線分を足すときは LiChaoTree を使う。傾きが単調なら、こちらが 3〜6 倍速い。
// procon-judge の atcoder-dp-z と yuki-409 で書き比べた核 (problems/_shared/lines/cht_deque.hpp) を写したもの。記録は
// algo-notes の notes/dp-optimization.md にある。
template <class T= long long, bool MAX= false> class MonotoneCHT {
 struct Line {
  T a, b;
 };
 std::vector<Line> h_;
 std::size_t head_= 0;  // query_monotone がこれより前を見ない
 // a1 > a2 > a3 のとき、l2 が包絡に出ないか。
 static bool hidden(const Line& l1, const Line& l2, const Line& l3) { return (__int128)(l3.b - l1.b) * (l1.a - l2.a) <= (__int128)(l2.b - l1.b) * (l1.a - l3.a); }
 static T eval(const Line& e, T x) { return e.a * x + e.b; }
public:
 void reserve(std::size_t n) { h_.reserve(n); }
 [[gnu::always_inline]] void add_line(T a, T b) {
  if constexpr(MAX) a= -a, b= -b;
  if(!h_.empty() && h_.back().a == a) {
   if(h_.back().b <= b) return;
   h_.pop_back();
  }
  const Line l{a, b};
  while(h_.size() >= 2 && hidden(h_[h_.size() - 2], h_.back(), l)) h_.pop_back();
  h_.push_back(l);
  // query_monotone が指していた直線が除かれたら、残った末尾の 1 つ前から見直す。
  if(head_ + 2 > h_.size()) head_= h_.size() >= 2 ? h_.size() - 2 : 0;
 }
 // 前の問いより小さくない x での最小 (MAX なら最大)。
 [[gnu::always_inline]] T query_monotone(T x) {
  while(head_ + 1 < h_.size() && eval(h_[head_ + 1], x) <= eval(h_[head_], x)) ++head_;
  const T v= eval(h_[head_], x);
  return MAX ? -v : v;
 }
 // 任意の x での最小 (MAX なら最大)。
 [[gnu::always_inline]] T query(T x) const {
  std::size_t lo= 0, hi= h_.size() - 1;
  while(lo < hi) {
   const std::size_t mid= (lo + hi) >> 1;
   if(eval(h_[mid], x) <= eval(h_[mid + 1], x)) hi= mid;
   else lo= mid + 1;
  }
  const T v= eval(h_[lo], x);
  return MAX ? -v : v;
 }
};
