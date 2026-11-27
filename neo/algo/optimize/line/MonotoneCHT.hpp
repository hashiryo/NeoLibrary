#pragma once
#include <cstddef>
#include <vector>
// 傾きが単調な順に直線 y = a x + b を足しながら、点 x での最小 (MAX が true なら最大) を引く Convex Hull Trick。
// - 足す順の約束: 最小なら傾き a が増えない順、最大なら減らない順 (同じ傾きが続いてよい)。逆の順のときは、x の符号を
//   反転して (傾きを -a にして) 使う。
// - query(x) は任意の x を二分探索で引く。query_monotone(x) は、前の問いより小さくない x を、包絡の先頭から進めて引く。
// 包絡に出る直線を配列に持ち、足すたびに末尾の要らなくなった直線を除く。要るかどうかは交点の位置の比較で、切片の差と
// 傾きの差の積を __int128 で比べる。傾きと切片は別の配列に持つ (問いと追加で読む値が 8 byte ずつ詰まって並ぶ)。
// T は 64 bit までの符号付き整数で、a x + b が T に収まること。直線を 1 本も足さないうちの問いは未定義。最大は、傾きと
// 切片の符号を反転して最小で持つ。
// 傾きの順を問わないときや線分を足すときは LiChaoTree を使う。傾きが単調なら、こちらが 3〜6 倍速い。
// procon-judge の atcoder-dp-z と yuki-409 で書き比べた核 (problems/_shared/lines/cht_deque.hpp) を写したもの。記録は
// algo-notes の notes/dp-optimization.md にある。
template <class T= long long, bool MAX= false> class MonotoneCHT {
 std::vector<T> as_, bs_;  // 包絡に出る直線の傾きと切片
 std::size_t head_= 0;     // query_monotone がこれより前を見ない
 // a1 > a2 > a3 のとき、l2 が包絡に出ないか。
 static bool hidden(T a1, T b1, T a2, T b2, T a3, T b3) { return (__int128)(b3 - b1) * (a1 - a2) <= (__int128)(b2 - b1) * (a1 - a3); }
public:
 void reserve(std::size_t n) { as_.reserve(n), bs_.reserve(n); }
 [[gnu::always_inline]] void add_line(T a, T b) {
  if constexpr(MAX) a= -a, b= -b;
  if(!as_.empty() && as_.back() == a) {
   if(bs_.back() <= b) return;
   as_.pop_back(), bs_.pop_back();
  }
  while(as_.size() >= 2 && hidden(as_[as_.size() - 2], bs_[bs_.size() - 2], as_.back(), bs_.back(), a, b)) as_.pop_back(), bs_.pop_back();
  as_.push_back(a), bs_.push_back(b);
  // query_monotone が指していた直線が除かれたら、残った末尾の 1 つ前から見直す。
  if(head_ + 2 > as_.size()) head_= as_.size() >= 2 ? as_.size() - 2 : 0;
 }
 // 前の問いより小さくない x での最小 (MAX なら最大)。
 [[gnu::always_inline]] T query_monotone(T x) {
  while(head_ + 1 < as_.size() && as_[head_ + 1] * x + bs_[head_ + 1] <= as_[head_] * x + bs_[head_]) ++head_;
  const T v= as_[head_] * x + bs_[head_];
  return MAX ? -v : v;
 }
 // 任意の x での最小 (MAX なら最大)。
 [[gnu::always_inline]] T query(T x) const {
  std::size_t lo= 0, hi= as_.size() - 1;
  while(lo < hi) {
   const std::size_t mid= (lo + hi) >> 1;
   if(as_[mid] * x + bs_[mid] <= as_[mid + 1] * x + bs_[mid + 1]) hi= mid;
   else lo= mid + 1;
  }
  const T v= as_[lo] * x + bs_[lo];
  return MAX ? -v : v;
 }
};
