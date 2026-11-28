---
title: 区分線形な凸関数 (slope trick)
documentation_of: ../../neo/algo/optimize/slope/PiecewiseLinearConvex.hpp
---

区分線形な凸関数 f を持ち、足し算、平行移動、累積最小、窓の最小、最小値と argmin の問いを載せます。

```cpp
template <class T = long long> class PiecewiseLinearConvex {
  PiecewiseLinearConvex();                                                   // f(x) = 0
  explicit PiecewiseLinearConvex(const std::vector<std::pair<T, T>>& ramps);  // f(x) = sum max(0, a(x-x0))
  void add_const(D c);                  // f(x) += c
  void add_linear(T a);                 // f(x) += a x
  void add_abs(T a, T x0);              // f(x) += a|x-x0|
  void add_ramp(T a, T x0);             // f(x) += max(0, a(x-x0))
  void add_max(T a, T b, T x0);         // f(x) += max(a(x-x0), b(x-x0)) (a < b)
  void add_inf(bool right, T x0);       // x < x0 (right=false) か x > x0 (right=true) で +inf
  void subtract_abs(T a, T x0);         // f(x) -= a|x-x0|
  void subtract_ramp(T a, T x0);        // f(x) -= max(0, a(x-x0))
  void subtract_max(T a, T b, T x0);    // f(x) -= max(a(x-x0), b(x-x0)) (a < b)
  void shift(T x0);                     // f(x) <- f(x-x0)
  void chmin_cum(bool rev = false);     // f(x) <- min_{y<=x} f(y)、rev=true で min_{x<=y}
  void chmin_slide_win(T lb, T ub);     // f(x) <- min_{lb<=y<=ub} f(x-y)
  std::optional<D> operator()(T x);     // f(x)。定義域の外は nullopt
  std::optional<D> min();               // 最小値。下に有界でなければ nullopt
  std::array<T, 2> argmin();            // 最小を取る x の区間。解が無ければ {1, 0}
  PiecewiseLinearConvex& operator+=(const PiecewiseLinearConvex& g);
  PiecewiseLinearConvex operator+(PiecewiseLinearConvex g) const;
};
```

`T` は 32 bit か 64 bit の符号付き整数です。値の計算は、`T` が 32 bit なら 64 bit、64 bit なら 128 bit (`D`) で行います。`add_inf` で定義域を狭めると、外側では `operator()` と `min()` が nullopt を返します。

## 使い方

```cpp
// 列を狭義単調増加にする最小の変更量: f(x) = min_{y<=x} f(y) + |x - a_i|
PiecewiseLinearConvex<long long> f;
for (long long v : a) {
  f.chmin_cum();
  f.add_abs(1, v);
}
long long ans = (long long)f.min().value();
```

## 中身

折れ点と、そこで傾きが増える量を treap に載せ、最小を取る折れ点の左の部分木と右の部分木に分けて持ちます。各部分木は x に足す遅延を持ちます。足し算は 2 つの treap の併合 (`unite`) です。

今の Library の PiecewiseLinearConvex から、永続 (`persistent`) とクラス static のノードプールを外したものです。ノードはインスタンスごとの配列に置くので、`reset` や `rebuild` は要りません。操作の範囲は、procon-judge の 33 問の提出が使うものに合わせています。

書き比べの記録は、algo-notes の notes/dp-optimization.md にあります。
