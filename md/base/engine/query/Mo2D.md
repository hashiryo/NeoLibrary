---
title: 二次元の Mo's algorithm
documentation_of: ../../neo/base/engine/query/Mo2D.hpp
---

部分行列のクエリをオフラインで処理する Mo's algorithm です。$xl \le x < xr$、$yl \le y < yr$ の長方形を 4 つの端点で表し、端点を動かすたびに行か列を足し引きします。

```cpp
struct Mo2D {
  void query(int xl, int xr, int yl, int yr);
  int size() const;
  bool empty() const;

  template <class AX, class EX, class AY, class EY, class Out>
  void run(const AX& add_x, const EX& erase_x, const AY& add_y, const EY& erase_y, const Out& out) const;
};
```

`add_x(x, yl, yr)` は行 `x` を列の範囲 `yl <= y < yr` で足します。`erase_x` は外します。`add_y(y, xl, xr)` は列 `y` を行の範囲 `xl <= x < xr` で足し、`erase_y` は外します。行や列を足すときは、そのときのもう一方の範囲が渡るので、その範囲を走査して状態を更新します。

```cpp
// 部分行列の和
std::vector<std::vector<int>> a(n, std::vector<int>(m));
Mo2D mo;
for (int i = 0; i < q; ++i) {
  int xl, xr, yl, yr;
  std::cin >> xl >> xr >> yl >> yr;
  mo.query(xl, xr, yl, yr);
}
long long sum = 0;
std::vector<long long> ans(q);
auto add_x = [&](int x, int yl, int yr) {
  for (int y = yl; y < yr; ++y) sum += a[x][y];
};
auto erase_x = [&](int x, int yl, int yr) {
  for (int y = yl; y < yr; ++y) sum -= a[x][y];
};
auto add_y = [&](int y, int xl, int xr) {
  for (int x = xl; x < xr; ++x) sum += a[x][y];
};
auto erase_y = [&](int y, int xl, int xr) {
  for (int x = xl; x < xr; ++x) sum -= a[x][y];
};
auto out = [&](int id) { ans[id] = sum; };
mo.run(add_x, erase_x, add_y, erase_y, out);
```

並べ方は $xl$ と $yl$ のブロックで、$xr$ と $yr$ を交互に並べます。移動回数は 4 端点の動きで決まり、1 回の移動で走査する範囲の長さが掛かります。行や列の足し引きが $O(1)$ で済む問題では、ブロックの幅を $N/\sqrt{Q}$ にして $O(N^2\sqrt{Q})$ 程度になります。
