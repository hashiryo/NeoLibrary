---
title: 半環
documentation_of: ../../neo/algebra/semiring.hpp
---

行列 (`neo/linalg/Matrix.hpp`、`neo/linalg/SquareMatrix.hpp`) の要素やセグメント木の値に使う、半環の値の型です。和 (`+`)、積 (`*`)、静的メンバ関数 `zero()` と `one()` を持ちます。

```cpp
template <class T> struct MinPlus;   // 和は min、積は +。0 は INF、1 は 0
template <class T> struct MaxPlus;   // 和は max、積は +。0 は -INF、1 は 0
template <class T> struct MaxMin;    // 和は max、積は min。0 は T の最小、1 は T の最大
struct OrAnd;                        // 和は or、積は and。0 は false、1 は true
```

どの型も値を公開のメンバ `x` に持ち、`T` (`OrAnd` は bool) から暗黙に作れます。`+=`、`*=`、`==`、`!=` があり、どれも constexpr です。既定の構築は 0 (和の単位元) です。

用途は、`MinPlus` が最短路、`MaxPlus` が最長路、`MaxMin` が道の幅の最大、`OrAnd` が到達です。たとえば `Matrix<MinPlus<long long>>` の $k$ 乗は、ちょうど $k$ 本の辺を通る最短路の長さになります。

整数の INF は `T` の最大の半分 (`MinPlus<T>::INF`) です。`MinPlus` の積は、片方が INF ならそのまま INF にし、和が INF 以上になれば INF に丸めます。負の重みがあっても INF は INF のままです。`MaxPlus` も符号を逆にして同じです。浮動小数点数の INF は infinity で、積はそのまま足します。

INF の扱いのほかに丸めはしないので、有限の値どうしの和が `T` の範囲を出ないように型を選んでください。

負の重みがあると、INF に丸めるせいで、積の結合法則は成り立たないことがあります。$(A B) C$ と $A (B C)$ で、途中で INF を超えた和に負の重みを足す順が変わるためです。

```cpp
static constexpr bool fast_range(const MinPlus* p, const MinPlus* q, std::size_t n);
static constexpr MinPlus fast_mul(MinPlus l, MinPlus r);   // MaxPlus にも同じもの
```

行列の積の速い道です。`fast_range` は、p と q の n 個の値が、INF の比べを省いた積 `fast_mul` でも同じ値になる範囲 (`MinPlus` は 0 以上 INF 以下、`MaxPlus` は -INF 以上 0 以下) にあるかを返します。`Matrix` と `SquareMatrix` の積は、両方の行列がこの範囲にあれば `fast_mul` で掛けます。

## 中身

整数の `MinPlus` の積は、$\max(x, y, x + y) \ge \mathrm{INF}$ なら INF、そうでなければ $x + y$ です。和は符号なしの算術で先に求め、INF にするかはマスクの算術で選びます。比べを条件式で書くと、x64 の GCC 15 はデータで向きの変わる分岐にしていました。`MaxPlus` も符号を逆にして同じ形です。

この形でも、x64 の GCC 15 は 3 × 3 の積をベクトル化しません。INF の比べを持たない `min(x + y, INF)` なら、i-k-j の列の向きを 2 個と 1 個に分けて vpaddd と vpminsd で掛けます。llvm-mca の Zen 3 の見積もりでは、3 × 3 の積 1 回が 27 サイクルと 70 サイクルの差でした。そこで、値がすべて 0 以上 INF 以下なら、`min(x + y, INF)` が比べのある積と同じ値になることを使います。INF は $2^k - 1$ の形なので、範囲に入るかは値の or 1 回で分かります。`fast_range` と `fast_mul` はこのための関数です。

procon-judge の cf-1609-e は、3 × 3 の min-plus の行列をセグメント木に載せる問題です。EPYC 7763 の x64-gcc で、比べを分岐で書いていた形の lib-neo は 86.67 ms (run 40306505343)、分岐なしにした形は 68.06 ms (run 40310861931) でした。INF の比べを持たない minplus_fixed は 39.32 ms です。速い道を足したあとの手元の M2 では、lib-neo と minplus_fixed が 38.8 ms から 40.2 ms と 37.3 ms で並びます。x64 では測り直していません。
