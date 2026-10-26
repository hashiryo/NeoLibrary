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

## 中身

x64 の GCC 15 は、`MinPlus` の積の INF の比べを分岐にします。マスクの算術で分岐をなくした形も試しましたが、cf-1609-e の形 (3 × 3 の行列のセグメント木) では手元の M2 の clang で 7% から 14% 遅くなりました。clang は元の形を分岐のない命令に直していたので、元の形のままにしています。x64 での差はまだ測っていません。
