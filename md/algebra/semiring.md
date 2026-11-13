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
static constexpr bool fast_one_sided(const MinPlus* p, std::size_t n);   // MaxPlus にも同じもの
static constexpr bool fast_one_sided(const MinPlus* p, const MinPlus* q, std::size_t n);
static constexpr bool fast_bounded(const MinPlus* p, std::size_t n);
static constexpr MinPlus fast_mul(MinPlus l, MinPlus r);
static constexpr MinPlus fast_fix(MinPlus v, bool bounded);
static void fast_fix_all(MinPlus* p, std::size_t n, bool bounded);
```

行列の積の速い道です。積は INF の比べのない和 `fast_mul` ($x + y$) とし、和 (min か max) で足します。足し終えた成分は `fast_fix` で INF へ直します。掛ける 2 つの行列の値がどれも `fast_one_sided` の範囲にあるか、どれも `fast_bounded` の範囲にあれば、厳密な積と同じ値になります。`fast_one_sided` の範囲は、`MinPlus` なら 0 以上 INF 以下、`MaxPlus` なら -INF 以上 0 以下です。`fast_bounded` の範囲は ±INF か、絶対値が INF / 4 以下の値で、負の重みと正の重みが混ざっていてもかまいません。`bounded` には、`fast_bounded` の範囲で掛けたかを渡します。`fast_fix_all` は、配列の n 個に `fast_fix` を当てます。2 つの配列を受ける `fast_one_sided` は、p と q の n 個ずつを 1 本のループで調べます (x86 で AVX2 が使えれば、先頭から 32 byte ずつの塊を vptest で調べます)。`Matrix` と `SquareMatrix` の積は、まず `fast_one_sided` を調べ、外れたときだけ `fast_bounded` を調べます。

## 中身

整数の `MinPlus` の積は、$\max(x, y, x + y) \ge \mathrm{INF}$ なら INF、そうでなければ $x + y$ です。和は符号なしの算術で先に求め、INF にするかはマスクの算術で選びます。比べを条件式で書くと、x64 の GCC 15 はデータで向きの変わる分岐にしていました。`MaxPlus` も符号を逆にして同じ形です。

この形でも、x64 の GCC 15 は 3 × 3 の積をベクトル化しません。INF の比べを持たない `min(x + y, INF)` なら、i-k-j の列の向きを 2 個と 1 個に分けて vpaddd と vpminsd で掛けます。llvm-mca の Zen 3 の見積もりでは、3 × 3 の積 1 回が 27 サイクルと 70 サイクルの差でした。そこで、値がすべて 0 以上 INF 以下なら、`min(x + y, INF)` が比べのある積と同じ値になることを使います。INF は $2^k - 1$ の形なので、範囲に入るかは値の or 1 回で分かります。min は単調なので、和 $x + y$ の min を取ってから INF 以上を INF にしても同じです。

負の重みがあるときは、値がどれも INF か、絶対値が $L = \mathrm{INF} / 4$ 以下なら速い道に載せます。有限の値どうしの和は絶対値が $2L$ 以下で、INF / 2 未満です。INF を含む和は $\mathrm{INF} - L$ 以上で、$2\,\mathrm{INF}$ 以下なので `T` の範囲を出ません。和の min を取ってから INF / 2 以上を INF に直せば、厳密な積の min と同じになります。手元の M2 では、300 次の `MaxPlus<int>` の積が、両側の重みで 2.81 ms、片側の重みで 2.79 ms、範囲を外れた値 (厳密な道) で 4.49 ms でした。

`Matrix` の積の 1 行分は、行と行が重ならないことを `__restrict` で伝えます。伝えないと、x64 の GCC 15 は 64 bit の値の min のループをベクトル化しません。`SquareMatrix` の片側の範囲では、和を最後に直さず、項ごとに `fast_fix(fast_mul(x, y), false)` ($\min(x + y, \mathrm{INF})$) を足します。最後に直す形は、手元 (arm の clang) の cf-1609-e で 15% 遅くなりました。

procon-judge の cf-1609-e は、3 × 3 の min-plus の行列をセグメント木に載せる問題です。EPYC 7763 の x64-gcc で、比べを分岐で書いていた形の lib-neo は 86.67 ms (run 40306505343)、分岐なしにした形は 68.06 ms (run 40310861931) でした。INF の比べを持たない minplus_fixed は 39.32 ms です。速い道を足したあとは 48.28 ms で (run 42352552289)、minplus_fixed (39.36 ms) との差の大半は積ごとの範囲の判定でした。判定を省いた写しは 39.97 ms です。

2 つの配列を受ける `fast_one_sided` は、x86 で AVX2 が使えれば、先頭から 32 byte ずつの塊を or して vptest で調べ、残りをスカラーの or で調べます。スカラーの or のループだけでは、x64 の GCC 15 は 32 byte の or を横に畳む 7 命令を付けて組み、3 × 3 の行列の判定が 14 命令ほどになっていました。vptest の形では 8 命令です。同じ回の写しで、x64-gcc は 48.28 ms から 43.32 ms に、x64-clang は 41.09 ms から 33.66 ms になりました。最後の塊を末尾から重ねて読む形も試しましたが、3 × 3 の行列では直前の積の 32 byte と 4 byte の書き込みにまたがり、書き込みからの受け渡しが効かずに 53.62 ms と遅くなりました。arm ではスカラーの or のままにします (SIMDe の形は arm の clang で 6% 遅くなりました)。

`fast_bounded` と `fast_fix_all` も、x86 で AVX2 が使えれば 8 個 (64 bit なら 4 個) ずつ調べて直します。スカラーのループのままでは、x64 の GCC 15 はどちらもベクトル化せず、`fast_bounded` が 1 個に 11 命令、`fast_fix` の一巡が 1 個に 8 命令かかっていました。procon-judge の cf-147-b では、辺 300 の積 1 回ごとに 2 つの行列の 18 万個を調べ、9 万個を直すので、写し元のタイルとの差の大半はここと見ています。
