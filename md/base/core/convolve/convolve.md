---
title: 畳み込みの共通の口
documentation_of: ../../neo/base/core/convolve/convolve.hpp
---

積を求める核をテンプレート引数で受けて、積の全体、積の一部の範囲、middle product を求めます。このヘッダは核を含みません。核のヘッダを別に include して使います。名前は仮です。

```cpp
template <class Conv, class K> std::vector<K> convolve(const std::vector<K>& a, const std::vector<K>& b);
template <class Conv, class K> std::vector<K> convolve_range(const std::vector<K>& a, const std::vector<K>& b, int lo, int hi);
template <class Conv, class K> std::vector<K> middle_product(const std::vector<K>& a, const std::vector<K>& b);
```

| 関数 | 返すもの |
| --- | --- |
| `convolve<Conv>(a, b)` | 積 $c_k = \sum_{i+j=k} a_i b_j$ の全体。どちらかが空なら空 |
| `convolve_range<Conv>(a, b, lo, hi)` | 積の `lo` 次から `hi - 1` 次の係数。積の次数を超える部分は 0 |
| `middle_product<Conv>(a, b)` | $c_i = \sum_j a_{i+j} b_j$ ($0 \le i \le n - m$、$n$ と $m$ は `a` と `b` の長さ)。`b` が空か、`a` が `b` より短ければ空 |

`Conv` は核の構造体です。今あるのは、素朴な積の `Naive` (`neo/base/core/convolve/naive.hpp`) です。

```cpp
#include "neo/base/core/convolve/convolve.hpp"
#include "neo/base/core/convolve/naive.hpp"
std::vector<ZMod<998244353>> c= convolve<Naive>(a, b);
std::vector<ZMod<998244353>> d= convolve_range<Naive>(a, b, 0, 10);  // 先頭 10 項
```

ポインタで受ける関数もあります。ほかの関数の中の繰り返しで、配列を確保せずに呼ぶためのものです。

```cpp
template <class Conv, class K> void convolve_range(const K* a, int n, const K* b, int m, int lo, int hi, K* out);
template <class Conv, class K> void middle_product(const K* a, int n, const K* b, int m, K* out);
```

`out` は、`a` や `b` と重ならない配列にします。`convolve_range` は、範囲に届かない項を落としてから核を呼びます。先頭 `len` 項だけを求めるときは、`a` と `b` の `len` 次以上の項を読みません。

## 核の選び方

`Conv` を省くと、係数の型 `K` の既定の核を使います。既定は、核のヘッダが宣言する `default_kernel(const K*)` の戻り値の型で決まります。

- include した核のうち、`K` を受け持つものが 1 つなら、それが既定になります。
- 1 つも無いか、2 つ以上あると、コンパイルエラーになります。そのときは `Conv` を書いて選びます。
- `Naive` は既定になりません。核を include し忘れたときに、$O(nm)$ の積で黙って動くのを避けるためです。

同じ係数の型で 2 つの核を使い分けるとき (核どうしを比べるときなど) は、どちらも include して `Conv` を書きます。

## 核を書くとき

核は、次の static 関数を持つ構造体です。

```cpp
struct MyKernel {
 // out[k - lo] = Σ_{i+j=k} a[i] b[j] (lo <= k < hi)
 template <class K> static void convolve_range(const K* a, int n, const K* b, int m, int lo, int hi, K* out);
};
MyKernel default_kernel(const MyType*);  // MyType の既定の核にするときだけ宣言する (定義は要らない)
```

この関数は `n >= 1`、`m >= 1`、`0 <= lo < hi <= n + m - 1` で呼ばれます。`a` と `b` が同じ配列 (`a == b` かつ `n == m`) のこともあります。
