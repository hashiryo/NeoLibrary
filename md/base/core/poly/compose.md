---
title: 合成
documentation_of: ../../neo/base/core/poly/compose.hpp
---

$f(g) \bmod x^n$ を求めます。

```cpp
template <class Conv, class K>
std::vector<K> compose(const std::vector<K>& f, const std::vector<K>& g, int n);
```

- `compose<Conv>(f, g, n)`: $f(g) \bmod x^n$ の $n$ 項の係数。手間は $O(M(n) \log n)$。
- $g$ の定数項は 0 であることが要ります。$g$ の定数項が 0 でなければ、先に $f$ を $g(0)$ だけずらして ($f(g_0 + h) = \sum_j f^{(j)}(g_0)/j! \cdot h^j$)、$h = g - g_0$ と合成します。
- 核はテンプレート引数 `Conv` で受け、畳み込みの共通の口の `convolve_range` に渡します。省くと、include した核のうち `K` を受け持つ既定の核になります。

## 中身

Kinoshita–Li の転置の power projection です。$g$ の分母

$$q_t(x, y) = \prod_{j < 2^t} (1 - y\, g(\omega^j x^{1/2^t}))$$

を $t = 0$ から順に作ります。$q_{t+1}$ は $q_t(\sqrt{x}, y) q_t(-\sqrt{x}, y)$ の偶数次の部分です。

y の次数が 2 倍になるので、y は Kronecker 置換で 1 本の配列に潰します。積は 1 回の畳み込みで取ります。

そのあと $t$ を大きい方から 0 まで下げながら、$q_t$ を逆向きに掛けて窓を取ります。これが転置の段で、$f$ の係数が答えの係数に移ります。

転置の段も分母を作る段も畳み込みだけで書けるので、核に変換 (`dft` など) が無くても動きます。
