---
title: 多項式の商と余り (fast)
documentation_of: ../../neo/base/core/poly/poly_divmod_fast.hpp
---

多項式の商と余りを $O(n \log n)$ で求めます。

```cpp
template <class Conv, class K>
std::pair<std::vector<K>, std::vector<K>> poly_divmod_fast(const std::vector<K>& a, const std::vector<K>& b);
```

$a = q b + r$ ($\deg r < \deg b$) の組 $(q, r)$ を返します。多項式は係数の `vector` ($k$ 次の係数が `[k]`) で、末尾の 0 は取り除いて返します (0 は空の `vector`)。`b` は 0 でないことが要ります。`K` は体で、`ZMod` なら法は素数です。素朴な $O(nm)$ の版は `poly_divmod` (`poly_divmod.hpp`) にあります。

核はテンプレート引数 `Conv` で受け、`fps::inv` と畳み込みの共通の口の `convolve_range` に渡します。省くと、include した核のうち `K` を受け持つ既定の核になります。

## 中身

`a` と `b` の係数を逆順にした多項式の割り算に直します。$q$ の係数を逆順にしたものは、$\mathrm{rev}(a)$ と $\mathrm{rev}(b)$ の商の先頭 $\deg a - \deg b + 1$ 項で、これは `fps::inv` と 1 回の畳み込みで作れます。余りは $a - q b$ の $\deg b$ 未満の部分です。
