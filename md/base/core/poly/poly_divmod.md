---
title: 多項式の商と余り
documentation_of: ../../neo/base/core/poly/poly_divmod.hpp
---

多項式の商と余りと gcd を、素朴な $O(nm)$ で求めます。

```cpp
template <class K> std::pair<std::vector<K>, std::vector<K>> poly_divmod(std::vector<K> a, std::vector<K> b);
template <class K> std::vector<K> poly_gcd(std::vector<K> a, std::vector<K> b);
```

多項式は係数の `vector` ($k$ 次の係数が `[k]`) で、末尾の 0 は取り除いて返します (0 は空の `vector`)。`poly_divmod(a, b)` は $a = qb + r$ ($\deg r < \deg b$) の組 $(q, r)$ を返し、$b$ は 0 でないことが要ります。`poly_gcd(a, b)` は monic な gcd を返し、両方 0 なら空を返します。`K` は体で、`ZMod` なら法は素数であることが要ります。

## 中身

割る式を monic にしてから、商の係数 $c$ ごとに余りの側へ $-c$ 倍の割る式を足します。`ZMod` では、掛ける数 $-c$ が決まっているので、Shoup の方法で 8 個ずつ掛けて足します。
