---
title: min-plus convolution
documentation_of: ../../neo/algo/optimize/monge/min_plus_convolution.hpp
---

凸な列どうし、および凸と任意の列の min-plus convolution を作ります。`c[k] = min_{i+j=k} a[i] + b[j]` です。

```cpp
template <class T> std::vector<T> min_plus_convolution_convex_convex(const std::vector<T>& a, const std::vector<T>& b);
template <class T> std::vector<T> min_plus_convolution_convex(const std::vector<T>& a, const std::vector<T>& b);
```

`min_plus_convolution_convex_convex` は、a と b の両方が凸 (差分が広義単調増加) のときに使えます。差分の列を併合するので O(n + m) です。

`min_plus_convolution_convex` は、a が凸で b が任意のときに使えます。行 k 列 j の値を a[k-j] + b[j] にすると行の最小の列が行について単調になるので、[SMAWK](smawk.md) で O(n + m) 回の評価で作ります。a[i] + b[j] が `T` に収まる必要があります。分割統治でも作れ、そちらが速いコンパイラもあるので、計測を見て中身を選び直すことがあります。

## 使い方

```cpp
std::vector<long long> c = min_plus_convolution_convex(a, b);            // a が凸、b は任意
std::vector<long long> d = min_plus_convolution_convex_convex(a2, b2);  // a2 も b2 も凸
```

procon-judge の min_plus_convolution_convex_convex で、併合は今の Library の monotone_minima より 7.4 倍速かったです。記録は algo-notes の notes/dp-optimization.md にあります。
