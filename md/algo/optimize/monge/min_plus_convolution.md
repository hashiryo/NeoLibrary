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

`min_plus_convolution_convex` は、a が凸で b が任意のときに使えます。行 k 列 j の値を a[k-j] + b[j] にすると行の最小の列が行について単調になるので、[分割統治](monotone_minima.md) で作ります。a[i] + b[j] が `T` に収まる必要があります。

中身は [SMAWK](smawk.md) でも作れます。評価回数は SMAWK のほうが少ないですが、計測では 1 候補 1 評価の分割統治のほうが x64 の gcc と clang の両方で 1.6 倍速かったため、分割統治にしています (記録は algo-notes の notes/dp-optimization.md)。

凹と任意の min-plus convolution は、凸と任意のようには作れません。逆順にしても凹は凹のままで、階段状の境界で行最小の単調性も崩れるためです (階段を扱う SMAWK などが要ります)。

## 使い方

```cpp
std::vector<long long> c = min_plus_convolution_convex(a, b);            // a が凸、b は任意
std::vector<long long> d = min_plus_convolution_convex_convex(a2, b2);  // a2 も b2 も凸
```

procon-judge の min_plus_convolution_convex_convex で、併合は今の Library の monotone_minima より 7.4 倍速かったです。記録は algo-notes の notes/dp-optimization.md にあります。
