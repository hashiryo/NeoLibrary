---
title: 凸の層
documentation_of: ../../neo/algo/compute/geo/convex_layers.hpp
---

凸包を剥がす操作を繰り返し、各点が何層目で消えるかを返します。

```cpp
template <class K> std::vector<int> convex_layers(const std::vector<Point<K>>& ps);
```

返り値の `i` 番目は、`ps[i]` が消える層です (1 始まり)。凸包の辺上にある点も、その層で消えます。点はすべて異なることを前提にします。

点を $y$ の順に並べ、区間ごとの左側の凸包の鎖を持つセグメント木を作ります。点を消すと、根までの節の鎖を共通接線で繋ぎ直します (Overmars–van Leeuwen の削除つき動的凸包)。右側の鎖は、点を反転して同じ木で持ちます。1 層ごとに両方の木の根の鎖を集め、その点を消します。計算量は $O(n \log^2 n)$ です。
