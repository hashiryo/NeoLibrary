---
title: 単調な行最小 (SMAWK)
documentation_of: ../../neo/algo/optimize/monge/smawk.hpp
---

行の最小の列が行について単調な行列の、各行の最小の列を SMAWK で求めます。`value` の評価は O(h + w) 回です。

```cpp
template <class T, class F> std::vector<int> smawk(int h, int w, const F& value);
```

`value(r, c)` は行 r 列 c の値を `T` で返します。有効でない (r, c) は `T` の最大値を返してください。行の最小の列が行について単調になるのは、行列が totally monotone のとき (Monge や逆 Monge を含む) です。使い方は [単調な行最小 (分割統治)](monotone_minima.md) と同じです。

## 中身

列を減らしながら奇数番目の行で再帰し、偶数番目の行は両隣の答えの列の間だけを見ます。列を減らす段では、同じ (行, 山の頂点) の値を覚えて 2 回目の評価を省きます。偶数番目の行の走査では、最良の値を覚えて 1 候補 1 評価にします。

評価回数は分割統治より少ないです。計測では、x64-gcc で SMAWK が分割統治より 25% 速く、x64-clang では分割統治が 22% 速かったです (比較の分岐が予測しにくいためと見ています)。記録は algo-notes の notes/dp-optimization.md にあります。
