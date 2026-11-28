---
title: 順列の common interval decomposition tree
documentation_of: ../../neo/algo/compute/permutation/permutation_tree.hpp
---

順列の common interval decomposition tree を作ります。connected interval は `[l, r)` で値の集合も区間になるもの、strong interval はどの connected interval とも交差しない (包含か非交) もので、strong interval 全体の包含のハッセ図がこの木になります。

```cpp
struct PermutationTreeNode {
 int parent;  // 根は -1
 int l, r;    // 半開 [l, r)
 bool linear;  // linear 節なら true、prime 節なら false
};
std::vector<PermutationTreeNode> permutation_tree(const std::vector<int>& p);
std::vector<PermutationTreeNode> renumber_preorder(std::vector<PermutationTreeNode> nodes);
```

`p` は 0 から n - 1 の順列です。節は構築順で、葉が 0 から n - 1、内部の節がそのあとに並びます。親は子より後の番号になることがあります。問題文のように親を子より前にしたいときは、`renumber_preorder` で行きがけ順 (子は左端の昇順) に振り直します。

Uno–Yagiura の common interval 列挙を Xuan–Habib–Paul が木の構築に直したスタック法で、頂点の数に比例する時間で動きます。「値域が左端より前の値を含むか」は、値域の双方向リストを後ろから抜いて作る前後の値で O(1) に見ます。linear 節は子を左端の順に並べたとき、どの連続部分列の和集合も connected interval になる節です。子が 2 個以下の節は必ず linear になります。
