---
title: 素集合の合併
documentation_of: ../../neo/base/ds/set/UnionFind.hpp
---

要素をグループに分け、2 つのグループの合併と、同じグループかの判定をほぼ定数時間で行います。

```cpp
struct UnionFind {
  explicit UnionFind(int n = 0);
  int find(int x);
  bool same(int a, int b);
  void unite(int a, int b);
  int size(int x);
};
```

`find(x)` は `x` の属するグループの代表を返し、経路を縮めます。`same(a, b)` は 2 つが同じグループかを返します。`unite(a, b)` は 2 つのグループを合併します。合併は大きいほうの代表に寄せるので、木の高さは $O(\log n)$ に収まります。`size(x)` は `x` のグループの要素数を返します。

グループの代表は `p[x] < 0` で表し、その絶対値がグループの大きさです。
