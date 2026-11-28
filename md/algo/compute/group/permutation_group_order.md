---
title: 置換群の位数
documentation_of: ../../neo/algo/compute/group/permutation_group_order.hpp
---

置換の列が生成する群の位数を、mod で割った余りで返します。

```cpp
unsigned permutation_group_order(int n, const std::vector<std::vector<int>>& generators, unsigned mod = 998244353);
```

置換 p は p[i] を i の行き先とする全単射で、合成は (p q)[i] = p[q[i]] です。生成元に恒等置換や重複があってよく、逆置換を明示しなくても構いません (有限の置換の半群は群になるので、逆元は正の冪で入ります)。

incremental Schreier–Sims (Knuth の Algorithm A/B) で、各レベルの軌道の大きさの積を掛けます。割り算をしないので、mod は素数でなくてよいです。最初に、2-サイクルがちょうど 1 つで他のサイクルがすべて奇数長の生成元から転置を作り、その共役でできる辺が全点を連結なら S_n として n! を返す近道を使います (推移的で転置を持っても S_n とは限らないので、辺の連結まで見ます)。

位数は最大で n! 程度になるので、必要なら mod を分けて CRT で戻します。n = 1 のときは 1 を返します。
