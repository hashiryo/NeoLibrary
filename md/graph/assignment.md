---
title: 割り当て問題
documentation_of: ../../neo/graph/assignment.hpp
---

行列の各行に違う列を 1 つずつ割り当て、割り当てた成分の和を最小 (または最大) にします。

```cpp
template <class Cost> struct Assignment { Cost cost; std::vector<int> col; std::vector<Cost> u, v; };
template <class Cost> Assignment<Cost> assignment(int n, int m, const std::vector<Cost>& a, bool maximize = false);
```

`a` は n 行 m 列の行列を行優先で並べたもので、大きさは nm です。n ≤ m なら各行に違う列を 1 つずつ、n > m なら各列に違う行を 1 つずつ割り当てます。`maximize` が true なら和を最大にします。

返す `cost` は割り当てた成分の和、`col[i]` は行 i に割り当てた列です。割り当てられない行は -1 になります。`u` と `v` は行と列の双対です。最小化なら、すべての組で $u_i + v_j \le a_{ij}$ が成り立ち、割り当てた組では等号になります。最大化なら不等号の向きが逆です。割り当てない列 (n < m のとき) の `v` と、割り当てない行 (n > m のとき) の `u` は 0 です。

`Cost` は整数の型か浮動小数点数の型です。整数の型では、費用の絶対値の最大を A として、$(n + m) A < 2^{56}$ である必要があります。計算量は $O(\min(n, m)^2 \max(n, m))$ です。

## 中身

Hungarian 法で、行を 1 つずつ足し、最短の増加路を Dijkstra で探します。1 行を足す間に動かしたポテンシャルの量の累計を持ち、使った列にはまとめて反映するので、列の側の走査は 1 歩につき 1 回で済みます。整数の費用ではこの走査を AVX2 で 4 列ずつ回し、浮動小数点数の費用では同じ形をスカラーで回します。正方行列のときは、列のポテンシャルを列ごとの最小から始め、縮約費用が同じ列のうちでは空いている列を選びます。どの割り当ても同じ費用になる行列などで、増加路が短くなります。行が列より多ければ、転置して解きます。

procon-judge の yosupo-assignment ($N \le 500$) で書き比べました。EPYC 7763 の x64-gcc で、最も遅いケースが 28.8 ms でした。同じ回に、今の Library の NetworkSimplex で解くと 225 ms でした。実数の費用は、ABC 373 G (ユークリッド距離の和が最小の完全マッチング) で確かめています。記録は algo-notes の notes/network-flow.md にあります。
