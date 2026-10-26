---
title: 階数
documentation_of: ../../neo/linalg/rank.hpp
---

行列の階数を求めます。

```cpp
template <class T> int rank(Matrix<T> A);
```

`T` は体です。`ZMod` なら法は素数でなければなりません (合成数の法では要の逆元を取れないことがあり、`inv()` の assert で止まります)。浮動小数点数では、列の中で絶対値が最大の値を要に選び、絶対値が $10^{-8}$ 以下なら 0 とみなします。$h \times w$ の行列に $O(hw \min(h, w))$ の時間がかかります。

`using namespace std;` のもとでは、型の特性 `std::rank` と名前がぶつかって `rank(A)` が曖昧になるので、`::rank(A)` と書きます。

## 中身

`ZMod` は、縦長なら転置して行の数を $\min(h, w)$ にし、行を u64 で持って、要を 8 本ずつ貯めて残りの各行を 1 回だけ読み書きします (遅らせた更新)。手順は `det` と同じで、要は列の中で値が 0 でない最初の行です。最後の組の更新は、残りの行をもう読まないので足しません。要の行を貯める領域は $\min(h, w, 8)$ 行分だけ取ります ($1 \times 250000$ のような行列で 8 行分を取ると 16 MB になるためです)。$\max(h, w) \le 16$ の `ZMod` とほかの体は、素朴な掃き出しにします。

procon-judge の yosupo-matrix-rank で書き比べた zmod_rank_blk8 を写しました。写し元は x64-gcc の Xeon 8370C で 7.05 ms と首位で、今の Library の 8 倍の速さです。lib-neo の提出は、手元の M2 で最大のケースが写し元と 1% の差で、全ケースの和は 5% 遅くなりました。差は $250000 \times 1$ のような縦長の行列で、`Matrix` に一度写す分です。記録は algo-notes の notes/linear-algebra.md にあります。
