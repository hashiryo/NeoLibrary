---
title: 階数
documentation_of: ../../neo/linalg/rank.hpp
---

行列の階数を求めます。

```cpp
template <class T> int matrix_rank(Matrix<T> A);
```

`T` は体です。`ZMod` なら法は素数でなければなりません (合成数の法では要の逆元を取れないことがあり、`inv()` の assert で止まります)。浮動小数点数では、列の中で絶対値が最大の値を要に選び、絶対値が $10^{-8}$ 以下なら 0 とみなします。$h \times w$ の行列に $O(hw \min(h, w))$ の時間がかかります。

名前を `rank` にしないのは、`using namespace std;` のもとで型の特性 `std::rank` とぶつかり、`rank(A)` が曖昧になるためです。

## 中身

`ZMod` は、縦長なら転置して行の数を $\min(h, w)$ にし、行を u64 で持って、要を 8 本ずつ貯めて残りの各行を 1 回だけ読み書きします (遅らせた更新)。手順は `det` と同じで、要は列の中で値が 0 でない最初の行です。最後の組の更新は、残りの行をもう読まないので足しません。要の行を貯める領域は $\min(h, w, 8)$ 行分だけ取ります ($1 \times 250000$ のような行列で 8 行分を取ると 16 MB になるためです)。$\max(h, w) \le 16$ の `ZMod` とほかの体は、素朴な掃き出しにします。

法が奇数なら、要を探す列、要の行、行ごとの掛ける数を AVX2 で 4 個ずつ求め、行列を huge page に置きます。形は `det` と同じです。

procon-judge の yosupo-matrix-rank で書き比べた zmod_rank_blk8 を写しました。法が奇数の道は、行列を huge page に置いて要を探す列をベクトル化した zmod_rank_hpv の形です。EPYC 9V45 の x64-gcc (run 41414562027) で、zmod_rank_hpv は 3.30 ms、zmod_rank_blk8 は 4.65 ms でした。Library Checker の首位の写しは 14.12 ms、今の Library は 43.94 ms です。記録は algo-notes の notes/linear-algebra.md にあります。
