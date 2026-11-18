---
title: Berlekamp–Massey
documentation_of: ../../neo/polynomial/berlekamp_massey.hpp
---

長さ $n$ の列 $a$ について、$a_i = \sum_{j=1}^{d} c_j a_{i-j}$ ($d \le i < n$) を満たす最短の係数の列 $c$ (長さ $d$) を求めます。

```cpp
template <class K> std::vector<K> berlekamp_massey(const std::vector<K>& a);
```

`K` は体で、`ZMod` なら法は素数であることが要ります。時間は $O(n^2)$ です。返す列は今の Library の `berlekamp_massey` と同じ並び (`c[0]` が $c_1$) で、第 k 項は `linear_recurrence(c, a, k)` で求まります。

## 中身

`ZMod` は値を $M$ 未満の u32 で持ちます。ずれの内積は、列を逆順に並べた配列と AVX2 で 8 個ずつ掛け、u64 にためて 8 回ごとに $2^{32}$ の位を畳みます。漸化式の更新 $C \leftarrow C - x z^m B$ は掛ける数 $x$ が決まっているので、$\lfloor x 2^{32} / M \rfloor$ を作って Shoup の方法で 8 個ずつ掛けます。長さが伸びる回は、新しい $C$ を $B$ の置き場に後ろから書いて置き場を入れ替え、写しを作りません。

procon-judge の yosupo-find-linear-recurrence ($n \le 10^4$) で書き比べた bm_avx2 を、法を定数にしない形で写しました。EPYC 7763 の x64-gcc (run 44256846647) で、この形の最大は 8.38 ms でした。写し元の bm_avx2 は 8.63 ms、Library Checker の最速の提出の写しは 12.10 ms、今の Library は 201.61 ms です。
