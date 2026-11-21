---
title: 更新付きの部分集合の和
documentation_of: ../../neo/base/ds/range/SubsetSum.hpp
---

長さ $2^n$ の列に、加算と部分集合の和の問いを好きな順で行います。1 点に足す形と、区間 ($a \subseteq c \subseteq b$ の $c$ 全部) に足す形の 2 つがあります。

```cpp
template <class T> struct SubsetSum {
  explicit SubsetSum(int n, int low = -1);                 // 0 で埋める。low は下の bit の数 (既定 n / 2)
  explicit SubsetSum(const std::vector<T>& a, int low = -1);
  void add(unsigned x, T v);   // a[x] += v
  T sum(unsigned x) const;     // x の部分集合 y の a[y] の和
};
template <class T> struct SubsetIntervalSum {
  explicit SubsetIntervalSum(int n);                       // 0 で埋める
  explicit SubsetIntervalSum(const std::vector<T>& a);
  void add(unsigned a, unsigned b, T v);  // a ⊆ c ⊆ b の c 全部に v を足す
  T sum(unsigned a, unsigned b);          // a ⊆ c ⊆ b の c 全部の和
};
```

## 1 点の形

bit を下の `low` 個と上の残りに分け、下の bit だけゼータ変換して持ちます。`add` は $2^{\mathrm{low} - |x \text{ の下}|}$、`sum` は $2^{|x \text{ の上}|}$ かかり、どちらも既定の分け方で $O(2^{n/2})$ です。`low` を大きくすると問いが軽く加算が重くなり、小さくするとその逆です。上位集合の和は、添字を補集合にして使います。`T` は `T()` が 0 で `+=` があればよいです。

```cpp
// CodeChef DANYANUM: m を含む要素の数 (上位集合の数) を、補集合の部分集合の数として問う
SubsetSum<int> cnt(k, 2 * k / 3);   // 問いが加算より多いので low を大きく
cnt.add(full ^ v, 1);               // v を足す (消すなら -1)
bool ok = cnt.sum(full ^ m) >= x;   // m を含む要素が x 個以上か
```

procon-judge の codechef-DANYANUM ($k \le 18$、問いは加算の 18 倍) で、`low` を $2k/3$ にした形は $k/2$ の 3.8 倍速くなりました (156.89 ms と 607.62 ms、EPYC 9V74 の x64-gcc)。

## 区間の形

操作を時間の向きにバケットに分け、バケットの初めの列とその部分集合の和を持ちます。

- バケットより前の分の問いは、自由な bit ($b \setminus a$) をたどるか、$a$ の 1 の bit を包除 ($[a \subseteq c \subseteq b] = \sum_{d \subseteq a} (-1)^{|d|} [c \subseteq b \setminus d]$) で外して部分集合の和を使うかの、少ない方で求めます ($2^{n/2}$ 以下)。
- バケットの中の加算 $(a', b', k)$ は、問いごとに、$a \cup a' \subseteq b \cap b'$ なら $k \cdot 2^{|b \cap b'| - |a \cup a'|}$ を足します。`T` が `unsigned int` なら AVX2 で 8 個ずつ足します。
- 加算は同じ分け方で「点に足す」か「部分集合全部に足す」の印にしておき、作り直すときに配ります。

作り直すのは、バケットの中の加算を問いで足した数の合計が、作り直しの手間の目安 (`unsigned int` は $4 n 2^n$、ほかの型は $n 2^n$) を超えたときです。1 回あたり均して $O(\sqrt{n 2^n})$ で、作り直す時期を使う側で決める必要はありません。$a \subseteq b$ でない区間は空で、`add` は何もせず `sum` は 0 を返します。`T` は `T()` が 0 で、`+`、`-`、`*` と整数からの変換があればよいです。

洛谷 P6562 ($n \le 16$、$q \le 2 \times 10^5$、法 $2^{32}$) の書き比べで決めた形です。EPYC 7763 の x64-gcc で、題解の方針 (バケットは $\sqrt{n 2^n}$ 個) の参照実装が 313.88 ms、この形が 60.35 ms でした。0 の bit も包除で外して上位集合の和を使う分け方 ($2^{n/3}$ 以下) も試しましたが、作り直す配列が倍になるぶん遅くなりました。記録は algo-notes の notes/set-power-series.md にあります。
