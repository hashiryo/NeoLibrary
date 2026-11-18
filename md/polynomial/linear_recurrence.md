---
title: 線形漸化式の第 k 項
documentation_of: ../../neo/polynomial/linear_recurrence.hpp
---

$[x^k] P(x)/Q(x)$ と、定数係数の線形漸化式の第 k 項を、Bostan–Mori の方法で求めます。

```cpp
template <class T> T div_at(std::vector<T> P, std::vector<T> Q, unsigned long long k);
template <class T> T linear_recurrence(const std::vector<T>& c, const std::vector<T>& a, unsigned long long k);
```

`div_at` には $Q(0)$ が可逆な $Q$ を渡します。`linear_recurrence` は、$a_n = \sum_{j=1}^{d} c_j a_{n-j}$ ($n \ge d$) の第 k 項を返します。`c[j - 1]` が $c_j$ で、`a` には初項の $d$ 個を渡します。`a` が長ければ先頭の $d$ 個を使い、k が `a` の範囲なら `a[k]` をそのまま返します。並びは今の Library の `linear_recurrence` と同じです。

`T` は可換環で、`ZMod` の法が合成数でもかまいません。割り算は、`div_at` が最初に $Q(0)$ の逆元を掛けるときだけです。時間は $O(d^2 \log k)$ ($d$ は $Q$ の次数) で、FFT を使いません。

## 中身

1 段で $Q(x) Q(-x) = V(x^2)$ と、$P(x) Q(-x)$ の偶数部か奇数部を作り、k を半分にします。$Q = E(x^2) + x O(x^2)$ と $P = P_e(x^2) + x P_o(x^2)$ に分けると、$V(y) = E(y)^2 - y O(y)^2$ です。$P(x) Q(-x)$ の偶数部と奇数部は、それぞれ $P_e E - y P_o O$ と $P_o E - P_e O$ になります。半分の大きさの積 4 回で 1 段進みます。

`ZMod` の積は `convolve_naive` と同じ核で、短い側が 16 以下ならスカラーの道を通ります。作業の配列は呼び出しをまたいで使い回すので、次数の小さい問いをたくさん聞いても確保が重なりません。

procon-judge の self-poly-kth-term (次数 3000 以下、mod 998244353) で書き比べた bmo を写しました。Xeon 8573C の x64-gcc (run 44231884007) で、bmo の最大は 200.61 ms でした。今の Library は 1331.06 ms です。
