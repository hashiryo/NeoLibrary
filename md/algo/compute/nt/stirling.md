---
title: スターリング数 (小さい p、大きい n)
documentation_of: ../../neo/algo/compute/nt/stirling.hpp
---

小さい素数 $p < 2^{15}$ を法とする第 1 種・第 2 種のスターリング数を、$n \le 10^{18}$ の範囲で返します。

```cpp
class StirlingNumber {
 explicit StirlingNumber(unsigned p, bool first = true, bool second = true);
 int nCk(unsigned long long n, unsigned long long k) const;  // C(n, k) mod p
 int S1(unsigned long long n, unsigned long long k) const;   // [n k] mod p
 int S2(unsigned long long n, unsigned long long k) const;   // {n k} mod p
};
```

$p$ は素数です。`first` と `second` で、作る三角形を選びます。片方だけ要る問題では、要らない方を作りません。

$n < p$ の範囲は、$C(n, k)$ の三角形と、スターリング数の三角形を $p$ で埋めておきます。$n \ge p$ は Lucas の定理と同じ形で、$n$ と $k$ を $p$ 進の桁に分けて求めます。今の Library の `StirlingNumber` と同じ算法で、三角形は `vector<vector<uint16_t>>` ではなく 1 本の配列 (行 $i$ は $i(i+1)/2$ から $i+1$ 個) に詰めます。

前計算は三角形の大きさで $O(p^2)$、1 問は桁の数だけなので $O(\log_p n)$ です。三角形は $p = 4999$ で 25 MB、$C$ と合わせて 50 MB です。

`yosupo-stirling-number-of-the-first-kind-small-p-large-n` と `yosupo-stirling-number-of-the-second-kind-small-p-large-n` ($T = 2 \times 10^5$、$p \le 4999$) の手元の Apple M2 Max で、第 1 種は最大 1308 ms、50.4 MB、今の Library は 1332 ms、83.9 MB でした。第 2 種は 1313 ms、50.3 MB と、1291 ms、83.8 MB です。メモリが 4 割ほど減っています。
