---
title: 多倍長整数
documentation_of: ../../neo/base/core/nt/BigInt.hpp
---

符号つきの多倍長整数です。桁は `u32` で $10^9$ 進法に持ち、積は 3 素数の NTT、割り算は Knuth の D と Newton 法で行います。

```cpp
struct BigInt {
  BigInt();
  template <class T> BigInt(T v);          // 整数型から
  explicit BigInt(const std::string& s);   // 10 進の文字列から
  std::string to_string() const;
  // 比較 (== != < <= > >=)、+ - * / %、単項 -、+= -= *= /= %=
};
friend void divmod(const BigInt& a, const BigInt& b, BigInt& q, BigInt& r);
```

除算は C++ と同じ切り捨てで、余りの符号は割られる数に合わせます。0 の符号は正で、`to_string` は "0" を返します。$10^9$ 進法にしたのは、10 進の入出力が桁数に比例する割り算だけで済むためです ($2^{32}$ 進法だと基数の変換にもう 1 回多倍長の割り算が要ります)。桁の値は $10^9$ 未満で、NTT の素数 (最小 880803841) の 2 倍未満なので、法を取るのは条件つきの引き算 1 回で足ります。

積は、小さい側が 64 桁以下なら素朴な二重ループ、それより大きければ 3 素数の NTT (`crt_engine<Ntt>` の `conv_at`) と Garner です。係数の上界 $\min(n, m)(10^9 - 1)^2$ は 3 素数の積 (約 $2^{89.2}$) 未満でなければなりませんが、実際の大きさでは問題になりません。

割り算は、除数か商が 24 桁以下のときと除数が 512 桁以下のときは Knuth の D、大きいときは Newton 法です。Newton 法は正規化した除数の逆数 $I \approx \lfloor \mathrm{BASE}^{m+K} / b \rfloor$ を 2 段で詰めます。作業の桁 ($m \times (p + 1)$ の積) で $p$ が半分まで倍々に詰めてから、答えの桁 (全体の積) で 2 回だけ詰めます。切り捨てを答えの桁で行うのが肝で、作業中の桁で切ると、切り捨ての 1 が答えの $\mathrm{BASE}^{K-p}$ 倍になり、誤差の漸化式 $e \leftarrow (b/\mathrm{BASE}^m) e^2 + 1$ に不動点が無くなって収束しません。

$10^9$ 進法の桁をそのまま持つので、yosupo の division-of-big-integers (200 万桁まで) はこの 1 本で通ります。今の Library の BigInt ($10^{15}$ 進、4 素数の NTT) と比べて 1.5 倍ほど遅いです。
