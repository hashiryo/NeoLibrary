---
title: 多素数の CRT の核
documentation_of: ../../neo/base/core/convolve/crt.hpp
---

任意の法 (素数でも合成数でもよい) の `ZMod<M>` の畳み込みを、いくつかの NTT の素数で求めて Garner の方法で戻す核です。

```cpp
template <class Ntt> struct crt_engine {
 static constexpr int max_primes;
 static constexpr u32 prime(int t);
 static int count_for(int bound_len, int sz, int k, u32 M);
 template <u32 P> static void conv_at(const u32* a, int n, const u32* b, int m, int sz, u32* out);
 static void conv_at(int t, const u32* a, int n, const u32* b, int m, int sz, u32* out);
 template <class F> static void each(int cnt, F&& f);
 static u32 garner(const u32* r, int cnt, u32 M);
 template <class K> static void convolve_range(const K* a, int n, const K* b, int m, int lo, int hi, K* out);
};
using Crt= crt_engine<Ntt>;
```

- 素数は 2^20 以上 2^23 までの 2 冪を $P - 1$ が持つ、$2^{30}$ 未満の 12 本です。2 の指数の大きい順に並べるので、長さ `len` の畳み込みに使える素数は先頭からの本数になります。
- `count_for(bound_len, sz, k, M)` は、値の上界 `bound_len^(1+k) M^2` を覆うのに要る素数の本数を返します (覆えなければ 0)。$k = 0$ は普通の畳み込み (上界 $\min(n, m) M^2$)、$k = |S|$ は Hurwitz 級数環の積です。
- `conv_at` は、素数 1 本の下での畳み込みです。`each(cnt, f)` は、先頭 `cnt` 本の素数それぞれで `f<T, P>()` を呼びます ($T$ は番号、$P$ は素数)。素数ごとの畳み込みの前後に掛け算を挟む道 (Hurwitz 級数環の積) が使います。
- `garner(r, cnt, M)` は、各素数の剰余 `r[t]` から、$x$ を法 $M$ の値に戻します。$x$ は素数の積未満の一意な値です。
- `convolve_range` は畳み込みの共通の口の形です。上界を覆えない大きさ (素因数の多い法の大きい形) は素朴な積に回します。

この核を include すると、任意の `ZMod<M>` について既定の核になります。`ntt_friendly` な素数については、単独の NTT の核 (`ntt.hpp`) の宣言が優先されます。

```cpp
#include "neo/base/core/convolve/convolve.hpp"
#include "neo/base/core/convolve/crt.hpp"
std::vector<ZMod<0>> c= convolve<Crt>(a, b);  // 法は ZMod<0>::set_mod で決めたもの
```
