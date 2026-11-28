---
title: 単独の NTT の核
documentation_of: ../../neo/base/core/convolve/ntt.hpp
---

`ZMod<P>` (P は NTT に使える素数) の畳み込みを、基数 4 の NTT で求める核です。実装は 2 つあり、include したときの環境で選びます。

- `__AVX2__` (x86 で AVX2 が使えるとき) か `SIMDE_ENABLE_NATIVE_ALIASES` (SIMDe を入れてあるとき) は、8 レーンの AVX2 の実装 (`ntt_avx2_impl.hpp`) を使います。x86 では intrinsics、arm などでは SIMDe で、同じソースが通ります。
- どちらも無いときは、スカラー基数 4 の実装 (`ntt_scalar_impl.hpp`) に落ちます。`NEO_NTT_SCALAR` を定義すると、速さを比べるためにスカラーの実装を選べます。

どちらの実装も口と評価点の並びは同じです。

```cpp
struct Ntt {
 template <class K> static void convolve_range(const K* a, int n, const K* b, int m, int lo, int hi, K* out);
 template <unsigned P> requires ntt_friendly<P> static void dft(ZMod<P>* a, int n);
 template <unsigned P> requires ntt_friendly<P> static void idft(ZMod<P>* a, int n);
 template <unsigned P> requires ntt_friendly<P> static void dft_doubling(ZMod<P>* a, int n);
 template <unsigned P> requires ntt_friendly<P> static void even_dft(ZMod<P>* a, int n);
 template <unsigned P> requires ntt_friendly<P> static void odd_dft(ZMod<P>* a, int n);
};
```

`convolve_range` は畳み込みの共通の口 (`convolve.hpp`) から呼ばれる形です。`out[k - lo] = Σ_{i+j=k} a[i] b[j]` を書き、`a` と `b` が同じ配列なら 2 乗として求めます。巡回の長さは `max(hi, n + m - 1 - lo)` 以上の 2 冪です。仕事の配列は長さごとに残すので、小さい畳み込みを何度も呼ぶ道具 (relaxed など) でも確保の費用がかかりません。

- `dft` は DIF で、bit 反転の並べ替えをせず、index の 4 進の桁を逆にした順に評価値を置きます (2 の指数が奇数のときは最後に基数 2 の段を 1 つ)。`idft` は DIT でその逆です。入出力はどちらもふつうの表現 (ZMod の値) です。
- `dft_doubling` は、長さ $n$ の変換から長さ $2n$ の変換を作ります (係数に戻して 0 を詰めて変換する素直な形)。
- `even_dft` と `odd_dft` は、長さ $n$ の変換から、偶数次の部分と奇数次の部分の長さ $n/2$ の変換を作ります。
- AVX2 の実装は値とひねりを Montgomery の表現で持って段を回し、段の入出力で直します。段のひねりの表は長さごとに 1 回だけ作って残します。

`ntt_friendly<P>` は、$P > 2$ が素数で、$P - 1$ が $2^{20}$ で割り切れることです (2^20 までの変換を約束します)。`ntt_friendly` でない素数、実行時の法 `ZMod<0>`、`ZMod` 以外の係数は、素朴な積 (`Naive`) に回します。

この核を include すると、`ntt_friendly` な素数について既定の核になります。多素数の CRT の核 (`crt.hpp`) と一緒に include したときは、`ntt_friendly` な素数ではこちらが (制約の付いた宣言が優先されて) 既定になります。
