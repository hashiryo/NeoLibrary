---
title: 疎な形式的冪級数
documentation_of: ../../neo/base/core/poly/sparse_fps.hpp
---

0 でない項が $K$ 個の形式的冪級数について、逆元、割り算、log、exp、冪、平方根の先頭 $n$ 項を $O(nK)$ で求めます。

```cpp
namespace sparse_fps {
template <class K> std::vector<K> inv(const std::vector<K>& f, int n);
template <class K> std::vector<K> div(const std::vector<K>& f, const std::vector<K>& g, int n);
template <class K> std::vector<K> log(const std::vector<K>& f, int n);
template <class K> std::vector<K> exp(const std::vector<K>& f, int n);
template <class K> std::vector<K> pow(const std::vector<K>& f, long long k, int n);
template <class K> std::vector<K> sqrt(const std::vector<K>& f, int n);
template <class K> std::vector<K> pow_mul_pow(const std::vector<K>& f, long long k, const std::vector<K>& g, long long l, int n);
template <class K> std::vector<K> pow_of_div(const std::vector<K>& f, const std::vector<K>& g, long long k, int n);
template <class K> std::vector<K> exp_of_div(const std::vector<K>& f, const std::vector<K>& g, int n);
}
```

冪級数は密な `vector` で受け (0 の項を読み飛ばします)、長さ $n$ の `vector` で答えを返します。どれも今の Library の `sfps` の関数に合わせた名前です。ただし今の Library は最後の添字 $N$ を受けて長さ $N + 1$ を返すのに対し、こちらは長さ $n$ を受けます。名前空間も、同じ翻訳単位で今の Library とぶつからないように `sparse_fps` にしました。

- `inv(f, n)` は $1/f$ で、$f(0) \ne 0$ を求めます。
- `div(f, g, n)` は $f/g$ を返します。$f$ は密でかまわず、$g$ が疎で $g(0) \ne 0$ であることを求めます。
- `log(f, n)` には $f(0) = 1$ の $f$ を、`exp(f, n)` には $f(0) = 0$ の $f$ を渡します。
- `pow(f, k, n)` は $f^k$ で、$0^0 = 1$ とします。$k$ が負なら、$f(0) \ne 0$ を求めます。
- `sqrt(f, n)` は $g^2 = f$ となる $g$ を 1 つ返し、無ければ空を返します。`K` は `ZMod` だけです。
- `pow_mul_pow(f, k, g, l, n)` は $f^k g^l$ を返します。
- `pow_of_div(f, g, k, n)` は $(f/g)^k$ を返します。
- `exp_of_div(f, g, n)` は $\exp(f/g)$ を返します。$g$ の最も低い次数より、$f$ の最も低い次数が大きいことを求めます。

`K` は体です。`log`、`exp`、冪、平方根は $1/i$ ($i < n$) を使うので、`ZMod` なら法は $n$ より大きい素数であることが要ります。

## 中身

逆元と割り算は $h_i = f_i/g_0 - \sum_k (g_k/g_0) h_{i-k}$ を前から回します。log は $h_i = i (\log f)_i$ が $h_i = i f_i - \sum_{k \ge 1} f_k h_{i-k}$ を満たすので、$h$ を求めてから $1/i$ を掛けます。exp は $i g_i = \sum_k k f_k g_{i-k}$ です。冪は $f = a x^e u$ ($u(0) = 1$) と分けます。$v = u^k$ は $v' u = k u' v$ を満たします。$u$ の 0 でない項の次数を $d_j$、係数を $u_j$ とすると、$v_i = (\sum_j u_j (k+1) d_j v_{i-d_j})/i - \sum_j u_j v_{i-d_j}$ を前から回せます。$u^k$ の係数は $k$ の多項式なので、`ZMod` なら $k$ の代わりに $k \bmod M$ を使えます。平方根は、この式で $k = 1/2$ としたものです。

`ZMod` は値を $M$ 未満の u32 で持ち、項の積を u64 に 15 個までためてから剰余を取ります ($15 M^2 < 2^{64}$)。$i$ が項の次数を越えるたびに使う項が 1 つ増えるので、区間ごとに項の数をテンプレートの引数にします。16 項以上は、実行時の数のループで 15 個ごとに剰余を取ります。いちばん次数の低い項の次数が 1 のときは、そこが前の値を待つ鎖になります。前の値はレジスタに持ち、その掛け算は Shoup の方法か 1 回の Barrett にします。$1/i$ の表は 8 本の鎖でまとめて求めます。$i \equiv b \pmod 8$ の数の累積積をレーン $b$ が受け持ち、AVX2 の Montgomery の積で前向きに置いてから、レーンごとの総積の逆元から後ろ向きに戻します。答えは返す `vector` の置き場に直接書き、$1/i$ の表は Linux では huge page に置きます。`pow_mul_pow`、`pow_of_div`、`exp_of_div` は、今の Library と同じ漸化式を値の型のまま回します。

書き比べの問題は procon-judge の yosupo の inv/exp/log/pow/sqrt-of-formal-power-series-sparse ($n \le 10^6$、$K \le 10$) です。そこで書き比べた sp1 を、法を定数にしない形で写しました。x64-gcc (run 44256846647) の最大は次のとおりです。

| 問題 | CPU | この形 | 写し元の sp1 | 今の Library |
| --- | --- | ---: | ---: | ---: |
| inv | Xeon 6975P-C | 11.45 ms | 5.59 ms | 28.69 ms |
| exp | EPYC 7763 | 26.52 ms | 14.12 ms | 91.94 ms |
| log | Xeon 6975P-C | 15.23 ms | 8.97 ms | 36.35 ms |
| pow | EPYC 7763 | 24.14 ms | 22.42 ms | 71.76 ms |
| sqrt | EPYC 7763 | 25.15 ms | 26.59 ms | 87.01 ms |

この形の提出は、入力の項から長さ $n$ の `vector<ZMod>` を作り、答えを u32 の列に直してハーネスへ返します。sp1 は長さ $n$ の u32 の列を 1 本だけ使います。inv、exp、log の差の大半は、この 4 MB の配列 2〜3 本ぶんです (CI のランナーで 4 KB のページに初めて触る手間が 1 本 2 ms ほど)。
