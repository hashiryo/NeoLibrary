---
title: 大きさを型で決める正方行列
documentation_of: ../../neo/linalg/SquareMatrix.hpp
---

大きさ $N$ を型で決める $N \times N$ の行列です。`std::array` に行優先で持つので、セグメント木の要素にしてもヒープを使いません。$N$ が定数なので、内側のループはコンパイラが展開できます。

```cpp
template <class T, int N> struct SquareMatrix {
 std::array<T, N * N> a;                       // a[i * N + j] が i 行 j 列
 constexpr SquareMatrix();                     // 0 で埋める
 constexpr SquareMatrix(const std::array<T, N * N>& a);
 static constexpr SquareMatrix identity();
 constexpr T* operator[](int i);               // i 行目の先頭
 constexpr SquareMatrix pow(uint64_t k) const;
};
SquareMatrix operator+(SquareMatrix l, const SquareMatrix& r);   // -、+=、-= も
SquareMatrix operator*(const SquareMatrix& l, const SquareMatrix& r);   // *= も
SquareMatrix operator*(SquareMatrix l, const T& c);              // c * A、A *= c も
std::array<T, N> operator*(const SquareMatrix& l, const std::array<T, N>& v);
bool operator==(const SquareMatrix& l, const SquareMatrix& r);
```

要素の型の約束は `Matrix` と同じです。和、積、0、1 を持つ半環であればよく、0 と 1 は `zero()` と `one()` があればそれを、なければ `T(0)` と `T(1)` を使います。引き算とスカラー倍は、`T` がそれを持つときだけ使えます。どの操作も constexpr で使えます。`pow(k)` は二分累乗です。

```cpp
using mint= ZMod<998244353>;
SquareMatrix<mint, 2> f;
f[0][0]= f[0][1]= f[1][0]= 1;
mint fib= f.pow(n)[0][1];
SquareMatrix<MinPlus<int>, 3> d;   // 全部 INF
```

## 中身

`T` が `fast_range` と `fast_mul` を持つとき (`neo/algebra/semiring.hpp` の `MinPlus` と `MaxPlus`) は、両方の行列の値がその範囲にあるかを見ます。範囲にあれば、積を `fast_mul` にした i-k-j で掛けます。GCC が列の向きをベクトル化します。それ以外の一般の `T` では、$N \le 8$ なら成分ごとの内積をレジスタにためて求め、それより大きければ i-k-j の順に要素の演算子で足します。i-k-j のままだと、積を 1 つ足すたびに和をメモリに読み書きします。

`T` が `ZMod` なら、積の成分ごとに $N$ 個の積を u64 の和にためてから、1 回だけ還元します。値は表現のまま ($2M$ 未満) 掛け、$N$ 個の和が u64 を溢れうるとき (法が $2^{30}$ に近ければ $N > 4$) だけ $M$ 未満に直してから掛けます。$N > 16$ なら、17 個目からは 8 段ごとに和の $2^{32}$ の位を $2^{32} \bmod M$ に掛けて畳みます。ベクトルとの積も同じです。

procon-judge の atcoder-abc246-h は、3 × 3 の行列をセグメント木に載せる問題です。EPYC 9V45 の x64-gcc で、この形の zmod_fixed_acc は 40.30 ms、積と和のたびに還元する zmod_fixed は 63.95 ms でした。速さの比は x64-gcc で 1.6 倍、x64-clang で 1.9 倍です。大きさを実行時に決める形より 2.7 倍から 5 倍速く、今の Library の Matrix の 4.9 倍です。`SquareMatrix<ZMod<998244353>, 3>` を載せた lib-neo の提出は、手元の M2 で 43.7 ms から 45.2 ms でした。zmod_fixed_acc (42.6 ms から 43.1 ms) と並びます。

cf-1609-e (min-plus の 3 × 3) の lib-neo の速さは、semiring.md の「中身」にあります。

C++20 以上が要ります。
