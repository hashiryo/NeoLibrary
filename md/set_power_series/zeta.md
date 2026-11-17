---
title: ゼータ変換と OR・AND・XOR 畳み込み
documentation_of: ../../neo/set_power_series/zeta.hpp
---

部分集合と上位集合のゼータ変換・メビウス変換、アダマール変換と、OR、AND、XOR の畳み込みです。和の代わりにモノイドの演算で畳むゼータ変換もあります。長さは $2^n$ で、変換はその場で書き換えます。どれも $O(n 2^n)$ です。

```cpp
template <class T> void sps::subset_zeta(std::vector<T>& f);    // f[S] <- sum_{T ⊆ S} f[T]
template <class T> void sps::subset_mobius(std::vector<T>& f);  // subset_zeta の逆
template <class T> void sps::supset_zeta(std::vector<T>& f);    // f[S] <- sum_{S ⊆ T} f[T]
template <class T> void sps::supset_mobius(std::vector<T>& f);  // supset_zeta の逆
template <class T> void sps::hadamard(std::vector<T>& f);       // f[S] <- sum_T (-1)^{|S ∩ T|} f[T]
template <class T> std::vector<T> sps::or_convolve(std::vector<T> f, std::vector<T> g);
template <class T> std::vector<T> sps::and_convolve(std::vector<T> f, std::vector<T> g);
template <class T> std::vector<T> sps::xor_convolve(std::vector<T> f, std::vector<T> g);
template <class M> void sps::subset_zeta(std::vector<typename M::T>& f);  // f[S] <- op_{T ⊆ S} f[T]
template <class M> void sps::supset_zeta(std::vector<typename M::T>& f);  // f[S] <- op_{S ⊆ T} f[T]
// k 進 (長さ k^n、k >= 2)
template <class T> void sps::subset_zeta(std::vector<T>& f, int k);     // f[x] <- sum_{y ≤ x (桁ごと)} f[y]
template <class T> void sps::subset_mobius(std::vector<T>& f, int k);
template <class T> void sps::supset_zeta(std::vector<T>& f, int k);     // f[x] <- sum_{y ≥ x (桁ごと)} f[y]
template <class T> void sps::supset_mobius(std::vector<T>& f, int k);
template <class M> void sps::subset_zeta(std::vector<typename M::T>& f, int k);
template <class M> void sps::supset_zeta(std::vector<typename M::T>& f, int k);
```

## 環の変換

値の型 `T` の演算子を使います。ゼータ変換は `+=` を、メビウス変換とアダマール変換は `-=` も使います。`xor_convolve` は最後に $2^n$ で割るので、法が偶数の `ZMod` では使えません。

- `ZMod` と、32 bit か 64 bit の整数は、AVX2 の核 ([zeta_simd.hpp](zeta_simd.md)) で行う。整数の和は $2^{32}$ か $2^{64}$ を法とする。
- ほかの型は、今の Library と同じ汎用の道 ([generic.hpp](generic.md)) で行う。`+=` を書いた自前の型 (和を min にした `MinPlus` など) もここで受ける。

`xor_convolve` は、$2^n$ を `T` に直してから割ります。今の Library の `xor_convolve` は `size_t` で割るので、`long long` で答えが負になると誤ります。

## モノイドのゼータ

`subset_zeta<M>` と `supset_zeta<M>` は、和の代わりにモノイドの演算で畳むゼータ変換です。`M` は [SegmentTree](../data_structure/SegmentTree.md) と同じ決まりの構造体で、`using T` と `static T op(T a, T b)` を使います (`ti` は使いません)。`op` は結合的で可換なものにします。メビウス変換は逆元が要るので、モノイドの形にはありません。

標準の構造体 `monoid::Sum<V>`、`monoid::Min<V>`、`monoid::Max<V>` は、neo/data_structure/segtree/Sum.hpp と Min.hpp にあります。`V` が 32 bit か 64 bit の整数なら、これらは AVX2 の核で行います。`monoid::Sum<ZMod<P>>` も同じです。ほかの構造体は、`op` を 1 つずつ呼ぶ汎用の道で行います。

```cpp
// ARC100 E (Or Plus Max): K ごとに、i | j <= K を満たす i < j の a_i + a_j の最大を出す
struct Top2 {
  using T = std::pair<int, int>;  // 大きい順の 2 個
  static T ti() { return {INT_MIN, INT_MIN}; }
  static T op(T x, T y) {
    if (x.first < y.first) std::swap(x, y);
    return {x.first, std::max(x.second, y.first)};
  }
};
std::vector<Top2::T> f(1 << n);
for (int i = 0; i < (1 << n); ++i) f[i] = {a[i], INT_MIN};
sps::subset_zeta<Top2>(f);  // f[K] は K の部分集合の添字の値の上位 2 個

std::vector<int> g = ...;
sps::subset_zeta<monoid::Max<int>>(g);  // 標準の構造体は AVX2 の核に回る
```

## k 進

`k` を渡す形は、長さ $k^n$ の列の添字を k 進の n 桁と見て、桁ごとに $y \le x$ (`subset`) か $y \ge x$ (`supset`) で畳みます。多次元の累積和で、$O(n k^n)$ です。各要素を $k - 1$ 個まで重ねてよい多重集合と見ると、部分多重集合の和になり、$k = 2$ なら上の形と同じ答えです。環の形とモノイドの形があり、どちらも上と同じ型で AVX2 の核に回ります (桁の重みが 32 byte 以上の段だけ。下の段は値を 1 つずつ合わせる)。

```cpp
// ARC136 D (Without Carry): 筆算で繰り上がらない組 (i < j) の数
std::vector<long long> c(1000000);  // c[x] = x の個数
for (int v : a) ++c[v];
sps::subset_zeta(c, 10);  // c[x] = 桁ごとに x 以下の値の個数
long long ans = 0;
for (int v : a) ans += c[999999 - v] - (どの桁も 4 以下 ? 1 : 0);
ans /= 2;
```

次の 2 つの畳み込みは置いていません。どちらも 1 の k 乗根や次数の変数が要り、NTT の側の話に近いからです。

- 桁ごとに mod k で足す畳み込み (XOR 畳み込みの k 進版)。
- 桁ごとに繰り上がりなしで足す畳み込み。$K[X]/(X_i^k)$ の積で、$k = 2$ なら部分集合畳み込みになる。

## 速さ

procon-judge の self-sps-zeta で、x64-gcc の EPYC 9V45 の回 (run 44075183890) の値です (n = 20)。汎用の道の列は、上の 3 行が標準の構造体と同じ演算の自前の構造体を渡したもので、下の 3 行が sps_generic (今の Library の写し) です。今の Library は、max と min を `+=` に書いた包みの型で呼んでいます。

| 変換 | このヘッダ | 汎用の道 | 今の Library |
| --- | ---: | ---: | ---: |
| 部分集合の向きの max のゼータ (int) | 0.49 ms | 3.25 ms | 2.79 ms |
| 上位集合の向きの min のゼータ (long long) | 1.11 ms | 4.67 ms | 3.79 ms |
| 部分集合の向きの和のゼータ (u64) | 0.96 ms | 3.01 ms | 2.69 ms |
| 部分集合の向きのゼータ (mod 998244353) | 0.52 ms | 5.16 ms | 4.88 ms |
| 上位集合の向きのメビウス変換 (mod 998244353) | 0.50 ms | 3.94 ms | 3.84 ms |
| アダマール変換 (mod 998244353) | 0.79 ms | 7.57 ms | 7.29 ms |
