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
