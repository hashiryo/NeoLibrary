---
title: ゼータ変換と OR・AND・XOR 畳み込み
documentation_of: ../../neo/set_power_series/zeta.hpp
---

部分集合と上位集合のゼータ変換・メビウス変換、アダマール変換と、OR、AND、XOR の畳み込みです。長さは $2^n$ で、変換はその場で書き換えます。

```cpp
template <class T> void sps::subset_zeta(std::vector<T>& f);    // f[S] <- sum_{T ⊆ S} f[T]
template <class T> void sps::subset_mobius(std::vector<T>& f);  // subset_zeta の逆
template <class T> void sps::supset_zeta(std::vector<T>& f);    // f[S] <- sum_{S ⊆ T} f[T]
template <class T> void sps::supset_mobius(std::vector<T>& f);  // supset_zeta の逆
template <class T> void sps::hadamard(std::vector<T>& f);       // f[S] <- sum_T (-1)^{|S ∩ T|} f[T]
template <class T> std::vector<T> sps::or_convolve(std::vector<T> f, std::vector<T> g);
template <class T> std::vector<T> sps::and_convolve(std::vector<T> f, std::vector<T> g);
template <class T> std::vector<T> sps::xor_convolve(std::vector<T> f, std::vector<T> g);
```

値の型 `T` の演算子だけを使います。ゼータ変換は `+=` だけを使うので、和を min にした `MinPlus` や、`+=` を書いた自前の型 (上位 2 個の max など) もそのまま使えます。メビウス変換とアダマール変換は `-=` も使います。`xor_convolve` は最後に $2^n$ で割るので、法が偶数の `ZMod` では使えません。どれも $O(n 2^n)$ です。今の Library と同じ汎用の道で、速い核はまだありません。
