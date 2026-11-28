#pragma once
// 多項式のテイラーシフト f(x + c)。f の全体 (次数 + 1 項) を返す。
// f(x+c) = Σ_j x^j/j! · Σ_{i>=j} a_i i! c^{i-j}/(i-j)! で、後ろの和は 1 回の畳み込み。
// 核はテンプレート引数 Conv で受け、畳み込みの共通の口 (convolve.hpp) の convolve に渡す。
// K は体で、|f| は法 (標数) より小さいこと (1/j! を使う)。
#include <vector>
#include "neo/base/core/convolve/convolve.hpp"

template <class Conv= convolve_internal::Auto, class K>
std::vector<K> taylor_shift(const std::vector<K>& f, K c) {
 const int n= int(f.size());
 if (n == 0 || c == K()) return f;
 std::vector<K> fact(n), ifact(n);
 fact[0]= K(1);
 for (int i= 1; i < n; ++i) fact[i]= fact[i - 1] * K(i);
 ifact[n - 1]= fact[n - 1].inv();
 for (int i= n - 1; i > 0; --i) ifact[i - 1]= ifact[i] * K(i);
 // A[i] = a_{n-1-i} (n-1-i)!、B[k] = c^k / k!。(A * B)[n-1-j] = Σ_{i>=j} a_i i! c^{i-j}/(i-j)!。
 std::vector<K> a(n), b(n);
 for (int i= 0; i < n; ++i) a[i]= f[n - 1 - i] * fact[n - 1 - i];
 b[0]= K(1);
 K cp= K(1);
 for (int k= 1; k < n; ++k) cp*= c, b[k]= cp * ifact[k];
 const std::vector<K> t= convolve<Conv>(a, b);
 std::vector<K> g(n);
 for (int j= 0; j < n; ++j) g[j]= t[n - 1 - j] * ifact[j];
 return g;
}
