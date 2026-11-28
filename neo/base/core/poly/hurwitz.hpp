#pragma once
// Hurwitz 級数環 (二項畳み込みの環) の道具。係数の列 (a_n) に、項ごとの和と
//   (a b)_n = Σ_k C(n, k) a_k b_{n-k}
// を入れた環です (Keigher 1997)。数え上げでは EGF の環と呼びます。微分は添字のずらし (∂a)_n = a_{n+1} で、
// inv、log、exp に整数の割り算が要りません。合成数の法でも 2^64 でも定義できます (今の実装の係数の型は ZMod だけです)。
// 係数の型 K は std::vector<K> で受け、先頭 n 項だけを求める。核はテンプレート引数で受け、
//   Conv: ふつうの畳み込みの核 (二項係数を階乗で分ける道で使う)
//   Crt: 多素数の核 (階乗が法で割り切れるときの道で使う。crt.hpp の default_crt_engine で既定が決まる)
// どちらも省くと、include した核の既定になる。Crt の既定が無いときは、階乗が法で割り切れる法だけを扱う。
// ふつうの冪級数 (fps.hpp) と共通のフレームワークは作らない (適用範囲が違い、共通にできるのは Newton 法の骨組みだけ)。
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <type_traits>
#include <vector>
#include "neo/base/core/algebra/ZMod.hpp"
#include "neo/base/core/convolve/convolve.hpp"

namespace hurwitz_internal {
using u32= unsigned;
using u64= unsigned long long;
using u128= unsigned __int128;

// 法 M の相異なる素因数 (昇順) の控え。M < 2^30 なので 9 個まで。
struct FactorCache {
 u32 m= 0;
 int cnt= 0;
 u32 p[10];
};
inline FactorCache& factor_cache() {
 static FactorCache c;
 return c;
}
// M の相異なる素因数。同じ M の続く呼び出しは控えを使う。
inline const FactorCache& factor_of(u32 M) {
 auto& c= factor_cache();
 if (c.m == M) return c;
 c.m= M, c.cnt= 0;
 u32 x= M;
 for (u32 p= 2; u64(p) * p <= x; ++p) {
  if (x % p) continue;
  c.p[c.cnt++]= p;
  while (x % p == 0) x/= p;
 }
 if (x > 1) c.p[c.cnt++]= x;
 return c;
}

template <u32 P> struct barrett_q {
 static constexpr u64 mu= (~0ull) / P;
 static u32 red(u64 x) {
  const u64 q= u64((u128(x) * mu) >> 64);
  const u64 r= x - q * P;
  return u32(r >= P ? r - P : r);
 }
 static u32 mul(u32 a, u32 b) { return red(u64(a) * b); }
 static u32 pow(u32 a, u64 e) {
  u32 r= 1;
  while (e) {
   if (e & 1) r= mul(r, a);
   a= mul(a, a), e>>= 1;
  }
  return r;
 }
 static u32 inv(u32 a) { return pow(a, u64(P) - 2); }
 // M 未満の値 x の mod P (P >= 880803841、M < 2^30 < 2P なので 1 回引けばよい)。
 static u32 red_small(u32 x) { return x >= P ? x - P : x; }
};

// 二項畳み込みの素朴な形。min(n, m) が小さいときと、素数の積で上界を覆えないときに使う。
// 二項係数は C(k, j) = C(k-1, j) + C(k-1, j-1) で足し算だけから作る (割り算が要らない)。
// 表は短い側の長さだけ持ち、長い側の添字 i は C(k, i) = C(k, k-i) で短い側の j = k - i に読み替える。
template <class K> std::vector<K> convolve_direct(const std::vector<K>& a, const std::vector<K>& b) {
 const int n= int(a.size()), m= int(b.size()), sz= n + m - 1;
 std::vector<K> c(sz);
 if (n <= m) {
  std::vector<K> bin(n);  // bin[i] = C(k, i)
  bin[0]= K(1);
  for (int k= 0; k < sz; ++k) {
   K acc= K();
   for (int i= std::max(0, k - m + 1), e= std::min(n - 1, k); i <= e; ++i) acc+= a[i] * b[k - i] * bin[i];
   c[k]= acc;
   for (int i= std::min(n - 1, k + 1); i > 0; --i) bin[i]+= bin[i - 1];
  }
 } else {
  std::vector<K> bin(m);  // bin[j] = C(k, j)
  bin[0]= K(1);
  for (int k= 0; k < sz; ++k) {
   K acc= K();
   for (int i= std::max(0, k - m + 1), e= std::min(n - 1, k); i <= e; ++i) acc+= a[i] * b[k - i] * bin[k - i];
   c[k]= acc;
   for (int j= std::min(m - 1, k + 1); j > 0; --j) bin[j]+= bin[j - 1];
  }
 }
 return c;
}

// 階乗で分ける道: c_n = n! [x^n] (Σ a_k x^k / k!) (Σ b_j x^j / j!)。1, ..., sz-1 が法で割り切れない前提。
template <class Conv, class K> std::vector<K> convolve_scaled(const std::vector<K>& a, const std::vector<K>& b) {
 const int n= int(a.size()), m= int(b.size()), sz= n + m - 1;
 std::vector<K> fact(sz), ifact(sz);
 fact[0]= K(1);
 for (int i= 0; i < sz - 1; ++i) fact[i + 1]= fact[i] * K(i + 1);
 ifact[sz - 1]= fact[sz - 1].inv();
 for (int i= sz - 1; i > 0; --i) ifact[i - 1]= ifact[i] * K(i);
 std::vector<K> fa(n), fb(m);
 for (int i= 0; i < n; ++i) fa[i]= a[i] * ifact[i];
 for (int i= 0; i < m; ++i) fb[i]= b[i] * ifact[i];
 // 同じ列のときは、核の 2 乗の道 (a == b の判定) に乗せる。
 std::vector<K> c= (&a == &b) ? convolve<Conv>(fa, fa) : convolve<Conv>(fa, fb);
 for (int i= 0; i < sz; ++i) c[i]*= fact[i];
 return c;
}

// 階乗が法で割り切れるときの道 (EntropyIncreaser「任意模数二项卷积」)。
// M の素因数のうち sz 未満のものを S とし、k! = u_k ∏_{p∈S} p^{v_p(k!)} と分ける。u_k は M と互いに素。
//   c_n = u_n Σ_k (∏_{p∈S} p^{d_p}) (a_k u_k^{-1}) (b_{n-k} u_{n-k}^{-1}),  d_p = v_p(n!) - v_p(k!) - v_p((n-k)!) <= log_p n
// の和 S_n は整数で N^{1+|S|} M^2 未満。NTT の素数 q の下で p は単元なので、
//   A_k = (a_k u_k^{-1}) ∏_p p^{-v_p(k!)},  B_j も同じ、を畳み込んで ∏_p p^{v_p(n!)} を掛けると S_n mod q になり、
// 上界を覆う本数の素数で Garner すれば S_n が戻る。最後に u_n を掛けて c_n mod M にする。
template <class Crt, class K> std::vector<K> convolve_composite(const std::vector<K>& a, const std::vector<K>& b) {
 const int n= int(a.size()), m= int(b.size()), sz= n + m - 1;
 const u32 M= K::mod();
 const FactorCache& f= factor_of(M);
 u32 S[10];
 int ns= 0;
 for (int i= 0; i < f.cnt; ++i)
  if (f.p[i] < u32(sz)) S[ns++]= f.p[i];
 const int cnt= Crt::count_for(sz, sz, ns, M);
 if (cnt == 0) return convolve_direct(a, b);  // 覆えない (素因数が多い大きい法)。遅い道に落ちる。
 // 1. ppart[k] = ∏_{p∈S} p^{v_p(k)}、strip[k] = k / ppart[k]。k = 0 はどちらも 1。
 std::vector<u32> ppart(sz), strip(sz);
 ppart[0]= strip[0]= 1;
 for (int k= 1; k < sz; ++k) {
  u32 x= u32(k), q= 1;
  for (int i= 0; i < ns; ++i)
   while (x % S[i] == 0) x/= S[i], q*= S[i];
  strip[k]= x, ppart[k]= q;
 }
 // 2. u[k] = k! の p を抜いた部分 mod M と、その逆元 (前積を 1 回だけ逆にして配る)。
 std::vector<K> u(sz), iu(sz);
 u[0]= K(1);
 for (int k= 1; k < sz; ++k) u[k]= u[k - 1] * K(strip[k]);
 K acc= K(1);
 for (int k= 0; k < sz; ++k) iu[k]= acc, acc*= u[k];
 acc= acc.inv();
 for (int k= sz - 1; k >= 0; --k) {
  iu[k]*= acc;
  acc*= u[k];
 }
 // 3. 各素数で、両側に ∏ p^{-v_p(k!)} を掛けて畳み込み、∏ p^{v_p(k!)} を掛ける。
 std::vector<u32> res(size_t(cnt) * sz);
 Crt::each(cnt, [&]<int T, u32 P>() {
  using B= barrett_q<P>;
  std::vector<u32> qv(sz), pv(sz), xa(n), xb(m), cv(sz);
  // qv[k] = ∏_{j<=k} ppart[j] = ∏ p^{v_p(k!)}、pv[k] はその逆 (前積の 1 回の逆から配る)。
  u32 a0= 1;
  for (int k= 0; k < sz; ++k) a0= B::mul(a0, ppart[k]), qv[k]= a0;
  a0= B::inv(qv[sz - 1]);
  for (int k= sz - 1; k >= 0; --k) pv[k]= a0, a0= B::mul(a0, ppart[k]);
  for (int k= 0; k < n; ++k) xa[k]= B::mul(B::red_small((a[k] * iu[k]).val()), pv[k]);
  for (int k= 0; k < m; ++k) xb[k]= B::mul(B::red_small((b[k] * iu[k]).val()), pv[k]);
  Crt::template conv_at<P>(xa.data(), n, xb.data(), m, sz, cv.data());
  for (int k= 0; k < sz; ++k) res[size_t(T) * sz + k]= B::mul(cv[k], qv[k]);
 });
 // 4. Garner で S_n を戻し、u_n を掛ける。
 u32 r[12];
 std::vector<K> c(sz);
 for (int k= 0; k < sz; ++k) {
  for (int t= 0; t < cnt; ++t) r[t]= res[size_t(t) * sz + k];
  c[k]= u[k] * K(Crt::garner(r, cnt, M));
 }
 return c;
}
}  // namespace hurwitz_internal

namespace hurwitz {
// 二項畳み込み c_n = Σ_k C(n, k) a_k b_{n-k}。どちらかが空なら空を返す。
template <class Conv= convolve_internal::Auto, class Crt= convolve_internal::Auto, class K> std::vector<K> convolve(const std::vector<K>& a, const std::vector<K>& b) {
 using namespace hurwitz_internal;
 const int n= int(a.size()), m= int(b.size());
 if (n == 0 || m == 0) return {};
 const int sz= n + m - 1;
 if (std::min(n, m) < 60) return convolve_direct(a, b);
 if constexpr (is_zmod_v<K>) {
  const u32 M= K::mod();
  if (M == 1) return std::vector<K>(sz);
  const FactorCache& f= factor_of(M);
  if (f.p[0] >= u32(sz)) return convolve_scaled<Conv>(a, b);
  if constexpr (requires(const K* p) { default_crt_engine(p); }) {
   return convolve_composite<decltype(default_crt_engine(static_cast<const K*>(nullptr))), K>(a, b);
  } else {
   assert(false && "階乗が法で割り切れる法には crt.hpp (の default_crt_engine) が要ります");
   return {};
  }
 } else {
  // ZMod 以外の型は、1, ..., sz-1 が逆元を持つ前提で階乗の道を使う。
  return convolve_scaled<Conv>(a, b);
 }
}

// f の先頭 n 項の逆数 g (二項畳み込みの環で f g = 1 mod x^n)。f[0] は逆元を持つこと。
// Newton 法: g <- g (2 - f g)。1 段は、t = (f g) の m 次以上の尾を次数 m に置いた列で g[m, m2) = -(g t) の [m, m2)。
// 尾を次数 m に置き直すのは、二項畳み込みでは x^m を掛けても添字がずれるだけにならないため (普通の冪級数の形は使えない)。
template <class Conv= convolve_internal::Auto, class Crt= convolve_internal::Auto, class K> std::vector<K> inv(const std::vector<K>& f, int n) {
 assert(n >= 1);
 assert(!f.empty() && f[0] != K());
 std::vector<K> g(1, f[0].inv());
 for (int m= 1; m < n;) {
  const int m2= std::min(2 * m, n), cnt= m2 - m;
  const int fl= std::min(int(f.size()), m2);
  std::vector<K> fg= convolve<Conv, Crt>(std::vector<K>(f.begin(), f.begin() + fl), g);
  fg.resize(std::max(int(fg.size()), m2));  // f が短いとき、足りない次数は 0
  std::vector<K> e(m2);
  for (int i= 0; i < cnt; ++i) e[m + i]= fg[m + i];
  std::vector<K> ge= convolve<Conv, Crt>(g, e);
  g.resize(m2);
  for (int i= 0; i < cnt; ++i) g[m + i]= -ge[m + i];
  m= m2;
 }
 return g;
}

// f の先頭 n 項の log (f[0] = 1)。∂(log f) = (∂f) / f を積分する (Hurwitz 級数環の積分は添字をずらすだけ)。
template <class Conv= convolve_internal::Auto, class Crt= convolve_internal::Auto, class K> std::vector<K> log(const std::vector<K>& f, int n) {
 assert(n >= 1);
 assert(!f.empty() && f[0] == K(1));
 std::vector<K> g(n);
 if (n == 1) return g;
 const int fl= std::min(int(f.size()), n), dl= fl - 1;
 std::vector<K> df(dl);
 for (int i= 0; i < dl; ++i) df[i]= f[i + 1];
 std::vector<K> iv= inv<Conv, Crt>(f, n - 1);
 std::vector<K> t= dl > 0 ? convolve<Conv, Crt>(df, iv) : std::vector<K>();
 t.resize(std::max(int(t.size()), n - 1));  // f が定数のとき、∂f = 0 なので log = 0
 for (int k= 1; k < n; ++k) g[k]= t[k - 1];
 return g;
}

// f の先頭 n 項の exp (f[0] = 0)。Newton 法 g <- g (1 + f - log g) を倍々に回す。
template <class Conv= convolve_internal::Auto, class Crt= convolve_internal::Auto, class K> std::vector<K> exp(const std::vector<K>& f, int n) {
 assert(n >= 1);
 assert(f.empty() || f[0] == K());
 std::vector<K> g(1, K(1));
 for (int m= 1; m < n;) {
  const int m2= std::min(2 * m, n);
  std::vector<K> gp(m2);
  std::copy(g.begin(), g.end(), gp.begin());
  std::vector<K> lg= log<Conv, Crt>(gp, m2);
  std::vector<K> w(m2);
  for (int i= 0; i < m2; ++i) w[i]= (i < int(f.size()) ? f[i] : K()) - lg[i];
  w[0]+= K(1);
  std::vector<K> t= convolve<Conv, Crt>(g, w);
  g.resize(m2);
  for (int i= 0; i < m2; ++i) g[i]= t[i];
  m= m2;
 }
 return g;
}
}  // namespace hurwitz
