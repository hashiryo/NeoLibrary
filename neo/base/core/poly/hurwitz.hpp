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

// 実行時に決まる法 m の Barrett 還元。mu = floor((2^64 - 1) / m)。x < 2^64 で商の誤差は 2 までなので、
// 2 回引けば m 未満になる。法がテンプレートの定数のときは % でよいが、ZMod<0> の法のような実行時の値では
// 割り算のままになるので、こちらを使う (Garner の 1 要素あたりの還元を軽くする)。
struct barrett_dyn {
 u64 m= 1, mu= ~0ull;
 barrett_dyn()= default;
 explicit barrett_dyn(u32 md): m(md), mu(md ? (~0ull) / md : 0) {}
 u32 red(u64 x) const {
  u64 r= x - u64((u128(x) * mu) >> 64) * m;
  r= r >= m ? r - m : r;
  return u32(r >= m ? r - m : r);
 }
};
// m を法とする a^e (Garner の W_t の逆元を作るだけに使う)。
inline u32 pow_mod(u64 a, u64 e, u64 m) {
 u64 r= 1;
 for (a%= m; e; e>>= 1, a= a * a % m)
  if (e & 1) r= r * a % m;
 return u32(r);
}

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
 // 割り算を使わず、素数 p の倍数 j をなめて作る。j が p で割り切れるとき、S の元を全部抜いた値は
 // strip[j] = strip[j/p]、抜いた積は ppart[j] = ppart[j/p] p になる (j/p < j なので同じ素数の列の中で前が埋まっている)。
 std::vector<u32> ppart(sz, 1), strip(sz);
 for (int k= 0; k < sz; ++k) strip[k]= u32(k);
 strip[0]= 1;
 for (int i= 0; i < ns; ++i) {
  const u32 p= S[i];
  for (u32 j= p; j < u32(sz); j+= p) strip[j]= strip[j / p], ppart[j]= ppart[j / p] * p;
 }
 // 2. u[k] = k! の p を抜いた部分 mod M。u[k] = u[k-1] strip[k] なので、逆元は 1 回の逆元から
 // u[k-1]^{-1} = u[k]^{-1} strip[k] で後ろへ配れる。どちらも 1 本の鎖にすると乗算の待ちが並ぶので、
 // 8 つずつの区画に分ける (区画の積は互いに独立、区画の前積は 1/8 の長さ、中は 8 つの鎖が並ぶ)。
 std::vector<K> u(sz);
 std::vector<u32> wa(n), wb(m);
 {
  constexpr int LU= 8;
  const size_t nb= size_t(sz + LU - 1) / LU;
  std::vector<K> bpr(nb), bst(nb);
  for (size_t bb= 0; bb < nb; ++bb) {
   K x= K(1);
   const int e= std::min(sz, int((bb + 1) * LU));
   for (int j= int(bb) * LU; j < e; ++j) x*= K(strip[j]);
   bpr[bb]= x;
  }
  K acc= K(1);
  for (size_t bb= 0; bb < nb; ++bb) {
   const K t= acc;
   acc*= bpr[bb], bst[bb]= t;
  }
  for (size_t bb= 0; bb < nb; ++bb) {
   K x= bst[bb];
   const int e= std::min(sz, int((bb + 1) * LU));
   for (int j= int(bb) * LU; j < e; ++j) x*= K(strip[j]), u[j]= x;
  }
  // 3. wa[k] = a_k u_k^{-1} mod M (b も同じ)。素数の下では、これに pv[k] = ∏ p^{-v_p(k!)} を
  // 掛けるだけで Â_k になる。逆元は後ろの区画から配る。
  K iv= u[sz - 1].inv();
  for (size_t bb= nb; bb-- > 0;) {
   const K t= iv;
   iv*= bpr[bb], bst[bb]= t;
  }
  for (size_t bb= nb; bb-- > 0;) {
   K x= bst[bb];
   const int e= std::min(sz, int((bb + 1) * LU));
   for (int j= e - 1; j >= int(bb) * LU; --j) {
    if (j < n) wa[j]= (a[j] * x).val();
    if (j < m) wb[j]= (b[j] * x).val();
    x*= K(strip[j]);
   }
  }
 }
 // 4. 各素数で、両側に p^{-v_p(k!)} を掛けて畳み込み、p^{v_p(k!)} を掛ける。
 // pv[k] と qv[k] は配列に持たない。pv は後ろから配りながら xa/xb に掛け、qv は畳み込みのあとに
 // 前から作って掛ける。配列は素数ごとに作らず、使い回す。
 std::vector<u32> res(size_t(cnt) * sz), xa(n), xb(m), bprod(size_t(sz) / 4 + 2);
 u32 qs[12];
 Crt::each(cnt, [&]<int T, u32 P>() {
  using B= barrett_q<P>;
  qs[T]= P;
  // xa/xb に要るのは k < max(n, m) までなので、pv[max(n, m) - 1] から始める。
  const int k0= std::max(n, m) - 1;
  u32 acc= 1;
  for (int i= 0; i < ns; ++i) {
   u32 e= 0;
   for (u32 x= u32(k0); x; x/= S[i]) e+= x / S[i];
   acc= B::mul(acc, B::pow(S[i], e));
  }
  acc= B::inv(acc);
  for (int k= k0; k >= 0; --k) {
   if (k < n) xa[k]= B::mul(B::red_small(wa[k]), acc);
   if (k < m) xb[k]= B::mul(B::red_small(wb[k]), acc);
   acc= B::mul(acc, ppart[k]);
  }
  u32* out= res.data() + size_t(T) * sz;
  Crt::template conv_at<P>(xa.data(), n, xb.data(), m, sz, out);
  // qv[k] = ∏_{j<=k} ppart[j] を掛ける。1 本の鎖で配ると Barrett の掛け算の待ちがそのまま並ぶので、
  // 4 つずつの区画に分ける (区画の積は互いに独立、区画の前積は 1/4 の長さ、中は 4 つの鎖が並ぶ)。
  constexpr int LB= 4;
  const int nb= (sz + LB - 1) / LB;
  for (int b= 0; b < nb; ++b) {
   u32 x= 1;
   const int e= std::min(sz, (b + 1) * LB);
   for (int j= b * LB; j < e; ++j) x= B::mul(x, ppart[j]);
   bprod[b]= x;
  }
  acc= 1;
  for (int b= 0; b < nb; ++b) {
   const u32 t= acc;
   acc= B::mul(acc, bprod[b]), bprod[b]= t;
  }
  for (int b= 0; b < nb; ++b) {
   u32 x= bprod[b];
   const int e= std::min(sz, (b + 1) * LB);
   for (int j= b * LB; j < e; ++j) x= B::mul(x, ppart[j]), out[j]= B::mul(out[j], x);
  }
 });
 // 5. Garner で S_n を戻し、u_n を掛ける。c_t は前の素数での剰余を引きずらない形で作り、
 // 要素ごとに t = 0, 1, ... と進めずに素数ごとのパスに分ける (前の素数での値に依存して待ちが続くのを避ける)。
 // パス t は res の行 t を c_t で上書きする (r_t はパス t でしか使わない)。
 // 法 q_t が実行時の値でも割り算にならないよう、還元は Barrett で行う。
 u32 wmt[12][12]= {}, iwt[12]= {}, wt[12]= {};
 barrett_dyn bq[12];
 for (int t= 0; t < cnt; ++t) {
  bq[t]= barrett_dyn(qs[t]);
  u64 w= 1;
  for (int s= 0; s < t; ++s) wmt[t][s]= u32(w % qs[t]), w= w * qs[s] % qs[t];
  iwt[t]= pow_mod(w, u64(qs[t]) - 2, qs[t]);
 }
 {
  u64 w= 1;
  for (int t= 0; t < cnt; ++t) wt[t]= u32(w % M), w= w * qs[t] % M;
 }
 for (int t= 1; t < cnt; ++t) {
  const u32 qt= qs[t];
  u32* ct= res.data() + size_t(t) * sz;
  for (int k= 0; k < sz; ++k) {
   u64 xq= 0;
   for (int s= 0; s < t; ++s) xq+= u64(res[size_t(s) * sz + k]) * wmt[t][s];
   const u32 d= ct[k] + qt - bq[t].red(xq);
   ct[k]= bq[t].red(u64(d >= qt ? d - qt : d) * iwt[t]);
  }
 }
 const barrett_dyn bm(M);
 std::vector<K> c(sz);
 for (int k= 0; k < sz; ++k) {
  u64 acc= 0;
  for (int t= 0; t < cnt; ++t) acc+= u64(res[size_t(t) * sz + k]) * wt[t];
  c[k]= u[k] * K(bm.red(acc));
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
