#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <cstdint>
#include <new>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
#include "neo/polynomial/convolve_naive.hpp"
// 多項式の補間。K は体で、ZMod なら法は素数であること。
// - interpolate_iota(y, x): y_i = f(i) (0 ≤ i < n) を満たす n-1 次以下の f の f(x)。O(n)。
// - interpolate(x, y): 点 (x_i, y_i) を通る n-1 次以下の多項式の係数 (x_i は相異なる)。O(n^2)。
// - sum_exp_poly(y, r, n): y_i = f(i) (0 ≤ i ≤ d) の d 次以下の f について Σ_{i=0}^{n-1} r^i f(i)。O(d)。
// - sum_exp_poly_limit(y, r): Σ_{i=0}^{∞} r^i f(i) (r ≠ 1、和は分母 (1 - r)^{d+1} の有理式として定める)。O(d)。
// interpolate_iota は f(x) = P/(n-1)! Σ_k y_k (-1)^{n-1-k} C(n-1, k) / (x - k) (P = Π_j (x - j)) で求める。ZMod では
// 1/k と 1/(x - k) を、i ≡ b (mod 8) の数の累積積をレーン b が受け持つ 8 本の鎖 (AVX2 の Montgomery の積) でまとめて求め、
// 1/(x - k) は後ろ向きの段でそのまま和にためて配列に置かない。C(n-1, k) は比 (n - k)/k の累積積で、8 個の区間の鎖に分けて
// から前の区間の積を掛ける。2 MB 以上の作業の配列は、Linux では 2 MB 境界の mmap に置いて huge page を頼む。
// sum_exp_poly は、S(x) = Σ_{i<x} r^i f(i) が c + r^x g(x) (g は d 次以下) と書けることを使う (min_25 の方法)。
// g の d+1 階差分が 0 なので c = r^{d+1} Σ_k C(d+1, k) (-1)^{d+1-k} r^{-k} S(k) / (1 - r)^{d+1}、g(k) = (S(k) - c) r^{-k}。
// procon-judge の yosupo-sum-of-exponential-times-polynomial と -limit で書き比べた se1 を写し、作業の配列を減らしたもの。
namespace lagrange_internal {
using u32= unsigned;
using u64= unsigned long long;
template <class T> struct HugeAlloc {
 using value_type= T;
 static constexpr std::size_t H= std::size_t(1) << 21, P= 4096;
 HugeAlloc()= default;
 template <class U> HugeAlloc(const HugeAlloc<U>&) {}
 T* allocate(std::size_t n) {
#ifdef __linux__
  if(n * sizeof(T) >= H) {
   const std::size_t bytes= (n * sizeof(T) + P - 1) & ~(P - 1);
   void* m= mmap(nullptr, bytes + H, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
   if(m == MAP_FAILED) throw std::bad_alloc();
   char *p= static_cast<char*>(m), *a= reinterpret_cast<char*>((reinterpret_cast<std::uintptr_t>(p) + H - 1) & ~(H - 1));
   if(a != p) munmap(p, a - p);
   munmap(a + bytes, p + H - a);
   madvise(a, bytes, MADV_HUGEPAGE);
   return reinterpret_cast<T*>(a);
  }
#endif
  return std::allocator<T>().allocate(n);
 }
 void deallocate(T* p, std::size_t n) {
#ifdef __linux__
  if(n * sizeof(T) >= H) {
   munmap(p, (n * sizeof(T) + P - 1) & ~(P - 1));
   return;
  }
#endif
  std::allocator<T>().deallocate(p, n);
 }
 template <class U> bool operator==(const HugeAlloc<U>&) const { return true; }
 template <class U> bool operator!=(const HugeAlloc<U>&) const { return false; }
};
using hvec= std::vector<u32, HugeAlloc<u32>>;
// x y 2^{-32} mod M を 8 個ずつ (M は奇数、x と y が 2M 未満なら結果も 2M 未満)。
[[gnu::always_inline]] inline __m256i mont8(__m256i x, __m256i y, __m256i mv, __m256i mi) {
 const __m256i te= _mm256_mul_epu32(x, y), to= _mm256_mul_epu32(_mm256_srli_epi64(x, 32), _mm256_srli_epi64(y, 32));
 const __m256i ue= _mm256_add_epi64(te, _mm256_mul_epu32(_mm256_mul_epu32(te, mi), mv));
 const __m256i uo= _mm256_add_epi64(to, _mm256_mul_epu32(_mm256_mul_epu32(to, mi), mv));
 return _mm256_blend_epi32(_mm256_srli_epi64(ue, 32), uo, 0b10101010);
}
template <class Z> struct Mont {
 u32 M;
 __m256i mv, mi;
 Mont(): M(Z::mod()) {
  u32 x= M;
  for(int i= 0; i < 5; ++i) x*= 2 - M * x;
  mv= _mm256_set1_epi32(M), mi= _mm256_set1_epi32(-x);
 }
};
// gen(t) が返す 8 個 (添字 8t..8t+7 の値、1 以上 M 未満) の逆元を、後ろ向きに back(t, 逆元の 8 個) で渡す。全部の値の積を返す。
// W は 8 steps 個以上の作業の置き場。レーンの総積は (Π 値) 2^{-32 steps} なので、全体の積は 8 本の積に 2^{256 steps} を掛ける。
template <class Z, class Gen, class Back> Z batch_inv(int steps, u32* W, const Mont<Z>& mt, Gen gen, Back back) {
 if(steps == 0) return Z(1);
 __m256i acc= _mm256_set1_epi32(1);
 for(int t= 0; t < steps; ++t) acc= mont8(acc, gen(t), mt.mv, mt.mi), _mm256_storeu_si256((__m256i*)(W + 8 * t), acc);
 alignas(32) u32 c[8];
 _mm256_store_si256((__m256i*)c, acc);
 Z prod(1);
 for(int b= 0; b < 8; ++b) {
  const Z v(c[b]);
  prod*= v, c[b]= v.inv().val();
 }
 prod*= Z(u64(1) << 32).pow(u64(8) * steps);
 __m256i cur= _mm256_load_si256((const __m256i*)c);
 for(int t= steps - 1; t >= 0; --t) {
  const __m256i pre= t ? _mm256_loadu_si256((const __m256i*)(W + 8 * (t - 1))) : _mm256_set1_epi32(1);
  const __m256i r= mont8(cur, pre, mt.mv, mt.mi);
  back(t, _mm256_min_epu32(r, _mm256_sub_epi32(r, mt.mv)));
  cur= mont8(cur, gen(t), mt.mv, mt.mi);
 }
 return prod;
}
// C(N, k) (0 ≤ k ≤ N) を c[0..N] に置き (c は 8 の倍数ぶん 0 で延ばす)、N! を返す。比 (N - k + 1)/k の累積積を 8 個の区間に
// 分けて区間ごとに求め、あとで前の区間の積を掛ける。
template <class Z> Z binom_row(int N, hvec& c, hvec& W, const Mont<Z>& mt) {
 const int steps= (N + 7) / 8;
 c.assign(8 * steps + 8, 0), W.resize(8 * steps + 8);
 c[0]= Z(1).val();
 if(N == 0) return Z(1);
 const __m256i one= _mm256_set1_epi32(1), lane= _mm256_setr_epi32(1, 2, 3, 4, 5, 6, 7, 8);
 // k = 8t + b + 1 (k > N の分は 1 にして積に効かせない)
 auto gen= [&](int t) {
  const __m256i k= _mm256_add_epi32(lane, _mm256_set1_epi32(8 * t));
  return _mm256_blendv_epi8(k, one, _mm256_cmpgt_epi32(k, _mm256_set1_epi32(N)));
 };
 const Z fact= batch_inv<Z>(steps, W.data(), mt, gen, [&](int t, __m256i iv) { _mm256_storeu_si256((__m256i*)(c.data() + 1 + 8 * t), iv); });
 // c[k] = (N - k + 1) / k の累積積
 constexpr int S= 8;
 const int len= (N + S - 1) / S;
 Z acc[S];
 for(int s= 0; s < S; ++s) acc[s]= Z(1);
 for(int t= 0; t < len; ++t)
  for(int s= 0; s < S; ++s) {
   const int k= 1 + s * len + t;
   if(k <= N) acc[s]*= Z::raw(c[k]) * Z(N - k + 1), c[k]= acc[s].val();
  }
 Z carry= acc[0];
 for(int s= 1; s < S && 1 + s * len <= N; ++s) {
  const int lo= 1 + s * len, hi= std::min(N, lo + len - 1);
  const auto cf= carry.fixed();
  for(int k= lo; k <= hi; ++k) c[k]= (Z::raw(c[k]) * cf).val();
  carry*= acc[s];
 }
 for(int k= N + 1; k < int(c.size()); ++k) c[k]= 0;
 return fact;
}
// y (長さ K) の 0..K-1 での値の多項式の x での値 (x mod M ≥ K)。c は C(K-1, k)、fact は (K-1)!。
template <class Z> Z iota_core(const std::vector<Z>& y, u32 xm, const hvec& c, Z fact, hvec& W, const Mont<Z>& mt) {
 const int K= int(y.size()), steps= (K + 7) / 8;
 // w[k] = y_k C(K-1, k)
 hvec w(8 * steps + 8, 0);
 for(int k= 0; k < K; ++k) w[k]= (y[k] * Z::raw(c[k])).val();
 W.resize(8 * steps + 8);
 const __m256i one= _mm256_set1_epi32(1), lane= _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
 auto gen= [&](int t) {
  const __m256i k= _mm256_add_epi32(lane, _mm256_set1_epi32(8 * t));
  const __m256i v= _mm256_sub_epi32(_mm256_set1_epi32(xm), k);  // x - k (k < K ≤ x mod M なので正)
  return _mm256_blendv_epi8(v, one, _mm256_cmpgt_epi32(k, _mm256_set1_epi32(K - 1)));
 };
 const __m256i neg= ((K - 1) & 1) ? _mm256_setr_epi32(-1, 0, -1, 0, -1, 0, -1, 0) : _mm256_setr_epi32(0, -1, 0, -1, 0, -1, 0, -1);
 const __m256i m2= _mm256_set1_epi32(2 * mt.M), lo32= _mm256_set1_epi64x(0xffffffffll);
 __m256i ae= _mm256_setzero_si256(), ao= _mm256_setzero_si256();
 const Z P= batch_inv<Z>(steps, W.data(), mt, gen, [&](int t, __m256i iv) {
  __m256i v= mont8(_mm256_loadu_si256((const __m256i*)(w.data() + 8 * t)), iv, mt.mv, mt.mi);
  v= _mm256_blendv_epi8(v, _mm256_sub_epi32(m2, v), neg);
  ae= _mm256_add_epi64(ae, _mm256_and_si256(v, lo32)), ao= _mm256_add_epi64(ao, _mm256_srli_epi64(v, 32));
 });
 alignas(32) u64 s[4];
 _mm256_store_si256((__m256i*)s, _mm256_add_epi64(ae, ao));
 // Montgomery の積 1 回の 2^{-32} を戻す
 const Z sum= (Z(s[0]) + Z(s[1]) + Z(s[2]) + Z(s[3])) * Z(u64(1) << 32);
 return P * fact.inv() * sum;
}
template <class Z> Z iota_zmod(const std::vector<Z>& y, u64 x) {
 const int K= int(y.size());
 if(x < u64(K)) return y[x];
 const u32 xm= u32(x % Z::mod());
 if(xm < u32(K)) return y[xm];
 const Mont<Z> mt;
 hvec c, W;
 const Z fact= binom_row<Z>(K - 1, c, W, mt);
 return iota_core<Z>(y, xm, c, fact, W, mt);
}
template <class K> K iota_generic(const std::vector<K>& y, K x) {
 const int n= int(y.size());
 std::vector<K> suf(n + 1, K(1));
 for(int i= n - 1; i >= 0; --i) suf[i]= suf[i + 1] * (x - K(i));
 std::vector<K> fact(n, K(1));
 for(int i= 1; i < n; ++i) fact[i]= fact[i - 1] * K(i);
 K pre(1), res(0);
 for(int i= 0; i < n; ++i) {
  K t= y[i] * pre * suf[i + 1] / (fact[i] * fact[n - 1 - i]);
  res+= ((n - 1 - i) & 1) ? -t : t;
  pre*= x - K(i);
 }
 return res;
}
// S(0..d+1) (S(k) = Σ_{i<k} r^i y_i)。
template <class Z> std::vector<Z> prefix_sums(const std::vector<Z>& y, Z r) {
 const int D= int(y.size()) - 1;
 std::vector<Z> S(D + 2);
 const Z r4= r.pow(4);
 const auto f4= r4.fixed();
 Z q[4]= {Z(1), r, r * r, r * r * r};
 S[0]= Z(0);
 for(int i= 0; i <= D; ++i) S[i + 1]= S[i] + q[i & 3] * y[i], q[i & 3]*= f4;
 return S;
}
// S(x) = c + r^x g(x) の c (bc は C(d, k) の行、d = |S| - 2)。g が要れば g(0..d) を置く。
template <class Z> Z const_term(const std::vector<Z>& S, Z r, const hvec& bc, std::vector<Z>* g) {
 const int D= int(S.size()) - 2;
 const Z ir= r.inv(), ir4= ir.pow(4);
 const auto f4= ir4.fixed();
 Z q[4]= {Z(1), ir, ir * ir, ir * ir * ir};
 Z sp(0), sn(0);
 for(int k= 0; k <= D + 1; ++k) {
  const Z cb= Z::raw(k <= D ? bc[k] : 0) + Z::raw(k ? bc[k - 1] : 0);
  const Z w= cb * q[k & 3] * S[k];
  if((D + 1 - k) & 1) sn+= w;
  else sp+= w;
  q[k & 3]*= f4;
 }
 const Z c= (sp - sn) * r.pow(D + 1) * (Z(1) - r).pow(D + 1).inv();
 if(g) {
  g->resize(D + 1);
  Z p[4]= {Z(1), ir, ir * ir, ir * ir * ir};
  for(int k= 0; k <= D; ++k) (*g)[k]= (S[k] - c) * p[k & 3], p[k & 3]*= f4;
 }
 return c;
}
}  // namespace lagrange_internal
// y_i = f(i) (0 ≤ i < n) を満たす n-1 次以下の f の f(x)。
template <class K> K interpolate_iota(const std::vector<K>& y, unsigned long long x) {
 if(y.empty()) return K();
 if constexpr(poly_internal::is_zmod_v<K>) return lagrange_internal::iota_zmod(y, x);
 else return x < y.size() ? y[x] : lagrange_internal::iota_generic(y, K(x));
}
// 点 (x_i, y_i) を通る n-1 次以下の多項式の係数。f = Σ_i y_i / P'(x_i) P(X)/(X - x_i) (P = Π (X - x_j)) を、組立除法で求める。
template <class K> std::vector<K> interpolate(const std::vector<K>& x, const std::vector<K>& y) {
 const int n= int(x.size());
 if(n == 0) return {};
 std::vector<K> P(n + 1, K());
 P[0]= K(1);
 for(int j= 0; j < n; ++j) {
  for(int k= j + 1; k > 0; --k) P[k]= P[k - 1] - x[j] * P[k];
  P[0]= -x[j] * P[0];
 }
 std::vector<K> f(n, K()), q(n);
 for(int i= 0; i < n; ++i) {
  // q = P / (X - x_i)、q_{n-1} = 1、q_{k-1} = Q_k + x_i q_k
  q[n - 1]= P[n];
  for(int k= n - 1; k > 0; --k) q[k - 1]= P[k] + x[i] * q[k];
  K d(0);
  for(int k= n - 1; k >= 0; --k) d= d * x[i] + q[k];  // P'(x_i) = q(x_i)
  const K w= y[i] / d;
  for(int k= 0; k < n; ++k) f[k]+= w * q[k];
 }
 return f;
}
// Σ_{i=0}^{n-1} r^i f(i) (y_i = f(i)、0 ≤ i ≤ d)。0^0 = 1。
template <class K> K sum_exp_poly(const std::vector<K>& y, K r, unsigned long long n) {
 using namespace lagrange_internal;
 if(n == 0 || y.empty()) return K();
 if(r == K()) return y[0];
 const int D= int(y.size()) - 1;
 if(r == K(1)) {
  std::vector<K> F(D + 2);
  F[0]= K();
  for(int k= 0; k <= D; ++k) F[k + 1]= F[k] + y[k];
  return interpolate_iota(F, n);
 }
 if constexpr(poly_internal::is_zmod_v<K>) {
  const auto S= prefix_sums(y, r);
  if(n <= u64(D) + 1) return S[n];
  const Mont<K> mt;
  hvec bc, W;
  const K fact= binom_row<K>(D, bc, W, mt);
  std::vector<K> g;
  const K c= const_term(S, r, bc, &g);
  const u32 xm= u32(n % K::mod());
  const K gn= xm <= u32(D) ? g[xm] : iota_core<K>(g, xm, bc, fact, W, mt);
  return c + r.pow(n) * gn;
 } else {
  std::vector<K> S(D + 2);
  S[0]= K();
  K q(1);
  for(int i= 0; i <= D; ++i) S[i + 1]= S[i] + q * y[i], q*= r;
  if(n <= u64(D) + 1) return S[n];
  // C(D+1, k) と r^{-k}
  std::vector<K> C(D + 2, K());
  C[0]= K(1);
  for(int k= 1; k <= D + 1; ++k) C[k]= C[k - 1] * K(D + 2 - k) / K(k);
  const K ir= K(1) / r;
  K num(0), rk(1);
  for(int k= 0; k <= D + 1; ++k) {
   const K w= C[k] * rk * S[k];
   num+= ((D + 1 - k) & 1) ? -w : w;
   rk*= ir;
  }
  K rd(1), om(1);
  for(int k= 0; k <= D; ++k) rd*= r, om*= K(1) - r;
  const K c= num * rd / om;
  std::vector<K> g(D + 1);
  rk= K(1);
  for(int k= 0; k <= D; ++k) g[k]= (S[k] - c) * rk, rk*= ir;
  K rn(1), b= r;
  for(auto e= n; e; e>>= 1, b*= b)
   if(e & 1) rn*= b;
  return c + rn * interpolate_iota(g, n);
 }
}
// Σ_{i=0}^{∞} r^i f(i) (r ≠ 1)。
template <class K> K sum_exp_poly_limit(const std::vector<K>& y, K r) {
 using namespace lagrange_internal;
 if(y.empty()) return K();
 if(r == K()) return y[0];
 const int D= int(y.size()) - 1;
 if constexpr(poly_internal::is_zmod_v<K>) {
  const auto S= prefix_sums(y, r);
  const Mont<K> mt;
  hvec bc, W;
  binom_row<K>(D, bc, W, mt);
  return const_term<K>(S, r, bc, nullptr);
 } else {
  std::vector<K> S(D + 2);
  S[0]= K();
  K q(1);
  for(int i= 0; i <= D; ++i) S[i + 1]= S[i] + q * y[i], q*= r;
  std::vector<K> C(D + 2, K());
  C[0]= K(1);
  for(int k= 1; k <= D + 1; ++k) C[k]= C[k - 1] * K(D + 2 - k) / K(k);
  const K ir= K(1) / r;
  K num(0), rk(1);
  for(int k= 0; k <= D + 1; ++k) {
   const K w= C[k] * rk * S[k];
   num+= ((D + 1 - k) & 1) ? -w : w;
   rk*= ir;
  }
  K rd(1), om(1);
  for(int k= 0; k <= D; ++k) rd*= r, om*= K(1) - r;
  return num * rd / om;
 }
}
