#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>
#include "neo/number_theory/sqrt_mod.hpp"
#include "neo/polynomial/convolve_naive.hpp"
// 疎な形式的冪級数の演算 (namespace sparse_fps)。0 でない項が K 個の f について、先頭 n 項を O(nK) で求める。
// f は密な vector で受け (0 の項は読み飛ばす)、答えは長さ n の vector で返す。K は体で、log、exp、冪、平方根は 1/i
// (i < n) を使うので、ZMod なら法は n より大きい素数であること。
// ZMod は値を M 未満の u32 で持ち、項の積を u64 に 15 個までためてから剰余を取る (15 M^2 < 2^64)。i が項の次数を越えるたびに
// 使う項が 1 つ増えるので、区間ごとに項の数をテンプレートの引数にする (16 項以上は実行時の数のループ)。いちばん次数の低い項の
// 次数が 1 のときは、そこが前の値を待つ鎖になるので、前の値はレジスタに持ち、掛け算は Shoup の方法か 1 回の Barrett にする。
// 1/i の表は、i ≡ b (mod 8) の数の累積積をレーン b が受け持つ 8 本の鎖を AVX2 の Montgomery の積で作り、レーンごとの総積の
// 逆元から後ろ向きに戻して求める。
// procon-judge の yosupo の inv/exp/log/pow/sqrt-of-formal-power-series-sparse で書き比べた sp1 を写したもの。f/g、f^a g^b、
// exp(f/g) は、今の Library の sparse_fps.hpp と同じ漸化式を、値の型のまま回す。
namespace sparse_fps_internal {
using u32= unsigned;
using u64= unsigned long long;
template <class Z> struct Ar {
 u32 M;
 u64 BX;  // floor((2^64 - 1) / M)
 Ar(): M(Z::mod()), BX(~u64(0) / Z::mod()) {}
 // u64 の t を M 未満に
 [[gnu::always_inline]] u32 red(u64 t) const {
  u64 r= t - u64((__uint128_t(t) * BX) >> 64) * M;
  r= r >= 2 * u64(M) ? r - 2 * u64(M) : r;
  return u32(r >= M ? r - M : r);
 }
 [[gnu::always_inline]] u32 quo(u32 x) const { return u32((u64(x) << 32) / M); }
 // a に x を掛けて 2M 未満で返す (q = quo(x))。
 [[gnu::always_inline]] u32 shoup(u32 a, u32 x, u32 q) const { return a * x - u32((u64(a) * q) >> 32) * M; }
 [[gnu::always_inline]] u32 sub1(u32 r) const { return r >= M ? r - M : r; }
};
// x y 2^{-32} mod M を 8 個ずつ (M は奇数、x と y が 2M 未満なら結果も 2M 未満)。
[[gnu::always_inline]] inline __m256i mont8(__m256i x, __m256i y, __m256i mv, __m256i mi) {
 const __m256i te= _mm256_mul_epu32(x, y), to= _mm256_mul_epu32(_mm256_srli_epi64(x, 32), _mm256_srli_epi64(y, 32));
 const __m256i ue= _mm256_add_epi64(te, _mm256_mul_epu32(_mm256_mul_epu32(te, mi), mv));
 const __m256i uo= _mm256_add_epi64(to, _mm256_mul_epu32(_mm256_mul_epu32(to, mi), mv));
 return _mm256_blend_epi32(_mm256_srli_epi64(ue, 32), uo, 0b10101010);
}
// iv[i] = 1/i mod M (1 ≤ i < n、iv[0] = 0、M は n より大きい奇素数)。前向きに累積積を置き、レーンごとの総積の逆元から
// 後ろ向きに戻す (Montgomery の 2^{-32} は前向きと後ろ向きで打ち消し合う)。
template <class Z> poly_internal::hvec inverses(int n) {
 if(n <= 1) return poly_internal::hvec(std::max(n, 0));
 const u32 M= Z::mod();
 u32 minv= M;
 for(int i= 0; i < 5; ++i) minv*= 2 - M * minv;
 const __m256i mv= _mm256_set1_epi32(M), mi= _mm256_set1_epi32(-minv), eight= _mm256_set1_epi32(8);
 const int steps= (n - 1 + 7) / 8;
 poly_internal::hvec iv(8 * steps + 8);
 __m256i acc= _mm256_set1_epi32(1), idx= _mm256_setr_epi32(1, 2, 3, 4, 5, 6, 7, 8);
 for(int t= 0; t < steps; ++t) {
  acc= mont8(acc, idx, mv, mi);
  _mm256_storeu_si256((__m256i*)(iv.data() + 1 + 8 * t), acc);
  idx= _mm256_add_epi32(idx, eight);
 }
 alignas(32) u32 a[8];
 _mm256_store_si256((__m256i*)a, acc);
 for(int b= 0; b < 8; ++b) a[b]= Z(a[b]).inv().val();
 __m256i cur= _mm256_load_si256((const __m256i*)a);
 idx= _mm256_sub_epi32(idx, eight);
 for(int t= steps - 1; t >= 0; --t) {
  const __m256i pre= t ? _mm256_loadu_si256((const __m256i*)(iv.data() + 1 + 8 * (t - 1))) : _mm256_set1_epi32(1);
  const __m256i r= mont8(cur, pre, mv, mi);
  _mm256_storeu_si256((__m256i*)(iv.data() + 1 + 8 * t), _mm256_min_epu32(r, _mm256_sub_epi32(r, mv)));
  cur= mont8(cur, idx, mv, mi);
  idx= _mm256_sub_epi32(idx, eight);
 }
 iv.resize(n), iv[0]= 0;
 return iv;
}
template <int... Ts, class F> void dispatch(int t, std::integer_sequence<int, Ts...>, F&& f) { ((t == Ts ? f(std::integral_constant<int, Ts>{}) : void()), ...); }
constexpr int TMAX= 15;
// h_i = (h_i の初期値) + Σ_{k<t} w_k h_{i - e_k} を [lo, hi) で求める (h は M 未満)。
template <int T, class Z> void seg_lin(u32* h, int lo, int hi, const int* e, const u32* w, int t, const Ar<Z>& ar) {
 if(lo >= hi) return;
 const int tt= T >= 0 ? T : t;
 if(e[0] == 1) {
  const u32 w0= w[0], q0= ar.quo(w0);
  u32 prev= h[lo - 1];
  for(int i= lo; i < hi; ++i) {
   u32 rest;
   if constexpr(T >= 0) {
    u64 s= h[i];
    for(int k= 1; k < T; ++k) s+= u64(w[k]) * h[i - e[k]];
    rest= ar.red(s);
   } else {
    u64 s= h[i];
    for(int k= 1, c= 0; k < tt; ++k) {
     s+= u64(w[k]) * h[i - e[k]];
     if(++c == TMAX - 1) s= ar.red(s), c= 0;
    }
    rest= ar.red(s);
   }
   u32 r= rest + ar.shoup(prev, w0, q0);  // 3M 未満
   r= r >= 2 * ar.M ? r - 2 * ar.M : r;
   h[i]= prev= ar.sub1(r);
  }
 } else {
  for(int i= lo; i < hi; ++i) {
   u64 s= h[i];
   if constexpr(T >= 0) {
    for(int k= 0; k < T; ++k) s+= u64(w[k]) * h[i - e[k]];
   } else {
    for(int k= 0, c= 0; k < tt; ++k) {
     s+= u64(w[k]) * h[i - e[k]];
     if(++c == TMAX - 1) s= ar.red(s), c= 0;
    }
   }
   h[i]= ar.red(s);
  }
 }
}
// i g_i = Σ_{k<t} b_k g_{i - e_k} を [lo, hi) で求める (g は M 未満、iv は 1/i)。
template <int T, class Z> void seg_exp(u32* g, const u32* iv, int lo, int hi, const int* e, const u32* b, int t, const Ar<Z>& ar) {
 if(lo >= hi) return;
 const int tt= T >= 0 ? T : t;
 auto sum= [&](int i, int k0) {
  u64 s= 0;
  for(int k= k0, c= 0; k < tt; ++k) {
   s+= u64(b[k]) * g[i - e[k]];
   if(T < 0 && ++c == TMAX) s= ar.red(s), c= 0;
  }
  return ar.red(s);
 };
 if(e[0] == 1) {
  // g_i = (b_0 / i) g_{i-1} + (Σ_{k≥1} b_k g_{i-e_k}) / i。鎖に乗るのは前の値との積 1 回だけ。
  const u32 b0= b[0], q0= ar.quo(b0);
  u32 prev= g[lo - 1];
  for(int i= lo; i < hi; ++i) {
   const u32 ivi= iv[i], c= ar.sub1(ar.shoup(ivi, b0, q0)), d= ar.red(u64(sum(i, 1)) * ivi);
   u32 r= ar.red(u64(prev) * c) + d;
   g[i]= prev= ar.sub1(r);
  }
 } else {
  for(int i= lo; i < hi; ++i) g[i]= ar.red(u64(sum(i, 0)) * iv[i]);
 }
}
// u = 1 + Σ_k c_k x^{d_k} の u^α を [lo, hi) で求める。v' u = α u' v から
// v_i = (Σ_k al_k v_{i-d_k}) / i - Σ_k c_k v_{i-d_k} (al_k = c_k (α + 1) d_k)。
template <int T, class Z> void seg_pow(u32* v, const u32* iv, int lo, int hi, const int* d, const u32* c, const u32* al, int t, const Ar<Z>& ar) {
 if(lo >= hi) return;
 const int tt= T >= 0 ? T : t;
 auto sums= [&](int i, int k0, u32& A, u32& B) {
  u64 sa= 0, sb= 0;
  for(int k= k0, cnt= 0; k < tt; ++k) {
   sa+= u64(al[k]) * v[i - d[k]], sb+= u64(c[k]) * v[i - d[k]];
   if(T < 0 && ++cnt == TMAX) sa= ar.red(sa), sb= ar.red(sb), cnt= 0;
  }
  A= ar.red(sa), B= ar.red(sb);
 };
 if(d[0] == 1) {
  // v_i = (al_0 / i - c_0) v_{i-1} + (Σ_{k≥1} al_k v_{i-d_k}) / i - Σ_{k≥1} c_k v_{i-d_k}
  const u32 a0= al[0], q0= ar.quo(a0), c0= ar.M - c[0];
  u32 prev= v[lo - 1];
  for(int i= lo; i < hi; ++i) {
   u32 A, B;
   sums(i, 1, A, B);
   const u32 ivi= iv[i];
   const u32 gm= ar.sub1(ar.sub1(ar.shoup(ivi, a0, q0)) + c0);
   const u32 dl= ar.sub1(ar.red(u64(A) * ivi) + ar.M - B);
   v[i]= prev= ar.sub1(ar.red(u64(prev) * gm) + dl);
  }
 } else {
  for(int i= lo; i < hi; ++i) {
   u32 A, B;
   sums(i, 0, A, B);
   v[i]= ar.sub1(ar.red(u64(A) * iv[i]) + ar.M - B);
  }
 }
}
// f の 0 でない項 (次数 < lim) を [from, ...) から集める。
template <class K> std::vector<std::pair<int, K>> terms(const std::vector<K>& f, int from, int lim) {
 std::vector<std::pair<int, K>> t;
 for(int i= from, e= std::min<int>(lim, int(f.size())); i < e; ++i)
  if(!(f[i] == K())) t.emplace_back(i, f[i]);
 return t;
}
template <class K> int min_deg(const std::vector<K>& f) {
 for(int i= 0; i < int(f.size()); ++i)
  if(!(f[i] == K())) return i;
 return -1;
}
// ZMod: h_i = init_i + Σ_k w_k h_{i - e_k} (0 ≤ i < n、h_0 = init_0)。
template <class Z> void run_lin(u32* h, int n, const std::vector<std::pair<int, Z>>& ts, const Ar<Z>& ar) {
 const int t= int(ts.size());
 std::vector<int> e(t);
 std::vector<u32> w(t);
 for(int k= 0; k < t; ++k) e[k]= ts[k].first, w[k]= ts[k].second.val();
 for(int j= 1; j <= t; ++j) {
  const int lo= e[j - 1], hi= j == t ? n : e[j];
  if(j <= TMAX) dispatch(j, std::make_integer_sequence<int, TMAX + 1>{}, [&](auto T) { seg_lin<T, Z>(h, lo, hi, e.data(), w.data(), j, ar); });
  else seg_lin<-1, Z>(h, lo, hi, e.data(), w.data(), j, ar);
 }
}
// 値の型のまま: u^a v^{b-1} の先頭 N+1 項 (今の Library の pw_mul_pw と同じ漸化式)。p と q は u - 1、v - 1 の項。
template <class K> void pw_mul_pw(K* f, const std::vector<std::pair<int, K>>& p, const std::vector<std::pair<int, K>>& q, K a, K b, int N) {
 std::vector<K> x(N + 1, K()), y(std::max(N, 1), K());
 f[0]= K(1);
 for(int i= 0, k= 1; i < N; i= k++) {
  for(auto&& [j, v]: q) {
   if(i < j) break;
   y[i]+= f[i - j] * v;
  }
  y[i]= (y[i] + f[i]) * a;
  for(auto&& [j, v]: p) {
   if(k < j) break;
   x[k]+= v * (y[k - j] * K(j) - x[k - j]);
  }
  f[k]= x[k];
  for(auto&& [j, v]: q) {
   if(k < j) break;
   f[k]+= v * (b * K(j) - K(k)) * f[k - j];
  }
  f[k]/= K(k);
 }
}
// K の値として整数 m を作る (負でもよい)。
template <class K> K from_ll(long long m) {
 if constexpr(poly_internal::is_zmod_v<K>) return K(m);
 else return m < 0 ? -K((unsigned long long)(-m)) : K((unsigned long long)m);
}
template <class K> K pow_ll(K a, long long k) {
 if(k < 0) a= K(1) / a, k= -k;
 K r= K(1);
 for(; k; k>>= 1, a*= a)
  if(k & 1) r*= a;
 return r;
}
}  // namespace sparse_fps_internal
namespace sparse_fps {
// f / g の先頭 n 項 (g は疎、g の最も低い次数の係数で割る。g(0) ≠ 0)。f は密でよい。
template <class K> std::vector<K> div(const std::vector<K>& f, const std::vector<K>& g, int n) {
 using namespace sparse_fps_internal;
 assert(!g.empty() && !(g[0] == K()));
 const K ig= K(1) / g[0];
 auto ts= terms(g, 1, n);
 for(auto& [j, v]: ts) v= -v * ig;
 if constexpr(poly_internal::is_zmod_v<K>) {
  const Ar<K> ar;
  std::vector<K> res(n);  // 答えの置き場に直接書く (u32 の列として)
  u32* h= poly_internal::raw(res);
  const auto x= ig.fixed();
  for(int i= 0; i < std::min<int>(n, int(f.size())); ++i) h[i]= (f[i] * x).val();
  if(!ts.empty()) run_lin<K>(h, n, ts, ar);
  return res;
 } else {
  std::vector<K> h(n, K());
  for(int i= 0; i < std::min<int>(n, int(f.size())); ++i) h[i]= f[i] * ig;
  for(int i= 1; i < n; ++i)
   for(auto&& [j, v]: ts) {
    if(i < j) break;
    h[i]+= v * h[i - j];
   }
  return h;
 }
}
// 1 / f の先頭 n 項 (f(0) ≠ 0)。
template <class K> std::vector<K> inv(const std::vector<K>& f, int n) { return div<K>(std::vector<K>{K(1)}, f, n); }
// log f の先頭 n 項 (f(0) = 1)。h_i = i (log f)_i が h_i = i f_i - Σ_{k≥1} f_k h_{i-k} を満たすので、h を求めてから 1/i を掛ける。
template <class K> std::vector<K> log(const std::vector<K>& f, int n) {
 using namespace sparse_fps_internal;
 assert(!f.empty() && f[0] == K(1));
 auto ts= terms(f, 1, n);
 if constexpr(poly_internal::is_zmod_v<K>) {
  const Ar<K> ar;
  std::vector<K> res(n);
  if(ts.empty()) return res;
  u32* h= poly_internal::raw(res);
  for(auto& [j, v]: ts) h[j]= (v * K(j)).val(), v= -v;
  run_lin<K>(h, n, ts, ar);
  const auto iv= inverses<K>(n);
  for(int i= 0; i < n; ++i) h[i]= ar.red(u64(h[i]) * iv[i]);
  return res;
 } else {
  std::vector<K> h(n, K());
  for(auto& [j, v]: ts) h[j]= v * K(j);
  for(int i= 1; i < n; ++i)
   for(auto&& [j, v]: ts) {
    if(i <= j) break;
    h[i]-= v * h[i - j];
   }
  for(int i= 1; i < n; ++i) h[i]/= K(i);
  if(n) h[0]= K();
  return h;
 }
}
// exp f の先頭 n 項 (f(0) = 0)。i g_i = Σ_k k f_k g_{i-k}。
template <class K> std::vector<K> exp(const std::vector<K>& f, int n) {
 using namespace sparse_fps_internal;
 assert(f.empty() || f[0] == K());
 if(n <= 0) return {};
 auto ts= terms(f, 1, n);
 for(auto& [j, v]: ts) v*= K(j);
 if constexpr(poly_internal::is_zmod_v<K>) {
  const Ar<K> ar;
  std::vector<K> res(n);
  res[0]= K(1);
  if(ts.empty()) return res;
  u32* g= poly_internal::raw(res);
  g[0]= K(1).val();
  const int t= int(ts.size());
  std::vector<int> e(t);
  std::vector<u32> b(t);
  for(int k= 0; k < t; ++k) e[k]= ts[k].first, b[k]= ts[k].second.val();
  const auto iv= inverses<K>(n);
  for(int j= 1; j <= t; ++j) {
   const int lo= e[j - 1], hi= j == t ? n : e[j];
   if(j <= TMAX) dispatch(j, std::make_integer_sequence<int, TMAX + 1>{}, [&](auto T) { seg_exp<T, K>(g, iv.data(), lo, hi, e.data(), b.data(), j, ar); });
   else seg_exp<-1, K>(g, iv.data(), lo, hi, e.data(), b.data(), j, ar);
  }
  return res;
 } else {
  std::vector<K> g(n, K());
  g[0]= K(1);
  for(int i= 1; i < n; ++i) {
   for(auto&& [j, v]: ts) {
    if(i < j) break;
    g[i]+= v * g[i - j];
   }
   g[i]/= K(i);
  }
  return g;
 }
}
// f^k の先頭 n 項 (k は負でもよく、そのときは f(0) ≠ 0)。0^0 = 1。f = a x^e u (u(0) = 1) として u^k を v' u = k u' v の
// 漸化式で求める。u^k の係数は k の多項式なので、ZMod では k を法で割った余りでよい。
template <class K> std::vector<K> pow(const std::vector<K>& f, long long k, int n) {
 using namespace sparse_fps_internal;
 std::vector<K> F(std::max(n, 0), K());
 if(n <= 0) return F;
 if(k == 0) return F[0]= K(1), F;
 const int e0= min_deg(f);
 assert(!(e0 < 0 && k < 0));
 if(e0 < 0) return F;
 assert(!(e0 > 0 && k < 0));
 if(e0 > 0 && (unsigned long long)k > (unsigned long long)(n - 1) / e0) return F;
 const int o= int((long long)e0 * k), L= n - o;
 const K a0= f[e0], ia= K(1) / a0, am= pow_ll(a0, k);
 auto ts= terms(f, e0 + 1, e0 + L);
 for(auto& [j, v]: ts) j-= e0, v*= ia;
 if constexpr(poly_internal::is_zmod_v<K>) {
  const Ar<K> ar;
  u32* v= poly_internal::raw(F) + o;  // 答えの置き場に直接書く
  v[0]= K(1).val();
  const int t= int(ts.size());
  if(t) {
   const u32 m1= (from_ll<K>(k) + K(1)).val();
   std::vector<int> d(t);
   std::vector<u32> c(t), al(t);
   for(int q= 0; q < t; ++q) d[q]= ts[q].first, c[q]= ts[q].second.val(), al[q]= (ts[q].second * K::raw(m1) * K(d[q])).val();
   const auto iv= inverses<K>(L);
   for(int j= 1; j <= t; ++j) {
    const int lo= d[j - 1], hi= j == t ? L : d[j];
    if(j <= TMAX) dispatch(j, std::make_integer_sequence<int, TMAX + 1>{}, [&](auto T) { seg_pow<T, K>(v, iv.data(), lo, hi, d.data(), c.data(), al.data(), j, ar); });
    else seg_pow<-1, K>(v, iv.data(), lo, hi, d.data(), c.data(), al.data(), j, ar);
   }
  }
  if(!(am == K(1))) {
   const auto fa= am.fixed();
   for(int i= o; i < n; ++i) F[i]*= fa;
  }
 } else {
  pw_mul_pw<K>(F.data() + o, ts, {}, from_ll<K>(k), K(), L - 1);
  for(int i= o; i < n; ++i) F[i]*= am;
 }
 return F;
}
// f^k g^l の先頭 n 項 (k、l は負でもよく、負なら対応する f(0) か g(0) は 0 でない)。
template <class K> std::vector<K> pow_mul_pow(const std::vector<K>& f, long long k, const std::vector<K>& g, long long l, int n) {
 using namespace sparse_fps_internal;
 if(k == 0) return pow(g, l, n);
 if(l == 0) return pow(f, k, n);
 std::vector<K> F(std::max(n, 0), K());
 const int x= min_deg(f), y= min_deg(g);
 assert(!(x < 0 && k < 0) && !(y < 0 && l < 0));
 if(n <= 0 || x < 0 || y < 0) return F;
 const __int128 o= __int128(x) * k + __int128(y) * l;
 assert(o >= 0);
 if(o >= n) return F;
 const int oo= int(o), L= n - oo;
 const K p0= f[x], q0= g[y], ip= K(1) / p0, iq= K(1) / q0, a= pow_ll(p0, k) * pow_ll(q0, l);
 auto p= terms(f, x + 1, x + L), q= terms(g, y + 1, y + L);
 for(auto& [j, v]: p) j-= x, v*= ip;
 for(auto& [j, v]: q) j-= y, v*= iq;
 pw_mul_pw<K>(F.data() + oo, p, q, from_ll<K>(k), from_ll<K>(l) + K(1), L - 1);
 for(int i= oo; i < n; ++i) F[i]*= a;
 return F;
}
// (f / g)^k の先頭 n 項。
template <class K> std::vector<K> pow_of_div(const std::vector<K>& f, const std::vector<K>& g, long long k, int n) { return pow_mul_pow(f, k, g, -k, n); }
// g^2 = f となる g の先頭 n 項を 1 つ。無ければ空の vector。K は ZMod (法は素数)。
template <class K> std::vector<K> sqrt(const std::vector<K>& f, int n) {
 static_assert(poly_internal::is_zmod_v<K>, "sparse_fps::sqrt は ZMod だけ");
 using namespace sparse_fps_internal;
 std::vector<K> F(std::max(n, 0), K());
 const int e0= min_deg(f);
 if(e0 < 0 || n <= 0) return F;
 if(e0 & 1) return {};
 const long long s= sqrt_mod(f[e0].val(), K::mod());
 if(s < 0) return {};
 const int o= e0 / 2;
 if(o >= n) return F;
 const int L= n - o;
 const K ia= K(1) / f[e0], i2= K(1) / K(2);
 auto ts= terms(f, e0 + 1, e0 + L);
 for(auto& [j, v]: ts) j-= e0, v*= ia;
 const Ar<K> ar;
 u32* v= poly_internal::raw(F) + o;  // 答えの置き場に直接書く
 v[0]= K(1).val();
 const int t= int(ts.size());
 if(t) {
  const u32 m1= (i2 + K(1)).val();
  std::vector<int> d(t);
  std::vector<u32> c(t), al(t);
  for(int q= 0; q < t; ++q) d[q]= ts[q].first, c[q]= ts[q].second.val(), al[q]= (ts[q].second * K::raw(m1) * K(d[q])).val();
  const auto iv= inverses<K>(L);
  for(int j= 1; j <= t; ++j) {
   const int lo= d[j - 1], hi= j == t ? L : d[j];
   if(j <= TMAX) dispatch(j, std::make_integer_sequence<int, TMAX + 1>{}, [&](auto T) { seg_pow<T, K>(v, iv.data(), lo, hi, d.data(), c.data(), al.data(), j, ar); });
   else seg_pow<-1, K>(v, iv.data(), lo, hi, d.data(), c.data(), al.data(), j, ar);
  }
 }
 const auto fs= K(s).fixed();
 for(int i= o; i < n; ++i) F[i]*= fs;
 return F;
}
// exp(f / g) の先頭 n 項 (g の最も低い次数を y として、f の最も低い次数は y より大きい)。今の Library の exp_of_div と同じ漸化式。
template <class K> std::vector<K> exp_of_div(const std::vector<K>& f, const std::vector<K>& g, int n) {
 using namespace sparse_fps_internal;
 const int x= min_deg(f), y= min_deg(g);
 assert(y >= 0 && (x < 0 || x > y));
 const int N= n - 1;
 std::vector<K> F(std::max(n, 0), K());
 if(n <= 0) return F;
 F[0]= K(1);
 if(N == 0) return F;
 const K iv= K(1) / g[y];
 auto p= terms(f, 0, N + y + 1), q= terms(g, y + 1, N + y + 1);
 for(auto& [j, v]: p) v*= iv, j-= y;
 for(auto& [j, v]: q) v*= iv, j-= y;
 std::vector<K> dF(N, K()), a(N, K()), b(N, K());
 for(int i= 0; i < N; ++i) {
  for(auto&& [j, v]: p) {
   if(i < j) break;
   b[i]+= v * F[i - j];
  }
  for(auto&& [j, v]: q) {
   if(i < j) break;
   a[i]-= v * (a[i - j] + b[i + 1 - j] * K(j));
   dF[i]-= v * dF[i - j];
  }
  for(auto&& [j, v]: p) {
   if(i + 1 < j) break;
   dF[i]+= v * K(j) * F[i + 1 - j];
  }
  F[i + 1]= (dF[i]+= a[i]) / K(i + 1);
 }
 return F;
}
}  // namespace sparse_fps
