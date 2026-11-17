#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>
#include "neo/set_power_series/generic.hpp"
// 引き算の無い半環の部分集合畳み込みと exp の AVX2 の核。sps の convolve.hpp と exp.hpp から呼ぶ内部の道具で、直接は使わない。
// 型は MinPlus、MaxPlus、MaxMin (neo/algebra/semiring.hpp) で、値が 32 bit か 64 bit の符号付き整数のもの。ほかの半環は
// sps_generic の 3^n の道で行う。
// h[S] = ⊕_{T ⊆ S} f[T] ⊗ g[S \ T] を、添字の下の B bit (値が 4 byte なら 3、8 byte なら 2) を 1 本の __m256i のレーンに置いて
// 求める。上の bit は、T の上の部分 t と、S の上の部分 t | v (v は t の補集合の部分集合) で回す。t ごとに、f の塊をレーンの
// xor で並べ替えた L 本 (レーン s に f[t, s ^ u] を置き、u ⊆ s でなければ和の単位元になる値にする) を作っておく。v ごとに
// g の塊の L 個を 1 つずつ配って掛け (足し算か min)、h の塊に畳む (min か max)。(t, v) の組は 3^{n-B} 個で、1 組で 3^B 項になる。
// MinPlus と MaxPlus の積は INF で止める足し算だが、核では止めずに足し、最後に INF の側へ外れた値を INF に直す (semiring.hpp
// の行列の積の速い道と同じ)。値がすべて 0 以上 INF 以下 (MaxPlus は -INF 以上 0 以下) か、すべて ±INF か絶対値が INF / 4 以下
// なら、厳密な積と同じ値になる。どちらでもなければ sps_generic の道で行う。MaxMin は min と max だけなので、いつも核で行う。
// exp は最大の要素を含む塊を選ぶ形で、h[2^i + S] = ⊕_{T ⊆ S} f[2^i + T] ⊗ h[S \ T] の部分集合畳み込みを i = 0, ..., n-1 と
// 続ける (sps_generic の semiring_exp は最小の要素を含む塊を選ぶ。どちらも集合の分割を 1 回ずつ数える)。
template <class T> struct MinPlus;
template <class T> struct MaxPlus;
template <class T> struct MaxMin;
namespace sps_semi {
using V= __m256i;
#define SS_AI [[gnu::always_inline]] inline
template <class X> constexpr bool word= std::is_integral_v<X> && std::is_signed_v<X> && (sizeof(X) == 4 || sizeof(X) == 8);
// kind は 0 が和 min と積の足し算 (MinPlus)、1 が和 max と積の足し算 (MaxPlus)、2 が和 max と積 min (MaxMin)。
template <class T> struct traits {
 static constexpr bool ok= false;
 static constexpr int kind= -1;
 using W= int;
};
template <class X> struct traits<MinPlus<X>> {
 static constexpr bool ok= word<X>;
 static constexpr int kind= 0;
 using W= X;
};
template <class X> struct traits<MaxPlus<X>> {
 static constexpr bool ok= word<X>;
 static constexpr int kind= 1;
 using W= X;
};
template <class X> struct traits<MaxMin<X>> {
 static constexpr bool ok= word<X>;
 static constexpr int kind= 2;
 using W= X;
};
template <class W> SS_AI V ld(const W* p) { return _mm256_loadu_si256((const V*)p); }
template <class W> SS_AI void st(W* p, V x) { _mm256_storeu_si256((V*)p, x); }
template <int E, class W> SS_AI V vset(W x) {
 if constexpr(E == 4) return _mm256_set1_epi32(int(x));
 else return _mm256_set1_epi64x((long long)x);
}
template <int E> SS_AI V vadd(V a, V b) {
 if constexpr(E == 4) return _mm256_add_epi32(a, b);
 else return _mm256_add_epi64(a, b);
}
// 64 bit の min と max は AVX2 に命令が無いので、cmpgt で比べて混ぜる。
template <int E> SS_AI V vmin(V a, V b) {
 if constexpr(E == 4) return _mm256_min_epi32(a, b);
 else return _mm256_blendv_epi8(a, b, _mm256_cmpgt_epi64(a, b));
}
template <int E> SS_AI V vmax(V a, V b) {
 if constexpr(E == 4) return _mm256_max_epi32(a, b);
 else return _mm256_blendv_epi8(a, b, _mm256_cmpgt_epi64(b, a));
}
template <int E, int K> SS_AI V oplus(V a, V b) {
 if constexpr(K == 0) return vmin<E>(a, b);
 else return vmax<E>(a, b);
}
template <int E, int K> SS_AI V otimes(V a, V b) {
 if constexpr(K == 2) return vmin<E>(a, b);
 else return vadd<E>(a, b);
}
// レーン s (s & u == u) の 32 bit のレーンの印。
constexpr int lane_mask(int E, int u) {
 int r= 0;
 for(int s= 0; s < 32 / E; ++s)
  if((s & u) == u) r|= E == 4 ? 1 << s : 3 << (2 * s);
 return r;
}
// レーン s に F のレーン s ^ u を置き、u ⊆ s でないレーンは H にする。
template <int E, int u> SS_AI V pxor(V F, V H) {
 if constexpr(E == 4) {
  const V idx= _mm256_setr_epi32(0 ^ u, 1 ^ u, 2 ^ u, 3 ^ u, 4 ^ u, 5 ^ u, 6 ^ u, 7 ^ u);
  return _mm256_blend_epi32(H, _mm256_permutevar8x32_epi32(F, idx), lane_mask(4, u));
 } else {
  constexpr int p= (0 ^ u) | (1 ^ u) << 2 | (2 ^ u) << 4 | (3 ^ u) << 6;
  return _mm256_blend_epi32(H, _mm256_permute4x64_epi64(F, p), lane_mask(8, u));
 }
}
// h[S] = ⊕_{T ⊆ S} f[T] ⊗' g[S \ T] (n >= B、⊗' は止めない足し算か min)。hv は和の単位元になる値 (MinPlus は INF、MaxPlus は
// -INF、MaxMin は最小の値)。h の値は直さずに返す。
template <int E, int K, class W> void conv_core(const W* f, const W* g, W* h, int n, W hv) {
 constexpr int L= 32 / E, B= E == 4 ? 3 : 2;
 const std::size_t M= std::size_t(1) << (n - B), mask= M - 1;
 const V ID= vset<E>(K == 0 ? std::numeric_limits<W>::max() : std::numeric_limits<W>::lowest()), H= vset<E>(hv);
 for(std::size_t i= 0; i < (M << B); i+= L) st(h + i, ID);
 for(std::size_t t= 0; t < M; ++t) {
  const V F= ld(f + t * L);
  const std::size_t c= mask ^ t;
  if constexpr(E == 4) {
   const V f0= pxor<4, 0>(F, H), f1= pxor<4, 1>(F, H), f2= pxor<4, 2>(F, H), f3= pxor<4, 3>(F, H);
   const V f4= pxor<4, 4>(F, H), f5= pxor<4, 5>(F, H), f6= pxor<4, 6>(F, H), f7= pxor<4, 7>(F, H);
   for(std::size_t v= c;; v= (v - 1) & c) {
    W* hp= h + (t | v) * L;
    const W* gp= g + v * L;
    const V p0= otimes<4, K>(f0, vset<4>(gp[0])), p1= otimes<4, K>(f1, vset<4>(gp[1]));
    const V p2= otimes<4, K>(f2, vset<4>(gp[2])), p3= otimes<4, K>(f3, vset<4>(gp[3]));
    const V p4= otimes<4, K>(f4, vset<4>(gp[4])), p5= otimes<4, K>(f5, vset<4>(gp[5]));
    const V p6= otimes<4, K>(f6, vset<4>(gp[6])), p7= otimes<4, K>(f7, vset<4>(gp[7]));
    const V a= oplus<4, K>(oplus<4, K>(oplus<4, K>(p0, p1), oplus<4, K>(p2, p3)), oplus<4, K>(oplus<4, K>(p4, p5), oplus<4, K>(p6, p7)));
    st(hp, oplus<4, K>(ld(hp), a));
    if(!v) break;
   }
  } else {
   const V f0= pxor<8, 0>(F, H), f1= pxor<8, 1>(F, H), f2= pxor<8, 2>(F, H), f3= pxor<8, 3>(F, H);
   for(std::size_t v= c;; v= (v - 1) & c) {
    W* hp= h + (t | v) * L;
    const W* gp= g + v * L;
    const V p0= otimes<8, K>(f0, vset<8>(gp[0])), p1= otimes<8, K>(f1, vset<8>(gp[1]));
    const V p2= otimes<8, K>(f2, vset<8>(gp[2])), p3= otimes<8, K>(f3, vset<8>(gp[3]));
    st(hp, oplus<8, K>(ld(hp), oplus<8, K>(oplus<8, K>(p0, p1), oplus<8, K>(p2, p3))));
    if(!v) break;
   }
  }
 }
}
// 厳密な積で 3^n (小さい n と、exp の始めの段)。h[S] = ⊕_{T ⊆ S} f[T] ⊗ g[S \ T]。
template <class T> void conv_scalar(const T* f, const T* g, T* h, std::size_t N) {
 for(std::size_t s= 0; s < N; ++s) {
  T acc= f[s] * g[0];
  for(std::size_t t= (s - 1) & s; t > 0; t= (t - 1) & s) acc+= f[t] * g[s ^ t];
  if(s) acc+= f[0] * g[s];
  h[s]= acc;
 }
}
template <class T> constexpr int lanes_bits() { return sizeof(typename traits<T>::W) == 4 ? 3 : 2; }
// 速い道で掛けてよいか。0 はよくない、1 は値が片側 (MinPlus は 0 以上 INF 以下)、2 は ±INF か絶対値が INF / 4 以下。MaxMin は 1。
template <class T> inline int mode(const T* p, const T* q, std::size_t n) {
 if constexpr(traits<T>::kind == 2) return 1;
 else {
  if(T::fast_one_sided(p, q, n)) return 1;
  if(T::fast_bounded(p, n) && T::fast_bounded(q, n)) return 2;
  return 0;
 }
}
// 速い道の部分集合畳み込み (n >= B)。h の値は直してから返す。
template <class T> inline void conv_fast(const T* f, const T* g, T* h, int n, int md) {
 using W= typename traits<T>::W;
 constexpr int K= traits<T>::kind, E= sizeof(W);
 static_assert(sizeof(T) == sizeof(W) && std::is_standard_layout_v<T>);
 W hv;
 if constexpr(K == 0) hv= T::INF;
 else if constexpr(K == 1) hv= -T::INF;
 else hv= std::numeric_limits<W>::lowest();
 conv_core<E, K>(reinterpret_cast<const W*>(f), reinterpret_cast<const W*>(g), reinterpret_cast<W*>(h), n, hv);
 if constexpr(K != 2) T::fast_fix_all(h, std::size_t(1) << n, md == 2);
}
template <class T> std::vector<T> convolve(const std::vector<T>& f, const std::vector<T>& g) {
 if constexpr(traits<T>::ok) {
  const std::size_t N= f.size();
  const int n= __builtin_ctzll(N);
  if(n >= lanes_bits<T>())
   if(const int md= mode(f.data(), g.data(), N)) {
    std::vector<T> h(N);
    conv_fast(f.data(), g.data(), h.data(), n, md);
    return h;
   }
 }
 return sps_generic::semiring_convolve(f, g);
}
// 速い道で掛けられない段が出たら、全体を sps_generic の道でやり直す (値が INF に近いときの止め方を、今までの形と同じにするため)。
template <class T> std::vector<T> exp(const std::vector<T>& f) {
 if constexpr(traits<T>::ok) {
  constexpr int B= lanes_bits<T>();
  const std::size_t N= f.size();
  const int n= __builtin_ctzll(N);
  if(n > B) {
   std::vector<T> h(N);
   h[0]= T::one();
   bool ok= true;
   for(int i= 0; i < n && ok; ++i) {
    const std::size_t L= std::size_t(1) << i;
    if(i < B) conv_scalar(f.data() + L, h.data(), h.data() + L, L);
    else if(const int md= mode(f.data() + L, h.data(), L)) conv_fast(f.data() + L, h.data(), h.data() + L, i, md);
    else ok= false;
   }
   if(ok) return h;
  }
 }
 return sps_generic::semiring_exp(f, T::one());
}
#undef SS_AI
}  // namespace sps_semi
