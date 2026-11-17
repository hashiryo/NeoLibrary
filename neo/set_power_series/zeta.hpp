#pragma once
#include <cassert>
#include <cstddef>
#include <type_traits>
#include <vector>
#include "neo/set_power_series/generic.hpp"
#include "neo/set_power_series/zeta_simd.hpp"
// 部分集合と上位集合のゼータ変換・メビウス変換、アダマール変換、OR と AND と XOR の畳み込み。長さは 2 冪で、その場で書き換える。
//   subset_zeta(f):    f[S] <- sum_{T ⊆ S} f[T]
//   subset_mobius(f):  subset_zeta の逆
//   supset_zeta(f):    f[S] <- sum_{S ⊆ T} f[T]
//   supset_mobius(f):  supset_zeta の逆
//   hadamard(f):       f[S] <- sum_T (-1)^{|S ∩ T|} f[T]
//   or_convolve(f, g)、and_convolve(f, g)、xor_convolve(f, g): h[S] = sum_{U op T = S} f[U] g[T]
//   subset_zeta<M>(f): f[S] <- op_{T ⊆ S} f[T] (モノイドのゼータ)
//   supset_zeta<M>(f): f[S] <- op_{S ⊆ T} f[T]
// 環の形は値の型 T の演算子を使う。ゼータは +=、メビウスとアダマールは -= も使う。xor_convolve は 2^n で割る (T(1) / 2^n が
// 0 でなければそれを掛け、0 なら T の割り算で 2^n で割る。法が偶数の ZMod では割れない)。ZMod (neo/algebra/ZMod.hpp) と
// 32 bit か 64 bit の整数は AVX2 の核で行う (整数の和は 2^32 か 2^64 を法とする)。ほかの型は今の Library の写しの汎用の道。
// モノイドの形の M は SegmentTree と同じ決まりの構造体で、using T と static T op(T a, T b) を使う (ti は使わない)。op は結合的で
// 可換なこと。標準の構造体 monoid::Sum<V>、monoid::Min<V>、monoid::Max<V> (neo/data_structure/segtree/Sum.hpp と Min.hpp) で、
// V が 32 bit か 64 bit の整数 (Sum は ZMod も) なら AVX2 の核、ほかの M は op を 1 つずつ呼ぶ汎用の道。どれも O(n 2^n)。
template <unsigned MOD> class ZMod;
namespace monoid {
template <class V> struct Sum;
template <class V> struct Min;
template <class V> struct Max;
}  // namespace monoid
namespace sps {
namespace zeta_internal {
template <class T> struct is_zmod: std::false_type {};
template <unsigned M> struct is_zmod<ZMod<M>>: std::true_type {};
template <class T> constexpr bool is_word= std::is_integral_v<T> && !std::is_same_v<T, bool> && (sizeof(T) == 4 || sizeof(T) == 8);
// 標準の構造体の見分け。kind は 0 がほか、1 が Sum、2 が Min、3 が Max。
template <class M> struct std_monoid {
 static constexpr int kind= 0;
 using V= char;  // is_word が偽になる型
};
template <class W> struct std_monoid<monoid::Sum<W>> {
 static constexpr int kind= 1;
 using V= W;
};
template <class W> struct std_monoid<monoid::Min<W>> {
 static constexpr int kind= 2;
 using V= W;
};
template <class W> struct std_monoid<monoid::Max<W>> {
 static constexpr int kind= 3;
 using V= W;
};
template <class T> inline __m256i zmod_m2() { return _mm256_set1_epi32(int(T::mod() << 1)); }
// 環の向きのある変換。INV ならメビウス。速い核で受けたら true。
template <bool UP, bool INV, class T> inline bool ring_fast(std::vector<T>& f) {
 if constexpr(is_zmod<T>::value) {
  static_assert(sizeof(T) == 4);
  if constexpr(INV) return sps_zeta::directed<4, UP>(f.data(), f.size(), sps_zeta::ZSub{zmod_m2<T>()});
  else return sps_zeta::directed<4, UP>(f.data(), f.size(), sps_zeta::ZAdd{zmod_m2<T>()});
 } else if constexpr(is_word<T>) {
  if constexpr(sizeof(T) == 4) {
   if constexpr(INV) return sps_zeta::directed<4, UP>(f.data(), f.size(), sps_zeta::Sub32{});
   else return sps_zeta::directed<4, UP>(f.data(), f.size(), sps_zeta::Add32{});
  } else {
   if constexpr(INV) return sps_zeta::directed<8, UP>(f.data(), f.size(), sps_zeta::Sub64{});
   else return sps_zeta::directed<8, UP>(f.data(), f.size(), sps_zeta::Add64{});
  }
 } else return false;
}
// モノイドの汎用の道。今の Library の汎用の道と同じく、半分ずつ再帰してから合わせる。
template <class M, bool UP> void monoid_rec(typename M::T* A, std::size_t l) {
 if(l > 128) {
  l>>= 1, monoid_rec<M, UP>(A, l), monoid_rec<M, UP>(A + l, l);
  for(std::size_t s= 0; s < l; ++s)
   if constexpr(UP) A[s + l]= M::op(A[s + l], A[s]);
   else A[s]= M::op(A[s], A[s + l]);
 } else
  for(std::size_t k= 1; k < l; k<<= 1)
   for(std::size_t i= 0; i < l; i+= 2 * k)
    for(std::size_t j= i; j < i + k; ++j)
     if constexpr(UP) A[j + k]= M::op(A[j + k], A[j]);
     else A[j]= M::op(A[j], A[j + k]);
}
template <class M, bool UP> inline bool monoid_fast(std::vector<typename M::T>& f) {
 using S= std_monoid<M>;
 using V= typename S::V;
 if constexpr(S::kind == 1 && is_zmod<V>::value) return sps_zeta::directed<4, UP>(f.data(), f.size(), sps_zeta::ZAdd{zmod_m2<V>()});
 else if constexpr(S::kind != 0 && is_word<V>) {
  static_assert(std::is_same_v<typename M::T, V>);
  constexpr bool SG= std::is_signed_v<V>, MX= S::kind == 3;
  if constexpr(S::kind == 1) {
   if constexpr(sizeof(V) == 4) return sps_zeta::directed<4, UP>(f.data(), f.size(), sps_zeta::Add32{});
   else return sps_zeta::directed<8, UP>(f.data(), f.size(), sps_zeta::Add64{});
  } else {
   if constexpr(sizeof(V) == 4) return sps_zeta::directed<4, UP>(f.data(), f.size(), sps_zeta::MinMax32<MX, SG>{});
   else return sps_zeta::directed<8, UP>(f.data(), f.size(), sps_zeta::MinMax64<MX, SG>{});
  }
 } else return false;
}
template <class M, bool UP> inline void monoid_zeta(std::vector<typename M::T>& f) {
 assert(!(f.size() & (f.size() - 1)));
 if(!monoid_fast<M, UP>(f)) monoid_rec<M, UP>(f.data(), f.size());
}
}  // namespace zeta_internal
template <class T> void subset_zeta(std::vector<T>& f) {
 if(!zeta_internal::ring_fast<true, false>(f)) sps_generic::subset_zeta(f);
}
template <class T> void subset_mobius(std::vector<T>& f) {
 if(!zeta_internal::ring_fast<true, true>(f)) sps_generic::subset_mobius(f);
}
template <class T> void supset_zeta(std::vector<T>& f) {
 if(!zeta_internal::ring_fast<false, false>(f)) sps_generic::supset_zeta(f);
}
template <class T> void supset_mobius(std::vector<T>& f) {
 if(!zeta_internal::ring_fast<false, true>(f)) sps_generic::supset_mobius(f);
}
template <class T> void hadamard(std::vector<T>& f) {
 bool done= false;
 if constexpr(zeta_internal::is_zmod<T>::value) {
  const __m256i m2= zeta_internal::zmod_m2<T>();
  done= sps_zeta::butterfly<4>(f.data(), f.size(), sps_zeta::Bfly<sps_zeta::ZAdd, sps_zeta::ZSub>{{m2}, {m2}});
 } else if constexpr(zeta_internal::is_word<T>) {
  if constexpr(sizeof(T) == 4) done= sps_zeta::butterfly<4>(f.data(), f.size(), sps_zeta::Bfly<sps_zeta::Add32, sps_zeta::Sub32>{});
  else done= sps_zeta::butterfly<8>(f.data(), f.size(), sps_zeta::Bfly<sps_zeta::Add64, sps_zeta::Sub64>{});
 }
 if(!done) sps_generic::hadamard(f);
}
template <class M> void subset_zeta(std::vector<typename M::T>& f) { zeta_internal::monoid_zeta<M, true>(f); }
template <class M> void supset_zeta(std::vector<typename M::T>& f) { zeta_internal::monoid_zeta<M, false>(f); }
template <class T> std::vector<T> or_convolve(std::vector<T> f, std::vector<T> g) {
 subset_zeta(f), subset_zeta(g);
 for(std::size_t s= 0; s < f.size(); ++s) f[s]*= g[s];
 return subset_mobius(f), f;
}
template <class T> std::vector<T> and_convolve(std::vector<T> f, std::vector<T> g) {
 supset_zeta(f), supset_zeta(g);
 for(std::size_t s= 0; s < f.size(); ++s) f[s]*= g[s];
 return supset_mobius(f), f;
}
template <class T> std::vector<T> xor_convolve(std::vector<T> f, std::vector<T> g) {
 hadamard(f), hadamard(g);
 for(std::size_t s= 0; s < f.size(); ++s) f[s]*= g[s];
 hadamard(f);
 // 2^n は T に直してから割る (今の Library の汎用の道は size_t で割るので、long long の負の値で誤る)。
 const T N= T(f.size());
 if(T iv= T(1) / N; iv == 0)
  for(std::size_t s= 0; s < f.size(); ++s) f[s]/= N;
 else
  for(std::size_t s= 0; s < f.size(); ++s) f[s]*= iv;
 return f;
}
}  // namespace sps
