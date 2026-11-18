#pragma once
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <utility>
#include <vector>
#include "neo/set_power_series/ranked_xor.hpp"
// 有向グラフの部分集合ごとの数え上げ (graph.hpp の strongly_connected_digraph と acyclic_digraph) の内部の核。直接は使わない。
// どちらも、頂点の集合 S ごとの値を
//   h[S] = c[S] - sum_{T ⊊ S} h[T] X(S \ T) F(T → S \ T)   (h[∅] = h0)
// で S の小さい順に求める 3^n の和に帰着する。F(T → U) は T から U へ向かう辺の (1 + w_e) の積、X(U) は U の中の辺 (自己ループ
// を含む) の (1 + w_e) の積 (WX が偽なら 1)。強連結は WX、c[S] = X(S)、h0 = 0 で g を求めて f = -log(1 - g)、DAG は辺を逆向きに
// したグラフで WX を偽、c = 0、h0 = 1 として b を求め、a[S] = (-1)^{|S|} b[S] (b[S] = -sum_{U⊊S} b[U] F(S\U → U) のため)。
// 汎用の道は、S ごとに S の部分集合 T の F(T → S \ T) を、S の bit を 1 つずつ決めて両側の表の値を掛ける形で 2^{|S|} 個作る。
// 法が奇数の ZMod (p < 2^30) で n >= 4 なら、添字の下の 3 bit を AVX2 のレーンに置く (procon-judge の luogu-P11714 の scc_simd)。
// 組 (T の上の部分 T_h、U の上の部分 U_h) ごとに、T_h の側と U_h の側でレーンの向きに作った 8 本の列の内積を取り、組ごとの値
// F(T_h → U_h) を掛ける (X(U_h) は U_h の側の列に掛けておく)。同じ塊の中の項 (U_h = ∅) は塊ごとにレーンの小さい順に引く。
// 重みが無いときは辺の因子がどれも 2 なので、組ごとの値を 2 の冪の肩 (辺の数) の足し算で求めて表から引く。
namespace sps_digraph_internal {
using u32= std::uint32_t;
using u64= std::uint64_t;
// 頂点の組ごとの辺の (1 + w) の積 Q[a * n + b] (a → b、辺の無い組は 1)
template <class T> std::vector<T> edge_factors(int n, const std::vector<std::pair<int, int>>& es, const std::vector<T>& w, bool rev) {
 std::vector<T> Q(std::size_t(n) * n, T(1));
 for(std::size_t i= 0; i < es.size(); ++i) {
  auto [a, b]= es[i];
  assert(0 <= a && a < n && 0 <= b && b < n);
  if(rev) std::swap(a, b);
  Q[std::size_t(a) * n + b]*= w.empty() ? T(2) : T(1) + w[i];
 }
 return Q;
}
// outq[v][X] = prod_{b ∈ X} Q[v][b]、inq[v][X] = prod_{a ∈ X} Q[a][v] (X は全体の部分集合)
template <class T> void vertex_tables(int n, const std::vector<T>& Q, std::vector<T>& outq, std::vector<T>& inq) {
 const std::size_t N= std::size_t(1) << n;
 outq.assign(std::size_t(n) * N, T(1)), inq.assign(std::size_t(n) * N, T(1));
 for(int v= 0; v < n; ++v)
  for(std::size_t X= 1; X < N; ++X) {
   const int u= __builtin_ctzll(X);
   outq[v * N + X]= outq[v * N + (X & (X - 1))] * Q[std::size_t(v) * n + u];
   inq[v * N + X]= inq[v * N + (X & (X - 1))] * Q[std::size_t(u) * n + v];
  }
}
// X(S) (S の中の辺と自己ループの積)
template <class T> std::vector<T> inner_table(int n, const std::vector<T>& Q, const std::vector<T>& outq, const std::vector<T>& inq) {
 const std::size_t N= std::size_t(1) << n;
 std::vector<T> X(N, T(1));
 for(std::size_t S= 1; S < N; ++S) {
  const int v= __builtin_ctzll(S);
  const std::size_t R= S & (S - 1);
  X[S]= X[R] * Q[std::size_t(v) * n + v] * outq[v * N + R] * inq[v * N + R];
 }
 return X;
}
// 汎用の道 (どの可換環でも)
template <class T> std::vector<T> twisted_generic(int n, const std::vector<T>& Q, bool WX, T h0) {
 const std::size_t N= std::size_t(1) << n;
 std::vector<T> outq, inq;
 vertex_tables(n, Q, outq, inq);
 std::vector<T> X= WX ? inner_table(n, Q, outq, inq) : std::vector<T>(N, T(1));
 std::vector<T> h(N, T()), arr(N);
 std::vector<u32> tm(N);
 h[0]= h0;
 int bits[32];
 for(std::size_t S= 1; S < N; ++S) {
  int k= 0;
  for(std::size_t t= S; t; t&= t - 1) bits[k++]= __builtin_ctzll(t);
  // arr[t] = F(T_t → (先に決めた bit のうち T_t でないもの))、tm[t] = T_t (頂点の集合)。bit を 1 つずつ決めて倍に伸ばす
  arr[0]= T(1), tm[0]= 0;
  u32 pre= 0;
  for(int i= 0; i < k; ++i) {
   const int v= bits[i];
   const std::size_t half= std::size_t(1) << i;
   for(std::size_t t= 0; t < half; ++t) {
    const u32 Tt= tm[t], Ut= pre ^ Tt;
    arr[t | half]= arr[t] * outq[v * N + Ut], tm[t | half]= Tt | (1u << v);
    arr[t]*= inq[v * N + Tt];
   }
   pre|= 1u << v;
  }
  T acc= WX ? X[S] : T();
  for(std::size_t t= 0; t + 1 < (std::size_t(1) << k); ++t) acc-= h[tm[t]] * X[S ^ tm[t]] * arr[t];
  h[S]= acc;
 }
 return h;
}
// ---- 法が奇数の ZMod の AVX2 の道 ----
using V= __m256i;
struct MM {
 u32 P, NP, R1, R2;
 explicit MM(u32 p): P(p), NP(0u - sps_xor::inv32(p)), R1(u32((u64(1) << 32) % p)), R2(u32(u64(R1) * R1 % p)) {}
 u32 mul(u32 a, u32 b) const {  // a b R^{-1} を [0, p) に
  const u64 t= u64(a) * b;
  const u64 r= (t + u64(u32(t) * NP) * P) >> 32;
  return u32(r >= P ? r - P : r);
 }
 u32 to(u32 x) const { return mul(x % P, R2); }
 u32 from(u32 x) const { return mul(x, 1); }
};
// Q は Montgomery の表現 (p 未満)。返す h も Montgomery の表現。cnt は、辺の因子がどれも 2 の冪 (重みが無い) のときの、頂点の組ごとの
// 辺の数 (Q[a * n + b] = 2^{cnt[a * n + b]})。あれば、組ごとの値 F(T_h → U_h) を 2 の冪の肩の足し算で求める。
inline std::vector<u32> twisted_simd(int n, const std::vector<u32>& Q, bool WX, u32 h0, const sps_xor::MK& k0, const std::vector<int>* cnt= nullptr) {
 const sps_xor::MK k= k0;
 const MM mm(k.P);
 const int H= n - 3;
 const std::size_t N= std::size_t(1) << n, NH= std::size_t(1) << H;
 const u32 one= mm.R1;
 auto q= [&](int a, int b) { return Q[std::size_t(a) * n + b]; };
 // 下の 3 bit の中の表: Fl[t][d] = F(t → d)、Xl[d] = X(d)
 u32 Fl[8][8], Xl[8], E[8][8];
 for(int t= 0; t < 8; ++t)
  for(int d= 0; d < 8; ++d) {
   u32 r= one;
   if(!(t & d))
    for(int a= 0; a < 3; ++a)
     if(t >> a & 1)
      for(int b= 0; b < 3; ++b)
       if(d >> b & 1) r= mm.mul(r, q(a, b));
   Fl[t][d]= r;
  }
 for(int d= 0; d < 8; ++d) {
  u32 r= one;
  for(int a= 0; a < 3; ++a)
   if(d >> a & 1)
    for(int b= 0; b < 3; ++b)
     if(d >> b & 1) r= mm.mul(r, q(a, b));
  Xl[d]= WX ? r : one;
 }
 for(int t= 0; t < 8; ++t)
  for(int d= 0; d < 8; ++d) E[t][d]= (t & d) ? 0 : mm.mul(Xl[d], Fl[t][d]);
 // 上の頂点 (番号 j + 3) と下の集合の表: Bt[T_h][d] = F(T_h → d)、Ct[t][U_h] = F(t → U_h)
 std::vector<u32> Bt(NH * 8, one), Ct(8 * NH, one);
 for(std::size_t th= 1; th < NH; ++th) {
  const int v= __builtin_ctzll(th) + 3;
  for(int d= 0; d < 8; ++d) {
   u32 r= Bt[(th & (th - 1)) * 8 + d];
   for(int b= 0; b < 3; ++b)
    if(d >> b & 1) r= mm.mul(r, q(v, b));
   Bt[th * 8 + d]= r;
  }
  for(int t= 0; t < 8; ++t) {
   u32 r= Ct[t * NH + (th & (th - 1))];
   for(int a= 0; a < 3; ++a)
    if(t >> a & 1) r= mm.mul(r, q(a, v));
   Ct[t * NH + th]= r;
  }
 }
 // 上の頂点どうしの表: oh[j][X_h] = prod_{b ∈ X_h} q(j+3, b+3)、ih[j][X_h] = prod_{a ∈ X_h} q(a+3, j+3)
 std::vector<u32> oh(std::size_t(H) * NH, one), ih(std::size_t(H) * NH, one);
 for(int j= 0; j < H; ++j)
  for(std::size_t X= 1; X < NH; ++X) {
   const int u= __builtin_ctzll(X);
   oh[j * NH + X]= mm.mul(oh[j * NH + (X & (X - 1))], q(j + 3, u + 3));
   ih[j * NH + X]= mm.mul(ih[j * NH + (X & (X - 1))], q(u + 3, j + 3));
  }
 // X(U_h) と、全体の X(S) (c[S] に使う)
 std::vector<u32> XH(NH, one), XS;
 if(WX) {
  for(std::size_t X= 1; X < NH; ++X) {
   const int v= __builtin_ctzll(X);
   const std::size_t R= X & (X - 1);
   XH[X]= mm.mul(mm.mul(mm.mul(XH[R], q(v + 3, v + 3)), oh[v * NH + R]), ih[v * NH + R]);
  }
  XS.resize(N);
  for(std::size_t S= 0; S < N; ++S) {
   const std::size_t sh= S >> 3, sl= S & 7;
   // X(S) = X(S_h) X(S_l) F(S_h → S_l) F(S_l → S_h)
   XS[S]= mm.mul(mm.mul(XH[sh], Xl[sl]), mm.mul(Bt[sh * 8 + sl], Ct[sl * NH + sh]));
  }
 }
 // Cd[U_h][d] (8 レーン): [d ⊆ S_l] F(S_l^d → U_h) X(U_h ↔ d) X(U_h)。X(U_h) は組ごとの値から外してここで掛ける
 V* Cd= static_cast<V*>(std::aligned_alloc(32, NH * 8 * sizeof(V)));
 V* Gd= static_cast<V*>(std::aligned_alloc(32, NH * 8 * sizeof(V)));
 for(std::size_t uh= 0; uh < NH; ++uh)
  for(int d= 0; d < 8; ++d) {
   alignas(32) u32 lane[8];
   const u32 x= WX ? mm.mul(mm.mul(Bt[uh * 8 + d], Ct[d * NH + uh]), XH[uh]) : one;
   for(int sl= 0; sl < 8; ++sl) lane[sl]= (d & ~sl) ? 0 : mm.mul(Ct[(sl ^ d) * NH + uh], x);
   Cd[uh * 8 + d]= _mm256_load_si256((const V*)lane);
  }
 std::vector<u32> h(N, 0), arr(NH), tmh(NH);
 // 重みが無いとき: co[j][X] = sum_{b ∈ X} cnt(j+3 → b+3)、ci[j][X] = sum_{a ∈ X} cnt(a+3 → j+3)、ex[t] は arr[t] の 2 の冪の肩
 const bool P2W= cnt != nullptr;
 std::vector<int> co, ci, ex;
 std::vector<u32> pw;
 if(P2W) {
  co.assign(std::size_t(H) * NH, 0), ci.assign(std::size_t(H) * NH, 0), ex.assign(NH, 0);
  int m= 0;
  for(int x: *cnt) m+= x;
  pw.resize(m + 1);
  pw[0]= one;
  const u32 two= mm.to(2);
  for(int i= 1; i <= m; ++i) pw[i]= mm.mul(pw[i - 1], two);
  for(int j= 0; j < H; ++j)
   for(std::size_t X= 1; X < NH; ++X) {
    const int u= __builtin_ctzll(X);
    co[j * NH + X]= co[j * NH + (X & (X - 1))] + (*cnt)[std::size_t(j + 3) * n + u + 3];
    ci[j * NH + X]= ci[j * NH + (X & (X - 1))] + (*cnt)[std::size_t(u + 3) * n + j + 3];
   }
 }
 int bits[32];
 const V P2= k.p2, P1= k.p;
 for(std::size_t sh= 0; sh < NH; ++sh) {
  alignas(32) u32 acc[8];
  for(int sl= 0; sl < 8; ++sl) acc[sl]= WX ? XS[(sh << 3) | sl] : 0;
  V a= _mm256_load_si256((const V*)acc);
  if(sh) {
   // arr[t] = F(T_h → U_h) (T_h = tmh[t]、U_h = S_h \ T_h)。S_h の bit を 1 つずつ決めて倍に伸ばす
   int kb= 0;
   for(std::size_t t= sh; t; t&= t - 1) bits[kb++]= __builtin_ctzll(t);
   arr[0]= one, tmh[0]= 0;
   if(P2W) ex[0]= 0;
   u32 pre= 0;
   for(int i= 0; i < kb; ++i) {
    const int j= bits[i];
    const std::size_t half= std::size_t(1) << i;
    for(std::size_t t= 0; t < half; ++t) {
     const u32 Tt= tmh[t], Ut= pre ^ Tt;
     tmh[t | half]= Tt | (1u << j);
     if(P2W) ex[t | half]= ex[t] + co[j * NH + Ut], ex[t]+= ci[j * NH + Tt];
     else arr[t | half]= mm.mul(arr[t], oh[j * NH + Ut]), arr[t]= mm.mul(arr[t], ih[j * NH + Tt]);
    }
    pre|= 1u << j;
   }
   for(std::size_t t= 0; t + 1 < (std::size_t(1) << kb); ++t) {
    const std::size_t th= tmh[t], uh= sh ^ th;
    const u32 A= P2W ? pw[ex[t]] : arr[t];
    const V* gp= Gd + th * 8;
    const V* cp= Cd + uh * 8;
    V se= _mm256_setzero_si256(), so= se;
    for(int d= 0; d < 8; ++d) {
     const V x= gp[d], y= cp[d];
     se= _mm256_add_epi64(se, _mm256_mul_epu32(x, y)), so= _mm256_add_epi64(so, _mm256_mul_epu32(_mm256_srli_epi64(x, 32), _mm256_srli_epi64(y, 32)));
    }
    V r= sps_xor::norm1(sps_xor::redc2(se, so, k), k);  // 8 項の和は 2^63 未満なので、還元は 3p 未満で、norm1 のあとは 2p 未満
    const V wv= _mm256_set1_epi32(int(A));
    r= sps_xor::redc2(_mm256_mul_epu32(r, wv), _mm256_mul_epu32(_mm256_srli_epi64(r, 32), wv), k);  // 2p 未満
    a= _mm256_sub_epi32(_mm256_add_epi32(a, P2), r);
    a= _mm256_min_epu32(a, _mm256_sub_epi32(a, P2));
    a= _mm256_min_epu32(a, _mm256_sub_epi32(a, P1));
   }
  }
  _mm256_store_si256((V*)acc, a);
  // 同じ塊の中 (U_h = ∅、T_l ⊊ S_l)。T = ∅ の項は h0
  u32* hs= &h[sh << 3];
  for(int sl= 0; sl < 8; ++sl) {
   if(sh == 0 && sl == 0) {
    hs[0]= h0;
    continue;
   }
   u64 v= acc[sl] >= k.P ? acc[sl] - k.P : acc[sl];
   for(int tl= (sl - 1) & sl;; tl= (tl - 1) & sl) {
    const int ul= sl ^ tl;
    v+= mm.mul(k.P - hs[tl], mm.mul(Bt[sh * 8 + ul], E[tl][ul]));
    if(!tl) break;
   }
   hs[sl]= u32(v % k.P);
  }
  for(int d= 0; d < 8; ++d) {
   alignas(32) u32 lane[8];
   for(int sl= 0; sl < 8; ++sl) lane[sl]= (d & ~sl) ? 0 : mm.mul(hs[sl ^ d], mm.mul(Bt[sh * 8 + d], E[sl ^ d][d]));
   Gd[sh * 8 + d]= _mm256_load_si256((const V*)lane);
  }
 }
 std::free(Cd), std::free(Gd);
 return h;
}
template <class T> std::vector<T> twisted(int n, const std::vector<std::pair<int, int>>& es, const std::vector<T>& w, bool rev, bool WX, T h0) {
 assert(w.empty() || w.size() == es.size());
 const std::vector<T> Q= edge_factors(n, es, w, rev);
 if constexpr(sps_xor::is_zmod<T>::value) {
  if((T::mod() & 1) && T::mod() >= 3 && n >= 4) {
   const sps_xor::MK k(T::mod());
   const MM mm(T::mod());
   std::vector<u32> Qm(Q.size());
   for(std::size_t i= 0; i < Q.size(); ++i) Qm[i]= mm.to(Q[i].val());
   std::vector<int> cnt;  // 重みが無いときの頂点の組ごとの辺の数
   if(w.empty()) {
    cnt.assign(std::size_t(n) * n, 0);
    for(auto [a, b]: es) {
     if(rev) std::swap(a, b);
     ++cnt[std::size_t(a) * n + b];
    }
   }
   const std::vector<u32> hm= twisted_simd(n, Qm, WX, mm.to(h0.val()), k, w.empty() ? &cnt : nullptr);
   std::vector<T> h(hm.size());
   for(std::size_t i= 0; i < hm.size(); ++i) h[i]= T::raw(mm.from(hm[i]));
   return h;
  }
 }
 return twisted_generic(n, Q, WX, h0);
}
}  // namespace sps_digraph_internal
