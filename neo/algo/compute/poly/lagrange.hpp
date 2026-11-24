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
#include "neo/base/core/convolve/naive.hpp"
// 多項式の補間。K は体で、ZMod なら法は素数であること。
// - interpolate_iota(y, x): y_i = f(i) (0 ≤ i < n) を満たす n-1 次以下の f の f(x)。O(n)。
// - interpolate(x, y): 点 (x_i, y_i) を通る n-1 次以下の多項式の係数 (x_i は相異なる)。O(n^2)。ZMod は冪の和と畳み込みの形。
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
using convolve_internal::hvec;
// x y 2^{-32} mod M を 8 個ずつ (M は奇数、x と y が 2M 未満なら結果も 2M 未満)。
[[gnu::always_inline]] inline __m256i mont8(__m256i x, __m256i y, __m256i mv, __m256i mi) {
 const __m256i te= _mm256_mul_epu32(x, y), to= _mm256_mul_epu32(_mm256_srli_epi64(x, 32), _mm256_srli_epi64(y, 32));
 const __m256i ue= _mm256_add_epi64(te, _mm256_mul_epu32(_mm256_mul_epu32(te, mi), mv));
 const __m256i uo= _mm256_add_epi64(to, _mm256_mul_epu32(_mm256_mul_epu32(to, mi), mv));
 return _mm256_blend_epi32(_mm256_srli_epi64(ue, 32), uo, 0b10101010);
}
// 法 M (奇数) の Montgomery の定数。R1 = 2^32 mod M (Montgomery の 1)、R2 = R1^2、R3 = R1^3。
template <class Z> struct Mont {
 u32 M, R1, R2, R3;
 __m256i mv, mi;
 Mont(): M(Z::mod()) {
  u32 x= M;
  for(int i= 0; i < 5; ++i) x*= 2 - M * x;
  mv= _mm256_set1_epi32(M), mi= _mm256_set1_epi32(-x);
  R1= Z(u64(1) << 32).val(), R2= (Z::raw(R1) * Z::raw(R1)).val(), R3= (Z::raw(R2) * Z::raw(R1)).val();
 }
 [[gnu::always_inline]] __m256i mul(__m256i x, __m256i y) const { return mont8(x, y, mv, mi); }
 [[gnu::always_inline]] __m256i red(__m256i x) const { return _mm256_min_epu32(x, _mm256_sub_epi32(x, mv)); }
 [[gnu::always_inline]] __m256i add(__m256i x, __m256i y) const { return red(_mm256_add_epi32(x, y)); }  // M 未満どうし
 [[gnu::always_inline]] __m256i neg(__m256i x) const { return red(_mm256_sub_epi32(mv, x)); }            // M 未満
};
// gen(j) が返す 8 個 (1 以上 M 未満) の逆元に sigma を掛けたものを、後ろ向きに back(j, sigma/値) で渡す。全部の値の積を返す。
// レーン b は添字 ≡ b (mod 8) の値の累積積を受け持つ (8 本の鎖)。レーンの総積は (Π 値) 2^{-32 steps} なので、全体の積は
// 8 本の積に 2^{256 steps} を掛けたもの。Montgomery の 2^{-32} は前向きと後ろ向きで打ち消し合う。
template <class Z, class Gen, class Back> Z batch_inv(int steps, u32* W, u32 sigma, const Mont<Z>& mt, Gen gen, Back back) {
 if(steps == 0) return Z(1);
 __m256i acc= _mm256_set1_epi32(1);
 for(int j= 0; j < steps; ++j) acc= mt.mul(acc, gen(j)), _mm256_storeu_si256((__m256i*)(W + 8 * j), acc);
 alignas(32) u32 c[8];
 _mm256_store_si256((__m256i*)c, acc);
 Z prod(1);
 for(int b= 0; b < 8; ++b) {
  const Z v(c[b]);
  prod*= v, c[b]= (v.inv() * Z::raw(sigma)).val();
 }
 prod*= Z::raw(mt.R1).pow(u64(8) * steps);
 __m256i cur= _mm256_load_si256((const __m256i*)c);
 for(int j= steps - 1; j >= 0; --j) {
  const __m256i pre= j ? _mm256_loadu_si256((const __m256i*)(W + 8 * (j - 1))) : _mm256_set1_epi32(1);
  back(j, mt.red(mt.mul(cur, pre)));
  cur= mt.mul(cur, gen(j));
 }
 return prod;
}
// 比 rho[i] (Montgomery の形、i < N の先は R1 で埋める) の累積積 P_i = Π_{j≤i} rho_j を区間 4 つの鎖で作り、
// A = Σ_{i<N} val[i] P_i と B = Σ_{i<N} P_i を返す (どちらも本当の値)。val は M 未満 (N 以降は読むが足さない)。
// ベクトルの中の累積積は 3 段の並べ替えと積で作り、前のベクトルからの繰り上げを掛ける。区間ごとの和は区間の始まりからの
// 累積積で作り、最後に前の区間の総積を掛けて足す。累積積は配列に置かない。
template <class Z> void scan_acc(const u32* rho, const u32* val, int N, const Mont<Z>& mt, Z& A, Z& B) {
 constexpr int SG= 4;
 const int nv= (N + 7) / 8, per= (nv + SG - 1) / SG;
 const __m256i one= _mm256_set1_epi32(mt.R1), lo32= _mm256_set1_epi64x(0xffffffffll), lane= _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
 const __m256i nn= _mm256_set1_epi32(N), i1= _mm256_setr_epi32(0, 0, 1, 2, 3, 4, 5, 6), i2= _mm256_setr_epi32(0, 0, 0, 1, 2, 3, 4, 5);
 const __m256i i4= _mm256_setr_epi32(0, 0, 0, 0, 0, 1, 2, 3), i7= _mm256_set1_epi32(7), u1= _mm256_set1_epi32(1);
 __m256i carry[SG], ae[SG], ao[SG], be[SG], bo[SG];
 for(int g= 0; g < SG; ++g) carry[g]= one, ae[g]= ao[g]= be[g]= bo[g]= _mm256_setzero_si256();
 auto step= [&](int g, int j) {
  __m256i x= _mm256_loadu_si256((const __m256i*)(rho + 8 * j));
  x= mt.mul(x, _mm256_blend_epi32(_mm256_permutevar8x32_epi32(x, i1), one, 0b00000001));
  x= mt.mul(x, _mm256_blend_epi32(_mm256_permutevar8x32_epi32(x, i2), one, 0b00000011));
  x= mt.mul(x, _mm256_blend_epi32(_mm256_permutevar8x32_epi32(x, i4), one, 0b00001111));
  const __m256i p= mt.mul(x, carry[g]);
  carry[g]= _mm256_permutevar8x32_epi32(p, i7);
  const __m256i ok= _mm256_cmpgt_epi32(nn, _mm256_add_epi32(lane, _mm256_set1_epi32(8 * j)));
  const __m256i v= _mm256_and_si256(mt.mul(_mm256_loadu_si256((const __m256i*)(val + 8 * j)), p), ok);
  const __m256i q= _mm256_and_si256(mt.mul(p, u1), ok);
  ae[g]= _mm256_add_epi64(ae[g], _mm256_and_si256(v, lo32)), ao[g]= _mm256_add_epi64(ao[g], _mm256_srli_epi64(v, 32));
  be[g]= _mm256_add_epi64(be[g], _mm256_and_si256(q, lo32)), bo[g]= _mm256_add_epi64(bo[g], _mm256_srli_epi64(q, 32));
 };
 for(int t= 0; t < per; ++t)
  for(int g= 0; g < SG; ++g)
   if(g * per + t < nv) step(g, g * per + t);
 Z start(1), sa(0), sb(0);
 const Z iR= Z::raw(mt.R1).inv();
 for(int g= 0; g < SG && g * per < nv; ++g) {
  alignas(32) u64 x[4], y[4];
  _mm256_store_si256((__m256i*)x, _mm256_add_epi64(ae[g], ao[g]));
  _mm256_store_si256((__m256i*)y, _mm256_add_epi64(be[g], bo[g]));
  sa+= start * (Z(x[0]) + Z(x[1]) + Z(x[2]) + Z(x[3])), sb+= start * (Z(y[0]) + Z(y[1]) + Z(y[2]) + Z(y[3]));
  alignas(32) u32 c[8];
  _mm256_store_si256((__m256i*)c, carry[g]);
  start*= Z(c[0]) * iR;  // 区間の総積 (Montgomery の形の R を外す)
 }
 A= sa, B= sb;
}
// t[k] (0 ≤ k ≤ D、M 未満) を S(k) = Σ_{i<k} r^i t_i に置き換え、S(D+1) を返す。r^i は 8 レーンの鎖 (Montgomery の形) で作り、
// ベクトルの中の累積和と前のベクトルからの繰り上げで足す。
template <class Z> Z prefix_sums(u32* t, int D, Z r, const Mont<Z>& mt) {
 const int n= D + 1, nv= n / 8;
 alignas(32) u32 a[8];
 {
  Z x(1);
  for(int b= 0; b < 8; ++b) a[b]= (x * Z::raw(mt.R1)).val(), x*= r;
 }
 __m256i rv= _mm256_load_si256((const __m256i*)a);
 const __m256i r8= _mm256_set1_epi32((r.pow(8) * Z::raw(mt.R1)).val());
 const __m256i i1= _mm256_setr_epi32(0, 0, 1, 2, 3, 4, 5, 6), i2= _mm256_setr_epi32(0, 0, 0, 1, 2, 3, 4, 5), i4= _mm256_setr_epi32(0, 0, 0, 0, 0, 1, 2, 3);
 const __m256i m1= _mm256_setr_epi32(0, -1, -1, -1, -1, -1, -1, -1), m2= _mm256_setr_epi32(0, 0, -1, -1, -1, -1, -1, -1), m4= _mm256_setr_epi32(0, 0, 0, 0, -1, -1, -1, -1);
 const __m256i i7= _mm256_set1_epi32(7);
 __m256i carry= _mm256_setzero_si256();
 for(int j= 0; j < nv; ++j) {
  const __m256i w= mt.red(mt.mul(_mm256_loadu_si256((const __m256i*)(t + 8 * j)), rv));  // r^k t_k
  rv= mt.mul(rv, r8);
  __m256i s= mt.add(w, _mm256_and_si256(_mm256_permutevar8x32_epi32(w, i1), m1));
  s= mt.add(s, _mm256_and_si256(_mm256_permutevar8x32_epi32(s, i2), m2));
  s= mt.add(s, _mm256_and_si256(_mm256_permutevar8x32_epi32(s, i4), m4));
  const __m256i incl= mt.add(s, carry);
  _mm256_storeu_si256((__m256i*)(t + 8 * j), mt.red(_mm256_add_epi32(_mm256_sub_epi32(incl, w), mt.mv)));  // 排他的な累積和
  carry= _mm256_permutevar8x32_epi32(incl, i7);
 }
 _mm256_store_si256((__m256i*)a, carry);
 Z S= Z::raw(a[0]);
 _mm256_store_si256((__m256i*)a, rv);
 Z rp= Z(a[0]) * Z::raw(mt.R1).inv();  // r^{8 nv}
 for(int k= 8 * nv; k < n; ++k) {
  const Z x= Z::raw(t[k]);
  t[k]= S.val(), S+= rp * x, rp*= r;
 }
 return S;
}
// 0..DD での値 val (M 未満) の多項式に重み rf^k を掛けた Σ_k val_k L_k(x) rf^k を、x mod M > DD として求める (rq = 1/rf)。
// W_k = L_k(x) rf^k の比は W_{k+1}/W_k = -(DD-k)(x-k) / ((k+1)(x-k-1) rq) で、q_k = (k+1)(x-k-1) rq の逆元をまとめて求め、
// 比の累積積を和にためる。W_0 = (-1)^DD P'/DD! (P' = Π_{j=1}^{DD} (x-j)) は、Q = Π q_k = DD! P' rq^DD から
// (-1)^DD P'^2 rq^DD / Q と求める。返すのは W_0 と、A = Σ_{k≥1} val_k W_k / W_0、B = Σ_{k≥1} W_k / W_0。
template <class Z> void lagrange_scan(const u32* val, int DD, u32 xm, Z rq, const Mont<Z>& mt, Z& W0, Z& A, Z& B) {
 if(DD == 0) return W0= Z(1), A= B= Z(0), void();
 const int steps= (DD + 7) / 8;
 hvec W(8 * steps + 8), rho(8 * steps + 8, mt.R1);
 const __m256i lane= _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7), one= _mm256_set1_epi32(1), last= _mm256_set1_epi32(DD - 1);
 const __m256i rqR2= _mm256_set1_epi32((rq * Z::raw(mt.R2)).val()), xv= _mm256_set1_epi32(xm), x1= _mm256_set1_epi32(xm - 1);
 // P' = Π_{j=1}^{DD} (x - j)
 __m256i pacc= one;
 for(int s= 0; s < steps; ++s) {
  const __m256i k= _mm256_add_epi32(lane, _mm256_set1_epi32(8 * s));
  pacc= mt.mul(pacc, _mm256_blendv_epi8(_mm256_sub_epi32(x1, k), one, _mm256_cmpgt_epi32(k, last)));
 }
 Z Pp(1);
 {
  alignas(32) u32 c[8];
  _mm256_store_si256((__m256i*)c, pacc);
  for(int b= 0; b < 8; ++b) Pp*= Z(c[b]);
  Pp*= Z::raw(mt.R1).pow(u64(8) * steps);
 }
 auto gen= [&](int s) {  // q_k = (k+1)(x-k-1) rq (k = 8s + b < DD、それ以外は 1)
  const __m256i k= _mm256_add_epi32(lane, _mm256_set1_epi32(8 * s));
  const __m256i q= mt.red(mt.mul(mt.mul(_mm256_add_epi32(k, one), _mm256_sub_epi32(x1, k)), rqR2));
  return _mm256_blendv_epi8(q, one, _mm256_cmpgt_epi32(k, last));
 };
 const __m256i dd= _mm256_set1_epi32(DD);
 const Z Q= batch_inv<Z>(steps, W.data(), mt.R3, mt, gen, [&](int s, __m256i iq) {
  // rho_k = -(DD-k)(x-k) / q_k を Montgomery の形で: iq = R^3/q、z = (DD-k)(x-k)/R なので mont(iq, z) = (DD-k)(x-k) R / q
  const __m256i k= _mm256_add_epi32(lane, _mm256_set1_epi32(8 * s));
  const __m256i z= mt.mul(_mm256_sub_epi32(dd, k), _mm256_sub_epi32(xv, k));
  const __m256i v= mt.neg(mt.red(mt.mul(iq, z)));
  _mm256_storeu_si256((__m256i*)(rho.data() + 8 * s), _mm256_blendv_epi8(v, _mm256_set1_epi32(mt.R1), _mm256_cmpgt_epi32(k, last)));
 });
 W0= Pp * Pp * rq.pow(DD) * Q.inv();
 if(DD & 1) W0= -W0;
 scan_acc<Z>(rho.data(), val + 1, DD, mt, A, B);
}
// S(x) = c + r^x g(x) の c (t[0..D+1] は S(0..D+1)、r ≠ 0, 1)。c = u_0 Σ_{j≥1} S(j) u_j/u_0 / (1 - r)^{D+1}、
// u_j = C(D+1, j) (-1)^{D+1-j} r^{D+1-j}、比 u_j/u_{j-1} = -(D+2-j) / (j r)。
template <class Z> Z const_term(const u32* t, int D, Z r, const Mont<Z>& mt) {
 const int N= D + 1, steps= (N + 7) / 8;
 hvec W(8 * steps + 8), rho(8 * steps + 8, mt.R1);
 const __m256i lane= _mm256_setr_epi32(1, 2, 3, 4, 5, 6, 7, 8), one= _mm256_set1_epi32(1), nn= _mm256_set1_epi32(N);
 const __m256i rR= _mm256_set1_epi32((r * Z::raw(mt.R1)).val()), d2= _mm256_set1_epi32(D + 2);
 auto gen= [&](int s) {  // j r (j = 8s + b + 1 ≤ N、それ以外は 1)
  const __m256i j= _mm256_add_epi32(lane, _mm256_set1_epi32(8 * s));
  return _mm256_blendv_epi8(mt.red(mt.mul(j, rR)), one, _mm256_cmpgt_epi32(j, nn));
 };
 batch_inv<Z>(steps, W.data(), mt.R3, mt, gen, [&](int s, __m256i iv) {
  const __m256i j= _mm256_add_epi32(lane, _mm256_set1_epi32(8 * s));
  const __m256i z= mt.mul(_mm256_sub_epi32(d2, j), one);  // (D+2-j)/R
  const __m256i v= mt.neg(mt.red(mt.mul(iv, z)));
  _mm256_storeu_si256((__m256i*)(rho.data() + 8 * s), _mm256_blendv_epi8(v, _mm256_set1_epi32(mt.R1), _mm256_cmpgt_epi32(j, nn)));
 });
 Z A, B;
 scan_acc<Z>(rho.data(), t + 1, N, mt, A, B);
 Z u0= r.pow(D + 1);
 if((D + 1) & 1) u0= -u0;
 return u0 * A * (Z(1) - r).pow(D + 1).inv();
}
// y を u32 の作業の配列 (huge page、8 の倍数ぶん 0 で延ばす) に写す。
template <class Z> hvec copy_vals(const std::vector<Z>& y, int extra) {
 hvec t(8 * ((y.size() + extra + 7) / 8) + 16, 0);
 for(std::size_t i= 0; i < y.size(); ++i) t[i]= y[i].val();
 return t;
}
template <class Z> Z iota_zmod(const std::vector<Z>& y, u64 x) {
 const int K= int(y.size());
 if(x < u64(K)) return y[x];
 const u32 xm= u32(x % Z::mod());
 if(xm < u32(K)) return y[xm];
 const Mont<Z> mt;
 const hvec t= copy_vals(y, 0);
 Z W0, A, B;
 lagrange_scan<Z>(t.data(), K - 1, xm, Z(1), mt, W0, A, B);
 return W0 * (y[0] + A);
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
// a x mod M を 8 個ずつ (a < M、x < M、q = floor(x 2^32 / M))。M 未満で返す。
[[gnu::always_inline]] inline __m256i shoup8(__m256i a, __m256i x, __m256i q, __m256i p) {
 const __m256i ev= _mm256_srli_epi64(_mm256_mul_epu32(a, q), 32), od= _mm256_mul_epu32(_mm256_srli_epi64(a, 32), q);
 const __m256i r= _mm256_sub_epi32(_mm256_mullo_epi32(a, x), _mm256_mullo_epi32(_mm256_blend_epi32(ev, od, 0b10101010), p));
 return _mm256_min_epu32(r, _mm256_sub_epi32(r, p));
}
template <class Z> std::vector<Z> interp_zmod(const std::vector<Z>& X, const std::vector<Z>& Y) {
 const int n= int(X.size());
 const Mont<Z> mt;
 const u32 M= mt.M;
 auto red8= [&](__m256i v) { return _mm256_min_epu32(v, _mm256_sub_epi32(v, mt.mv)); };
 const u32 R1= Z(u64(1) << 32).val(), R2= (Z::raw(R1) * Z::raw(R1)).val();
 // 1. P (長さ n + 1)。P_new[k] = P[k-1] + c P[k] (c = -x_j) を上から 8 個ずつ (読みを先に済ませる)。
 std::vector<u32> P(n + 1 + 16, 0);
 P[0]= Z(1).val();
 for(int j= 0; j < n; ++j) {
  const u32 c= (-X[j]).val(), q= u32((u64(c) << 32) / M);
  const __m256i vc= _mm256_set1_epi32(c), vq= _mm256_set1_epi32(q);
  int k= j + 2;  // P_new の添字は 0..j+1
  for(; k - 8 >= 0; k-= 8) {
   const int lo= k - 8;
   const __m256i cur= _mm256_loadu_si256((const __m256i*)(P.data() + lo));
   const __m256i prev= lo ? _mm256_loadu_si256((const __m256i*)(P.data() + lo - 1)) : _mm256_permutevar8x32_epi32(_mm256_and_si256(cur, _mm256_setr_epi32(-1, -1, -1, -1, -1, -1, -1, 0)), _mm256_setr_epi32(7, 0, 1, 2, 3, 4, 5, 6));
   _mm256_storeu_si256((__m256i*)(P.data() + lo), red8(_mm256_add_epi32(prev, shoup8(cur, vc, vq, mt.mv))));
  }
  for(--k; k >= 0; --k) P[k]= ((k ? Z::raw(P[k - 1]) : Z(0)) + Z::raw(c) * Z::raw(P[k])).val();
 }
 // 2. P'(x_i) を Horner の方法で (係数と x は Montgomery の形)。
 std::vector<u32> dPm(n);
 for(int k= 0; k < n; ++k) dPm[k]= (Z(k + 1) * Z::raw(P[k + 1]) * Z::raw(R1)).val();
 const int ng= (n + 7) / 8;
 std::vector<u32> xs(8 * ng + 32, 0), ds(8 * ng + 32, 1);
 for(int i= 0; i < n; ++i) xs[i]= X[i].val();
 const __m256i r2= _mm256_set1_epi32(R2), one= _mm256_set1_epi32(1);
 for(int g= 0; g < ng; g+= 4) {
  __m256i xm[4], d[4];
  for(int t= 0; t < 4; ++t) xm[t]= mont8(_mm256_loadu_si256((const __m256i*)(xs.data() + 8 * (g + t))), r2, mt.mv, mt.mi), d[t]= _mm256_setzero_si256();
  for(int k= n - 1; k >= 0; --k) {
   const __m256i c= _mm256_set1_epi32(dPm[k]);
   for(int t= 0; t < 4; ++t) d[t]= red8(_mm256_add_epi32(red8(mont8(d[t], xm[t], mt.mv, mt.mi)), c));
  }
  for(int t= 0; t < 4; ++t) _mm256_storeu_si256((__m256i*)(ds.data() + 8 * (g + t)), red8(mont8(d[t], one, mt.mv, mt.mi)));
 }
 // 3. w_i = y_i / P'(x_i) (逆元はまとめて)
 std::vector<u32> w(8 * ng + 32, 0);
 {
  std::vector<Z> pre(n);
  Z acc(1);
  for(int i= 0; i < n; ++i) pre[i]= acc, acc*= Z::raw(ds[i]);
  Z iv= acc.inv();
  for(int i= n - 1; i >= 0; --i) w[i]= (iv * pre[i] * Y[i]).val(), iv*= Z::raw(ds[i]);
 }
 // 4. s_m = Σ_i w_i x_i^m (m < n)。点の組 4 つずつ、m ごとに 8 レーンの和を sv[m] に足す。
 std::vector<u32> sv(8 * std::size_t(n), 0);
 for(int g= 0; g < ng; g+= 4) {
  __m256i xm[4], v[4];
  for(int t= 0; t < 4; ++t) {
   xm[t]= mont8(_mm256_loadu_si256((const __m256i*)(xs.data() + 8 * (g + t))), r2, mt.mv, mt.mi);
   v[t]= _mm256_loadu_si256((const __m256i*)(w.data() + 8 * (g + t)));
  }
  for(int m= 0; m < n; ++m) {
   __m256i acc= _mm256_loadu_si256((const __m256i*)(sv.data() + 8 * m));
   for(int t= 0; t < 4; ++t) acc= red8(_mm256_add_epi32(acc, v[t])), v[t]= red8(mont8(v[t], xm[t], mt.mv, mt.mi));
   _mm256_storeu_si256((__m256i*)(sv.data() + 8 * m), acc);
  }
 }
 std::vector<u32> s(n);
 for(int m= 0; m < n; ++m) {
  u64 t= 0;
  for(int b= 0; b < 8; ++b) t+= sv[8 * m + b];
  s[m]= Z(t).val();
 }
 // 5. f_k = Σ_m P_{k+1+m} s_m。Prev[j] = P_{n-j} として f_k = (Prev s)[n-1-k]。
 std::vector<u32> Prev(n), c(n);
 for(int j= 0; j < n; ++j) Prev[j]= P[n - j];
 convolve_internal::conv<Z>(Prev.data(), n, s.data(), n, c.data(), n);
 std::vector<Z> f(n);
 for(int k= 0; k < n; ++k) f[k]= Z::raw(c[n - 1 - k]);
 return f;
}
}  // namespace lagrange_internal
// y_i = f(i) (0 ≤ i < n) を満たす n-1 次以下の f の f(x)。
template <class K> K interpolate_iota(const std::vector<K>& y, unsigned long long x) {
 if(y.empty()) return K();
 if constexpr(convolve_internal::is_zmod_v<K>) return lagrange_internal::iota_zmod(y, x);
 else return x < y.size() ? y[x] : lagrange_internal::iota_generic(y, K(x));
}
// 点 (x_i, y_i) を通る n-1 次以下の多項式の係数。f = Σ_i w_i P(X)/(X - x_i) (P = Π (X - x_j)、w_i = y_i / P'(x_i))。
// ZMod (法は 2 より大きい素数) では、P/(X - x_i) の k 次の係数が Σ_{t>k} P_t x_i^{t-k-1} なので
// f_k = Σ_m P_{k+1+m} s_m (s_m = Σ_i w_i x_i^m) と書き、P は Shoup の方法で 8 個ずつ、P'(x_i) と s_m は点を 8 個ずつ
// レーンに載せた Montgomery の積 (4 組を並べる) で作り、f は逆に並べた P と s の畳み込みの下半分を素朴な積の核で求める。
// ほかの型は、点ごとに組立除法で P/(X - x_i) を作って足す。
template <class K> std::vector<K> interpolate(const std::vector<K>& x, const std::vector<K>& y) {
 const int n= int(x.size());
 if(n == 0) return {};
 if constexpr(convolve_internal::is_zmod_v<K>) {
  if(n > 8 && (K::mod() & 1)) return lagrange_internal::interp_zmod(x, y);
 }
 std::vector<K> P(n + 1, K());
 P[0]= K(1);
 for(int j= 0; j < n; ++j) {
  for(int k= j + 1; k > 0; --k) P[k]= P[k - 1] - x[j] * P[k];
  P[0]= -x[j] * P[0];
 }
 std::vector<K> f(n, K()), q(n);
 for(int i= 0; i < n; ++i) {
  // q = P / (X - x_i)、q_{n-1} = 1、q_{k-1} = P_k + x_i q_k
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
 if(r == K(1)) {  // F(x) = Σ_{i<x} f(i) は d+1 次の多項式で、F(0..d+1) から補間する
  if constexpr(convolve_internal::is_zmod_v<K>) {
   if(K::mod() & 1) {  // 作業の配列の上で累積和を取り、そのまま補間する
    const Mont<K> mt;
    hvec t= copy_vals(y, 1);
    t[D + 1]= prefix_sums(t.data(), D, K(1), mt).val();
    if(n <= u64(D) + 1) return K::raw(t[n]);
    const u32 xm= u32(n % K::mod());
    if(xm <= u32(D) + 1) return K::raw(t[xm]);
    K W0, A, B;
    lagrange_scan<K>(t.data(), D + 1, xm, K(1), mt, W0, A, B);
    return W0 * (K::raw(t[0]) + A);
   }
  }
  std::vector<K> F(D + 2);
  F[0]= K();
  for(int k= 0; k <= D; ++k) F[k + 1]= F[k] + y[k];
  return interpolate_iota(F, n);
 }
 if constexpr(convolve_internal::is_zmod_v<K>) {
  if(K::mod() & 1) {
   const Mont<K> mt;
   hvec t= copy_vals(y, 1);
   t[D + 1]= prefix_sums(t.data(), D, r, mt).val();
   if(n <= u64(D) + 1) return K::raw(t[n]);
   const K c= const_term(t.data(), D, r, mt);
   const u32 xm= u32(n % K::mod());
   K gn;
   if(xm <= u32(D)) gn= (K::raw(t[xm]) - c) * r.pow(xm).inv();
   else {  // g(n) = Σ_k (S_k - c) r^{-k} L_k(x) = W_0 ((S_0 - c) + A - c B)
    K W0, A, B;
    lagrange_scan<K>(t.data(), D, xm, r, mt, W0, A, B);
    gn= W0 * (K::raw(t[0]) - c + A - c * B);
   }
   return c + r.pow(n) * gn;
  }
 }
 {  // 汎用の道 (ZMod で法が偶数のときも)
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
 if constexpr(convolve_internal::is_zmod_v<K>) {
  if(K::mod() & 1) {
   const Mont<K> mt;
   hvec t= copy_vals(y, 1);
   t[D + 1]= prefix_sums(t.data(), D, r, mt).val();
   return const_term(t.data(), D, r, mt);
  }
 }
 {  // 汎用の道 (ZMod で法が偶数のときも)
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
