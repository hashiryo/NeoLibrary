#pragma once
// スカラーの実装 (基数 4)。ntt.hpp が、AVX2 も SIMDe も無い処理系でだけ include する。
// 段のひねりは段ごとの表 (w1、w2、w3) に持ち、内側の繰り返しで ω を掛け直さない。
namespace ntt_internal {

// 1 の原始 4 乗根 ω^{(P-1)/4} とその逆。基数 4 の段の i の掛け算に使う。
template <unsigned P> const ZMod<P>& root4() {
 static const ZMod<P> v= pow_mod(primitive_root<P>(), u64(P - 1) / 4);
 return v;
}
template <unsigned P> const ZMod<P>& iroot4() {
 static const ZMod<P> v= root4<P>().inv();
 return v;
}

// 段のひねりの表。長さ len の段では h = len/4 として、tw[j] = ω_len^j、tw[h + j] = ω_len^{2j}、tw[2h + j] = ω_len^{3j}。
// 同じ段のすべてのブロックで使い回す。大きさは 3 * (n/4) まで。
template <unsigned P> struct scratch {
 static inline thread_local std::vector<ZMod<P>> tw;
};
template <unsigned P> ZMod<P>* get_tw(int h) {
 auto& v= scratch<P>::tw;
 if ((int)v.size() < 3 * h) v.resize(3 * h);
 return v.data();
}

// tw に 3 本の等比数列を詰める。wl は ω_len (逆変換では ω_len^{-1})。
template <unsigned P> void fill_tw(ZMod<P>* tw, int h, ZMod<P> wl) {
 ZMod<P>* w1= tw;
 ZMod<P>* w2= tw + h;
 ZMod<P>* w3= tw + 2 * h;
 const ZMod<P> wl2= wl * wl, wl3= wl2 * wl;
 w1[0]= w2[0]= w3[0]= ZMod<P>(1);
 for (int j= 1; j < h; ++j) w1[j]= w1[j - 1] * wl, w2[j]= w2[j - 1] * wl2, w3[j]= w3[j - 1] * wl3;
}

// DIF。a[i] に Σ_j a[j] ω^{j k} が入る。k は位置 i の 4 進の桁を逆にした数 (2 の指数が奇数なら、最後の基数 2 の段で最上位の桁が最下位に来る)。
template <unsigned P> void dft_impl(ZMod<P>* a, int n) {
 if (n <= 1) return;
 assert((n & (n - 1)) == 0 && (P - 1) % unsigned(n) == 0);
 const ZMod<P> g= primitive_root<P>();
 const auto fi= root4<P>().fixed();
 int len= n;
 for (; len >= 4; len>>= 2) {
  const int h= len>> 2;
  ZMod<P>* tw= get_tw<P>(h);
  fill_tw<P>(tw, h, pow_mod(g, u64(P - 1) / unsigned(len)));
  for (int i= 0; i < n; i+= len) {
   ZMod<P>* x0= a + i;
   ZMod<P>* x1= x0 + h;
   ZMod<P>* x2= x1 + h;
   ZMod<P>* x3= x2 + h;
   for (int j= 0; j < h; ++j) {
    const ZMod<P> u0= x0[j], u1= x1[j], u2= x2[j], u3= x3[j];
    const ZMod<P> t0= u0 + u2, t1= u0 - u2, t2= u1 + u3, t3= u1 - u3;
    const ZMod<P> y= t3 * fi;
    x0[j]= t0 + t2;
    x1[j]= (t1 + y) * tw[j];
    x2[j]= (t0 - t2) * tw[h + j];
    x3[j]= (t1 - y) * tw[2 * h + j];
   }
  }
 }
 if (len == 2) {
  // 長さ 2 の段 (最後の 1 段)。ひねりは 1。
  for (int i= 0; i < n; i+= 2) {
   const ZMod<P> u= a[i], v= a[i + 1];
   a[i]= u + v, a[i + 1]= u - v;
  }
 }
}

// DIT。dft の逆。段ごとの 1/4 (最後の基数 2 の段の 1/2) は掛けず、最後に 1/n を掛ける。
template <unsigned P> void idft_impl(ZMod<P>* a, int n) {
 if (n <= 1) return;
 assert((n & (n - 1)) == 0 && (P - 1) % unsigned(n) == 0);
 const ZMod<P> g= primitive_root<P>();
 const auto fii= iroot4<P>().fixed();
 if (__builtin_ctz(unsigned(n)) & 1) {
  // dft の最後の基数 2 の段を先に戻す。
  for (int i= 0; i < n; i+= 2) {
   const ZMod<P> u= a[i], v= a[i + 1];
   a[i]= u + v, a[i + 1]= u - v;
  }
 }
 // dft の基数 4 の段の長さは、2 の指数が偶数なら 4, 16, ..., n、奇数なら 8, 32, ..., n。逆はその昇順。
 for (int len= (__builtin_ctz(unsigned(n)) & 1) ? 8 : 4; len <= n; len<<= 2) {
  const int h= len>> 2;
  ZMod<P>* tw= get_tw<P>(h);
  fill_tw<P>(tw, h, pow_mod(g, u64(P - 1) - u64(P - 1) / unsigned(len)));
  for (int i= 0; i < n; i+= len) {
   ZMod<P>* x0= a + i;
   ZMod<P>* x1= x0 + h;
   ZMod<P>* x2= x1 + h;
   ZMod<P>* x3= x2 + h;
   for (int j= 0; j < h; ++j) {
    const ZMod<P> u0= x0[j];
    const ZMod<P> y1= x1[j] * tw[j];
    const ZMod<P> y2= x2[j] * tw[h + j];
    const ZMod<P> y3= x3[j] * tw[2 * h + j];
    const ZMod<P> t0= u0 + y2, t2= u0 - y2, t1= y1 + y3, t3= (y1 - y3) * fii;
    x0[j]= t0 + t1;
    x1[j]= t2 + t3;
    x2[j]= t0 - t1;
    x3[j]= t2 - t3;
   }
  }
 }
 const ZMod<P> iv= ZMod<P>(unsigned(n)).inv();
 for (int i= 0; i < n; ++i) a[i]*= iv;
}

// a[0..n) は長さ n の変換。係数に戻して 0 を詰め、長さ 2n の変換を a[0..2n) に作る (素直な形。変換の使い回しはまだ作っていない)。
template <unsigned P> void dft_doubling_impl(ZMod<P>* a, int n) {
 std::vector<ZMod<P>> c(a, a + n);
 idft_impl<P>(c.data(), n);
 std::fill(a, a + 2 * n, ZMod<P>());
 std::copy(c.begin(), c.end(), a);
 dft_impl<P>(a, 2 * n);
}

// 畳み込みの作業の配列。長さごとに残す (relaxed のように小さい畳み込みを何度も呼ぶ道具のため)。
// 呼び出しのたびに確保すると、その費用がそのまま遅さになる。
template <unsigned P> struct conv_work {
 static inline thread_local std::vector<ZMod<P>> fa, fb;
};

// 畳み込み。Ntt::convolve_range が、最小の長さと法の条件を確かめてから呼ぶ。
template <unsigned P> void conv_kernel(const ZMod<P>* a, int n, const ZMod<P>* b, int m, int lo, int hi, ZMod<P>* out) {
 const int need= std::max(hi, n + m - 1 - lo);
 int l= 1;
 while (l < need) l<<= 1;
 auto& fa= conv_work<P>::fa;
 auto& fb= conv_work<P>::fb;
 if ((int)fa.size() < l) fa.resize(l);
 std::copy(a, a + n, fa.begin());
 std::fill(fa.begin() + n, fa.begin() + l, ZMod<P>());
 if (a == b && n == m) {
  dft_impl<P>(fa.data(), l);
  for (int i= 0; i < l; ++i) fa[i]*= fa[i];
 } else {
  if ((int)fb.size() < l) fb.resize(l);
  std::copy(b, b + m, fb.begin());
  std::fill(fb.begin() + m, fb.begin() + l, ZMod<P>());
  dft_impl<P>(fa.data(), l), dft_impl<P>(fb.data(), l);
  for (int i= 0; i < l; ++i) fa[i]*= fb[i];
 }
 idft_impl<P>(fa.data(), l);
 for (int k= lo; k < hi; ++k) out[k - lo]= fa[k];
}
}  // namespace ntt_internal
