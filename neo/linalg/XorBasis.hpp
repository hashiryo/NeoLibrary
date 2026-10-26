#pragma once
#include <algorithm>
#include <bit>
#include <cstdint>
#include <vector>
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
// GF(2) 上の線形空間の基底 (xor 基底)。U は 128 bit までの符号なし整数で、U の値を W bit (W は U の bit 数) のベクトルとみる。
// 基底の元は入れた順に持ち、あとから入った元は前の元の最上位の bit (要) をどれも持たない (入れるときに消してある)。z を基底で
// 消すときは、入れた順に z = min(z, z ^ b) を当てる。b の要が z に立っているときだけ z ^ b < z なのでその bit が消え、あとの元は
// その bit を持たないので、消した bit がまた立つことはない。残るのは z + (張る空間) のうち要がどれも 0 の唯一の元で、それが
// z + (張る空間) の最小になる。z = max(z, z ^ b) なら要がどれも 1 の元、つまり最大が残る。どちらも分岐がなく、要の位置を clz で
// 探す形より速い。procon-judge の yosupo-intersection-of-f2-vector-spaces で書き比べた zassenhaus_min と同じ形。
namespace xor_basis_internal {
using u64= std::uint64_t;
// x の最上位の bit の位置。x > 0。
template <class U> constexpr int top(U x) {
 if constexpr(sizeof(U) > 8) {
  const u64 h= u64(x >> 64);
  return h ? 127 - std::countl_zero(h) : 63 - std::countl_zero(u64(x));
 } else return 63 - std::countl_zero(u64(x));
}
// x の最下位の bit の位置。x > 0。
template <class U> constexpr int low(U x) {
 if constexpr(sizeof(U) > 8) {
  const u64 l= u64(x);
  return l ? std::countr_zero(l) : 64 + std::countr_zero(u64(x >> 64));
 } else return std::countr_zero(u64(x));
}
// k の下から i 番目の bit を、m の下から i 番目に立っている bit の位置に置く (BMI2 の pdep)。
inline u64 deposit64(u64 k, u64 m) {
#if defined(__x86_64__) && defined(__BMI2__)
 return _pdep_u64(k, m);
#else
 u64 r= 0;
 for(; m; m&= m - 1, k>>= 1) r|= m & -m & -(k & 1);
 return r;
#endif
}
template <class U> U deposit(U k, U m) {
 if constexpr(sizeof(U) > 8) {
  const u64 lo= u64(m);
  return U(deposit64(u64(k), lo)) | U(deposit64(u64(k >> std::popcount(lo)), u64(m >> 64))) << 64;
 } else return U(deposit64(k, m));
}
}
template <class U= std::uint64_t> class XorBasis {
 static_assert(U(-1) > U(0), "U は符号なし整数");
public:
 static constexpr int W= int(sizeof(U)) * 8;
private:
 U v[W]{}, piv= 0;  // v[0..r) が基底の元 (入れた順)、piv は要の集合
 int r= 0;
 template <class V> friend XorBasis<V> intersection(const XorBasis<V>&, const XorBasis<V>&);
 template <class V> friend class RangeXorBasis;
public:
 // 基底の元の数 (張る空間の次元)。
 int size() const { return r; }
 // i 番目に入った基底の元 (入れたときに前の元で消したもの)。
 U operator[](int i) const { return v[i]; }
 const U* begin() const { return v; }
 const U* end() const { return v + r; }
 // 張る空間の元の数 2^size()。W < 64 なら u64、W ≥ 64 なら unsigned __int128 で返す (W = 128 で size() = 128 なら 0)。
 auto count() const {
  if constexpr(W < 64) return std::uint64_t(1) << r;
  else {
   using u128= unsigned __int128;
   return r < 128 ? u128(1) << r : u128(0);
  }
 }
 // x を入れる。張る空間が大きくなったら true。
 bool insert(U x) {
  for(int i= 0; i < r; ++i) x= std::min<U>(x, x ^ v[i]);
  if(!x) return false;
  v[r++]= x, piv|= U(1) << xor_basis_internal::top(x);
  return true;
 }
 // o の張る空間を足す。
 void merge(const XorBasis& o) {
  for(int i= 0; i < o.r; ++i) insert(o.v[i]);
 }
 // x を基底で消した余り。x + (張る空間) の最小で、要の bit はどれも 0。
 U reduce(U x) const {
  for(int i= 0; i < r; ++i) x= std::min<U>(x, x ^ v[i]);
  return x;
 }
 // x が張る空間に入るか。
 bool contains(U x) const { return !reduce(x); }
 // x + (張る空間) の最大。
 U max_xor(U x= 0) const {
  for(int i= 0; i < r; ++i) x= std::max<U>(x, x ^ v[i]);
  return x;
 }
 // x + (張る空間) の元のうち、小さいほうから k 番目 (0 始まり)。0 <= k < 2^size()。
 // 2 つの元の大小は要の bit の並びの大小で決まるので、k の bit を要の位置に配った K について、要の bit が K になる元を返す。
 // K ^ reduce(K) は張る空間の元で要の bit が K のもの、reduce(x) は x + (張る空間) の元で要の bit が 0 のもの。
 U kth(U k, U x= 0) const {
  const U K= xor_basis_internal::deposit(k, piv);
  return reduce(x) ^ K ^ reduce(K);
 }
};
// 2 つの基底の張る空間の共通部分の基底。Zassenhaus の方法で、元の多いほうを A、少ないほうを B とし、2W bit の行 (x, t) を、
// A の元 a から (a, 0)、B の元 b から (b, b) として並べて x の側で消していく。x が 0 まで消えた B の行の t は共通部分に入る
// (t は B の元の和で、x - t はいつも A の元の和なので)。B の元は要の小さい順に流す。行の t は先に流した B の元の和なので、
// x が 0 まで消えた行の t の要は流した元の要と同じになり、得られる t は要の異なる階段の形になる (要の大きい順に並べれば
// そのまま基底の形)。2W ≤ 64 なら行を 1 語 (x << W) | t に詰めて min の形で消す。
template <class U> XorBasis<U> intersection(const XorBasis<U>& a, const XorBasis<U>& b) {
 using namespace xor_basis_internal;
 constexpr int W= XorBasis<U>::W;
 const bool sw= a.r < b.r;
 const XorBasis<U>&A= sw ? b : a, &B= sw ? a : b;
 const int ra= A.r, rb= B.r;
 if(!rb) return XorBasis<U>();
 // A が [0, 2^ra) をすべて張り、B の元がどれも 2^ra 未満なら、B は A に含まれる。
 if(ra == W || (A.piv == U((U(1) << ra) - 1) && !(B.piv >> ra))) return B;
 U ys[W], out[W];
 {
  std::uint8_t at[W];
  for(int i= 0; i < rb; ++i) at[top(B.v[i])]= std::uint8_t(i);
  int k= 0;
  for(U m= B.piv; m; m&= m - 1) ys[k++]= B.v[at[low(m)]];
 }
 int k= 0;
 if constexpr(2 * W <= 64) {
  using P= std::uint64_t;
  P h[W];
  int nh= 0;
  for(int i= 0; i < ra; ++i) h[nh++]= P(A.v[i]) << W;
  for(int i= 0; i < rb; ++i) {
   P z= P(ys[i]) << W | ys[i];
   for(int j= 0; j < nh; ++j) z= std::min(z, z ^ h[j]);
   if(z >> W) h[nh++]= z;
   else out[k++]= U(z);
  }
 } else {
  U hx[W], ht[W];
  int nh= 0;
  for(int i= 0; i < ra; ++i) hx[nh]= A.v[i], ht[nh++]= 0;
  for(int i= 0; i < rb; ++i) {
   U x= ys[i], t= ys[i];
   for(int j= 0; j < nh; ++j) {
    const U y= x ^ hx[j], m= U(0) - U(y < x);
    x= std::min(x, y), t^= ht[j] & m;
   }
   if(x) hx[nh]= x, ht[nh++]= t;
   else out[k++]= t;
  }
 }
 XorBasis<U> res;
 for(int i= 0; i < k; ++i) res.v[i]= out[k - 1 - i], res.piv|= U(1) << top(out[i]);
 res.r= k;
 return res;
}
// 列 a の区間 [l, r) の元の張る空間を問う。接頭辞 a[0..r) ごとに、要 h ごとに最も右の元を残す基底 (欄 h の元 b[h] と、それを
// 作った元の位置 p[h]) を考える。新しい元を入れるときは、上の bit から見て、その欄の元より新しいほうを欄に残し、2 つの xor を
// 古いほうの位置で下へ流す。こうすると、位置が l 以上の欄の元が a[l..r) の張る空間の基底になる。
// 写しは 16 個おきにだけ持つ。問い合わせでは r を 16 の倍数 kB に切り下げた写しの、位置が l 以上の欄 (a[l..kB) の基底) で、
// x と残りの元 a[max(kB, l)..r) (15 個以下) を上の bit から min の形で消す。消したあとの x と残りの元は写しの要をどれも
// 持たないので、残りの元を入れた順に互いに消してから x に当てれば、x が区間の元の和で書けるかが分かる。写しは欄が
// L = (a の元の最上位の bit の位置の最大) + 1 個で、(N / 16 + 1) L 個の (U, int) と a の写しを持つ (U = u64、N = 4 × 10^5、
// L = 60 で 21 MB ほど)。procon-judge の atcoder-abc223-h で書き比べた snap16_red の形。
template <class U= std::uint64_t> class RangeXorBasis {
 static constexpr int W= int(sizeof(U)) * 8, B= 16;
 int n= 0, L= 0;
 std::vector<U> a, sb;  // sb[k L + h] は a[0..kB) の基底の欄 h の元 (無ければ 0)
 std::vector<int> sp;   // sp[k L + h] はその位置 (無ければ -1)
 // e[0..16) に a[max(kB, l)..r) を置き、f(m) を写し k の位置が l 以上の欄の元 m (上の bit から順、ほかの欄は 0) で呼び
 // ながら e を消し、e を入れた順に互いに消す。e の数を返す。U が 64 bit なら、e を AVX2 の 4 本に持って欄ごとに消す。欄 h の
 // 元は最上位の bit が h なので、e の bit h を符号の位置へずらし、それが立っている lane だけ m との xor に取り替える。
 template <class F> int run(int l, int r, U* e, F&& f) const {
  const int k= r / B, s= std::max(k * B, l), ne= std::max(r - s, 0);
  for(int t= 0; t < ne; ++t) e[t]= a[s + t];
  if(l < k * B) {
   const U* b= sb.data() + std::size_t(k) * L;
   const int* p= sp.data() + std::size_t(k) * L;
   bool done= false;
   if constexpr(sizeof(U) == 8)
    if(ne) {
     done= true;
     for(int t= ne; t < B; ++t) e[t]= 0;
     const int nv= (ne + 3) >> 2;
     __m256i E0= _mm256_loadu_si256((const __m256i*)e), E1= _mm256_loadu_si256((const __m256i*)(e + 4)), E2= _mm256_loadu_si256((const __m256i*)(e + 8)), E3= _mm256_loadu_si256((const __m256i*)(e + 12));
     __m256i C= _mm256_set1_epi64x(63 - L);
     const __m256i one= _mm256_set1_epi64x(1);
     auto step= [&](__m256i& E, __m256i M) { E= _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(E), _mm256_castsi256_pd(_mm256_xor_si256(E, M)), _mm256_castsi256_pd(_mm256_sllv_epi64(E, C)))); };
     for(int h= L; h--;) {
      const U m= b[h] & (U(0) - U(p[h] >= l));
      f(m);
      const __m256i M= _mm256_set1_epi64x((long long)m);
      C= _mm256_add_epi64(C, one);
      switch(nv) {
       case 4: step(E3, M); [[fallthrough]];
       case 3: step(E2, M); [[fallthrough]];
       case 2: step(E1, M); [[fallthrough]];
       case 1: step(E0, M); [[fallthrough]];
       default: break;
      }
     }
     _mm256_storeu_si256((__m256i*)e, E0), _mm256_storeu_si256((__m256i*)(e + 4), E1), _mm256_storeu_si256((__m256i*)(e + 8), E2), _mm256_storeu_si256((__m256i*)(e + 12), E3);
    }
   if(!done) {
    for(int h= L; h--;) {
     const U m= b[h] & (U(0) - U(p[h] >= l));
     f(m);
     for(int t= 0; t < ne; ++t) e[t]= std::min<U>(e[t], e[t] ^ m);
    }
   }
  }
  for(int t= 1; t < ne; ++t) {
   U y= e[t];
   for(int q= 0; q < t; ++q) y= std::min<U>(y, y ^ e[q]);
   e[t]= y;
  }
  return ne;
 }
public:
 RangeXorBasis()= default;
 explicit RangeXorBasis(const std::vector<U>& a_): n(a_.size()), a(a_) {
  U all= 0;
  for(U x: a) all|= x;
  L= all ? xor_basis_internal::top(all) + 1 : 0;
  const std::size_t ns= n / B + 1;
  sb.assign(ns * L, 0), sp.assign(ns * L, -1);
  U b[W]{};
  int p[W];
  std::fill(p, p + W, -1);
  for(int i= 0; i < n; ++i) {
   // 空の欄 (位置 -1) に入れると、流す位置が -1 になり、以後どの欄も変えない。
   U x= a[i];
   int pos= i;
   // x の最上位の bit から下へ見て、x が 0 になるか空の欄に入ったら止める。選びはマスクの算術で書く (三項演算子だと gcc が
   // 分岐にする)。欄 h の元は最上位の bit が h なので、x の更新は min の形。
   for(int h= x ? xor_basis_internal::top(x) : 0; x && pos >= 0; --h) {
    const U bh= b[h];
    const int ph= p[h], hm= -int(x >> h & 1), sm= hm & -int(pos > ph);
    b[h]= bh ^ ((x ^ bh) & (U(0) - U(sm & 1))), p[h]= ph ^ ((pos ^ ph) & sm);
    pos= std::min(pos, ph ^ ((ph ^ 0x7fffffff) & ~hm)), x= std::min<U>(x, x ^ bh);
   }
   if((i + 1) % B == 0) {
    const std::size_t o= std::size_t(i + 1) / B * L;
    std::copy(b, b + L, sb.begin() + o), std::copy(p, p + L, sp.begin() + o);
   }
  }
 }
 int size() const { return n; }
 // x が a[l..r) の元の和で書けるか (空の和 0 も含む)。0 <= l <= r <= size()。
 bool contains(int l, int r, U x) const {
  U e[B];
  const int ne= run(l, r, e, [&](U m) { x= std::min<U>(x, x ^ m); });
  for(int t= 0; t < ne; ++t) x= std::min<U>(x, x ^ e[t]);
  return !x;
 }
 // x ^ (a[l..r) の元の和) の最大。
 U max_xor(int l, int r, U x= 0) const {
  U e[B];
  const int ne= run(l, r, e, [&](U m) { x= std::max<U>(x, x ^ m); });
  for(int t= 0; t < ne; ++t) x= std::max<U>(x, x ^ e[t]);
  return x;
 }
 // a[l..r) の元の張る空間の次元。
 int rank(int l, int r) const {
  U e[B];
  int c= 0;
  const int ne= run(l, r, e, [&](U m) { c+= m != 0; });
  for(int t= 0; t < ne; ++t) c+= e[t] != 0;
  return c;
 }
 // a[l..r) の元の張る空間の基底。写しの欄の元 (上の bit から順) と、消したあとの残りの元を並べたもので、XorBasis の形になる。
 XorBasis<U> basis(int l, int r) const {
  U e[B];
  XorBasis<U> res;
  const int ne= run(l, r, e, [&](U m) {
   if(m) res.v[res.r++]= m, res.piv|= U(1) << xor_basis_internal::top(m);
  });
  for(int t= 0; t < ne; ++t)
   if(e[t]) res.v[res.r++]= e[t], res.piv|= U(1) << xor_basis_internal::top(e[t]);
  return res;
 }
};
