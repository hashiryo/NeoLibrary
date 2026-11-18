#pragma once
#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <new>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include "neo/algebra/ZMod.hpp"
#include "neo/data_structure/segtree/Affine.hpp"
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#ifdef __linux__
#include <sys/mman.h>
#endif
// 列の木 (塊の形)。節が連続する値の配列 (塊、最大 CAP 個) を持つ treap で、木が扱う節の数を値の数の CAP 分の 1 ほどにして、木を
// 下りる段数を減らす。構造体 M の決まりは SegmentTree と同じで、モノイドは T、ti、op、作用は E、mp、cp (mp は区間の長さを受けて
// もよい)。M の代わりに int などの型を渡すと値を並べるだけの列になる。REV で反転を選び、可換でないモノイドを反転するときは逆順の
// 積も持つ (M に using commute = void を書くと省く)。永続にはできない (永続は WeightBalancedTree)。
// - 同じ型の木は節と塊の置き場を共有するので、別の木どうしで切り貼りできる。木の値を写すと同じ木を指すだけなので、切ったり
//   つなげたりしたあとは元の値を使わない。
// - 併合で境目の塊のどちらかが CAP / 4 より短ければ、2 つを 1 つにまとめる (各節に部分木の左端と右端の塊の長さを持つ)。満杯の塊
//   への挿入は塊を半分ずつに割り、空になった節は外す。
// - 区間の和と作用は形を変えずに根から下り、端の塊だけを触る。塊の遅延 (自分の塊へ当てる作用と反転の印) は子へ配る遅延と別に
//   持ち、塊の反転の印は添字の読み替えで済ませる。区間の反転は、両端の塊のうち区間の外にある部分を逆順にして互いの塊へ移して
//   から、塊の境目で割って印を付ける (塊を割らない、Library Checker の nandhagk の 402977 の方法)。
// - monoid::AffineSum<ZMod<P>> (法が定数) では、塊の大きさを 256 にし、塊の中の和、作用、遅延の当て直しを AVX2 で 8 個ずつ行う
//   (値は 2P 未満のそのままの値、掛け算は Shoup の方法)。塊の一部への作用は残っている塊の遅延と合わせて 1 回の走査で書き、
//   部分の和は短い側を足して求める。
// procon-judge の dynamic-sequence-range-affine-range-sum などで書き比べた (記録は algo-notes の notes/sequence-tree.md)。
namespace sequence_internal {
template <class M, class= void> struct monoid: std::false_type {};
template <class M> struct monoid<M, std::void_t<typename M::T, decltype(&M::op), decltype(&M::ti)>>: std::true_type {};
template <class M, class= void> struct action: std::false_type {};
template <class M> struct action<M, std::void_t<typename M::T, typename M::E, decltype(&M::mp), decltype(&M::cp)>>: std::true_type {};
template <class M, class= void> struct commute: std::false_type {};
template <class M> struct commute<M, std::void_t<typename M::commute>>: std::true_type {};
template <class M, class= void> struct value_of {
 using type= M;
};
template <class M> struct value_of<M, std::void_t<typename M::T>> {
 using type= typename M::T;
};
template <class M, class= void> struct action_of {
 using type= char;
};
template <class M> struct action_of<M, std::void_t<typename M::E>> {
 using type= typename M::E;
};
// 1 MB 以上の確保は 2 MB 境界に寄せ、Linux では huge page を頼む。小さい確保も 64 byte 境界にする (塊を AVX2 で読むため)。
template <class T> struct Alloc {
 using value_type= T;
 Alloc()= default;
 template <class U> Alloc(const Alloc<U>&) {}
 T* allocate(std::size_t n) {
  const std::size_t want= n * sizeof(T);
  if(want < (std::size_t(1) << 20)) return static_cast<T*>(::operator new(want, std::align_val_t(64)));
  const std::size_t H= std::size_t(1) << 21, bytes= (want + H - 1) & ~(H - 1);
  T* p= static_cast<T*>(::operator new(bytes, std::align_val_t(H)));
#ifdef __linux__
  madvise(p, bytes, MADV_HUGEPAGE);
#endif
  return p;
 }
 void deallocate(T* p, std::size_t n) {
  const std::size_t want= n * sizeof(T);
  if(want < (std::size_t(1) << 20)) ::operator delete(p, std::align_val_t(64));
  else ::operator delete(p, std::align_val_t(std::size_t(1) << 21));
 }
 template <class U> bool operator==(const Alloc<U>&) const { return true; }
};
// ::monoid::AffineSum<ZMod<P>> (P は定数) かどうか (この名前空間の monoid は、モノイドかを見る判定)。
template <class M> struct fast_affine: std::false_type {};
template <unsigned P> struct fast_affine<::monoid::AffineSum<ZMod<P>>>: std::bool_constant<P != 0> {};
using u32= unsigned;
using u64= unsigned long long;
// 塊の中の和と 1 次関数の作用を AVX2 で 8 個ずつ行う核。値はそのままの値を 2P 未満で持つ (ZMod と同じ)。
template <unsigned P> struct Kern {
 static constexpr u32 P2= 2 * P;
 // x -> m x + c。m と c は P 未満で、q は Shoup の商の近似 floor(m 2^32 / P)。
 struct Act {
  __m256i m, q, c;
 };
 static Act act(const std::array<ZMod<P>, 2>& f) {
  const u32 w= f[0].val();
  return {_mm256_set1_epi32(int(w)), _mm256_set1_epi32(int(u32((u64(w) << 32) / P))), _mm256_set1_epi32(int(f[1].val()))};
 }
 // x m mod P を 2P 未満で。m と q はどのレーンも同じ値であること (奇数のレーンの積にも偶数のレーンの q を使う)。
 [[gnu::always_inline]] static __m256i mul(__m256i x, __m256i m, __m256i q) {
  const __m256i e= _mm256_mul_epu32(x, q), o= _mm256_mul_epu32(_mm256_srli_epi64(x, 32), q);
  const __m256i h= _mm256_blend_epi32(_mm256_srli_epi64(e, 32), o, 0xaa);
  return _mm256_sub_epi32(_mm256_mullo_epi32(x, m), _mm256_mullo_epi32(h, _mm256_set1_epi32(int(P))));
 }
 // mul で、m と q がレーンごとに違ってよい形。
 [[gnu::always_inline]] static __m256i mul_lanes(__m256i x, __m256i m, __m256i q) {
  const __m256i e= _mm256_mul_epu32(x, q), o= _mm256_mul_epu32(_mm256_srli_epi64(x, 32), _mm256_srli_epi64(q, 32));
  const __m256i h= _mm256_blend_epi32(_mm256_srli_epi64(e, 32), o, 0xaa);
  return _mm256_sub_epi32(_mm256_mullo_epi32(x, m), _mm256_mullo_epi32(h, _mm256_set1_epi32(int(P))));
 }
 [[gnu::always_inline]] static __m256i add(__m256i a, __m256i b) {
  const __m256i s= _mm256_add_epi32(a, b);
  return _mm256_min_epu32(s, _mm256_sub_epi32(s, _mm256_set1_epi32(int(P2))));
 }
 static u32 hsum(__m256i s) {
  const __m128i p2= _mm_set1_epi32(int(P2));
  __m128i t= _mm_add_epi32(_mm256_castsi256_si128(s), _mm256_extracti128_si256(s, 1));
  t= _mm_min_epu32(t, _mm_sub_epi32(t, p2));
  t= _mm_add_epi32(t, _mm_shuffle_epi32(t, 0x4e));
  t= _mm_min_epu32(t, _mm_sub_epi32(t, p2));
  t= _mm_add_epi32(t, _mm_shuffle_epi32(t, 0xb1));
  t= _mm_min_epu32(t, _mm_sub_epi32(t, p2));
  return u32(_mm_cvtsi128_si32(t));
 }
 [[gnu::always_inline]] static __m256i iota8() { return _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7); }
 [[gnu::always_inline]] static __m256i ld(const u32* p) { return _mm256_load_si256((const __m256i*)p); }
 [[gnu::always_inline]] static void st(u32* p, __m256i x) { _mm256_store_si256((__m256i*)p, x); }
 // a[lo, hi) の和を 2P 未満で。
 static u32 sum_range(const u32* a, int lo, int hi) {
  if(lo >= hi) return 0;
  const __m256i L= _mm256_set1_epi32(lo - 1), H= _mm256_set1_epi32(hi), step= _mm256_set1_epi32(8);
  int k= lo >> 3;
  const int ke= (hi + 7) >> 3;
  __m256i idx= _mm256_add_epi32(iota8(), _mm256_set1_epi32(k << 3)), s0= _mm256_setzero_si256(), s1= s0;
  for(; k + 1 < ke; k+= 2) {
   const __m256i i1= _mm256_add_epi32(idx, step);
   const __m256i k0= _mm256_and_si256(_mm256_cmpgt_epi32(idx, L), _mm256_cmpgt_epi32(H, idx));
   const __m256i k1= _mm256_and_si256(_mm256_cmpgt_epi32(i1, L), _mm256_cmpgt_epi32(H, i1));
   s0= add(s0, _mm256_and_si256(ld(a + (k << 3)), k0));
   s1= add(s1, _mm256_and_si256(ld(a + (k << 3) + 8), k1));
   idx= _mm256_add_epi32(i1, step);
  }
  if(k < ke) s0= add(s0, _mm256_and_si256(ld(a + (k << 3)), _mm256_and_si256(_mm256_cmpgt_epi32(idx, L), _mm256_cmpgt_epi32(H, idx))));
  return hsum(add(s0, s1));
 }
 // a[0, len) に f を当てる。
 static void transform_all(u32* a, int len, const Act& f) {
  for(int k= 0; k < len; k+= 8) st(a + k, add(mul(ld(a + k), f.m, f.q), f.c));
 }
 // a[lo, hi) に f を当て、当てる前の和を 2P 未満で返す。
 static u32 transform_range(u32* a, int lo, int hi, const Act& f) {
  const __m256i L= _mm256_set1_epi32(lo - 1), H= _mm256_set1_epi32(hi), step= _mm256_set1_epi32(8);
  int k= lo >> 3;
  const int ke= (hi + 7) >> 3;
  __m256i idx= _mm256_add_epi32(iota8(), _mm256_set1_epi32(k << 3)), s= _mm256_setzero_si256();
  for(; k < ke; ++k) {
   const __m256i x= ld(a + (k << 3));
   const __m256i msk= _mm256_and_si256(_mm256_cmpgt_epi32(idx, L), _mm256_cmpgt_epi32(H, idx));
   st(a + (k << 3), _mm256_blendv_epi8(x, add(mul(x, f.m, f.q), f.c), msk));
   s= add(s, _mm256_and_si256(x, msk));
   idx= _mm256_add_epi32(idx, step);
  }
  return hsum(s);
 }
 // a[0, len) のうち [lo, hi) に in を、ほかに out を当て、当てたあとの和を 2P 未満で返す。
 static u32 transform2(u32* a, int len, int lo, int hi, const Act& in, const Act& out) {
  const __m256i L= _mm256_set1_epi32(lo - 1), H= _mm256_set1_epi32(hi), N= _mm256_set1_epi32(len), step= _mm256_set1_epi32(8);
  __m256i idx= iota8(), s= _mm256_setzero_si256();
  for(int k= 0; k < len; k+= 8) {
   const __m256i msk= _mm256_and_si256(_mm256_cmpgt_epi32(idx, L), _mm256_cmpgt_epi32(H, idx));
   const __m256i m= _mm256_blendv_epi8(out.m, in.m, msk), q= _mm256_blendv_epi8(out.q, in.q, msk), c= _mm256_blendv_epi8(out.c, in.c, msk);
   const __m256i y= add(mul_lanes(ld(a + k), m, q), c);
   st(a + k, y);
   s= add(s, _mm256_and_si256(y, _mm256_cmpgt_epi32(N, idx)));
   idx= _mm256_add_epi32(idx, step);
  }
  return hsum(s);
 }
};
template <class M, bool FA= fast_affine<M>::value> struct kern_of {
 using type= void;
};
template <class M> struct kern_of<M, true> {
 using type= Kern<M::T::mod()>;
};
}  // namespace sequence_internal
template <class M, bool REV= false, int B= 0> class Sequence {
 using T= typename sequence_internal::value_of<M>::type;
 using E= typename sequence_internal::action_of<M>::type;
 static constexpr bool MO= sequence_internal::monoid<M>::value, AC= sequence_internal::action<M>::value;
 static constexpr bool FA= sequence_internal::fast_affine<M>::value;  // 塊の中を AVX2 の核で回す
 using K= typename sequence_internal::kern_of<M>::type;
 static constexpr bool RS= MO && REV && !(sequence_internal::commute<M>::value || FA);  // 逆順の積も持つ
 static constexpr int CAP= B ? B : FA ? 256 : std::clamp<int>(int(512 / sizeof(T)) / 16 * 16, 16, 128), MIN= CAP / 4;
 static_assert(!FA || CAP % 8 == 0, "AVX2 の核は塊の大きさが 8 の倍数であること");
 using S= std::conditional_t<MO, T, char>;
 using SR= std::conditional_t<RS, T, char>;
 using EA= std::conditional_t<AC, E, char>;
 struct Node {
  int l, r, cnt, len, fr, bk;  // fr と bk は部分木の左端と右端の塊の長さ
  unsigned pri;
  bool rev, lz, crev, clz;  // rev と lz は子へ配る分、crev と clz は自分の塊の分
  S sum, cf;  // sum は部分木の積、cf は塊を配列の順に畳んだ積 (どちらも作用を当てたあとの値)
  SR rsum, cb;  // rsum は部分木の逆順の積、cb は塊を配列の逆順に畳んだ積
  EA e, ce;  // e は子へ配る作用、ce は自分の塊へ当てる作用
 };
 static Node blank(unsigned pri) {
  Node n{};
  n.pri= pri;
  if constexpr(MO) n.sum= n.cf= M::ti();
  if constexpr(RS) n.rsum= n.cb= M::ti();
  return n;
 }
 template <class U> using Vec= std::vector<U, sequence_internal::Alloc<U>>;
 static inline Vec<Node> nd= Vec<Node>(1, blank(0));  // 0 は空の木
 static inline Vec<T> dat;  // 節 x の塊は dat[x CAP, x CAP + len)
 static inline std::vector<int> freed;
 // 乱数の種は起動の時刻から取る (決め打ちの種は、公開したライブラリでは hack で狙える)。
 static unsigned seed() {
  const auto t= (unsigned long long)std::chrono::steady_clock::now().time_since_epoch().count();
  const unsigned v= unsigned(t ^ (t >> 32)) ^ 2463534242u;
  return v ? v : 2463534242u;
 }
 static inline unsigned rs= seed();
 static inline std::vector<int> path, path2;
 static inline std::vector<T> tu, tv;  // 反転で入れ替える塊の外の部分
 static unsigned rnd() { return rs^= rs << 13, rs^= rs >> 17, rs^= rs << 5; }
 static T* arr(int x) { return dat.data() + std::size_t(x) * CAP; }
 // ZMod は u32 を 1 つ持つ標準レイアウトの型なので、塊を u32 の列として読み書きする。
 static sequence_internal::u32* raw(T* a) { return reinterpret_cast<sequence_internal::u32*>(a); }
 static const sequence_internal::u32* raw(const T* a) { return reinterpret_cast<const sequence_internal::u32*>(a); }
 static int make(unsigned pri) {
  if(!freed.empty()) {
   const int x= freed.back();
   freed.pop_back();
   nd[x]= blank(pri);
   return x;
  }
  nd.push_back(blank(pri));
  if(dat.size() < nd.size() * CAP) dat.resize(std::max(nd.size() * CAP, dat.size() * 2));
  return (int)nd.size() - 1;
 }
 static void release(int x) { freed.push_back(x); }
 static void map1(T& v, const E& e, int len) {
  if constexpr(std::is_invocable_v<decltype(M::mp), T&, const E&, int>) M::mp(v, e, len);
  else M::mp(v, e);
 }
 // 配列の [lo, hi) を順に (rv なら逆順に) 畳む。
 static T fold(const T* a, int lo, int hi, bool rv) {
  if constexpr(FA) return T::raw(K::sum_range(raw(a), lo, hi));
  T s= M::ti();
  if(RS && rv)
   for(int i= hi; i-- > lo;) s= M::op(s, a[i]);
  else
   for(int i= lo; i < hi; ++i) s= M::op(s, a[i]);
  return s;
 }
 // 塊の列の順の積と、その逆順の積。
 static const S& lsum(const Node& n) {
  if constexpr(RS) return n.crev ? n.cb : n.cf;
  else return n.cf;
 }
 static const SR& lrsum(const Node& n) {
  if constexpr(RS) return n.crev ? n.cf : n.cb;
  else return n.cb;
 }
 // 塊の積を求め直す。塊の遅延が残っていれば、畳んだ積に当てる。
 static void recalc(int x) {
  if constexpr(MO) {
   Node& n= nd[x];
   const T* a= arr(x);
   n.cf= fold(a, 0, n.len, false);
   if constexpr(RS) n.cb= fold(a, 0, n.len, true);
   if constexpr(AC)
    if(n.clz) {
     map1(n.cf, n.ce, n.len);
     if constexpr(RS) map1(n.cb, n.ce, n.len);
    }
  }
 }
 static void pull(int x) {
  Node& n= nd[x];
  const Node &L= nd[n.l], &R= nd[n.r];
  n.cnt= L.cnt + n.len + R.cnt;
  n.fr= n.l ? L.fr : n.len, n.bk= n.r ? R.bk : n.len;
  if constexpr(MO) n.sum= M::op(M::op(L.sum, lsum(n)), R.sum);
  if constexpr(RS) n.rsum= M::op(M::op(R.rsum, lrsum(n)), L.rsum);
 }
 static void aff(int x, const E& e) {
  Node& n= nd[x];
  if constexpr(MO) {
   map1(n.sum, e, n.cnt), map1(n.cf, e, n.len);
   if constexpr(RS) map1(n.rsum, e, n.cnt), map1(n.cb, e, n.len);
  }
  if(n.lz) M::cp(n.e, e);
  else n.e= e, n.lz= true;
  if(n.clz) M::cp(n.ce, e);
  else n.ce= e, n.clz= true;
 }
 static void tog(int x) {
  Node& n= nd[x];
  std::swap(n.l, n.r), std::swap(n.fr, n.bk);
  if constexpr(RS) std::swap(n.sum, n.rsum);
  n.rev^= 1, n.crev^= 1;
 }
 static void push(int x) {
  Node& n= nd[x];
  if constexpr(REV)
   if(n.rev) {
    if(n.l) tog(n.l);
    if(n.r) tog(n.r);
    n.rev= false;
   }
  if constexpr(AC)
   if(n.lz) {
    const E e= n.e;
    if(n.l) aff(n.l, e);
    if(nd[x].r) aff(nd[x].r, e);
    nd[x].lz= false;
   }
 }
 // 塊の遅延の作用を配列に当てる (反転の印は残す)。
 static void materialize(int x) {
  if constexpr(AC) {
   Node& n= nd[x];
   if(!n.clz) return;
   T* a= arr(x);
   if constexpr(FA) K::transform_all(raw(a), n.len, K::act(n.ce));
   else
    for(int i= 0; i < n.len; ++i) map1(a[i], n.ce, 1);
   n.clz= false;
  }
 }
 // 作用を当てて、配列を列の順に並べる。
 static void forward(int x) {
  materialize(x);
  if constexpr(REV) {
   Node& n= nd[x];
   if(n.crev) {
    std::reverse(arr(x), arr(x) + n.len), n.crev= false;
    if constexpr(RS) std::swap(n.cf, n.cb);
   }
  }
 }
 static int phys(const Node& n, int o) { return REV && n.crev ? n.len - 1 - o : o; }
 // 塊の中の列の順の [cl, cr) の積。
 static T part_sum(int x, int cl, int cr) {
  const Node& n= nd[x];
  const bool rv= REV && n.crev;
  const int pl= rv ? n.len - cr : cl, pr= rv ? n.len - cl : cr;
  if constexpr(FA) {
   // 部分が塊の半分より長ければ、残りの部分の和を塊の和から引く (和は引ける)。
   const int w= pr - pl;
   const auto* a= raw(arr(x));
   if(2 * w <= n.len) {
    const T s= T::raw(K::sum_range(a, pl, pr));
    return n.clz ? n.ce[0] * s + n.ce[1] * T(w) : s;
   }
   const T s= T::raw(K::sum_range(a, 0, pl)) + T::raw(K::sum_range(a, pr, n.len));
   return n.cf - (n.clz ? n.ce[0] * s + n.ce[1] * T(n.len - w) : s);
  }
  T s= fold(arr(x), pl, pr, rv);
  if constexpr(AC)
   if(n.clz) map1(s, n.ce, cr - cl);
  return s;
 }
 // 塊の中の列の順の [cl, cr) (塊の全体ではない) に e を当てる。
 static void part_apply(int x, int cl, int cr, const E& e) {
  if constexpr(FA) {
   Node& n= nd[x];
   const bool rv= REV && n.crev;
   const int len= n.len, pl= rv ? len - cr : cl, pr= rv ? len - cl : cr;
   auto* a= raw(arr(x));
   if(n.clz) {  // 残っている塊の遅延と合わせて 1 回で書く (区間の内は遅延のあとに e、外は遅延だけ)
    const E in{e[0] * n.ce[0], e[0] * n.ce[1] + e[1]};
    n.cf= T::raw(K::transform2(a, len, pl, pr, K::act(in), K::act(n.ce)));
    n.clz= false;
   } else {
    const T old= T::raw(K::transform_range(a, pl, pr, K::act(e)));
    n.cf+= (e[0] - T(1)) * old + e[1] * T(pr - pl);
   }
   return;
  }
  materialize(x);
  const Node& n= nd[x];
  const bool rv= REV && n.crev;
  const int pl= rv ? n.len - cr : cl, pr= rv ? n.len - cl : cr;
  T* a= arr(x);
  for(int i= pl; i < pr; ++i) map1(a[i], e, 1);
  recalc(x);
 }
 // 先頭の k 個を a に、残りを b に分ける。k が塊の中に落ちたら、塊を 2 つの節に割る。
 static void split(int x, int k, int& a, int& b) {
  if(!x) return void(a= b= 0);
  push(x);
  const int ls= nd[nd[x].l].cnt, len= nd[x].len;
  if(k <= ls) {
   int p, q;
   split(nd[x].l, k, p, q);
   nd[x].l= q, pull(x), a= p, b= x;
  } else if(k >= ls + len) {
   int p, q;
   split(nd[x].r, k - ls - len, p, q);
   nd[x].r= p, pull(x), a= x, b= q;
  } else {
   const int o= k - ls;
   forward(x);
   const int y= make(rnd());
   std::copy(arr(x) + o, arr(x) + len, arr(y));
   nd[y].len= len - o, nd[x].len= o;
   recalc(x), recalc(y), pull(y);
   const int r= nd[x].r;
   nd[x].r= 0, pull(x);
   a= x, b= merge_(y, r);
  }
 }
 static int merge_(int a, int b) {
  if(!a || !b) return a | b;
  if(nd[a].pri > nd[b].pri) {
   push(a);
   const int r= merge_(nd[a].r, b);
   nd[a].r= r, pull(a);
   return a;
  }
  push(b);
  const int l= merge_(a, nd[b].l);
  nd[b].l= l, pull(b);
  return b;
 }
 // b の先頭の塊の節を外して返す (塊は列の順に並べてある)。b には残りを入れる。
 static int take_first(int& b) {
  path.clear();
  int x= b;
  for(;;) {
   push(x);
   if(!nd[x].l) break;
   path.push_back(x), x= nd[x].l;
  }
  const int r= nd[x].r;
  if(path.empty()) b= r;
  else nd[path.back()].l= r;
  for(int k= (int)path.size(); k--;) pull(path[k]);
  forward(x);
  return x;
 }
 // a の末尾の塊に、節 f の塊 (列の順に並べてある) を足して、f を返す。
 static void absorb_last(int a, int f) {
  path.clear();
  int x= a;
  for(;;) {
   push(x), path.push_back(x);
   if(!nd[x].r) break;
   x= nd[x].r;
  }
  forward(x);
  std::copy(arr(f), arr(f) + nd[f].len, arr(x) + nd[x].len);
  nd[x].len+= nd[f].len;
  recalc(x);
  for(int k= (int)path.size(); k--;) pull(path[k]);
  release(f);
 }
 // 併合。境目の塊のどちらかが短く、2 つが 1 つの塊に収まれば、まとめてからつなぐ。
 static int merge(int a, int b) {
  if(!a || !b) return a | b;
  if((nd[a].bk < MIN || nd[b].fr < MIN) && nd[a].bk + nd[b].fr <= CAP) {
   const int f= take_first(b);
   absorb_last(a, f);
   if(!b) return a;
  }
  return merge_(a, b);
 }
 static T prod(int x, int l, int r) {
  if(!x || l >= r) return M::ti();
  if(l == 0 && r == nd[x].cnt) return nd[x].sum;
  push(x);
  const int ls= nd[nd[x].l].cnt, len= nd[x].len;
  T s= M::ti();
  if(l < ls) s= prod(nd[x].l, l, std::min(r, ls));
  const int cl= std::max(l - ls, 0), cr= std::min(r - ls, len);
  if(cl < cr) s= M::op(s, cl == 0 && cr == len ? lsum(nd[x]) : part_sum(x, cl, cr));
  if(r > ls + len) s= M::op(s, prod(nd[x].r, std::max(l - ls - len, 0), r - ls - len));
  return s;
 }
 static void apply(int x, int l, int r, const E& e) {
  if(!x || l >= r) return;
  if(l == 0 && r == nd[x].cnt) return aff(x, e);
  push(x);
  const int ls= nd[nd[x].l].cnt, len= nd[x].len;
  if(l < ls) apply(nd[x].l, l, std::min(r, ls), e);
  const int cl= std::max(l - ls, 0), cr= std::min(r - ls, len);
  if(cl < cr) {
   if(cl == 0 && cr == len) {
    Node& n= nd[x];
    if constexpr(MO) {
     map1(n.cf, e, len);
     if constexpr(RS) map1(n.cb, e, len);
    }
    if(n.clz) M::cp(n.ce, e);
    else n.ce= e, n.clz= true;
   } else part_apply(x, cl, cr, e);
  }
  if(r > ls + len) apply(nd[x].r, std::max(l - ls - len, 0), r - ls - len, e);
  pull(x);
 }
 // 根 rt から位置 i を含む塊まで、遅延を配りながら下りる。end が真なら塊の後ろの端も含める。o に塊の中の位置を入れる。
 static int descend(int rt, int i, bool end, std::vector<int>& ps, int& o) {
  ps.clear();
  int x= rt;
  o= i;
  for(;;) {
   push(x), ps.push_back(x);
   const int ls= nd[nd[x].l].cnt;
   if(o < ls) x= nd[x].l;
   else if(o < ls + nd[x].len + (end ? 1 : 0)) {
    o-= ls;
    return x;
   } else o-= ls + nd[x].len, x= nd[x].r;
  }
 }
 static void pull_path(const std::vector<int>& ps) {
  for(int k= (int)ps.size(); k--;) pull(ps[k]);
 }
 // 塊の中の列の順の位置 o に v を入れる (空きがあること)。
 static void put(int x, int o, const T& v) {
  materialize(x);
  Node& n= nd[x];
  T* a= arr(x);
  const int p= REV && n.crev ? n.len - o : o;
  std::move_backward(a + p, a + n.len, a + n.len + 1);
  a[p]= v, ++n.len;
  if constexpr(MO) {
   if constexpr(RS) recalc(x);
   else n.cf= M::op(n.cf, v);  // 可換なので順は問わない
  }
 }
 static void collect(int x, std::vector<T>& out) {
  if(!x) return;
  push(x), forward(x);
  collect(nd[x].l, out);
  out.insert(out.end(), arr(x), arr(x) + nd[x].len);
  collect(nd[x].r, out);
 }
 template <class F> static int build(const F& at, int n, int lo, int hi, int depth) {  // 塊の番号の [lo, hi)
  if(lo >= hi) return 0;
  constexpr int H= CAP * 3 / 4;
  const int m= (lo + hi) / 2;
  const int x= make(~0u - (unsigned(depth) << 24) - (rnd() & 0xffffffu));  // 浅いほど優先度が高い
  const int s= m * H, e= std::min(s + H, n);
  T* a= arr(x);
  for(int i= s; i < e; ++i) a[i - s]= at(i);
  nd[x].len= e - s, recalc(x);
  const int l= build(at, n, lo, m, depth + 1), r= build(at, n, m + 1, hi, depth + 1);
  nd[x].l= l, nd[x].r= r, pull(x);
  return x;
 }
 template <class F> static int build(const F& at, int n) {
  constexpr int H= CAP * 3 / 4;
  return build(at, n, 0, (n + H - 1) / H, 0);
 }
 int root= 0;
 static Sequence of(int rt) {
  Sequence s;
  s.root= rt;
  return s;
 }
public:
 Sequence()= default;
 explicit Sequence(int n, const T& v= T()) { root= build([&](int) { return v; }, n); }
 Sequence(const std::vector<T>& v) { root= build([&](int i) { return v[i]; }, (int)v.size()); }
 template <class It, std::enable_if_t<!std::is_integral_v<It>, int> = 0> Sequence(It bg, It ed) {
  const std::vector<T> v(bg, ed);
  root= build([&](int i) { return v[i]; }, (int)v.size());
 }
 int size() const { return nd[root].cnt; }
 bool empty() const { return !root; }
 void clear() { root= 0; }
 std::pair<Sequence, Sequence> split(int k) {
  int a, b;
  split(root, k, a, b);
  return {of(a), of(b)};
 }
 std::tuple<Sequence, Sequence, Sequence> split3(int l, int r) {
  int a, b, c;
  split(root, r, b, c);
  split(b, l, a, b);
  return {of(a), of(b), of(c)};
 }
 Sequence& operator+=(Sequence rhs) { return root= merge(root, rhs.root), *this; }
 Sequence operator+(Sequence rhs) const { return of(merge(root, rhs.root)); }
 void insert(int k, const T& v) {
  if(!root) {
   root= make(rnd());
   arr(root)[0]= v, nd[root].len= 1, recalc(root), pull(root);
   return;
  }
  int o;
  const int x= descend(root, k, true, path2, o);
  if(nd[x].len < CAP) {
   put(x, o, v);
   pull_path(path2);
   return;
  }
  // 満杯なら、前半を x に残して後半を新しい節 y に移し、y を x のすぐ後ろへつなぐ。
  forward(x);
  constexpr int H= CAP / 2;
  const int y= make(rnd());
  std::copy(arr(x) + H, arr(x) + CAP, arr(y));
  nd[y].len= CAP - H, nd[x].len= H;
  recalc(x), recalc(y);
  if(o <= H) put(x, o, v);
  else put(y, o - H, v);
  pull(y), pull_path(path2);
  int a, b;
  split(root, k - o + nd[x].len, a, b);
  root= merge_(merge_(a, y), b);
 }
 void push_back(const T& v) { insert(size(), v); }
 void push_front(const T& v) { insert(0, v); }
 T erase(int k) {
  int o;
  const int x= descend(root, k, false, path2, o);
  Node& n= nd[x];
  T* a= arr(x);
  const int p= phys(n, o);
  T v= a[p];
  if constexpr(AC)
   if(n.clz) map1(v, n.ce, 1);
  std::move(a + p + 1, a + n.len, a + p);
  --n.len;
  if(!n.len) {  // 空になった節を外す (x の遅延は descend で配ってある)
   const int m= merge_(n.l, n.r);
   path2.pop_back();
   if(path2.empty()) root= m;
   else {
    const int par= path2.back();
    (nd[par].l == x ? nd[par].l : nd[par].r)= m;
   }
   release(x);
  } else recalc(x);
  pull_path(path2);
  return v;
 }
 T pop_back() { return erase(size() - 1); }
 T pop_front() { return erase(0); }
 T get(int k) {
  int o;
  const int x= descend(root, k, false, path2, o);
  T v= arr(x)[phys(nd[x], o)];
  if constexpr(AC)
   if(nd[x].clz) map1(v, nd[x].ce, 1);
  return v;
 }
 T operator[](int k) { return get(k); }
 void set(int k, const T& v) {
  int o;
  const int x= descend(root, k, false, path2, o);
  materialize(x);
  arr(x)[phys(nd[x], o)]= v;
  recalc(x);
  pull_path(path2);
 }
 void mul(int k, const T& v) {
  static_assert(MO, "mul はモノイドがあるときだけ");
  set(k, M::op(get(k), v));
 }
 T prod(int l, int r) {
  static_assert(MO, "prod はモノイドがあるときだけ");
  return prod(root, l, r);
 }
 T all_prod() const {
  static_assert(MO, "all_prod はモノイドがあるときだけ");
  return nd[root].sum;
 }
 void apply(int l, int r, const E& e) {
  static_assert(AC, "apply は作用があるときだけ");
  apply(root, l, r, e);
 }
 void reverse() {
  static_assert(REV, "reverse は REV のときだけ");
  if(root) tog(root);
 }
 void reverse(int l, int r) {
  static_assert(REV, "reverse は REV のときだけ");
  if(r - l < 2) return;
  int i, jj;
  const int u= descend(root, l, false, path, i);
  const int v= descend(root, r - 1, false, path2, jj);
  const int j= jj + 1;
  if(u == v) {  // 塊の中だけの反転
   const Node& n= nd[u];
   const int pl= n.crev ? n.len - j : i, pr= n.crev ? n.len - i : j;
   std::reverse(arr(u) + pl, arr(u) + pr);
   if constexpr(RS) recalc(u), pull_path(path);
   return;
  }
  const int nu= nd[u].len, nv= nd[v].len, tail= nv - j;
  if(tail + nu - i > CAP || j + i > CAP) {
   int a, b, c;
   split(root, l, a, b);
   split(b, r - l, b, c);
   tog(b);
   root= merge(merge(a, b), c);
   return;
  }
  if(i || tail) {
   // u の [0, i) と v の [j, nv) を逆順にして入れ替える: u = rev(v[j, nv)) + u[i, nu)、v = v[0, j) + rev(u[0, i))。
   forward(u), forward(v);
   T *U= arr(u), *V= arr(v);
   tu.assign(U, U + i), tv.assign(V + j, V + nv);
   if(tail > i) std::move_backward(U + i, U + nu, U + tail + nu - i);
   else std::move(U + i, U + nu, U + tail);
   for(int k= 0; k < tail; ++k) U[k]= tv[tail - 1 - k];
   for(int k= 0; k < i; ++k) V[j + k]= tu[i - 1 - k];
   nd[u].len= tail + nu - i, nd[v].len= j + i;
   recalc(u), recalc(v);
   pull_path(path2), pull_path(path);
  }
  // 塊の境目 (u の先頭と v の末尾) で割って、間の木に反転の印を付ける。
  int a, b, c;
  split(root, l - i, a, b);
  split(b, r + tail - (l - i), b, c);
  tog(b);
  root= merge_(merge_(a, b), c);
 }
 std::vector<T> dump() {
  std::vector<T> out;
  out.reserve(size());
  collect(root, out);
  return out;
 }
};
