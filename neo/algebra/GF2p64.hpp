#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#include <simde/x86/clmul.h>
#endif
#include <tuple>
#include <iostream>
#include <cassert>
namespace gf2p64_internal {
using u64= unsigned long long;
using u32= unsigned;
using u16= unsigned short;
using u8= unsigned char;
template <class T, int N> struct Arr {
 T t[N];
};
template <class U, int N> constexpr void lin(U (&t)[N][256], const U* b) {
 for(int i= 0; i < N * 8; ++i)
  for(u32 h= 1 << (i & 7), j= h; j--;) t[i >> 3][h | j]= t[i >> 3][j] ^ b[i];
}
template <class U, int N> struct Lin {
 U t[N][256];
 constexpr Lin(const U* b): t{} { lin(t, b); }
 constexpr U operator()(u64 x) const {
  if constexpr(N == 2) return t[0][u8(x)] ^ t[1][u8(x >> 8)];
  else if constexpr(N == 4) return t[0][u8(x)] ^ t[1][u8(x >> 8)] ^ t[2][u8(x >> 16)] ^ t[3][u8(x >> 24)];
  else return t[0][u8(x)] ^ t[1][u8(x >> 8)] ^ t[2][u8(x >> 16)] ^ t[3][u8(x >> 24)] ^ t[4][u8(x >> 32)] ^ t[5][u8(x >> 40)] ^ t[6][u8(x >> 48)] ^ t[7][u8(x >> 56)];
 }
 constexpr Lin operator*(const Lin& r) const {
  U b[N * 8]{};
  for(int i= 0; i < N * 8; ++i) b[i]= (*this)(r.t[i >> 3][1 << (i & 7)]);
  return Lin(b);
 }
};
using LinMap= Lin<u64, 8>;
constexpr u64 cmul(u64 a, u64 b) {
 u64 r= 0;
 for(int i= 64; i--;) r= r << 1 ^ (0x1b & -(r >> 63)) ^ (a & -(b >> i & 1));
 return r;
}
constexpr u64 cpow(u64 a, u64 e) {
 u64 r= 1;
 for(; e; e>>= 1, a= cmul(a, a)) r= e & 1 ? cmul(r, a) : r;
 return r;
}
template <class T, int N> constexpr Arr<T, N> inverse(const T* c) {
 Arr<T, N> a{}, r{};
 for(int i= 0; i < N; ++i) a.t[i]= c[i], r.t[i]= T(1) << i;
 for(int j= 0; j < N; ++j) {
  int p= j;
  while(!(a.t[p] >> j & 1)) ++p;
  if(p != j) a.t[j]^= a.t[p], r.t[j]^= r.t[p];
  for(int i= 0; i < N; ++i)
   if(i != j && a.t[i] >> j & 1) a.t[i]^= a.t[j], r.t[i]^= r.t[j];
 }
 return r;
}
constexpr u64 TO_NIM_B[64]= {0x0000000000000001, 0x5211145c804b6109, 0x7c8bc2cad259879f, 0x565854b4c60c1e0b, 0x4068acf7104c20c3, 0x662d2bd0f2739155, 0x7a90c83701fa8323, 0x21cfa750247e8755, 0x67d1044e545abf47, 0x4d9d3b5a8568f839, 0x567a9d7331b6b3c6, 0x1ca54bfdd6d1ae59, 0x454fa483275db25c, 0x6766df6fec4e9d44, 0x35cb621cec1fe7f9, 0x4c606d3e52faf263, 0x57640dc825a57954, 0x7aca87838b7f6315, 0x6d53c884ebf2b0ed, 0x3721d998bb50164b, 0x7aa7c62fd6cd53ab, 0x47cbb2c51f7c040f, 0x132063b7f5e42489, 0x0c1b36c8b2993f8a, 0x60119ecff680497a, 0x5175da444cc11791, 0x5792ff4554765b09, 0x0c9fdb8a01334e82, 0x2be0a763a68a4725, 0x3c2dc8260ad051f6, 0x6c4c9fed8816bb9c, 0x630062753ffaf766,
                             0x7b37d31b5d519225, 0x2364f7f79705691c, 0x453eb8a83e2fec71, 0x7c0121b37e828666, 0x59190d3250e66011, 0x103207f9dda18cae, 0x28233dce01c69b76, 0x4fa519899227a5e7, 0x4567ba46ee7bc6cd, 0x0a284773d021afd5, 0x63894079bbe3a824, 0x11013c7fdfaaa5c2, 0x1aa984f18574f3b0, 0x0cbaba126fd0c4db, 0x0b8797719e6dc725, 0x4a2845680aefaa72, 0x536d2535f6934e15, 0x01db7a57effcd689, 0x7e1ed0ad01e2a5ad, 0x0aedc9b3cee826f6, 0x7ba716eccf9f68e1, 0x5d5e23bc0f3dc38f, 0x0b5f2a3b88674d83, 0x2de9bafc2f00f8d4, 0x3b56712ad419c7e0, 0x3ab4be8c30c19253, 0x2708522ffaa654b0, 0x2b8bca57bf643598, 0x588825d1a5fa8e1c, 0x86adf8bf4d45962f, 0x51b4c15d8719dd73, 0xe4a2b3b59783d0aa};
constexpr LinMap TO_NIM= LinMap(TO_NIM_B), FROM_NIM= LinMap(inverse<u64, 64>(TO_NIM_B).t);
template <int K, int D= 0> constexpr LinMap F= F<K / 2, D> * F<K - K / 2, D>;
template <int D>
constexpr LinMap F<1, D> = []() {
 u64 g[64]{};
 for(int i= 0; i < 64; ++i) g[i]= cmul(u64(1) << i, u64(1) << i);
 return LinMap(g);
}();
inline u64 mul(u64 a, u64 b) {
 static constexpr u8 RED[]= {0, 27, 45, 54, 90, 65, 119, 108};
 __m128i v= _mm_clmulepi64_si128(_mm_cvtsi64_si128(a), _mm_cvtsi64_si128(b), 0);
 u64 h= v[1], d= h ^ (h << 1);
 return v[0] ^ RED[h >> 60] ^ d ^ (d << 3);
}
inline u64 sq(u64 a) {
 static constexpr u8 RED_SQ[4]= {0, 27, 90, 65};
 const __m128i MASK_LO= _mm_set1_epi8(0x0f);
 const __m128i SPR= _mm_set_epi64x(0xfffcf3f0cfccc3c0, 0x3f3c33300f0c0300);
 __m128i v= _mm_set_epi64x(0, a);
 __m128i x= _mm_shuffle_epi8(SPR, _mm_and_si128(_mm_unpacklo_epi8(v, _mm_srli_epi16(v, 4)), MASK_LO));
 u64 d= x[1];
 return (x[0] & 0x5555555555555555) ^ RED_SQ[a >> 62] ^ d ^ (d << 3);
}
inline __m256i sq2(const __m256i& v) {
 const __m256i MASK_LO= _mm256_set1_epi8(0x0f), EVEN= _mm256_set1_epi64x(0x5555555555555555);
 const __m256i SPR= _mm256_setr_epi64x(0x3f3c33300f0c0300, 0xfffcf3f0cfccc3c0, 0x3f3c33300f0c0300, 0xfffcf3f0cfccc3c0);
 const __m256i RED_SQ= _mm256_setr_epi64x(0x415a1b00, 0, 0x415a1b00, 0);
 const __m256i x= _mm256_shuffle_epi8(SPR, _mm256_and_si256(_mm256_unpacklo_epi8(v, _mm256_srli_epi16(v, 4)), MASK_LO)), d= _mm256_srli_si256(x, 8);
 return _mm256_xor_si256(_mm256_xor_si256(_mm256_and_si256(x, EVEN), _mm256_shuffle_epi8(RED_SQ, _mm256_srli_epi64(v, 62))), _mm256_xor_si256(d, _mm256_slli_epi64(d, 3)));
}
template <bool V> inline __m256i mul2(const __m256i& a_vec, const __m256i& b_vec) {
 const __m256i RED256= _mm256_setr_epi64x(0x6c77415a362d1b00, 0, 0x6c77415a362d1b00, 0);
 __m256i prod;
 if constexpr(V) prod= _mm256_clmulepi64_epi128(a_vec, b_vec, 0);
 else prod= _mm256_setr_m128i(_mm_clmulepi64_si128(_mm256_castsi256_si128(a_vec), _mm256_castsi256_si128(b_vec), 0), _mm_clmulepi64_si128(_mm256_extracti128_si256(a_vec, 1), _mm256_extracti128_si256(b_vec, 1), 0));
 __m256i h= _mm256_srli_si256(prod, 8);
 __m256i d= _mm256_xor_si256(h, _mm256_slli_epi64(h, 1));
 return _mm256_xor_si256(_mm256_xor_si256(prod, _mm256_shuffle_epi8(RED256, _mm256_srli_epi64(h, 60))), _mm256_xor_si256(d, _mm256_slli_epi64(d, 3)));
}
inline std::pair<u64, u64> unpack(const __m256i& vec) { return std::make_pair(u64(_mm256_extract_epi64(vec, 0)), u64(_mm256_extract_epi64(vec, 2))); }
inline __m256i pk(u64 h, u64 l) { return _mm256_set_epi64x(0, h, 0, l); }
#define B2(i) mul2<V>(pk(F<32>(t[e >> 4 * (i + 9) & 15]), F<32>(t[e >> 4 * (i + 8) & 15])), pk(t[e >> 4 * (i + 1) & 15], t[e >> 4 * (i) & 15]))
template <bool V> inline u64 pw(u64 a, u64 e) {
 u64 t[16]= {1, a, sq(a)};
 __m256i t12= pk(t[2], a);
 __m256i t34= mul2<V>(t12, _mm256_set1_epi64x(t[2]));
 std::tie(t[3], t[4])= unpack(t34);
 __m256i t4= _mm256_set1_epi64x(t[4]);
 __m256i t56= mul2<V>(t4, t12);
 std::tie(t[5], t[6])= unpack(t56);
 std::tie(t[7], t[8])= unpack(mul2<V>(t4, t34));
 __m256i t8= _mm256_set1_epi64x(t[8]);
 std::tie(t[9], t[10])= unpack(mul2<V>(t8, t12));
 std::tie(t[11], t[12])= unpack(mul2<V>(t8, t34));
 std::tie(t[13], t[14])= unpack(mul2<V>(t8, t56));
 t[15]= mul(t[7], t[8]);
 auto [b6, b7]= unpack(B2(6));
 auto [b4, b5]= unpack(B2(4));
 __m256i b23= B2(2), b01= B2(0);
 auto [b2, b3]= unpack(mul2<V>(pk(F<16>(b7), F<16>(b6)), b23));
 auto [b0, b1]= unpack(mul2<V>(pk(F<8>(b3), F<8>(b2)), mul2<V>(pk(F<16>(b5), F<16>(b4)), b01)));
 return mul(F<4>(b1), b0);
}
#undef B2
constexpr u64 EMB_B[]= {0x0000000000000001, 0x5fbfaec6aeac0002, 0xb06c601895640004, 0xb013b5277b7c0008, 0xb5ebb915248a0010, 0x109bb25b2c600020, 0xbf3bd95bd4190040, 0x0fc66342279b0080, 0xb6418f5e57c50100, 0xaa194bd4b83f0200, 0x1b5217b4dcc70400, 0xbb06fa73867a0800, 0x006fd55b23331000, 0x4ae8fb39198c2000, 0xfbd141b29b4f4000, 0x1d9ce1776be78000};
constexpr Lin<u64, 2> EMB= Lin<u64, 2>(EMB_B);
constexpr u32 MH_B[16]= {42619, 34034, 37264, 59687, 13661, 58726, 9805, 26873, 8763, 63546, 2437, 49325, 17957, 37424, 41924, 9918};
constexpr Arr<u32, 16> MHI= inverse<u32, 16>(MH_B);
constexpr u16 h641(u64 a) { return u16((a * 0xffef5fb99f1bf6e7) >> 50); }
constexpr u16 PHI_B[16]= {49349, 60640, 60091, 52204, 8753, 26688, 50952, 24030, 14026, 41051, 57150, 31936, 39252, 22252, 63476, 55223};
constexpr u32 MC_B[32]= {0xca137f44, 0x02f9ac22, 0x24119ddf, 0x677fa964, 0x1c3c90b8, 0x61acd330, 0x087e6d0e, 0x98f43405, 0x17ef3800, 0x46a70e74, 0xfdd52d61, 0x9767f2ed, 0xa06bb110, 0xf0ef2346, 0x88d7f773, 0x3bdf87f2, 0xb557b556, 0xaedbaed9, 0xb9ceb9ca, 0xce1bce13, 0x8c848c94, 0xb29cb2bc, 0x65706530, 0xacf1ac71, 0x2fef2eef, 0x48d34ad3, 0xd0b4d4b4, 0x658a6d8a, 0x117b017b, 0xd3a9f3a9, 0x7fa43fa4, 0xbc2d3c2d};
struct ClassTable65537 {
 u16 t[65535];
 u32 K0, s;
};
constexpr u32 OR_L= 6700417, OR_N= 58900, OR_S= 4337141, OR_I= 13, OR_BB= 14, OR_SLOTS= 8u << OR_BB, OR_PAD= 64;
inline u32 eqm(const __m256i& a, const __m256i& b) { return _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpeq_epi32(a, b))); }
constexpr u64 orbit_canon(u32 w) {
 if(!__builtin_is_constant_evaluated()) {
  const __m256i W= _mm256_set1_epi32(int(w));
  const __m256i r0= _mm256_or_si256(_mm256_sllv_epi32(W, _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7)), _mm256_srlv_epi32(W, _mm256_setr_epi32(32, 31, 30, 29, 28, 27, 26, 25)));
  const __m256i r1= _mm256_or_si256(_mm256_slli_epi32(r0, 8), _mm256_srli_epi32(r0, 24));
  const __m256i r2= _mm256_or_si256(_mm256_slli_epi32(r0, 16), _mm256_srli_epi32(r0, 16));
  const __m256i r3= _mm256_or_si256(_mm256_slli_epi32(r0, 24), _mm256_srli_epi32(r0, 8));
  __m256i m= _mm256_min_epu32(_mm256_min_epu32(r0, r1), _mm256_min_epu32(r2, r3));
  m= _mm256_min_epu32(m, _mm256_permute2x128_si256(m, m, 1));
  m= _mm256_min_epu32(m, _mm256_shuffle_epi32(m, 0x4e));
  m= _mm256_min_epu32(m, _mm256_shuffle_epi32(m, 0xb1));
  const u32 f= eqm(r0, m) | eqm(r1, m) << 8 | eqm(r2, m) << 16 | eqm(r3, m) << 24;
  return u64(__builtin_ctz(f)) << 32 | u32(_mm256_cvtsi256_si32(m));
 }
 const u32 nw= ~w;
 u32 z= nw;
 for(u32 t= 1, nz; t < 32 && (nz= z & (nw << t | nw >> (32 - t))); ++t) z= nz;
 u32 j= (31 - __builtin_ctz(z)) & 31, best= j ? w << j | w >> (32 - j) : w, k= j;
 for(z&= z - 1; z; z&= z - 1) {
  j= (31 - __builtin_ctz(z)) & 31;
  if(const u32 v= j ? w << j | w >> (32 - j) : w; v < best) best= v, k= j;
 }
 return u64(k) << 32 | best;
}
struct alignas(64) OrbitTab {
 u64 t[OR_SLOTS + OR_PAD];
 u8 cnt[(OR_SLOTS + OR_PAD) / 8];
 u32 p, c, n;
};
struct OrbitPow {
 u64 lo[256], mid[256], hi[(OR_L >> 16) + 1], hp[OR_I];
};
inline u32 orbit_home(u32 key) { return ((key * 0x9e3779b1u) >> (32 - OR_BB)) * 8; }
struct Ord16 {
 u16 t[65536];
 u32 x, y;
};
inline u64 fold(__uint128_t a) {
 const u64 lo= u64(a), t= lo + u64(a >> 64);
 return t + (t < lo);
}
template <int D> struct Log {
 static constexpr Arr<u16, 65536> INV16= []() {
  u32 f[2][256]{}, b[2][256]{};
  lin(f, MH_B), lin(b, MHI.t);
  Arr<u16, 65536> r{};
  r.t[1]= 1;
  for(u32 k= 32767, x= 1, y= 1; k--;) x= f[0][x & 255] ^ f[1][x >> 8], y= b[0][y & 255] ^ b[1][y >> 8], r.t[x]= y, r.t[y]= x;
  return r;
 }();
 template <bool V> static u64 iv(u64 a) {
  assert(a);
  u64 a32= F<32>(a), b= mul(a, a32);
  auto [g, c]= unpack(mul2<V>(pk(b, a32), _mm256_set1_epi64x(F<16>(b))));
  return mul(EMB(INV16.t[u16(c)]), g);
 }
 static constexpr Arr<u32, 65536> IL16= []() {
  u32 f[2][256]{}, b[2][256]{};
  lin(f, MH_B), lin(b, MHI.t);
  Arr<u32, 65536> r{};
  r.t[1]= 1 << 16;
  for(u32 l= 1, x= 1, y= 1; l < 32768; ++l) x= f[0][x & 255] ^ f[1][x >> 8], y= b[0][y & 255] ^ b[1][y >> 8], r.t[x]= y << 16 | l, r.t[y]= x << 16 | (65535 - l);
  return r;
 }();
 static constexpr Ord16 ord16(Ord16 r, u32 k, u32 e) {
  u32 f[2][256]{}, b[2][256]{};
  lin(f, MH_B), lin(b, MHI.t);
  u32 x= r.x, y= r.y;
  for(u32 l= k; l < e; ++l) x= f[0][x & 255] ^ f[1][x >> 8], y= b[0][y & 255] ^ b[1][y >> 8], r.t[x]= r.t[y]= (l % 3 ? 3 : 1) * (l % 5 ? 5 : 1) * (l % 17 ? 17 : 1) * (l % 257 ? 257 : 1);
  r.x= x, r.y= y;
  return r;
 }
 static constexpr Ord16 ORD16_0= ord16({{0, 1}, 1, 1}, 1, 16385), ORD16= ord16(ORD16_0, 16385, 32768);
 static constexpr Arr<u16, 1 << 14> LN641= []() {
  Arr<u16, 1 << 14> h{};
  for(u64 k= 0, cur= 1, q= cpow(2, ~0ull / 641); k < 641; ++k, cur= cmul(cur, q)) h.t[h641(cur)]= k * 613 % 641;
  return h;
 }();
 static u16 ln641(u64 a) { return LN641.t[h641(a)]; }
 static constexpr Lin<u16, 2> PHI= []() { return Lin<u16, 2>(PHI_B); }();
 static constexpr ClassTable65537 cls65537(ClassTable65537 r, u32 k, u32 e) {
  u32 m[4][256]{};
  lin(m, MC_B);
  u32 s= r.s, b0, b1, l1;
  for(; k < e; ++k) {
   s= m[0][s & 255] ^ m[1][s >> 8 & 255] ^ m[2][s >> 16 & 255] ^ m[3][s >> 24], b0= s & 65535, b1= s >> 16, l1= 65535 - (IL16.t[b1] & 65535);
   if(b0) r.t[((IL16.t[b0] & 65535) + l1) % 65535]= u16(k - 1);
   else r.K0= k;
   if(b0 ^ b1) r.t[((IL16.t[b0 ^ b1] & 65535) + l1) % 65535]= u16(65536 - k);
   else r.K0= 65537 - k;
  }
  r.s= s;
  return r;
 }
 static constexpr ClassTable65537 CLS65537_0= cls65537({{}, 0, 1}, 1, 16385), CLS65537= cls65537(CLS65537_0, 16385, 32769);
 static u32 log_65537(u64 n, u64 fn, const ClassTable65537& c= CLS65537) {
  const u16 b1= n ^ fn;
  if(!b1) return 0;
  const u16 b0= n ^ PHI(b1);
  if(!b0) return c.K0;
  u32 idx= (IL16.t[b0] & 65535) + 65535 - (IL16.t[b1] & 65535);
  if(idx >= 65535) idx-= 65535;
  return c.t[idx] + 1u;
 }
 static constexpr Lin<u32, 8> OR_LAM= []() {
  u32 b[64]{};
  for(int i= 0; i < 64; ++i)
   for(u64 j= 0, w= u64(1) << i; j < 32; ++j, w= F<63>(w)) b[i]|= u32(__builtin_parityll(w & 0x4206046c1dd56086)) << j;
  return Lin<u32, 8>(b);
 }();
 static constexpr u64 OR_G= []() { return cpow(2, ~0ull / OR_L); }(), OR_H1= []() { return cpow(OR_G, OR_L - OR_S); }();
 static constexpr Arr<u32, 32> or_mul(u64 a) {
  Arr<u32, 32> r{};
  for(u64 j= 0, w= 0x20000001aa3c9c8c; j < 32; ++j, w= cmul(w, w)) r.t[j]= OR_LAM(cmul(a, w));
  return r;
 }
 static constexpr Lin<u32, 4> OR_M1= []() { return Lin<u32, 4>(or_mul(OR_H1 ^ F<32>(OR_H1)).t); }();
 static constexpr Arr<u32, 32> OR_MG2= []() { return or_mul(cmul(OR_G, OR_G) ^ F<32>(cmul(OR_G, OR_G))); }();
 static constexpr OrbitTab orbit_part(OrbitTab r, u32 e) {
  u32 m[4][256]{};
  lin(m, OR_MG2.t);
  u32 p= r.p, c= r.c, n= r.n;
  for(; n <= e; ++n) {
   const u64 ck= orbit_canon(c);
   const u32 best= u32(ck), k= u32(ck >> 32);
   u32 g= (best * 0x9e3779b1u) >> (32 - OR_BB);
   while(r.cnt[g] == 8) ++g;
   r.t[g * 8 + r.cnt[g]++]= u64(best) << 32 | u32((u64(2 * n - 1) << k) % OR_L);
   const u32 x= m[0][c & 255] ^ m[1][c >> 8 & 255] ^ m[2][c >> 16 & 255] ^ m[3][c >> 24] ^ p;
   p= c, c= x;
  }
  r.p= p, r.c= c, r.n= n;
  return r;
 }
 static constexpr OrbitTab OR_TAB_0= orbit_part(OrbitTab{{}, {}, OR_LAM(OR_G), OR_LAM(OR_G), 1}, 5890), OR_TAB_1= orbit_part(OR_TAB_0, 11780), OR_TAB_2= orbit_part(OR_TAB_1, 17670), OR_TAB_3= orbit_part(OR_TAB_2, 23560), OR_TAB_4= orbit_part(OR_TAB_3, 29450), OR_TAB_5= orbit_part(OR_TAB_4, 35340), OR_TAB_6= orbit_part(OR_TAB_5, 41230), OR_TAB_7= orbit_part(OR_TAB_6, 47120), OR_TAB_8= orbit_part(OR_TAB_7, 53010), OR_TAB= orbit_part(OR_TAB_8, OR_N);
 static constexpr OrbitPow OR_POW= []() {
  OrbitPow r{};
  r.lo[0]= r.mid[0]= r.hi[0]= r.hp[0]= 1;
  const u64 g8= cpow(OR_G, 256), g16= cpow(g8, 256);
  for(int i= 1; i < 256; ++i) r.lo[i]= cmul(r.lo[i - 1], OR_G), r.mid[i]= cmul(r.mid[i - 1], g8);
  for(u32 i= 1; i <= (OR_L >> 16); ++i) r.hi[i]= cmul(r.hi[i - 1], g16);
  for(u32 i= 1; i < OR_I; ++i) r.hp[i]= cmul(r.hp[i - 1], OR_H1);
  return r;
 }();
 static constexpr Lin<u32, 8> OR_LAMH= []() {
  u32 b[64]{};
  for(int i= 0; i < 64; ++i) b[i]= OR_LAM(cmul(u64(1) << i, OR_H1));
  return Lin<u32, 8>(b);
 }();
 static constexpr Arr<u32, 32> OR_P2INV= []() {
  Arr<u32, 32> r{{1}};
  for(int k= 1; k < 32; ++k) r.t[k]= u64(r.t[k - 1]) * (OR_L / 2 + 1) % OR_L;
  return r;
 }();
 static u32 orbit_lookup(u32 key) {
  const __m256i K= _mm256_set1_epi64x((long long)(u64(key) << 32));
  for(u32 h= orbit_home(key);; h+= 8) {
   const __m256i* q= (const __m256i*)&OR_TAB.t[h];
   const u32 f= eqm(_mm256_loadu_si256(q), K) | eqm(_mm256_loadu_si256(q + 1), K) << 8;
   if(f & 0xaaaa) return u32(OR_TAB.t[h + (__builtin_ctz(f & 0xaaaa) >> 1)]);
   if(f & 0x5555) return 0;
  }
 }
 static u32 orbit_find(u64 y, u32 c, u64 ck, u32 i) {
  constexpr u32 L= OR_L;
  if(!c) return u32(u64(OR_S) * i % L);
  const u32 e= orbit_lookup(u32(ck));
  if(!e) return L;
  const u32 e1= u32(u64(e) * OR_P2INV.t[ck >> 32] % L);
  const u64 z= i ? mul(y, OR_POW.hp[i]) : y, g= mul(mul(OR_POW.lo[e1 & 255], OR_POW.mid[e1 >> 8 & 255]), OR_POW.hi[e1 >> 16]);
  return u32((u64(OR_S) * i + (z == g ? e1 : L - e1)) % L);
 }
 static void orbit_prefetch(u64 ck) { _mm_prefetch((const char*)&OR_TAB.t[orbit_home(u32(ck))], _MM_HINT_T0); }
 static u32 log_6700417(u64 y) {
  u32 c0= OR_LAM(y), c1= OR_LAMH(y);
  for(u32 i= 0; i < OR_I; i+= 2) {
   const u64 k0= c0 ? orbit_canon(c0) : 0, k1= c1 ? orbit_canon(c1) : 0;
   orbit_prefetch(k0), orbit_prefetch(k1);
   if(const u32 e= orbit_find(y, c0, k0, i); e != OR_L) return e;
   if(i + 1 < OR_I)
    if(const u32 e= orbit_find(y, c1, k1, i + 1); e != OR_L) return e;
   c0= OR_M1(c1) ^ c0, c1= OR_M1(c0) ^ c1;
  }
  return OR_L;
 }
 template <bool V> static u64 ln(u64 x) {
  assert(x);
  const u64 x32= F<32>(x), n= mul(x, x32), fn= F<16>(n);
  auto [x_f16, w]= unpack(mul2<V>(pk(sq(x32), n), _mm256_set1_epi64x(fn)));
  const u32 il= IL16.t[u16(x_f16)];
  const u64 s= mul(EMB(u16(il >> 16)), w), s2= sq(s), ye= mul(mul(s, s2), F<57, D>(s)), y= mul(ye, s2);
  u32 r3= log_6700417(y);
  u64 r0= ln641(mul(mul(mul(s, sq(ye)), F<51, D>(ye)), F<29, D>(y))), r2= log_65537(n, fn);
  return fold(0x663d80ff99c27f * r0 + __uint128_t(0x2f205935d0dfa6ca) * r3 + 0x1000100010001ull * (il & 65535) + 0xffff0000ffff * r2);
 }
 static u64 ord(u64 a) {
  assert(a);
  const u64 n= mul(a, F<32>(a)), fn= F<16>(n), a2= sq(a), x3= mul(a, a2), y= mul(x3, F<57, D>(a)), c= mul(y, a2);
  const u64 d= mul(mul(mul(y, F<21, D>(a)), F<60, D>(x3)), F<50, D>(y));
  return u64(ORD16.t[u16(mul(n, fn))]) * (fn != n ? 65537 : 1) * (F<32>(d) != d ? 641 : 1) * (F<32>(c) != c ? 6700417 : 1);
 }
 static constexpr u32 FP= 6700417, FN= 128, FB= FN * FN, FA= FP / FN;
 template <class T, u32 P, int N= P> static constexpr Arr<T, N> inv_tab(Arr<T, N> r, u32 lo, u32 hi) {
  for(u32 i= lo; i < hi; ++i) r.t[i]= T(P - u64(P / i) * r.t[P % i] % P);
  return r;
 }
 static constexpr Arr<u32, 255> T255= []() {
  Arr<u32, 255> r{};
  for(u32 x= 0; x < 255; ++x) {
   u32 y= x;
   for(int i= 1; i < 15; ++i) y= y * x % 255;
   r.t[x]= 8224 * y % 65535 | u32(x % 3 == 0) << 16 | u32(x % 5 == 0) << 17 | u32(x % 17 == 0) << 18;
  }
  return r;
 }();
 static constexpr Arr<u16, 257> T257= []() {
  Arr<u16, 257> r= inv_tab<u16, 257>({{0, 1}}, 2, 257);
  for(u32 x= 1; x < 257; ++x) r.t[x]= u16(8160 * r.t[x] % 65535);
  return r;
 }();
 static constexpr Arr<u16, 641> INV641= inv_tab<u16, 641>({{0, 590}}, 2, 641);
 static constexpr Arr<u32, 65537> INV65537_0= inv_tab<u32, 65537>({{0, 16384}}, 2, 32769), INV65537= inv_tab<u32, 65537>(INV65537_0, 32769, 65537);
 static constexpr ClassTable65537 ICLS65537= []() {
  ClassTable65537 r{};
  for(u32 i= 0; i < 65535; ++i) r.t[i]= INV65537.t[CLS65537.t[i] + 1] - 1;
  r.K0= INV65537.t[CLS65537.K0];
  return r;
 }();
 static constexpr Arr<u32, FB> FR= []() {
  Arr<u32, FB> r{};
  for(u32 a= 0, b= 1, c= 1, d= FN;;) {
   for(u32 i= a * FB / b + 1, e= c == d ? FB - 1 : c * FB / d; i <= e; ++i) r.t[i]= ((i + 1) * d >= c * FB ? c | d << 8 : a | b << 8) | (c | d << 8) << 16;
   if(c == 1 && d == 1) break;
   const u32 k= (FN + b) / d, e= k * c - a, f= k * d - b;
   a= c, b= d, c= e, d= f;
  }
  r.t[0]= 1 << 8 | 1 << 24;
  return r;
 }();
 static constexpr Arr<u32, FA + 1> FI= inv_tab<u32, FP, FA + 1>({{0, 3883315}}, 2, FA + 1);
 static u32 inv_6700417(u32 x) {
  const u32 f= FR.t[u64(x) * FB / FP], c1= f & 255, b1= f >> 8 & 255, c2= f >> 16 & 255, b2= f >> 24;
  const long long a1= (long long)b1 * x - (long long)c1 * FP, a2= (long long)b2 * x - (long long)c2 * FP, m1= a1 < 0 ? -a1 : a1, m2= a2 < 0 ? -a2 : a2;
  const bool s= m1 <= m2;
  const u32 r= u32(u64(s ? b1 : b2) * FI.t[s ? m1 : m2] % FP);
  return (s ? a1 : a2) < 0 && r ? FP - r : r;
 }
 template <bool V> static u64 ln(u64 a, u64 b) {
  constexpr u64 M= ~0ull;
  if(b == 1) return 0;
  if(a == b) return 1;
  if(a <= 1 || !b) return M;
  const u64 a32= F<32>(a), b32= F<32>(b);
  const __m256i x32= pk(b32, a32);
  const auto [na, nb]= unpack(mul2<V>(pk(b, a), x32));
  const u64 fna= F<16>(na), fnb= F<16>(nb);
  const __m256i fn2= pk(fnb, fna);
  const auto [ma, mb]= unpack(mul2<V>(pk(nb, na), fn2));
  const u32 ila= IL16.t[u16(ma)], ilb= IL16.t[u16(mb)], la= ila & 65535, lb= ilb & 65535, a255= la % 255, b255= lb % 255, a257= la % 257, b257= lb % 257;
  const auto [wa, wb]= unpack(mul2<V>(sq2(x32), fn2));
  const __m256i s2= mul2<V>(pk(EMB(u16(ilb >> 16)), EMB(u16(ila >> 16))), pk(wb, wa));
  const auto [sa, sb]= unpack(s2);
  const __m256i q2= sq2(s2), ye2= mul2<V>(mul2<V>(s2, q2), pk(F<57, D>(sb), F<57, D>(sa))), y2= mul2<V>(ye2, q2);
  const auto [yea, yeb]= unpack(ye2);
  const auto [ya, yb]= unpack(y2);
  const __m256i w2= mul2<V>(mul2<V>(s2, sq2(ye2)), pk(F<51, D>(yeb), F<51, D>(yea)));
  const auto [za, zb]= unpack(mul2<V>(w2, pk(F<29, D>(yb), F<29, D>(ya))));
  const u32 v0a= ln641(za), v0b= ln641(zb), i2a= log_65537(na, fna, ICLS65537), v2b= log_65537(nb, fnb);
  const u32 ra= log_6700417(ya), rb= log_6700417(yb), t255= T255.t[a255];
  if((t255 >> 16 & ~(T255.t[b255] >> 16)) || (!a257 && b257) || (!v0a && v0b) || (!i2a && v2b) || (!ra && rb)) return M;
  const u64 v1= b255 * (t255 & 65535) + b257 * T257.t[a257], v0= v0b * INV641.t[v0a], v2= u64(v2b) * i2a, v3= u64(rb) * inv_6700417(ra);
  return fold(__uint128_t(0x1000100010001) * v1 + __uint128_t(0x663d80ff99c27f) * v0 + __uint128_t(0xffff0000ffff) * v2 + __uint128_t(0x280fffffd7f) * v3) % (u64(t255 >> 16 & 1 ? 1 : 3) * (t255 >> 17 & 1 ? 1 : 5) * (t255 >> 18 & 1 ? 1 : 17) * (a257 ? 257 : 1) * (v0a ? 641 : 1) * (i2a ? 65537 : 1) * (ra ? 6700417 : 1));
 }
};
#ifdef __x86_64__
#define VP(f, ...) \
 if(__builtin_cpu_supports("vpclmulqdq")) return f<1>(__VA_ARGS__); \
 return f<0>(__VA_ARGS__)
#else
#define VP(f, ...) return f<0>(__VA_ARGS__)
#endif
template <int D> class GF2p64_ {
 u64 x;
 u64 iv_() const { VP(Log<D>::template iv, x); }
public:
 GF2p64_(u64 y= 0): x(y) {}
 static GF2p64_ from_nimber(u64 n) { return GF2p64_(FROM_NIM(n)); }
 static GF2p64_ generator() { return GF2p64_(2); }
 auto operator<=>(const GF2p64_&) const= default;
 GF2p64_& operator+=(GF2p64_ r) { return x^= r.x, *this; }
 GF2p64_& operator-=(GF2p64_ r) { return x^= r.x, *this; }
 GF2p64_& operator*=(GF2p64_ r) { return x= mul(x, r.x), *this; }
 GF2p64_& operator/=(GF2p64_ r) { return x= mul(x, r.iv_()), *this; }
 GF2p64_ operator+(GF2p64_ r) const { return x ^ r.x; }
 GF2p64_ operator-(GF2p64_ r) const { return x ^ r.x; }
 GF2p64_ operator*(GF2p64_ r) const { return mul(x, r.x); }
 GF2p64_ operator/(GF2p64_ r) const { return mul(x, r.iv_()); }
 GF2p64_ square() const { return sq(x); }
 GF2p64_ sqrt() const { return F<63>(x); }
 GF2p64_ inv() const { return iv_(); }
 GF2p64_ pow(u64 e) const { VP(pw, x, e); }
 u64 log() const { VP(Log<D>::template ln, x); }
 u64 log(GF2p64_ base) const { VP(Log<D>::template ln, base.x, x); }
 u64 ord() const { return Log<D>::ord(x); }
 u64 to_nimber() const { return TO_NIM(x); }
 explicit operator u64() const { return x; }
 explicit operator bool() const { return x != 0; }
 friend std::ostream& operator<<(std::ostream& os, const GF2p64_& r) { return os << r.x; }
 friend std::istream& operator>>(std::istream& is, GF2p64_& r) { return is >> r.x; }
};
#undef VP
}
using GF2p64= gf2p64_internal::GF2p64_<0>;
