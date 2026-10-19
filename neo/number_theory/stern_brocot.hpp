#pragma once
#include <array>
#include <utility>
#ifdef __x86_64__
#include <immintrin.h>
#endif
// Stern-Brocot 木の問を、1 ≤ a, b < 2^32 の既約分数について解く。道は連の長さを R、L、R、... の順に並べて表し、
// 最初の R の連だけ 0 になりうる (a/b < 1 のとき)。どの問も a/b の連分数展開 [q_0; q_1, ..., q_m] (互除法の商) を使う。
// a/b への道は R^{q_0} L^{q_1} R^{q_2} ... で、最後の連だけ q_m - 1。道の途中の節点は近似分数 h_i / k_i
// (h_i = q_i h_{i-1} + h_{i-2}、h_{-1} / k_{-1} = 1/0、h_{-2} / k_{-2} = 0/1) から作る。連 t_0, ..., t_j の道の行き先は
// (h_j + h_{j-1}) / (k_j + k_{j-1}) で、h_j は q_j を t_j に替えて作る。答えは 2^32 未満なので、途中は u32 の桁あふれのまま
// 計算してよい。互除法は 32 bit の割り算で回す (x86 の div は商と余りを一緒に出すので、1 段の鎖は割り算 1 回ぶん)。
// procon-judge の yosupo-stern-brocot-tree の 20 回目の sb20 の形。
namespace stern_brocot_internal {
using u32= unsigned;
using u64= unsigned long long;
using u128= unsigned __int128;
// r > h なら s を、そうでなければ r を返し、f = (r > h) とする。gcc は三項演算子の選びを分岐に戻すので、x86 では asm で書く
[[gnu::always_inline]] inline u32 pick(u32 r, u32 s, u32 h, u32& f) {
#ifdef __x86_64__
 u32 ff;
 asm("xorl %k1, %k1\n\t"
     "cmpl %k3, %k0\n\t"
     "seta %b1\n\t"
     "cmova %k2, %k0"
     : "+r"(r), "=&r"(ff)
     : "r"(s), "r"(h)
     : "cc");
 f= ff;
 return r;
#else
 f= r > h;
 return f ? s : r;
#endif
}
inline u64 inv64(u64 n) {
 u64 x= (3 * n) ^ 2, y= 1 - n * x;
 x*= 1 + y, y*= y;
 x*= 1 + y, y*= y;
 x*= 1 + y, y*= y;
 return x * (1 + y);
}
// v 2^{-c} mod M を [0, M] で返す (0 <= c <= 64、v <= M、M は奇数、Minv = M^{-1} mod 2^64)
inline u64 redc(u64 v, int c, u64 M, u64 Minv) {
 u128 t= (u128)v << (64 - c);
 u64 hi= u64(t >> 64), mh= u64((u128)(u64(t) * Minv) * M >> 64);
 return hi >= mh ? hi - mh : hi - mh + M;
}
// A x ≡ 1 (mod B) の x を [1, B] で返す。B は奇数、gcd(A, B) = 1、Binv = B^{-1} mod 2^64。
// neo/number_theory/inv_gcd.hpp と同じ Kaliski の形の Stein の方法を、u32 と gcd = 1 に絞ったもの。k < 64 なので還元は 1 回。
inline u32 inv_odd(u32 A, u32 B, u64 Binv) {
 constexpr int GAP= 8;
 if(A > B && __builtin_clz(B) - __builtin_clz(A) > GAP) A%= B;
 if(A == 0) return B;  // B = 1
 int k= __builtin_ctz(A);
 u32 a= B, b= A >> k, ca= 0, cb= 1, neg= 0;
 if(__builtin_clz(b) - __builtin_clz(a) > GAP) {
  const u32 q= a / b, r= a - q * b;
  if(r == 0) a= b, ca= q - 1;
  else {
   const int c= __builtin_ctz(r);
   a= r >> c, ca= q, cb<<= c, k+= c;
  }
 }
 while(a != b) {
  const u32 d= a - b, e= b - a;
  const int c= __builtin_ctz(d);
  const u32 mn= a < b ? a : b, t= a > b ? d : e, cm= a > b ? cb : ca;
  neg+= a < b;
  ca+= cb, cb= cm << c, k+= c;
  a= t >> c, b= mn;
 }
 u64 x= redc(cb, k, B, Binv);
 x= x >= B ? x - B : x;
 return u32((neg & 1) ? B - x : (x ? x : B));
}
}
// a/b への道の連の長さを out に書き、書き終わりを返す。連の数は 46 以下で、out は 46 個ぶんあれば足りる。
// 互除法の 1 段 a = q b + r で r > b/2 なら次の段の商は 1 なので、次の段を (b, b - r) にして 2 段を 1 回の割り算で済ませる
// (b を b - r で割った余りは r を b - r で割った余りと同じで、商は 1 大きい)。割り算は一様な分数で 7 割ほどに減る。
// 1 段で商と続く 1 の 2 つを書いてから 1 つか 2 つ進むので、仮に書いた 1 は次の段が上書きする。
inline unsigned* sb_encode(unsigned a, unsigned b, unsigned* out) {
 using namespace stern_brocot_internal;
 u32 adj= 0;
 while(true) {
  const u32 q= a / b, r= a % b, Q= q - adj;
  if(r == 0) {
   *out++= Q - 1;
   return out;
  }
  u32 f;
  const u32 nb= pick(r, b - r, b >> 1, f);
  out[0]= Q, out[1]= 1;
  out+= 1 + f;
  a= b, b= nb, adj= f;
 }
}
// 連の長さ [first, last) (R から) の道の行き先 (a, b)。答えは 2^32 未満であること。
// x86 では (lx, ly) と (rx, ry) を xmm の 64 bit の枠に置いて vpmuludq で掛ける (コンパイラは u32 のままだと vpmulld に
// ベクトル化し、Intel ではそれが 10 サイクル待つ)。
inline std::pair<unsigned, unsigned> sb_decode(const unsigned* first, const unsigned* last) {
 using namespace stern_brocot_internal;
 const u32 len= u32(last - first);
 u32 i= 0;
#ifdef __x86_64__
 __m128i L= _mm_set_epi64x(1, 0), R= _mm_set_epi64x(0, 1);
 for(; i + 1 < len; i+= 2) {
  L= _mm_add_epi64(L, _mm_mul_epu32(_mm_set1_epi64x(first[i]), R));
  R= _mm_add_epi64(R, _mm_mul_epu32(_mm_set1_epi64x(first[i + 1]), L));
 }
 if(i < len) L= _mm_add_epi64(L, _mm_mul_epu32(_mm_set1_epi64x(first[i]), R));
 const __m128i s= _mm_add_epi64(L, R);
 return {u32(_mm_cvtsi128_si64(s)), u32(_mm_extract_epi64(s, 1))};
#else
 u32 lx= 0, ly= 1, rx= 1, ry= 0;
 for(; i + 1 < len; i+= 2) {
  lx+= first[i] * rx, ly+= first[i] * ry;
  rx+= first[i + 1] * lx, ry+= first[i + 1] * ly;
 }
 if(i < len) lx+= first[i] * rx, ly+= first[i] * ry;
 return {lx + rx, ly + ry};
#endif
}
// a/b と c/d の最も深い共通の祖先。2 つの互除法を同じ段で進め、連が初めて食い違う段か、どちらかの道が終わる段で、
// 短いほうの連で止まる。
inline std::pair<unsigned, unsigned> sb_lca(unsigned a, unsigned b, unsigned c, unsigned d) {
 using namespace stern_brocot_internal;
 u32 h1= 1, k1= 0, h0= 0, k0= 1;
 while(true) {
  const u32 q= a / b, r= a % b;
  const u32 s= c / d, t= c % d;
  const u32 T= q - (r == 0), S= s - (t == 0);
  if(T != S || r == 0 || t == 0) {
   const u32 m= T < S ? T : S;
   return {m * h1 + h0 + h1, m * k1 + k0 + k1};
  }
  const u32 h= q * h1 + h0, k= q * k1 + k0;
  h0= h1, k0= k1, h1= h, k1= k;
  a= b, b= r, c= d, d= t;
 }
}
// a/b の祖先のうち深さ k のもの (1/1 の深さは 0)。無ければ (0, 0)。互除法の組 (x, y) (x ≥ y) から先の深さは x - 1 以下なので
// (x = q y + r として q + Σ(y, r) ≤ q + y ≤ x)、残りの深さが x 以上になったら無いと分かる。
inline std::pair<unsigned, unsigned> sb_ancestor(unsigned long long k, unsigned a, unsigned b) {
 using namespace stern_brocot_internal;
 if(k >= (a > b ? a : b)) return {0, 0};
 u32 kk= u32(k), h1= 1, k1= 0, h0= 0, k0= 1;
 while(true) {
  const u32 q= a / b, r= a % b;
  const u32 T= q - (r == 0);
  if(kk <= T) return {kk * h1 + h0 + h1, kk * k1 + k0 + k1};
  if(r == 0) return {0, 0};
  kk-= q;
  const u32 h= q * h1 + h0, k2= q * k1 + k0;
  h0= h1, k0= k1, h1= h, k1= k2;
  a= b, b= r;
  if(kk >= a) return {0, 0};
 }
}
// a/b の子孫の集まりである開区間の下端 f/g と上端 h/k を {f, g, h, k} で返す (0 は 0/1、+∞ は 1/0)。
// a g - b f = 1 なので、b が奇数なら g = a^{-1} mod b、f = (a g - 1) / b。b が偶数なら a が奇数なので h = b^{-1} mod a、
// k = (b h - 1) / a。逆元は割り算を使わない Stein の方法で、もう一方は奇数の法で割り切れる割り算なので法の逆元を掛ける。
inline std::array<unsigned, 4> sb_range(unsigned a, unsigned b) {
 using namespace stern_brocot_internal;
 u32 m= 0u - (~b & 1);  // b が偶数なら全ビット 1
 asm("" : "+r"(m));
 const u32 B= (a & m) | (b & ~m), A= a ^ b ^ B;
 const u64 Binv= inv64(B);
 const u32 x= inv_odd(A, B, Binv), y= u32((u64(A) * x - 1) * Binv);
 const u32 g= (x & ~m) | ((b - y) & m), f= (y & ~m) | ((a - x) & m);
 return {f, g, a - f, b - g};
}
