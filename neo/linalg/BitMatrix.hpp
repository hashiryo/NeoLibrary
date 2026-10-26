#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <optional>
#include <string>
#include <utility>
#include <vector>
// GF(2) の行列。行を 64 bit の語の列に詰めて持つ (行 i の j 列目は、行 i の語 j / 64 の bit j % 64)。行の語の数は ⌈w / 64⌉
// で、w 列目より後ろの bit はいつも 0 にしておく。
// 積 A B は Four Russians の方法で求める。C の列を 1024 bit の帯に分け、帯ごとに、B の 64 行分 (A の 1 語分) を 8 行ずつ
// 256 通りに XOR した表 8 枚を作ってから、A の各行の C の帯 (ymm 4 本) に表の行を 8 本 XOR する。C の帯の読み書きは、B の
// 64 行につき 1 回になる。A の行が 64 本以下なら、A の 1 の bit ごとに B の行を XOR するだけにする。
// 行列式、階数、連立方程式、核、逆行列は、要 (ピボット) を 16 本ずつまとめる掃き出し (M4RI の形) で求める。組の中では、残りの
// 行ごとに「どの要の行を足すか」の印と、要を探す列の語に貯めた更新を足した今の値だけを持つ。列の bit はその値の bit で、要を
// 足すのはその値に要の行の語を XOR するだけなので、4 行ずつ AVX2 で調べて書く。要のない列では、残りの行の値の OR から次に
// 1 のある列まで飛ぶ (階数の低い行列で、列ごとに全部の行をなめない)。16 本たまったら、要の行 8 本ずつの 256 通りの XOR の
// 表を 2 枚作り、残りの各行に表の行を 1 本ずつ XOR する (行が 64 本以下なら表を作らずに要の行を足す)。行列の読み書きは要
// 16 本につき 1 回になる。核と逆行列は、要の行を下から 16 本ずつ組にし、組の表で上の要の行を消して RREF に直してから読む。
// 核の基底は、要の行を要の列の位置に、自由な列 g には単位ベクトル e_g を置いた正方行列の列 g なので、その行列の 0 でない
// 64 × 64 の小行列だけを転置して取り出す。逆行列は [A | I] の I の列を、行を要にした順に付け替えて置くので、前進消去で読み書き
// する右の部分が要の数の分だけになる。
// 64 列以下の階数と 64 次以下の行列式は、1 語の行のまま掃き出す。行が 64 本以下の階数は、行を 1 本ずつ前の要で消して入れ、
// 行が列の 2 倍より多い行列の階数は転置してから消す。
// procon-judge の yosupo の mod 2 の 5 問で書き比べた gf2_mul_m4r8 と gf2_*_blk16 を写し、掃き出しの要を探す手順を書き直した
// もの。
namespace bit_matrix_internal {
using u32= unsigned;
using u64= unsigned long long;
constexpr int K= 16;  // 1 組の要の数
inline int round4(int x) { return (x + 3) & ~3; }
// x[w, W) ^= y[w, W)、x[w, W) ^= y ^ z、x[w, W) = y ^ z、x[w, W) = y。W は 4 の倍数で、w は 4 の倍数とは限らない。
inline void xor_in(u64* x, const u64* y, int w, int W) {
 for(; w & 3; ++w) x[w]^= y[w];
 for(; w < W; w+= 4) _mm256_storeu_si256((__m256i*)(x + w), _mm256_xor_si256(_mm256_loadu_si256((const __m256i*)(x + w)), _mm256_loadu_si256((const __m256i*)(y + w))));
}
inline void xor2_in(u64* x, const u64* y, const u64* z, int w, int W) {
 for(; w & 3; ++w) x[w]^= y[w] ^ z[w];
 for(; w < W; w+= 4) _mm256_storeu_si256((__m256i*)(x + w), _mm256_xor_si256(_mm256_loadu_si256((const __m256i*)(x + w)), _mm256_xor_si256(_mm256_loadu_si256((const __m256i*)(y + w)), _mm256_loadu_si256((const __m256i*)(z + w)))));
}
inline void xor_to(u64* x, const u64* y, const u64* z, int w, int W) {
 for(; w & 3; ++w) x[w]= y[w] ^ z[w];
 for(; w < W; w+= 4) _mm256_storeu_si256((__m256i*)(x + w), _mm256_xor_si256(_mm256_loadu_si256((const __m256i*)(y + w)), _mm256_loadu_si256((const __m256i*)(z + w))));
}
inline void copy_to(u64* x, const u64* y, int w, int W) {
 for(; w & 3; ++w) x[w]= y[w];
 for(; w < W; w+= 4) _mm256_storeu_si256((__m256i*)(x + w), _mm256_loadu_si256((const __m256i*)(y + w)));
}
// x[0, n) ^= y[0, n)。n は 4 の倍数とは限らない。
inline void xor_n(u64* x, const u64* y, int n) {
 int i= 0;
 for(; i + 4 <= n; i+= 4) _mm256_storeu_si256((__m256i*)(x + i), _mm256_xor_si256(_mm256_loadu_si256((const __m256i*)(x + i)), _mm256_loadu_si256((const __m256i*)(y + i))));
 for(; i < n; ++i) x[i]^= y[i];
}
inline bool bit(const u64* x, int j) { return x[j >> 6] >> (j & 63) & 1; }
// 文字 s[0, m) の下位 bit を語 dst に詰める (dst は 0 で埋めてあること)。32 文字ずつ、bit 0 を bit 7 へずらして movemask で取る。
inline void pack_row(const char* s, int m, u64* dst) {
 int j= 0;
 for(; j + 32 <= m; j+= 32) {
  const __m256i c= _mm256_loadu_si256((const __m256i*)(s + j));
  dst[j >> 6]|= u64(u32(_mm256_movemask_epi8(_mm256_slli_epi16(c, 7)))) << (j & 63);
 }
 for(; j < m; ++j) dst[j >> 6]|= u64(s[j] & 1) << (j & 63);
}
// 語 src の bit [0, m) を '0' と '1' の m 文字にして dst に書く。32 bit ずつ、各 byte に配ってから bit を調べる。
inline void unpack_row(const u64* src, int m, char* dst) {
 const __m256i sh= _mm256_setr_epi8(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3);
 const __m256i bm= _mm256_set1_epi64x(0x8040201008040201ll), zero= _mm256_set1_epi8('0');
 int j= 0;
 for(; j + 32 <= m; j+= 32) {
  const __m256i v= _mm256_shuffle_epi8(_mm256_set1_epi32(int(u32(src[j >> 6] >> (j & 63)))), sh);
  _mm256_storeu_si256((__m256i*)(dst + j), _mm256_sub_epi8(zero, _mm256_cmpeq_epi8(_mm256_and_si256(v, bm), bm)));
 }
 for(; j < m; ++j) dst[j]= char('0' + (src[j >> 6] >> (j & 63) & 1));
}
// 64 × 64 の bit 行列 x (x[r] の bit c が r 行 c 列) をその場で転置する。対角でない 2 つの小行列を入れ替えながら半分ずつにする。
inline void transpose64(u64* x) {
 u64 m= 0x00000000ffffffffull;
 for(int j= 32; j; j>>= 1, m^= m << j)
  for(int k= 0; k < 64; k= ((k | j) + 1) & ~j) {
   const u64 t= ((x[k] >> j) ^ x[k | j]) & m;
   x[k]^= t << j, x[k | j]^= t;
  }
}
// n 行、語の数 W (4 の倍数) の行列の前進消去。要は列の順に探し、その列で bit が立っている最初の行を選ぶ。
struct Elim {
 int n, W;
 std::vector<u64> a;           // n × W
 std::vector<int> rem;         // まだ要にしていない行
 std::vector<int> prow, pcol;  // 要にした行と列 (要にした順)
 std::vector<u64> U;           // 要の行 (要にした順に W 語ずつ)。前の要で消したあとの値で、要の列の語より左は 0
 // aug >= 0 なら、列 aug から右に単位行列を付けたものとして消す。単位行列の列は置かずに、行を t 番目の要にしたとき、その行の
 // 単位行列の 1 を列 aug + t に置く。すると要が t 本の間、右の部分の値は列 aug + t より左にしかなく、そこまでを読み書きすれば
 // よい。最後に右の部分の列 aug + t は元の行 prow[t] を表す。
 int aug= -1;
 // 語の数が sw の行 n 本 src を写す。
 Elim(int n, int W, const u64* src, int sw): n(n), W(W), rem(n) {
  if(sw == W) a.assign(src, src + (size_t)n * W);
  else {
   a.assign((size_t)n * W, 0);
   if(sw)
    for(int i= 0; i < n; ++i) std::memcpy(row(i), src + (size_t)i * sw, sw * sizeof(u64));
  }
  for(int i= 0; i < n; ++i) rem[i]= i;
 }
 u64* row(int r) { return a.data() + (size_t)r * W; }
 u64* urow(int t) { return U.data() + (size_t)t * W; }
 const u64* urow(int t) const { return U.data() + (size_t)t * W; }
 int rank() const { return prow.size(); }
 // 要が t 本のときに読み書きする語の数。
 int end_of(int t) const { return aug < 0 ? W : std::min(W, round4((aug + t + 63) >> 6)); }
 // 列 [0, lim) で要を探しながら前進消去する。要のない列があれば、stop なら false を返して止め、そうでなければ飛ばす。
 // 残りの行 rem[i] について、要を探す列の語に貯めた更新を足した今の値を v[i] に、足す要の印を f[i] に持つ。v があれば、列の
 // bit は v[i] の bit で、要を足すのは v[i] に要の行の語を XOR するだけになり、4 行ずつ AVX2 で調べて書ける。要のない列では、
 // v の OR から、この語の中で次に 1 のある列まで飛ぶ。
 bool run(int lim, bool stop) {
  U.assign((size_t)std::min(n, lim) * W, 0);
  std::vector<u64> v(n), f(n), T1, T2;
  const u64* Ub[K];
  int nr= n, col= 0, cwi= -1, bt= 0;
  while(col < lim && nr) {
   const int wi= col >> 6, b= col & 63;
   if(wi != cwi) gather(wi, nr, v.data(), f.data(), Ub, bt), cwi= wi;
   // 列 col が 1 の最初の行を探す。
   const __m128i sh= _mm_cvtsi32_si128(63 - b);
   int pi= -1, i= 0;
   __m256i acc= _mm256_setzero_si256();
   for(; i + 4 <= nr; i+= 4) {
    const __m256i x= _mm256_loadu_si256((const __m256i*)(v.data() + i));
    if(const int mk= _mm256_movemask_pd(_mm256_castsi256_pd(_mm256_sll_epi64(x, sh)))) {
     pi= i + __builtin_ctz(mk);
     break;
    }
    acc= _mm256_or_si256(acc, x);
   }
   if(pi < 0) {
    u64 any= 0;
    for(; i < nr; ++i) {
     if(v[i] >> b & 1) {
      pi= i;
      break;
     }
     any|= v[i];
    }
    if(pi < 0) {  // 要のない列
     if(stop) return false;
     alignas(32) u64 t[4];
     _mm256_store_si256((__m256i*)t, acc);
     any|= t[0] | t[1] | t[2] | t[3];
     any= b == 63 ? 0 : any >> (b + 1) << (b + 1);
     col= any ? (wi << 6) + __builtin_ctzll(any) : (wi + 1) << 6;
     continue;
    }
   }
   const int p= rem[pi];
   u64* up= urow(prow.size());
   const int t= prow.size(), e1= end_of(t + 1);
   copy_to(up, row(p), wi, e1);
   for(u64 m= f[pi]; m; m&= m - 1) xor_in(up, Ub[__builtin_ctzll(m)], wi, e1);
   if(aug >= 0) up[(aug + t) >> 6]|= u64(1) << ((aug + t) & 63);
   // pi より後ろの行のうち列 col が 1 のものに、この要の印を付けて v に要の行の語を足す (pi より前の行は、この列が 0)。
   const u64 uw= v[pi], ub= u64(1) << bt;
   i= pi + 1;
   for(; i < nr && (i & 3); ++i) {
    const u64 m= -(v[i] >> b & 1);
    v[i]^= uw & m, f[i]|= ub & m;
   }
   const __m256i z= _mm256_setzero_si256(), U4= _mm256_set1_epi64x(uw), B4= _mm256_set1_epi64x(ub);
   for(; i + 4 <= nr; i+= 4) {
    const __m256i x= _mm256_loadu_si256((const __m256i*)(v.data() + i)), m= _mm256_cmpgt_epi64(z, _mm256_sll_epi64(x, sh));
    _mm256_storeu_si256((__m256i*)(v.data() + i), _mm256_xor_si256(x, _mm256_and_si256(m, U4)));
    _mm256_storeu_si256((__m256i*)(f.data() + i), _mm256_or_si256(_mm256_loadu_si256((const __m256i*)(f.data() + i)), _mm256_and_si256(m, B4)));
   }
   for(; i < nr; ++i) {
    const u64 m= -(v[i] >> b & 1);
    v[i]^= uw & m, f[i]|= ub & m;
   }
   prow.push_back(p), pcol.push_back(col);
   Ub[bt++]= up;
   --nr, rem[pi]= rem[nr], v[pi]= v[nr], f[pi]= f[nr];
   // 16 本たまったら残りの行に足す。足したあとの行の語 col >> 6 は v と同じなので、同じ語なら v はそのまま使える。
   if(++col; bt == K) flush(Ub, bt, col >> 6, nr, f.data(), T1, T2), bt= 0;
  }
  if(bt) flush(Ub, bt, std::min(col, lim) >> 6, nr, f.data(), T1, T2);
  rem.resize(nr);
  return true;
 }
 // 残りの行の語 wi を v に集める。貯めている要 bt 本の分は、要の行の語 wi の 8 本ずつの XOR の表 2 枚で足す。
 void gather(int wi, int nr, u64* v, const u64* f, const u64* const* Ub, int bt) const {
  if(!bt) {
   for(int i= 0; i < nr; ++i) v[i]= a[(size_t)rem[i] * W + wi];
   return;
  }
  u64 t1[256], t2[256];
  t1[0]= t2[0]= 0;
  for(int m= 1; m < (1 << std::min(bt, 8)); ++m) t1[m]= t1[m & (m - 1)] ^ Ub[__builtin_ctz(m)][wi];
  for(int m= 1; m < (1 << std::max(bt - 8, 0)); ++m) t2[m]= t2[m & (m - 1)] ^ Ub[8 + __builtin_ctz(m)][wi];
  for(int i= 0; i < nr; ++i) v[i]= a[(size_t)rem[i] * W + wi] ^ t1[f[i] & 255] ^ t2[f[i] >> 8];
 }
 // T の行 m (1 <= m < 2^cnt) を、Ub[s0 + i] (m の bit i が立つもの) の XOR にする。語 [w0, we) だけ。
 void build(std::vector<u64>& T, const u64* const* Ub, int s0, int cnt, int w0, int we) const {
  std::fill(T.begin() + w0, T.begin() + we, 0);
  for(int m= 1; m < (1 << cnt); ++m) xor_to(&T[(size_t)m * W], &T[(size_t)(m & (m - 1)) * W], Ub[s0 + __builtin_ctz(m)], w0, we);
 }
 // 行 rows(i) = {x, m} (i < cnt) の語 [w0, we) に、要の行 Ub[s] (m の bit s が立つもの) を足す。行が 64 本以下なら、表を
 // 作らずに要の行を 1 本ずつ足す。
 template <class R> void apply(int cnt, R&& rows, const u64* const* Ub, int bt, int w0, int we, std::vector<u64>& T1, std::vector<u64>& T2) {
  if(cnt <= 64) {
   for(int i= 0; i < cnt; ++i) {
    auto [x, m]= rows(i);
    for(; m; m&= m - 1) xor_in(x, Ub[__builtin_ctzll(m)], w0, we);
   }
   return;
  }
  if(T1.empty()) T1.resize((size_t)W << 8), T2.resize((size_t)W << 8);
  build(T1, Ub, 0, std::min(bt, 8), w0, we);
  if(bt > 8) build(T2, Ub, 8, bt - 8, w0, we);
  for(int i= 0; i < cnt; ++i) {
   auto [x, m]= rows(i);
   if(!m) continue;
   if(m >> 8) xor2_in(x, &T1[(size_t)(m & 255) * W], &T2[(size_t)(m >> 8) * W], w0, we);
   else xor_in(x, &T1[(size_t)(m & 255) * W], w0, we);
  }
 }
 // 貯めた bt 本の要の更新を、残りの行の語 [wb, end_of(rank())) に足す。
 void flush(const u64* const* Ub, int bt, int wb, int nr, u64* f, std::vector<u64>& T1, std::vector<u64>& T2) {
  apply(
      nr,
      [&](int i) {
       const u64 m= f[i];
       f[i]= 0;
       return std::pair{row(rem[i]), m};
      },
      Ub, bt, wb, end_of(rank()), T1, T2);
 }
 // 前進消去のあと、要の行を RREF に直す (U の行 t が、要の列 pcol[t] で 1、ほかの要の列で 0 になる)。要の行を下から 16 本ずつ
 // 組にし、組の中を下から消してから、組の要の行の 2 枚の表で、それより上の要の行の組の要の列を消す。
 void back() {
  const int r= prow.size();
  std::vector<u64> T1, T2;
  for(int t1= r; t1 > 0; t1-= K) {
   // 下の組を足したあとの上の要の行は、右の部分が要の数 r の分まで広がるので、語は end_of(r) まで読み書きする。
   const int t0= std::max(0, t1 - K), bt= t1 - t0, w0= pcol[t0] >> 6, e1= end_of(r);
   for(int t= t1 - 2; t >= t0; --t) {
    u64* x= urow(t);
    for(int s= t + 1; s < t1; ++s)
     if(bit(x, pcol[s])) xor_in(x, urow(s), pcol[s] >> 6, e1);
   }
   const u64* Ub[K];
   for(int s= 0; s < bt; ++s) Ub[s]= urow(t0 + s);
   const bool consecutive= pcol[t1 - 1] - pcol[t0] == bt - 1;
   apply(
       t0,
       [&](int t) {
        u64* x= urow(t);
        u64 m= 0;
        if(consecutive) {  // 組の要の列が続いているなら、2 語から 16 bit をまとめて取る
         const int j= pcol[t0], wj= j >> 6, b= j & 63;
         u64 y= x[wj] >> b;
         if(b && wj + 1 < W) y|= x[wj + 1] << (64 - b);
         m= y & ((u64(1) << bt) - 1);
        } else
         for(int s= 0; s < bt; ++s) m|= u64(bit(x, pcol[t0 + s])) << s;
        return std::pair{x, m};
       },
       Ub, bt, w0, e1, T1, T2);
  }
 }
};
// C の語 [p, p + 4 NV) に A B の分を足す。B の行は KW 語、C の行は CW 語おき。T は 8 × 256 行 × 16 語の作業領域。
// bw は帯の中の B の語の数 (4 NV より小さければ、B の行を 0 で埋めた写しから表を作る)。
template <int NV> void mul_band(const u64* A, int MW, int n, int m, const u64* B, int KW, int p, int bw, u64* C, int CW, u64* T) {
 constexpr int PW= 16;
 for(int kw= 0; kw < MW; ++kw) {
  for(int g= 0; g < 8; ++g) {
   const int r0= kw * 64 + g * 8, cnt= std::clamp(m - r0, 0, 8);  // m 行目より後ろの B の行は 0 (A のその bit も 0)
   u64* tg= T + (size_t)g * 256 * PW;
   for(int v= 0; v < NV; ++v) _mm256_storeu_si256((__m256i*)(tg + 4 * v), _mm256_setzero_si256());
   alignas(32) u64 tmp[8][4 * NV];
   const u64* br[8];
   for(int j= 0; j < cnt; ++j) {
    const u64* s= B + (size_t)(r0 + j) * KW + p;
    if(bw == 4 * NV) br[j]= s;
    else {
     std::memcpy(tmp[j], s, bw * sizeof(u64));
     std::memset(tmp[j] + bw, 0, (4 * NV - bw) * sizeof(u64));
     br[j]= tmp[j];
    }
   }
   for(int s= 1; s < (1 << cnt); ++s) {
    u64* t= tg + (size_t)s * PW;
    const u64 *q= tg + (size_t)(s & (s - 1)) * PW, *b= br[__builtin_ctz(s)];
    for(int v= 0; v < NV; ++v) _mm256_storeu_si256((__m256i*)(t + 4 * v), _mm256_xor_si256(_mm256_loadu_si256((const __m256i*)(q + 4 * v)), _mm256_loadu_si256((const __m256i*)(b + 4 * v))));
   }
  }
  for(int i= 0; i < n; ++i) {
   const u64 x= A[(size_t)i * MW + kw];
   if(!x) continue;
   u64* c= C + (size_t)i * CW + p;
   __m256i c0= _mm256_loadu_si256((const __m256i*)c), c1, c2, c3;
   if constexpr(NV > 1) c1= _mm256_loadu_si256((const __m256i*)(c + 4));
   if constexpr(NV > 2) c2= _mm256_loadu_si256((const __m256i*)(c + 8));
   if constexpr(NV > 3) c3= _mm256_loadu_si256((const __m256i*)(c + 12));
   for(int g= 0; g < 8; ++g) {
    const u64* t= T + ((size_t)g * 256 + (x >> (8 * g) & 255)) * PW;
    c0= _mm256_xor_si256(c0, _mm256_loadu_si256((const __m256i*)t));
    if constexpr(NV > 1) c1= _mm256_xor_si256(c1, _mm256_loadu_si256((const __m256i*)(t + 4)));
    if constexpr(NV > 2) c2= _mm256_xor_si256(c2, _mm256_loadu_si256((const __m256i*)(t + 8)));
    if constexpr(NV > 3) c3= _mm256_xor_si256(c3, _mm256_loadu_si256((const __m256i*)(t + 12)));
   }
   _mm256_storeu_si256((__m256i*)c, c0);
   if constexpr(NV > 1) _mm256_storeu_si256((__m256i*)(c + 4), c1);
   if constexpr(NV > 2) _mm256_storeu_si256((__m256i*)(c + 8), c2);
   if constexpr(NV > 3) _mm256_storeu_si256((__m256i*)(c + 12), c3);
  }
 }
}
}
class BitMatrix {
 using u64= unsigned long long;
 int h= 0, w= 0, W= 0;  // W は行の語の数
 std::vector<u64> a;
 // 前進消去と RREF のあとの e から核の基底を作る。要の行 t を行 pcol[t] に、自由な列 g の行に e_g を置いた w × w の行列 Z の
 // 列 g が、g に対する基底のベクトル。Z を 64 × 64 の小行列ごとに転置して行で取り出す。Z の行 j は列 j より左が 0 なので、
 // 行の塊が列の塊より下の小行列は 0 で、要の行のない行の塊は対角の小行列のほかは 0 なので、それらは飛ばす。
 static BitMatrix kernel_of(const bit_matrix_internal::Elim& e, int w) {
  const int r= e.rank(), Wd= (w + 63) >> 6;
  BitMatrix k(w - r, w);
  std::vector<int> pt(w, -1), id(w, -1);  // 行 j の要の番号と、自由な列 g の基底の番号
  std::vector<char> hasp(Wd), hasf(Wd);
  for(int t= 0; t < r; ++t) pt[e.pcol[t]]= t, hasp[e.pcol[t] >> 6]= 1;
  for(int g= 0, d= 0; g < w; ++g)
   if(pt[g] < 0) id[g]= d++, hasf[g >> 6]= 1;
  const u64 last= (w & 63) ? (u64(1) << (w & 63)) - 1 : ~u64(0);  // 列 w より右 (b の列) を落とす
  u64 x[64];
  for(int q= 0; q < Wd; ++q) {
   if(!hasf[q]) continue;
   const u64 mq= q == Wd - 1 ? last : ~u64(0);
   for(int bi= 0; bi <= q; ++bi) {
    if(!hasp[bi] && bi != q) continue;
    const int rn= std::min(64, w - bi * 64);
    for(int c= 0; c < rn; ++c) {
     const int j= bi * 64 + c;
     x[c]= pt[j] >= 0 ? e.urow(pt[j])[q] & mq : bi == q ? u64(1) << c : 0;
    }
    std::fill(x + rn, x + 64, 0);
    bit_matrix_internal::transpose64(x);
    const int cn= std::min(64, w - q * 64);
    for(int c= 0; c < cn; ++c)
     if(id[q * 64 + c] >= 0) k[id[q * 64 + c]][bi]= x[c];
   }
  }
  return k;
 }
public:
 BitMatrix()= default;
 BitMatrix(int h, int w): h(h), w(w), W((w + 63) >> 6), a((size_t)h * W) {}
 // 各行を '0' と '1' の文字列で与える (文字の下位 bit を見る)。行の長さはそろえること。
 explicit BitMatrix(const std::vector<std::string>& s): BitMatrix(s.size(), s.empty() ? 0 : s[0].size()) {
  for(int i= 0; i < h; ++i) {
   assert((int)s[i].size() == w);
   bit_matrix_internal::pack_row(s[i].data(), w, (*this)[i]);
  }
 }
 static BitMatrix identity(int n) {
  BitMatrix r(n, n);
  for(int i= 0; i < n; ++i) r.set(i, i);
  return r;
 }
 int height() const { return h; }
 int width() const { return w; }
 // 行の語の数 (⌈w / 64⌉)。
 int words() const { return W; }
 // 行 i の語の列。j 列目は語 j / 64 の bit j % 64 で、w 列目より後ろの bit は 0 のままにすること。
 u64* operator[](int i) { return a.data() + (size_t)i * W; }
 const u64* operator[](int i) const { return a.data() + (size_t)i * W; }
 bool get(int i, int j) const { return (*this)[i][j >> 6] >> (j & 63) & 1; }
 void set(int i, int j, bool v= true) {
  u64& x= (*this)[i][j >> 6];
  const u64 b= u64(1) << (j & 63);
  x= v ? x | b : x & ~b;
 }
 void flip(int i, int j) { (*this)[i][j >> 6]^= u64(1) << (j & 63); }
 // 各行を '0' と '1' の文字列にする。
 std::vector<std::string> to_strings() const {
  std::vector<std::string> s(h);
  for(int i= 0; i < h; ++i) {
#if __cpp_lib_string_resize_and_overwrite >= 202110L
   s[i].resize_and_overwrite(w, [&](char* p, size_t n) {
    bit_matrix_internal::unpack_row((*this)[i], w, p);
    return n;
   });
#else
   s[i].resize(w);
   bit_matrix_internal::unpack_row((*this)[i], w, s[i].data());
#endif
  }
  return s;
 }
 // 64 × 64 の小行列ごとに転置する。
 BitMatrix transpose() const {
  BitMatrix t(w, h);
  u64 x[64];
  for(int bi= 0; bi < h; bi+= 64) {
   const int rn= std::min(64, h - bi);
   for(int q= 0; q < W; ++q) {
    for(int r= 0; r < rn; ++r) x[r]= a[(size_t)(bi + r) * W + q];
    std::fill(x + rn, x + 64, 0);
    bit_matrix_internal::transpose64(x);
    const int cn= std::min(64, w - q * 64);
    for(int c= 0; c < cn; ++c) t.a[(size_t)(q * 64 + c) * t.W + (bi >> 6)]= x[c];
   }
  }
  return t;
 }
 BitMatrix& operator+=(const BitMatrix& r) {
  assert(h == r.h && w == r.w);
  for(size_t i= 0; i < a.size(); ++i) a[i]^= r.a[i];
  return *this;
 }
 friend BitMatrix operator+(BitMatrix l, const BitMatrix& r) { return l+= r; }
 friend BitMatrix operator*(const BitMatrix& l, const BitMatrix& r) {
  using namespace bit_matrix_internal;
  assert(l.w == r.h);
  const int n= l.h, m= l.w, KW= r.W;
  BitMatrix c(n, r.w);
  if(!n || !m || !KW) return c;
  if(n <= 64) {  // A の 1 の bit ごとに B の行を足す
   for(int i= 0; i < n; ++i)
    for(int q= 0; q < l.W; ++q)
     for(u64 x= l[i][q]; x; x&= x - 1) xor_n(c[i], r[q * 64 + __builtin_ctzll(x)], KW);
   return c;
  }
  // C の行の語の数を 4 の倍数にした作業領域で求め、違えば最後に詰める。
  const int CW= round4(KW);
  std::vector<u64> buf(CW == KW ? 0 : (size_t)n * CW), T((size_t)8 * 256 * 16);
  u64* C= CW == KW ? c.a.data() : buf.data();
  for(int p= 0; p < KW; p+= 16) {
   const int bw= std::min(16, KW - p);
   switch((bw + 3) >> 2) {
    case 1: mul_band<1>(l.a.data(), l.W, n, m, r.a.data(), KW, p, bw, C, CW, T.data()); break;
    case 2: mul_band<2>(l.a.data(), l.W, n, m, r.a.data(), KW, p, bw, C, CW, T.data()); break;
    case 3: mul_band<3>(l.a.data(), l.W, n, m, r.a.data(), KW, p, bw, C, CW, T.data()); break;
    default: mul_band<4>(l.a.data(), l.W, n, m, r.a.data(), KW, p, bw, C, CW, T.data());
   }
  }
  if(CW != KW)
   for(int i= 0; i < n; ++i) std::memcpy(c[i], C + (size_t)i * CW, KW * sizeof(u64));
  return c;
 }
 BitMatrix& operator*=(const BitMatrix& r) { return *this= *this * r; }
 // 列ベクトルとの積 A v。v の各要素は 0 か 1。
 friend std::vector<int> operator*(const BitMatrix& l, const std::vector<int>& v) {
  assert(l.w == (int)v.size());
  std::vector<u64> x(l.W);
  for(int j= 0; j < l.w; ++j) x[j >> 6]|= u64(v[j] & 1) << (j & 63);
  std::vector<int> r(l.h);
  for(int i= 0; i < l.h; ++i) {
   u64 s= 0;
   for(int q= 0; q < l.W; ++q) s^= l[i][q] & x[q];
   r[i]= __builtin_parityll(s);
  }
  return r;
 }
 // A^k。二分累乗。
 BitMatrix pow(uint64_t k) const {
  assert(h == w);
  BitMatrix r= identity(h), b= *this;
  for(; k; k>>= 1) {
   if(k & 1) r= r * b;
   if(k > 1) b= b * b;
  }
  return r;
 }
 bool det() const {
  assert(h == w);
  const int n= h;
  if(n <= 64) {  // 1 語の行のまま掃き出す
   u64 r[64];
   for(int i= 0; i < n; ++i) r[i]= a[i];
   for(int c= 0; c < n; ++c) {
    int p= c;
    while(p < n && !(r[p] >> c & 1)) ++p;
    if(p == n) return false;
    std::swap(r[p], r[c]);
    for(int i= c + 1; i < n; ++i) r[i]^= r[c] & -(r[i] >> c & 1);
   }
   return true;
  }
  bit_matrix_internal::Elim e(n, bit_matrix_internal::round4(W), a.data(), W);
  return e.run(n, true);
 }
 int rank() const {
  using namespace bit_matrix_internal;
  if(!h || !w) return 0;
  if(w <= 64) {  // 行を 1 語ずつ、最上位の bit ごとの基底に入れる
   u64 basis[64]= {};
   int r= 0;
   for(int i= 0; i < h; ++i)
    for(u64 x= a[i]; x;) {
     const int b= 63 - __builtin_clzll(x);
     if(!basis[b]) {
      basis[b]= x;
      if(++r == w) return r;
      break;
     }
     x^= basis[b];
    }
   return r;
  }
  if(h <= 64) {  // 行を 1 本ずつ、前の要の列で消してから最も左の 1 を要にする
   std::vector<u64> buf((size_t)h * W);
   int pc[64], r= 0;
   for(int i= 0; i < h; ++i) {
    u64* x= buf.data() + (size_t)r * W;
    std::memcpy(x, (*this)[i], W * sizeof(u64));
    for(int t= 0; t < r; ++t)
     if(bit(x, pc[t])) xor_n(x + (pc[t] >> 6), buf.data() + (size_t)t * W + (pc[t] >> 6), W - (pc[t] >> 6));
    int q= 0;
    while(q < W && !x[q]) ++q;
    if(q < W) pc[r++]= q * 64 + __builtin_ctzll(x[q]);
   }
   return r;
  }
  // 行が列の 2 倍より多ければ転置して消す (要ごとに残りの行をなめる分が、行の少ない向きで減る)。
  if(h > 2 * w) return transpose().rank();
  Elim e(h, round4(W), a.data(), W);
  e.run(w, false);
  return e.rank();
 }
 // A x = b の解を 1 つ返す (自由な変数は 0)。解がなければ nullopt。b の各要素は 0 か 1。
 std::optional<std::vector<int>> solve(const std::vector<int>& b) const {
  using namespace bit_matrix_internal;
  assert((int)b.size() == h);
  Elim e(h, round4((w + 64) >> 6), a.data(), W);  // b を列 w に置く
  for(int i= 0; i < h; ++i) e.row(i)[w >> 6]|= u64(b[i] & 1) << (w & 63);
  e.run(w, false);
  for(int i: e.rem)
   if(bit(e.row(i), w)) return std::nullopt;
  // 後退代入。要の行 t は要の列より左が 0 なので、x と重なるのは後ろの要の列だけ。
  std::vector<u64> x(W);
  for(int t= e.rank(); t--;) {
   const u64* u= e.urow(t);
   const int c= e.pcol[t];
   u64 s= 0;
   for(int q= c >> 6; q < W; ++q) s^= u[q] & x[q];
   x[c >> 6]|= u64(bit(u, w) ^ __builtin_parityll(s)) << (c & 63);
  }
  std::vector<int> res(w);
  for(int j= 0; j < w; ++j) res[j]= bit(x.data(), j);
  return res;
 }
 // A x = b の解の 1 つ x と、A v = 0 の解空間の基底 (行列の各行) を返す。解がなければ nullopt。
 std::optional<std::pair<std::vector<int>, BitMatrix>> solve_with_kernel(const std::vector<int>& b) const {
  using namespace bit_matrix_internal;
  assert((int)b.size() == h);
  Elim e(h, round4((w + 64) >> 6), a.data(), W);
  for(int i= 0; i < h; ++i) e.row(i)[w >> 6]|= u64(b[i] & 1) << (w & 63);
  e.run(w, false);
  for(int i: e.rem)
   if(bit(e.row(i), w)) return std::nullopt;
  e.back();
  std::vector<int> x(w);
  for(int t= 0; t < e.rank(); ++t) x[e.pcol[t]]= bit(e.urow(t), w);
  return std::pair{std::move(x), kernel_of(e, w)};
 }
 // A v = 0 の解空間の基底を、行列の各行にして返す (w - rank 行)。
 BitMatrix kernel() const {
  using namespace bit_matrix_internal;
  Elim e(h, round4(W), a.data(), W);
  e.run(w, false);
  e.back();
  return kernel_of(e, w);
 }
 // 逆行列。正則でなければ nullopt。
 std::optional<BitMatrix> inverse() const {
  using namespace bit_matrix_internal;
  assert(h == w);
  const int n= h;
  // [A | I] を前進消去してから RREF に直すと、要の行 t (要の列は t) の右半分 R が、逆行列の行 t になる。I は語 W から置き、
  // 列を付け替える (Elim::aug) ので、R の列 s は元の行 prow[s] を表す。逆行列の列 prow[s] が R の列 s なので、R を転置して行を
  // 並べ替え、転置し直す。
  Elim e(n, round4(2 * W), a.data(), W);
  e.aug= 64 * W;
  if(!e.run(n, true)) return std::nullopt;
  e.back();
  BitMatrix r(n, n);
  for(int t= 0; t < n; ++t) std::memcpy(r[t], e.urow(t) + W, W * sizeof(u64));
  const BitMatrix rt= r.transpose();
  for(int s= 0; s < n; ++s) std::memcpy(r[e.prow[s]], rt[s], W * sizeof(u64));
  return r.transpose();
 }
 friend bool operator==(const BitMatrix&, const BitMatrix&)= default;
};
