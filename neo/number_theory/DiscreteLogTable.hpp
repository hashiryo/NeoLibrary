#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <utility>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 素数 p と原始根 g を固定して、g^x ≡ a (mod p) となる 0 ≤ x < p - 1 を表で引く (maspy の 1 問 O(1) の方法)。
// a ≡ num / den (mod p) で |num| ≤ L、den < 2^16 となる分数を、a / p を NB 個に区切った区画ごとの Farey の表 (Stern-Brocot
// の中間数で埋める) で引き、log a = log |num| + [num < 0] (p - 1) / 2 - log den とする。区画 floor(a NB / p) は
// floor(2^64 NB / p) との掛け算の上位 64 bit で求める (誤差は a 2^{-64} 未満で、0 < a < p では a NB / p の小数部が 1/p
// 以上なので割り算と一致する)。log の表は L までの整数について次の順に作る。
// - 最初の K1 個の素数: p - 1 の素因数 r ごとの Pohlig-Hellman。位数 r の baby-step giant-step は、z^j を 4 本の掛け算の鎖で
//   作って (鍵 << 32 | j) を 4 個ずつの籠に置き (籠の中は分岐なしで比べる)、K1 個の的を並べて引く。B = √(K1 r / 2)。
// - √p までの残りの素数 q: q c^{-t} (c = g^K) を Farey の表で分数にし、分子と分母の log が分かるまで t を進める。SM 以下の
//   素数の冪を掛けた smooth な数の log もここで広げる。
// - √p までの合成数は最小の素因数から求める。√p より大きい i は、2 か 3 で割り切れれば log 2 + log(i / 2) か
//   log 3 + log(i / 3) (順に読める)、そうでなければ p = (p / i) i + (p mod i) から log i = log(p mod i) + log(-1) - log(p / i)。
//   i が 1 増えると p / i は高々 1 しか減らないので、割り算なしで進める。
// p < M0 なら g の冪を 4 本の鎖でたどって全部の log を表にする。表は合わせて 9 MB ほどで、Linux では 2 MB 境界の mmap に
// 置いて huge page を頼む。剰余の掛け算は R = 2^32 の Montgomery。表を作るときは定数と表の場所を局所変数に写す (メンバの
// ままだと log の表への書き込みと別名になりうるので、gcc が輪の中で読み直す)。procon-judge の
// yosupo-discrete-logarithm-fixed-mod の 5 回目の dlf5 を写したもの。
namespace discrete_log_table_internal {
using u8= unsigned char;
using u32= unsigned;
using i32= int;
using u64= unsigned long long;
using u128= unsigned __int128;
// 0 で埋めた bytes の領域。2 MB 以上なら Linux では 2 MB 境界の mmap に置いて huge page を頼む。
class Zeroed {
 void* p= nullptr;
#ifdef __linux__
 size_t mapped= 0;
#endif
public:
 explicit Zeroed(size_t bytes) {
#ifdef __linux__
  constexpr size_t H= size_t(1) << 21;
  if(bytes >= H) {
   const size_t b= (bytes + 4095) & ~size_t(4095);
   void* m= mmap(nullptr, b + H, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
   if(m != MAP_FAILED) {
    char *q= static_cast<char*>(m), *a= reinterpret_cast<char*>((reinterpret_cast<uintptr_t>(q) + H - 1) & ~(H - 1));
    if(a != q) munmap(q, a - q);
    munmap(a + b, q + H - a);
    madvise(a, b, MADV_HUGEPAGE);
    p= a, mapped= b;
    return;
   }
  }
#endif
  p= calloc(bytes ? bytes : 1, 1);
 }
 ~Zeroed() {
#ifdef __linux__
  if(mapped) {
   munmap(p, mapped);
   return;
  }
#endif
  free(p);
 }
 Zeroed(const Zeroed&)= delete;
 Zeroed& operator=(const Zeroed&)= delete;
 u32* get() const { return static_cast<u32*>(p); }
};
// R = 2^32 の Montgomery。値は [0, 2p) のまま持つ (p < 2^30)
struct Mont {
 u32 p, pinv, r1, r2;
 explicit Mont(u32 p_): p(p_) {
  u32 x= (3 * p) ^ 2;  // 下位 5 ビットが正しい
  for(int i= 0; i < 3; ++i) x*= 2 - p * x;
  pinv= x, r1= (0u - p) % p, r2= (u32)((u64)r1 * r1 % p);
 }
 u32 mul(u32 a, u32 b) const {
  const u64 t= (u64)a * b;
  return (u32)(t >> 32) + p - (u32)((u64)((u32)t * pinv) * p >> 32);
 }
 u32 to(u32 a) const { return mul(a, r2); }
 u32 norm(u32 a) const { return a >= p ? a - p : a; }
 u32 from(u32 a) const { return norm(mul(a, 1)); }
 u32 pow(u32 b, u64 e) const {  // b は Montgomery 形
  u32 r= r1;
  for(; e; e>>= 1, b= mul(b, b))
   if(e & 1) r= mul(r, b);
  return r;
 }
 u32 powi(u32 b, u64 e) const { return from(pow(to(b), e)); }  // 普通の形で
};
inline u32 inv_mod(u32 a, u32 m) {
 long long x= 1, y= 0;
 for(u32 b= m; b;) {
  const u32 t= a / b;
  a-= t * b, std::swap(a, b);
  x-= (long long)t * y, std::swap(x, y);
 }
 return (u32)(x < 0 ? x + m : x);
}
// 位数 n の元 z の離散対数を引く baby-step giant-step。表は z^j (j < B) を (鍵 << 32 | j) で 4 個ずつの籠に前から詰めて
// 持ち、鍵は Montgomery 形を [0, p) に直した値 (0 にならないので 0 を空きに使う)
struct Bsgs {
 const Mont& M;
 u32 n, B, bmask, step;  // step = z^{-B} (Montgomery 形)
 std::vector<u64> tab;
 static u32 hash(u32 k) { return (k * 0x9E3779B1u) >> 8; }
 Bsgs(const Mont& M_, u32 z, u32 n_, u32 B_): M(M_), n(n_), B(std::min(B_, n_)) {
  u32 nb= 4;
  while(nb * 6 < B * 5) nb*= 2;  // 籠あたり 1.2 個以下
  bmask= nb - 1, tab.assign(4 * (size_t)nb, 0);
  const u32 z1= M.to(z), z2= M.mul(z1, z1), z3= M.mul(z2, z1), z4= M.mul(z3, z1);
  u32 x[4]= {M.r1, z1, z2, z3};
  for(u32 j= 0; j < B; j+= 4)
   for(u32 d= 0; d < 4; ++d) {
    if(j + d < B) {
     const u32 k= M.norm(x[d]);
     for(u32 b= hash(k) & bmask;; b= (b + 1) & bmask) {
      u64* e= &tab[4 * (size_t)b];
      const int t= (e[0] != 0) + (e[1] != 0) + (e[2] != 0) + (e[3] != 0);  // 前から詰めるので空きの先頭
      if(t < 4) {
       e[t]= (u64)k << 32 | (j + d);
       break;
      }
     }
    }
    x[d]= M.mul(x[d], z4);
   }
  step= M.pow(z1, n - B);
 }
 u32 look(u32 k) const {
  for(u32 b= hash(k) & bmask;; b= (b + 1) & bmask) {
   const u64* e= &tab[4 * (size_t)b];
   u32 hit= ~0u;
   for(int t= 0; t < 4; ++t) hit= (u32)(e[t] >> 32) == k ? (u32)e[t] : hit;
   if(hit != ~0u || !e[3]) return hit;
  }
 }
 // z^{out[i]} = c[i] (c は普通の形) を cnt 個 (64 以下) 並べて求める
 void solve(const u32* c, u32 cnt, u32* out) const {
  u32 x[64], r[64], act[64], na= cnt;
  for(u32 i= 0; i < cnt; ++i) x[i]= M.to(c[i]), r[i]= 0, act[i]= i;
  while(na)
   for(u32 k= 0; k < na;) {
    const u32 i= act[k], v= look(M.norm(x[i]));
    if(v != ~0u) {
     out[i]= r[i] + v, act[k]= act[--na];
     continue;
    }
    x[i]= M.mul(x[i], step), r[i]+= B, ++k;
   }
 }
};
constexpr u32 M0= 2000000, NB= 1000000, K1= 10, SM= 200, WK= 0x2545F491u;  // WK は 2 段目の歩みの指数
}
class DiscreteLogTable {
 using u32= discrete_log_table_internal::u32;
 using i32= discrete_log_table_internal::i32;
 using u64= discrete_log_table_internal::u64;
 using u128= discrete_log_table_internal::u128;
 u32 p, pm1, halfp, L;  // L は log の表の上限
 bool small;
 u64 bk= 0;  // floor(2^64 NB / p)
 discrete_log_table_internal::Zeroed dmem, fmem;
 u32 *dl, *fa;
 static u32 lim(u32 p) { return std::max<u32>(1300000, (u32)((u64)p * 13 / 10000 + 1)); }
 void build(u32 g);
public:
 // 2 ≤ p < 2^30 の素数 p と、その原始根 g。p が素数であることは確かめない
 DiscreteLogTable(unsigned p_, unsigned g): p(p_), pm1(p_ - 1), halfp(p_ / 2), L(lim(p_)), small(p_ < discrete_log_table_internal::M0), dmem(4 * (size_t)(small ? p_ : L + 1)), fmem(small ? 0 : 4 * (size_t)discrete_log_table_internal::NB) {
  assert(2 <= p && p < (1u << 30));
  dl= dmem.get(), fa= fmem.get();
  build(g % p);
 }
 // g^x ≡ a (mod p) となる 0 ≤ x < p - 1。a ≡ 0 (mod p) なら -1
 long long operator()(unsigned long long a) const {
  if(a >= p) a%= p;
  if(small) return a ? (long long)dl[a] : -1;
  const u32 b= (u32)a, f= fa[(u32)((u128)b * bk >> 64)], den= f & 0xffff;
  const i32 num= (i32)(den * b - (f >> 16) * p);
  const u32 s= (u32)(num >> 31), an= ((u32)num ^ s) - s;
  u32 r= dl[an] + (halfp & s) + pm1 - dl[den];  // < 3 (p - 1)
  u32 t= r - pm1;
  r= t + (pm1 & (u32)((i32)t >> 31));
  t= r - pm1;
  r= t + (pm1 & (u32)((i32)t >> 31));
  return a ? (long long)r : -1;
 }
};
inline void DiscreteLogTable::build(u32 g) {
 using namespace discrete_log_table_internal;
 // メンバを局所変数に写す (log の表への書き込みと別名にならない)
 const u32 p= this->p, pm1= this->pm1, halfp= this->halfp, L= this->L;
 u32* const dl= this->dl;
 u32* const fa= this->fa;
 auto addm= [pm1](u32 a, u32 b) {  // (a + b) mod (p - 1)、a, b < p - 1
  const u32 t= a + b - pm1;
  return t + (pm1 & (u32)((i32)t >> 31));
 };
 auto negm= [pm1](u32 a) { return a ? pm1 - a : 0; };
 const Mont M(p);
 if(small) {  // g の冪を 4 本の鎖でたどる
  const u32 g1= M.to(g), g2= M.mul(g1, g1), g3= M.mul(g2, g1), g4= M.mul(g3, g1);
  u32 x[4]= {g1, g2, g3, g4};
  for(u32 i= 1; i < pm1; i+= 4)
   for(u32 d= 0; d < 4; ++d) {
    if(i + d < pm1) dl[M.from(x[d])]= i + d;
    x[d]= M.mul(x[d], g4);
   }
  return;
 }
 const u64 bk= this->bk= (u64)(((u128)NB << 64) / p);
 auto bucket= [bk](u32 a) { return (u32)((u128)a * bk >> 64); };
 // Farey の表。区画 [x, y) を f1 と f2 の間の分数で埋める (分数は分子 << 16 | 分母)
 auto rec= [&](auto& self, u32 f1, u32 f2, u32 x, u32 y) -> void {
  const u32 f3= f1 + f2;
  u32 l= (u32)((((u64)p * (f3 >> 16) - L) * NB - 1) / ((u64)p * (f3 & 0xffff)) + 1);
  u32 r= (u32)((((u64)p * (f3 >> 16) + L) * NB) / ((u64)p * (f3 & 0xffff)));
  l= std::max(l, x), r= std::min(r, y);
  if(x < l) self(self, f1, f3, x, l);
  for(u32 i= l; i < r; ++i) fa[i]= f3;
  if(r < y) self(self, f3, f2, r, y);
 };
 {
  // [0, fx) は 0/1、[fy, NB/2) は 1/2 で足りる。後ろ半分は a → p - a で写す
  const u32 fx= (u32)std::min<u64>((u64)L * NB / p, NB / 2), fy= (u32)(((u64)(p - L) * NB - 1) / (p * 2ull) + 1);
  for(u32 i= 0; i < fx; ++i) fa[i]= 1;
  if(fx < fy) rec(rec, 1u, 0x10002u, fx, fy);
  for(u32 i= fy; i < NB / 2; ++i) fa[i]= 0x10002u;
  for(u32 i= NB / 2; i < NB; ++i) fa[i]= (fa[NB - 1 - i] * 0xffff0001u ^ 0xffff0000u) + 0x10000u;
 }
 // SQ (SQ^2 ≥ p) までの最小の素因数
 u32 SQ= (u32)std::sqrt((double)p);
 while((u64)SQ * SQ < p) ++SQ;
 std::vector<u32> spf(SQ + 1, 0), primes;
 primes.reserve(3600);
 for(u32 i= 2; i <= SQ; ++i) {
  if(!spf[i]) spf[i]= i, primes.push_back(i);
  for(u32 q: primes) {
   if(q * i > SQ) break;
   spf[q * i]= q;
   if(q == spf[i]) break;
  }
 }
 // 1 段目: 最初の K1 個の素数の log を Pohlig-Hellman で
 {
  std::vector<std::pair<u32, u32>> fac;
  u32 m= pm1;
  for(u32 q: primes) {
   if((u64)q * q > m) break;
   if(m % q == 0) {
    u32 e= 0;
    do ++e, m/= q;
    while(m % q == 0);
    fac.emplace_back(q, e);
   }
  }
  if(m > 1) fac.emplace_back(m, 1);
  const u32 ginv= inv_mod(g, p);
  std::array<u32, K1> ans{};
  u32 mod= 1;
  for(auto [q, e]: fac) {
   assert(M.powi(g, pm1 / q) != 1);  // g は原始根
   u32 qe= 1;
   for(u32 i= 0; i < e; ++i) qe*= q;
   const u32 ee= pm1 / qe, gg= M.powi(ginv, ee);  // gg は位数 qe の部分群の生成元の逆元
   const Bsgs bs(M, M.powi(g, pm1 / q), q, q < 100 ? q : (u32)std::sqrt((double)q * K1 / 2));
   const u32 iv= inv_mod(mod % qe, qe);
   std::array<u32, K1> h, x, c, d;
   for(u32 i= 0; i < K1; ++i) h[i]= M.powi(primes[i], ee), x[i]= 0;
   u32 pe= qe, mul= 1;
   for(u32 j= 0; j < e; ++j) {  // q 進の桁を下から。桁ごとに K1 個を並べて引く
    pe/= q;
    for(u32 i= 0; i < K1; ++i) c[i]= M.powi((u32)((u64)h[i] * M.powi(gg, x[i]) % p), pe);
    bs.solve(c.data(), K1, d.data());
    for(u32 i= 0; i < K1; ++i) x[i]+= d[i] * mul;
    mul*= q;
   }
   for(u32 i= 0; i < K1; ++i) ans[i]+= mod * (u32)((u64)((x[i] + qe - ans[i] % qe) % qe) * iv % qe);  // CRT
   mod*= qe;
  }
  for(u32 i= 0; i < K1; ++i) dl[primes[i]]= ans[i];
 }
 // 2 段目: SQ までの残りの素数と、SM 以下の素数の冪を掛けた smooth な数
 {
  const u32 cinvm= M.pow(M.to(inv_mod(g, p)), WK), wk= WK % pm1;  // c^{-1} (Montgomery 形)
  std::vector<u32> smooth{1};
  smooth.reserve(1 << 14);
  for(u32 q: primes) {
   if(!dl[q])
    for(u32 h= q, t= 0;; ++t, h= M.norm(M.mul(h, cinvm))) {  // h c^{-t} (h は普通の形)
     const u32 f= fa[bucket(h)], den= f & 0xffff;
     const i32 num= (i32)(den * h - (f >> 16) * p);
     const u32 an= (u32)std::abs(num);
     if((den != 1 && !dl[den]) || (an != 1 && !dl[an])) continue;
     dl[q]= addm(addm(addm(dl[an], num < 0 ? halfp : 0), negm(dl[den])), (u32)((u64)t * wk % pm1));
     break;
    }
   if(q > SM) continue;
   smooth.erase(std::remove_if(smooth.begin(), smooth.end(), [&](u32 v) { return v > L / q; }), smooth.end());
   for(size_t j= 0, nn= smooth.size(); j < nn; ++j) {
    u32 prev= smooth[j], ns= prev * q;
    do {
     if(ns <= L / q) smooth.push_back(ns);
     dl[ns]= addm(dl[prev], dl[q]);
     prev= ns, ns*= q;
    } while(prev <= L / q);
   }
  }
 }
 // 3 段目: SQ までの合成数
 for(u32 i= 2; i <= SQ; ++i)
  if(!dl[i]) dl[i]= addm(dl[i / spf[i]], dl[spf[i]]);
 // 4 段目: SQ より大きい i。p = q i + r (r = p mod i < i、q = p / i < i)。i が 1 増えると q は高々 1 減る (i (i + 1) > p)
 u32 i= SQ + 1, q= p / i, r= p - q * i, c6= i % 6;
 const u32 l2= dl[2], l3= dl[3];
 for(; i <= L; ++i) {
  if(c6 % 2 == 0) dl[i]= addm(l2, dl[i >> 1]);
  else if(c6 == 3) dl[i]= addm(l3, dl[(u32)((u64)i * 0xAAAAAAABu >> 33)]);
  else dl[i]= addm(addm(dl[r], halfp), negm(dl[q]));
  c6= c6 == 5 ? 0 : c6 + 1;
  if(r >= q) r-= q;
  else r+= i + 1 - q, --q;
 }
}
