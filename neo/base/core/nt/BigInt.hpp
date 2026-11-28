#pragma once
// 多倍長整数。桁は u32 で 10^9 進法に持ち (下位の桁から並べる)、符号は別に持つ。
// 積は 3 素数の NTT (crt.hpp の crt_engine<Ntt>::conv_at) で求め、Garner で u128 に戻してから桁上げする。
// 割り算は、除数が短いか商が短いときは Knuth の D、大きいときは Newton 法:
// 正規化した除数の、BASE^(m+K) / b の近似 I を倍々に詰め、q = floor(a I / BASE^(m+K)) を最後に 2、3 回直す。
// 近似は常に真の値を下から押さえるので、直す回数は定数で済む。
// 10^9 進法にしたのは、10 進の入出力が桁数に比例する割り算だけで済むため (2^32 進法だと基数の変換に
// もう 1 回多倍長の割り算が要る)。桁の値は 10^9 未満で、NTT の素数 (最小 880803841) の 2 倍未満なので、
// 法を取るのは 1 回の条件付きの引き算で足りる。
// 積の係数の上界 min(n, m) (10^9 - 1)^2 は 3 素数の積 (約 2^89.2) 未満でなければならないが、
// 10^9 進法では min(n, m) が 2^56 桁まで覆うので、実際の大きさでは気にしなくてよい。
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>
#include "neo/base/core/convolve/crt.hpp"
#include "neo/base/core/convolve/ntt.hpp"

namespace bigint_internal {
using u32= unsigned;
using u64= unsigned long long;
using u128= unsigned __int128;
constexpr u32 BASE= 1000000000;
constexpr int DIG= 9;

// 積に使う 3 素数 (crt_internal::PRIMES の先頭 3 本と同じ)。
constexpr u32 Q0= 998244353, Q1= 897581057, Q2= 880803841;

constexpr u64 pow_mod(u64 a, u64 e, u64 p) {
 u64 r= 1;
 while (e) {
  if (e & 1) r= r * a % p;
  a= a * a % p, e>>= 1;
 }
 return r;
}
// Garner の定数: Q0^{-1} mod Q1 と (Q0 Q1)^{-1} mod Q2。
constexpr u32 I01= u32(pow_mod(Q0, Q1 - 2, Q1));
constexpr u64 Q0Q1= u64(Q0) * Q1;
constexpr u32 I012= u32(pow_mod(Q0Q1 % Q2, Q2 - 2, Q2));

// 3 つの素数の下の値から、Q0 Q1 Q2 未満の真の値 (積の係数) を戻す。
inline u128 garner3(u32 r0, u32 r1, u32 r2) {
 const u32 c0= r0;
 const u32 c1= u32(u64((r1 + Q1 - c0 % Q1) % Q1) * I01 % Q1);
 const u64 x01= u64(c0) + u64(c1) * Q0;  // < Q0 Q1 < 2^63
 const u32 c2= u32(u64((r2 + Q2 - x01 % Q2) % Q2) * I012 % Q2);
 return u128(x01) + u128(c2) * Q0Q1;
}

inline void trim(std::vector<u32>& v) {
 while (!v.empty() && !v.back()) v.pop_back();
}

inline int cmp_abs(const std::vector<u32>& a, const std::vector<u32>& b) {
 if (a.size() != b.size()) return a.size() < b.size() ? -1 : 1;
 for (int i= int(a.size()); i--;)
  if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
 return 0;
}

inline std::vector<u32> add_abs(const std::vector<u32>& a, const std::vector<u32>& b) {
 const int n= int(a.size()), m= int(b.size()), e= std::max(n, m);
 std::vector<u32> r(e + 1);
 u64 c= 0;
 for (int i= 0; i < e; ++i) {
  const u64 t= c + (i < n ? a[i] : 0) + (i < m ? b[i] : 0);
  r[i]= u32(t >= BASE ? t - BASE : t), c= t >= BASE;
 }
 r[e]= u32(c);
 trim(r);
 return r;
}

// a >= b のとき、a - b を a に書く。
inline void sub_abs_inplace(std::vector<u32>& a, const std::vector<u32>& b) {
 assert(cmp_abs(a, b) >= 0);
 u64 borrow= 0;
 for (size_t i= 0; i < b.size(); ++i) {
  const long long t= (long long)a[i] - b[i] - (long long)borrow;
  if (t < 0) a[i]= u32(t + BASE), borrow= 1;
  else a[i]= u32(t), borrow= 0;
 }
 for (size_t i= b.size(); i < a.size() && borrow; ++i) {
  if (a[i]) a[i]--, borrow= 0;
  else a[i]= BASE - 1;
 }
 assert(!borrow);
 trim(a);
}

inline std::vector<u32> mul_small(const std::vector<u32>& a, u32 m) {
 if (a.empty() || m == 0) return {};
 std::vector<u32> r(a.size() + 1);
 u64 c= 0;
 for (size_t i= 0; i < a.size(); ++i) {
  const u64 t= u64(a[i]) * m + c;
  r[i]= u32(t % BASE), c= t / BASE;
 }
 r[a.size()]= u32(c);
 trim(r);
 return r;
}

// a を m で割った商。余りは rem に入れる。
inline std::vector<u32> divmod_small(const std::vector<u32>& a, u32 m, u32& rem) {
 std::vector<u32> q(a.size());
 u64 r= 0;
 for (size_t i= a.size(); i--;) {
  const u64 t= r * BASE + a[i];
  q[i]= u32(t / m), r= t % m;
 }
 rem= u32(r);
 trim(q);
 return q;
}

inline std::vector<u32> mul(const std::vector<u32>& a, const std::vector<u32>& b);

// a (n 桁) を b (m 桁) で割った商 q と余り r。n >= m >= 1、b の末尾は 0 でない (Knuth の D)。
inline void knuth_divmod(const u32* a, int n, const u32* b, int m, std::vector<u32>& q, std::vector<u32>& r) {
 if (m == 1) {
  q.assign(n, 0);
  u64 rem= 0;
  for (int i= n; i--;) {
   const u64 t= rem * BASE + a[i];
   q[i]= u32(t / b[0]), rem= t % b[0];
  }
  trim(q), r.assign(1, u32(rem)), trim(r);
  return;
 }
 // 正規化: d = BASE / (b の末尾 + 1)。b の末尾は BASE / 2 以上になる。
 const u32 d= u32(BASE / (u64(b[m - 1]) + 1));
 std::vector<u32> ua(n + 1), ub(m);
 u64 c= 0;
 for (int i= 0; i < n; ++i) {
  const u64 t= u64(a[i]) * d + c;
  ua[i]= u32(t % BASE), c= t / BASE;
 }
 ua[n]= u32(c);
 c= 0;
 for (int i= 0; i < m; ++i) {
  const u64 t= u64(b[i]) * d + c;
  ub[i]= u32(t % BASE), c= t / BASE;
 }
 assert(c == 0 && ub[m - 1] >= BASE / 2);
 q.assign(n - m + 1, 0);
 for (int j= n - m; j >= 0; --j) {
  const u64 num= u64(ua[j + m]) * BASE + ua[j + m - 1];
  u64 qhat= num / ub[m - 1], rhat= num % ub[m - 1];
  while (qhat >= BASE || u128(qhat) * ub[m - 2] > u128(rhat) * BASE + ua[j + m - 2]) {
   --qhat, rhat+= ub[m - 1];
   if (rhat >= BASE) break;
  }
  // ua[j .. j + m] から qhat * ub を引く。
  u64 carry= 0, borrow= 0;
  for (int i= 0; i < m; ++i) {
   const u64 t= qhat * ub[i] + carry;
   carry= t / BASE;
   const u64 lo= t % BASE + borrow;
   if (ua[j + i] < lo) ua[j + i]= u32(u64(ua[j + i]) + BASE - lo), borrow= 1;
   else ua[j + i]-= u32(lo), borrow= 0;
  }
  const u64 top= carry + borrow;
  if (u64(ua[j + m]) < top) {
   // 引きすぎた: qhat を 1 減らして ub を足し戻す (借りで回り込んだ分が戻る)。
   --qhat;
   u64 cc= 0;
   for (int i= 0; i < m; ++i) {
    const u64 t= u64(ua[j + i]) + ub[i] + cc;
    ua[j + i]= u32(t >= BASE ? t - BASE : t), cc= t >= BASE;
   }
   ua[j + m]= u32((u64(ua[j + m]) + cc) % BASE);
  } else {
   ua[j + m]= u32(u64(ua[j + m]) - top);
  }
  q[j]= u32(qhat);
 }
 // 余りは ua を d で割ったもの (割り切れる)。
 r.assign(m, 0);
 u64 rem= 0;
 for (int i= m; i--;) {
  const u64 t= rem * BASE + ua[i];
  r[i]= u32(t / d), rem= t % d;
 }
 assert(rem == 0);
 trim(q), trim(r);
}

// b は正規化済み (m 桁、末尾は BASE/2 以上)。floor(BASE^(m+K) / b) を K + 1 桁で返す。
// 2 段に分ける。
// 1. 作業の桁で倍々に詰める。x ≈ BASE^(m+p)/b を、E = BASE^(m+p) - b x (>= 0) として
//      x <- x BASE^p + floor(E x / BASE^m)
//    で p = 1, 2, 4, ... と倍にする (掛けるのは m x (p + 1) と (m + p) x (p + 1) で軽い)。
//    この形の誤差 e は e <- (b/BASE^m) e^2 + 1 で、b/BASE^m が 1/2 以上だと不動点が無く、段を重ねると
//    少しずつ膨らむ。p が半分まで来たら 2 に移る (絶対の誤差は p 桁の尺度で 1 桁分に収まっている)。
// 2. 答えの桁で 2 回詰める。I = x BASE^(K-p) から、E = BASE^(m+K) - b I、I += floor(E I / BASE^(m+K))。
//    切り捨ての 1 が答えの 1 桁分にしかならないので、誤差 e は e <- e^2 / (BASE^(m+K)/b) + 1 になり、
//    1 段で 1、2 に収まる (1 の相対誤差は BASE^(-p) 以下なので、最初の 1 段でほぼ 0 になる)。
//    こちらは全体の積を 2 回払うが、2 段で済む。
inline std::vector<u32> reciprocal(const u32* b, int m, int K) {
 std::vector<u32> dividend(m + 2, 0), x, rr;
 dividend[m + 1]= 1;  // BASE^(m+1)
 knuth_divmod(dividend.data(), m + 2, b, m, x, rr);  // x = floor(BASE^(m+1) / b) (2 桁)
 const std::vector<u32> bb(b, b + m);
 int p= 1;
 while (2 * p < K + 1) {
  std::vector<u32> E(m + p + 1, 0), prod= mul(bb, x);
  E[m + p]= 1;
  sub_abs_inplace(E, prod);  // E = BASE^(m+p) - b x >= 0
  std::vector<u32> T= mul(E, x);
  T.erase(T.begin(), T.begin() + std::min<int>(m, int(T.size())));
  x.insert(x.begin(), p, 0);
  x= add_abs(x, T);
  p*= 2;
 }
 std::vector<u32> I= x;
 I.insert(I.begin(), K - p, 0);  // I = x BASE^(K-p) (K >= p)
 std::vector<u32> target(m + K + 1, 0);
 target[m + K]= 1;
 for (int it= 0; it < 2; ++it) {
  std::vector<u32> prod= mul(bb, I);
  assert(cmp_abs(target, prod) >= 0);  // I は真の値を超えない
  std::vector<u32> E= target;
  sub_abs_inplace(E, prod);
  std::vector<u32> W= mul(E, I);
  W.erase(W.begin(), W.begin() + std::min<int>(m + K, int(W.size())));
  I= add_abs(I, W);
 }
 trim(I);
 return I;
}

inline std::vector<u32> mul(const std::vector<u32>& a, const std::vector<u32>& b) {
 const int n= int(a.size()), m= int(b.size());
 if (n == 0 || m == 0) return {};
 if (std::min(n, m) <= 64) {
  std::vector<u128> c(n + m);
  for (int i= 0; i < n; ++i)
   for (int j= 0; j < m; ++j) c[i + j]+= u128(a[i]) * b[j];
  std::vector<u32> r(n + m);
  u128 carry= 0;
  for (int i= 0; i < n + m; ++i) {
   const u128 t= c[i] + carry;
   r[i]= u32(t % BASE), carry= t / BASE;
  }
  assert(carry == 0);
  trim(r);
  return r;
 }
 const int sz= n + m - 1;
 std::vector<u32> xa(n), xb(m), rs(std::size_t(3) * sz);
 auto conv= [&]<u32 P>(u32* out) {
  for (int i= 0; i < n; ++i) xa[i]= a[i] >= P ? a[i] - P : a[i];
  for (int i= 0; i < m; ++i) xb[i]= b[i] >= P ? b[i] - P : b[i];
  crt_engine<Ntt>::template conv_at<P>(xa.data(), n, xb.data(), m, sz, out);
 };
 conv.template operator()<Q0>(rs.data());
 conv.template operator()<Q1>(rs.data() + sz);
 conv.template operator()<Q2>(rs.data() + 2 * sz);
 std::vector<u32> r(n + m);
 u128 carry= 0;
 for (int i= 0; i < sz; ++i) {
  const u128 t= garner3(rs[i], rs[sz + i], rs[2 * sz + i]) + carry;
  r[i]= u32(t % BASE), carry= t / BASE;
 }
 for (int i= sz; carry; ++i) r[i]= u32(carry % BASE), carry/= BASE;
 trim(r);
 return r;
}

// a / b の商と余り (どちらも非負)。a >= b >= 1 で、b の末尾は 0 でない。
inline void divmod_pos(const std::vector<u32>& a, const std::vector<u32>& b, std::vector<u32>& q, std::vector<u32>& r) {
 const int n0= int(a.size()), m0= int(b.size());
 if (m0 <= 512 || n0 - m0 + 1 <= 24) {
  knuth_divmod(a.data(), n0, b.data(), m0, q, r);
  return;
 }
 // Newton 法。正規化して逆数の近似を作り、q を求めてから最後に直す。
 const u32 d= u32(BASE / (u64(b[m0 - 1]) + 1));
 const std::vector<u32> A= mul_small(a, d), B= mul_small(b, d);
 const int n= int(A.size()), m= int(B.size()), K= n - m + 1;
 if (m <= 24 || K <= 24) {
  knuth_divmod(A.data(), n, B.data(), m, q, r);
  // 余りは d 倍されている。
  u32 z;
  r= divmod_small(r, d, z);
  assert(z == 0);
  return;
 }
 const std::vector<u32> I= reciprocal(B.data(), m, K);
 std::vector<u32> P= mul(A, I);
 P.erase(P.begin(), P.begin() + std::min<int>(m + K, int(P.size())));  // floor(A I / BASE^(m+K))
 q= std::move(P);
 trim(q);
 std::vector<u32> qb= mul(q, B);
 const std::vector<u32> one(1, 1);
 while (cmp_abs(A, qb) < 0) sub_abs_inplace(q, one), sub_abs_inplace(qb, B);
 std::vector<u32> rem= A;
 sub_abs_inplace(rem, qb);
 while (cmp_abs(rem, B) >= 0) q= add_abs(q, one), sub_abs_inplace(rem, B);
 // rem は d 倍されているので、d で割って戻す。
 u32 z;
 r= divmod_small(rem, d, z);
 assert(z == 0);
}
}  // namespace bigint_internal

// 多倍長整数。符号と 10^9 進の桁 (下位から、末尾に 0 を置かない) を持つ。0 の符号は正。
// 除算は C++ と同じ切り捨てで、余りの符号は割られる数に合わせる。
struct BigInt {
 using u32= bigint_internal::u32;
 bool neg= false;
 std::vector<u32> d;

 BigInt()= default;
 template <class T> requires std::is_integral_v<T> BigInt(T v) {
  unsigned long long x;
  if constexpr (std::is_signed_v<T>) {
   if (v < 0) neg= true, x= 0ull - (unsigned long long)v;
   else x= (unsigned long long)v;
  } else {
   x= (unsigned long long)v;
  }
  for (; x; x/= bigint_internal::BASE) d.push_back(u32(x % bigint_internal::BASE));
 }
 explicit BigInt(const std::string& s) {
  int p= 0, e= int(s.size());
  for (; p < e && (s[p] == '-' || s[p] == '+'); ++p)
   if (s[p] == '-') neg= !neg;
  while (p < e && s[p] == '0') ++p;
  const int L= e - p;
  d.resize((L + bigint_internal::DIG - 1) / bigint_internal::DIG);
  for (int i= 0; i < int(d.size()); ++i) {
   const int hi= e - i * bigint_internal::DIG, lo= std::max(p, hi - bigint_internal::DIG);
   u32 x= 0;
   for (int j= lo; j < hi; ++j) x= x * 10 + u32(s[j] - '0');
   d[i]= x;
  }
  bigint_internal::trim(d);
  if (d.empty()) neg= false;
 }
 std::string to_string() const {
  if (d.empty()) return "0";
  std::string s;
  if (neg) s.push_back('-');
  s+= std::to_string(d.back());
  char buf[bigint_internal::DIG];
  for (int i= int(d.size()) - 1; i-- > 0;) {
   u32 x= d[i];
   for (int j= bigint_internal::DIG; j--; x/= 10) buf[j]= char('0' + x % 10);
   s.append(buf, bigint_internal::DIG);
  }
  return s;
 }
 bool is_zero() const { return d.empty(); }
 void normalize() {
  bigint_internal::trim(d);
  if (d.empty()) neg= false;
 }
 BigInt abs() const {
  BigInt r= *this;
  r.neg= false;
  return r;
 }
 BigInt operator-() const {
  BigInt r= *this;
  if (!r.d.empty()) r.neg= !r.neg;
  return r;
 }
 friend int cmp(const BigInt& a, const BigInt& b) {
  if (a.neg != b.neg) return a.neg ? -1 : 1;
  const int c= bigint_internal::cmp_abs(a.d, b.d);
  return a.neg ? -c : c;
 }
 friend bool operator==(const BigInt& a, const BigInt& b) { return a.neg == b.neg && a.d == b.d; }
 friend bool operator!=(const BigInt& a, const BigInt& b) { return !(a == b); }
 friend bool operator<(const BigInt& a, const BigInt& b) { return cmp(a, b) < 0; }
 friend bool operator>(const BigInt& a, const BigInt& b) { return cmp(a, b) > 0; }
 friend bool operator<=(const BigInt& a, const BigInt& b) { return cmp(a, b) <= 0; }
 friend bool operator>=(const BigInt& a, const BigInt& b) { return cmp(a, b) >= 0; }
 friend BigInt operator+(const BigInt& a, const BigInt& b) {
  using namespace bigint_internal;
  BigInt r;
  if (a.neg == b.neg) {
   r.d= add_abs(a.d, b.d), r.neg= a.neg;
  } else if (cmp_abs(a.d, b.d) >= 0) {
   r.d= a.d, sub_abs_inplace(r.d, b.d), r.neg= a.neg;
  } else {
   r.d= b.d, sub_abs_inplace(r.d, a.d), r.neg= b.neg;
  }
  if (r.d.empty()) r.neg= false;
  return r;
 }
 friend BigInt operator-(const BigInt& a, const BigInt& b) { return a + (-b); }
 friend BigInt operator*(const BigInt& a, const BigInt& b) {
  BigInt r;
  r.d= bigint_internal::mul(a.d, b.d);
  r.neg= !r.d.empty() && (a.neg != b.neg);
  return r;
 }
 // 商と余り (切り捨て)。b は 0 でないこと。
 friend void divmod(const BigInt& a, const BigInt& b, BigInt& q, BigInt& r) {
  using namespace bigint_internal;
  assert(!b.d.empty());
  if (a.d.empty()) {
   q= BigInt(), r= BigInt();
   return;
  }
  if (cmp_abs(a.d, b.d) < 0) {
   q= BigInt(), r= a;
   return;
  }
  divmod_pos(a.d, b.d, q.d, r.d);
  q.neg= !q.d.empty() && (a.neg != b.neg), r.neg= !r.d.empty() && a.neg;
 }
 friend BigInt operator/(const BigInt& a, const BigInt& b) {
  BigInt q, r;
  divmod(a, b, q, r);
  return q;
 }
 friend BigInt operator%(const BigInt& a, const BigInt& b) {
  BigInt q, r;
  divmod(a, b, q, r);
  return r;
 }
 BigInt& operator+=(const BigInt& b) { return *this= *this + b; }
 BigInt& operator-=(const BigInt& b) { return *this= *this - b; }
 BigInt& operator*=(const BigInt& b) { return *this= *this * b; }
 BigInt& operator/=(const BigInt& b) { return *this= *this / b; }
 BigInt& operator%=(const BigInt& b) { return *this= *this % b; }
 friend std::ostream& operator<<(std::ostream& os, const BigInt& v) { return os << v.to_string(); }
};
