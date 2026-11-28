#pragma once
#include <cmath>
#include <cstdint>
#include <type_traits>
#include <utility>
#include <vector>
// モノイド作用の離散対数。X をモノイド、S を集合、X が S に左作用するとする。
//   act: X × S -> S (作用)、op: X × X -> X (合成)、hash: S -> 整数 (等しい値は等しいハッシュ)
// x ∈ X、s, t ∈ S、N に対し、x^n s = t かつ 0 ≤ n < N を満たす最小の n を返す。無ければ -1。
// 逆元は仮定しない。maspy の記事 (モノイド作用に関する離散対数問題) の方法。
//
// T = {x^i t | 1 ≤ i ≤ m} をハッシュ表に入れ、s に x^m を作用させて x^{km} s ∈ T となる最小の k を探す。
// 解 n があれば k = ceil((n+1)/m) で x^{km} s = x^{km-n} t ∈ T になるので、そのブロック [(k-1)m, km) を
// 全探索すれば最小の解が得られる。ハッシュの衝突で空振りすることがあるので、全探索が 2 回空振りしたら
// 解なしとする (x^{km} s が t の軌道に入ったあとは、サイクルの上を回って必ずもう一度当たる)。
//
// m は √N 以下の最大の 2 の冪。2 の冪にすると、単位元を要求せずに x^m を 2 乗の繰り返しで作れる。
// 作用は O(√N) 回、合成は O(log N) 回、ハッシュは O(√N) 回。メモリは O(√N)。
// 値の列 x^i t は前の値から作る (余分なコピーをしない)。s = t と N が小さいときは表を作らない。
//
//   auto act= [](Aff f, Mint x) { return f.a * x + f.b; };
//   auto op= [](Aff f, Aff g) { return Aff{f.a * g.a, f.a * g.b + f.b}; };
//   auto hash= [](Mint x) { return x.val(); };
//   DiscreteLogAction dlog(act, op, hash, 1000000000);
//   long long n= dlog(x, s, t);        // 0 <= n < 1000000000
//   long long n= dlog(x, s, t, N);     // 0 <= n < N
namespace discrete_log_action_internal {
using u32= unsigned;
using u64= unsigned long long;
using i64= long long;
// 関数オブジェクト (ラムダ) と関数ポインタから、戻り値の型と第 1 引数の型を取る。
// 呼び出し側の型 (Mint のところへ int を渡す、{A, B} と書く) を関数の型に合わせて変換させるために要る。
template <class F> struct fn_type {
 using type= decltype(&F::operator());
};
template <class R, class... A> struct fn_type<R (*)(A...)> {
 using type= R (*)(A...);
};
template <class F> using fn_type_t= typename fn_type<F>::type;
template <class T> struct fn_result;
template <class R, class... A> struct fn_result<R (*)(A...)> {
 using type= R;
};
template <class C, class R, class... A> struct fn_result<R (C::*)(A...)> {
 using type= R;
};
template <class C, class R, class... A> struct fn_result<R (C::*)(A...) const> {
 using type= R;
};
template <class C, class R, class... A> struct fn_result<R (C::*)(A...) noexcept> {
 using type= R;
};
template <class C, class R, class... A> struct fn_result<R (C::*)(A...) const noexcept> {
 using type= R;
};
template <class T> struct fn_first_arg;
template <class R, class A0, class... A> struct fn_first_arg<R (*)(A0, A...)> {
 using type= std::decay_t<A0>;
};
template <class C, class R, class A0, class... A> struct fn_first_arg<R (C::*)(A0, A...)> {
 using type= std::decay_t<A0>;
};
template <class C, class R, class A0, class... A> struct fn_first_arg<R (C::*)(A0, A...) const> {
 using type= std::decay_t<A0>;
};
template <class C, class R, class A0, class... A> struct fn_first_arg<R (C::*)(A0, A...) noexcept> {
 using type= std::decay_t<A0>;
};
template <class C, class R, class A0, class... A> struct fn_first_arg<R (C::*)(A0, A...) const noexcept> {
 using type= std::decay_t<A0>;
};
}  // namespace discrete_log_action_internal

template <class Act, class Op, class Hash> struct DiscreteLogAction {
 using X= typename discrete_log_action_internal::fn_first_arg<discrete_log_action_internal::fn_type_t<Act>>::type;
 using S= std::decay_t<typename discrete_log_action_internal::fn_result<discrete_log_action_internal::fn_type_t<Act>>::type>;
 Act act;
 Op op;
 Hash hash;
 discrete_log_action_internal::i64 lim;
 explicit DiscreteLogAction(Act act_, Op op_, Hash hash_, discrete_log_action_internal::i64 lim_= 1ll << 50): act(act_), op(op_), hash(hash_), lim(lim_) {}
 discrete_log_action_internal::i64 operator()(const X& x, S s, const S& t, discrete_log_action_internal::i64 N= -1) const {
  using discrete_log_action_internal::u32;
  using discrete_log_action_internal::u64;
  using discrete_log_action_internal::i64;
  if(N < 0) N= lim;
  if(N <= 0) return -1;
  if(s == t) return 0;
  if(N <= 64) {
   for(i64 n= 1; n < N; ++n) {
    s= act(x, s);
    if(s == t) return n;
   }
   return -1;
  }
  const u64 r= (u64)std::sqrt((double)N) + 1;
  const int m= (int)(1ull << (63 - __builtin_clzll(r))), mask= m - 1;
  std::vector<S> val(m), vs(m);
  std::vector<int> os(m + 1), so(m);
  val[0]= act(x, t);
  for(int i= 1; i < m; ++i) val[i]= act(x, val[i - 1]);
  for(int i= 0; i < m; ++i) ++os[(so[i]= (u32)hash(val[i]) & mask)];
  for(int i= 0; i < m; ++i) os[i + 1]+= os[i];
  for(int i= 0; i < m; ++i) vs[--os[so[i]]]= std::move(val[i]);
  X y= x;
  for(int k= m; k>>= 1;) y= op(y, y);
  bool failed= false;
  S s1;
  for(i64 n= 0;; s= s1) {
   const int a= (u32)hash(s1= act(y, s)) & mask;
   bool hit= false;
   for(int j= os[a]; j < os[a + 1]; ++j)
    if(s1 == vs[j]) {
     hit= true;
     break;
    }
   if(hit) {
    const i64 len= m < N - n ? (i64)m : N - n;
    S w= s;
    for(i64 i= 0; i < len; ++i) {
     if(w == t) return n + i;
     w= act(x, w);
    }
    if(failed) return -1;
    failed= true;
   }
   if((n+= m) >= N) break;
  }
  return -1;
 }
};
