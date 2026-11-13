#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <deque>
#include <functional>
#include <iterator>
#include <map>
#include <memory>
#include <ranges>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>
#include "neo/data_structure/HashMap.hpp"
// 有限オートマトン (DFA) と、受理される長さ n の記号の列の上の DP。今の Library の Automaton (mylib/counting/Automaton.hpp)
// の作り直しで、口 (遷移の関数から作る、NFA と ε 遷移、&、DP の op と f(値, 記号, 位置)) は同じにしてある。
// 構築は、初期状態から幅優先にたどって状態に番号を振り、状態の数 × 記号の数の表にする。番号を引く表は、状態の型が
// 64 bit に詰められる (整数、bool、それらの pair や tuple や array) なら NeoLibrary の HashMap、ほかは要素のハッシュを
// 混ぜる unordered_map、ハッシュを作れない型だけ std::map にする。NFA は、NFA の状態に先に番号を振ってから部分集合構成を
// する (状態が 64 個以下なら部分集合を u64 の bit で持つ)。
// 積 (&) は表を作らず、成分の表を共有して並べるだけにする。DP は 1 歩ずつ進め、成分ごとにその歩で来られる状態の集合を
// その成分だけで求めて 0 から番号を振り、成分ごとの番号を 2 の冪に切り上げた幅の bit に並べた値を、密な配列 (来た組の
// 置き場所を引く表) の添字にする。ハッシュも待ち行列も使わない。成分の集合が前の歩と同じになったら、以後は求め直さない。
// 集合が止まった成分を bit の下位に並べ、その表の和を前もって 1 つの表にまとめ、変わる成分どうしも毎歩まとめるので、
// 移る先を 1 つ求めるのに引く表はたいてい 2 つになる。
// 「N 以下」の DFA (位置を状態に持つ) との積では、その成分の集合が毎歩 2 つになるので、手で書く桁 DP の
// dp[i][tight][...] と同じ大きさの配列で進む。添字の bit が 24 を超えたら (成分どうしが連動していて組の大半には来られない
// とき)、bit の多い 2 つの成分を、来られる組だけの表にした積に替えて最初からやり直す。
namespace dfa_internal {
using u64= unsigned long long;
template <class T> struct is_pair: std::false_type {};
template <class A, class B> struct is_pair<std::pair<A, B>>: std::true_type {};
template <class T> struct is_tuple: std::false_type {};
template <class... A> struct is_tuple<std::tuple<A...>>: std::true_type {};
template <class T> struct is_array: std::false_type {};
template <class A, std::size_t N> struct is_array<std::array<A, N>>: std::true_type {};
template <class T> constexpr int pack_bits();
template <class T, std::size_t... I> constexpr int tuple_bits(std::index_sequence<I...>) {
 int b[]= {pack_bits<std::tuple_element_t<I, T>>()..., 0}, s= 0;
 for(int x: b) {
  if(x < 0) return -1;
  s+= x;
 }
 return s;
}
// T を u64 に単射で詰めるときの bit 数。詰められない型は -1 (64 を超えるものも呼ぶ側で除く)。
template <class T> constexpr int pack_bits() {
 if constexpr(std::is_same_v<T, bool>) return 1;
 else if constexpr(std::is_enum_v<T>) return pack_bits<std::underlying_type_t<T>>();
 else if constexpr(std::is_integral_v<T>) return sizeof(T) <= 8 ? int(sizeof(T)) * 8 : -1;
 else if constexpr(is_pair<T>::value) {
  constexpr int a= pack_bits<typename T::first_type>(), b= pack_bits<typename T::second_type>();
  return a < 0 || b < 0 ? -1 : a + b;
 } else if constexpr(is_array<T>::value) {
  constexpr int a= pack_bits<typename T::value_type>();
  return a < 0 ? -1 : a * int(std::tuple_size_v<T>);
 } else if constexpr(is_tuple<T>::value) return tuple_bits<T>(std::make_index_sequence<std::tuple_size_v<T>>());
 else return -1;
}
template <class T> constexpr bool packable= pack_bits<T>() >= 0 && pack_bits<T>() <= 64;
inline u64 shl(u64 x, int s) { return s >= 64 ? 0 : x << s; }
template <class T> u64 pack(const T& x) {
 if constexpr(std::is_same_v<T, bool>) return x;
 else if constexpr(std::is_enum_v<T>) return pack(static_cast<std::underlying_type_t<T>>(x));
 else if constexpr(std::is_integral_v<T>) return static_cast<std::make_unsigned_t<T>>(x);
 else if constexpr(is_pair<T>::value) return pack(x.first) | shl(pack(x.second), pack_bits<typename T::first_type>());
 else if constexpr(is_array<T>::value) {
  constexpr int b= pack_bits<typename T::value_type>();
  u64 r= 0;
  for(std::size_t i= 0; i < x.size(); ++i) r|= shl(pack(x[i]), b * int(i));
  return r;
 } else {
  u64 r= 0;
  int s= 0;
  std::apply([&](const auto&... a) { ((r|= shl(pack(a), s), s+= pack_bits<std::decay_t<decltype(a)>>()), ...); }, x);
  return r;
 }
}
inline u64 mix(u64 x) {
 x+= 0x9e3779b97f4a7c15ull;
 x= (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
 x= (x ^ (x >> 27)) * 0x94d049bb133111ebull;
 return x ^ (x >> 31);
}
template <class T> constexpr bool hashable() {
 if constexpr(packable<T> || std::is_integral_v<T> || std::is_floating_point_v<T>) return true;
 else if constexpr(is_pair<T>::value) return hashable<typename T::first_type>() && hashable<typename T::second_type>();
 else if constexpr(is_tuple<T>::value) return []<std::size_t... I>(std::index_sequence<I...>) { return (hashable<std::tuple_element_t<I, T>>() && ...); }(std::make_index_sequence<std::tuple_size_v<T>>());
 else if constexpr(is_array<T>::value) return hashable<typename T::value_type>();
 else if constexpr(requires(const T& x) { std::hash<T>{}(x); }) return true;
 else if constexpr(requires(const T& x) { std::begin(x) != std::end(x); }) return hashable<std::decay_t<decltype(*std::begin(std::declval<const T&>()))>>();
 else return false;
}
template <class T> u64 hash_of(const T& x) {
 if constexpr(packable<T>) return mix(pack(x));
 else if constexpr(std::is_integral_v<T>) {
  const std::make_unsigned_t<T> v= x;
  u64 h= 0;
  for(std::size_t i= 0; i < sizeof(T); i+= 8) h= mix(h ^ u64(v >> (8 * i)));
  return h;
 } else if constexpr(is_pair<T>::value) return mix(hash_of(x.first) ^ mix(hash_of(x.second)));
 else if constexpr(is_tuple<T>::value || is_array<T>::value) {
  u64 h= 0;
  std::apply([&](const auto&... a) { ((h= mix(h ^ hash_of(a))), ...); }, x);
  return h;
 } else if constexpr(requires { std::hash<T>{}(x); }) return mix(std::hash<T>{}(x));
 else {
  // 列。要素が詰められる型なら、要素ごとに掛け算 1 回で混ぜ、最後に splitmix64 で混ぜる (長い vector の状態で
  // 要素ごとに splitmix64 を連ねると遅い)。中身が連続した列は 4 本に分けて、掛け算の待ちを重ねる。
  using E= std::decay_t<decltype(*std::begin(x))>;
  constexpr u64 K= 0x9e3779b97f4a7c15ull;
  if constexpr(packable<E> && std::ranges::contiguous_range<const T&> && std::ranges::sized_range<const T&>) {
   const E* p= std::ranges::data(x);
   const std::size_t n= std::ranges::size(x);
   u64 a= 1, b= 2, c= 3, d= 4;
   std::size_t i= 0;
   for(; i + 4 <= n; i+= 4) a= (a ^ pack(p[i])) * K, b= (b ^ pack(p[i + 1])) * K, c= (c ^ pack(p[i + 2])) * K, d= (d ^ pack(p[i + 3])) * K;
   for(; i < n; ++i) a= (a ^ pack(p[i])) * K;
   return mix(a ^ (b << 16 | b >> 48) ^ (c << 32 | c >> 32) ^ (d << 48 | d >> 16) ^ n);
  } else if constexpr(packable<E>) {
   u64 h= 1, n= 0;
   for(const auto& e: x) h= (h ^ pack(e)) * K, ++n;
   return mix(h ^ n);
  } else {
   u64 h= 0;
   for(const auto& e: x) h= mix(h ^ hash_of(e));
   return h;
  }
 }
}
struct Hasher {
 template <class T> std::size_t operator()(const T& x) const { return hash_of(x); }
};
// 状態に 0 から番号を振る表。insert は (番号, 新しく入れたか) を返す。
template <class K> struct NumPack {
 HashMap<int> h;
 int n= 0;
 std::pair<int, bool> insert(const K& k) {
  int& v= h[pack(k)];
  if(v) return {v - 1, false};
  return v= ++n, std::pair{n - 1, true};
 }
};
template <class K> struct NumHash {
 std::unordered_map<K, int, Hasher> h;
 std::pair<int, bool> insert(const K& k) {
  auto [it, ins]= h.try_emplace(k, int(h.size()));
  return {it->second, ins};
 }
};
template <class K> struct NumTree {
 std::map<K, int> h;
 std::pair<int, bool> insert(const K& k) {
  auto [it, ins]= h.try_emplace(k, int(h.size()));
  return {it->second, ins};
 }
};
template <class K> using Numbering= std::conditional_t<packable<K>, NumPack<K>, std::conditional_t<hashable<K>(), NumHash<K>, NumTree<K>>>;
// 表にした DFA。状態 0 が初期状態で、tab[s * m + c] は状態 s から c 番目の記号で移る先 (捨てるなら -1)。
struct Comp {
 int n= 0;
 std::vector<int> tab;
 std::vector<char> acc;
};
using CompPtr= std::shared_ptr<const Comp>;
template <class R, class S> concept state_range= !std::is_convertible_v<R, S> && requires(const R& r) {
 std::begin(r) != std::end(r);
 requires std::is_convertible_v<decltype(*std::begin(r)), const S&>;
};
template <class F, class S, class A> concept nfa_transition= std::is_invocable_v<const F&, const S&, const A&> && state_range<std::invoke_result_t<const F&, const S&, const A&>, S>;
template <class F, class S, class A> concept dfa_transition= std::is_invocable_r_v<S, const F&, const S&, const A&> && !nfa_transition<F, S, A>;
template <class H, class S> concept eps_transition= std::is_invocable_v<const H&, const S&> && state_range<std::invoke_result_t<const H&, const S&>, S>;
// tr(s, j) は状態 s から j 番目の記号で移る先。rej(t) が true の先は捨てる。
template <class S, class Tr, class Ac, class Rej> CompPtr build(int m, const S& init, const Tr& tr, const Ac& ac, const Rej& rej) {
 Numbering<S> num;
 std::deque<S> dec;  // 足しても前の要素への参照が保たれる
 std::vector<int> tab;
 num.insert(init), dec.push_back(init);
 for(std::size_t i= 0; i < dec.size(); ++i) {
  const S& s= dec[i];
  tab.resize(tab.size() + m);
  int* row= tab.data() + i * m;
  for(int j= 0; j < m; ++j) {
   S t= tr(s, j);
   if(rej(t)) {
    row[j]= -1;
    continue;
   }
   auto [id, ins]= num.insert(t);
   if(ins) dec.push_back(std::move(t));
   row[j]= id;
  }
 }
 auto C= std::make_shared<Comp>();
 C->n= int(dec.size()), C->tab= std::move(tab), C->acc.resize(C->n);
 for(int i= 0; i < C->n; ++i) C->acc[i]= bool(ac(dec[i]));
 return C;
}
inline void sort_unique_tail(std::vector<int>& v, std::size_t b) {
 std::sort(v.begin() + b, v.end());
 v.erase(std::unique(v.begin() + b, v.end()), v.end());
}
// NFA の部分集合構成。tr(s, j) と eps(s) は状態の列 (std::set など) を返す。部分集合は ε 遷移で閉じたものを持つ。
template <bool EPS, class S, class Tr, class Ac, class Eps> CompPtr build_nfa(int m, const S& init, const Tr& tr, const Ac& ac, const Eps& eps) {
 Numbering<S> num;
 std::deque<S> dec;
 auto id= [&](const S& y) {
  auto [v, ins]= num.insert(y);
  if(ins) dec.push_back(y);
  return v;
 };
 // NFA の状態に番号を振り、移り先を (状態, 記号) ごとに並べる
 std::vector<int> off{0}, lst, eoff{0}, elst;
 id(init);
 for(std::size_t x= 0; x < dec.size(); ++x) {
  const S& s= dec[x];
  for(int j= 0; j < m; ++j) {
   const std::size_t b= lst.size();
   for(const auto& y: tr(s, j)) lst.push_back(id(y));
   sort_unique_tail(lst, b), off.push_back(int(lst.size()));
  }
  if constexpr(EPS) {
   const std::size_t b= elst.size();
   for(const auto& y: eps(s)) elst.push_back(id(y));
   sort_unique_tail(elst, b), eoff.push_back(int(elst.size()));
  }
 }
 const int N= int(dec.size());
 std::vector<char> nacc(N);
 for(int x= 0; x < N; ++x) nacc[x]= bool(ac(dec[x]));
 // ε 遷移で閉じた移り先 (toff, tlst) と初期状態の閉包
 std::vector<int> toff, tlst, start;
 if constexpr(EPS) {
  std::vector<int> coff{0}, clst, mk(N, -1);
  for(int x= 0; x < N; ++x) {
   const std::size_t b= clst.size();
   clst.push_back(x), mk[x]= x;
   for(std::size_t q= b; q < clst.size(); ++q)
    for(int e= eoff[clst[q]]; e < eoff[clst[q] + 1]; ++e)
     if(const int z= elst[e]; mk[z] != x) mk[z]= x, clst.push_back(z);
   std::sort(clst.begin() + b, clst.end()), coff.push_back(int(clst.size()));
  }
  toff.push_back(0);
  for(int x= 0, tick= N; x < N; ++x)
   for(int j= 0; j < m; ++j, ++tick) {
    const std::size_t b= tlst.size();
    for(int e= off[x * m + j]; e < off[x * m + j + 1]; ++e)
     for(int q= coff[lst[e]]; q < coff[lst[e] + 1]; ++q)
      if(const int z= clst[q]; mk[z] != tick) mk[z]= tick, tlst.push_back(z);
    std::sort(tlst.begin() + b, tlst.end()), toff.push_back(int(tlst.size()));
   }
  start.assign(clst.begin(), clst.begin() + coff[1]);
 } else toff= std::move(off), tlst= std::move(lst), start= {0};
 auto C= std::make_shared<Comp>();
 std::vector<int> tab;
 if(N <= 64) {
  std::vector<u64> tm(std::size_t(N) * m);
  u64 am= 0, s0= 0;
  for(int x= 0; x < N; ++x) {
   if(nacc[x]) am|= u64(1) << x;
   for(int j= 0; j < m; ++j)
    for(int e= toff[x * m + j]; e < toff[x * m + j + 1]; ++e) tm[x * m + j]|= u64(1) << tlst[e];
  }
  for(int x: start) s0|= u64(1) << x;
  HashMap<int> h;
  std::vector<u64> sub{s0};
  h[s0]= 1;
  for(std::size_t i= 0; i < sub.size(); ++i) {
   const u64 cur= sub[i];
   for(int j= 0; j < m; ++j) {
    u64 t= 0;
    for(u64 b= cur; b; b&= b - 1) t|= tm[__builtin_ctzll(b) * m + j];
    if(!t) {
     tab.push_back(-1);
     continue;
    }
    int& v= h[t];
    if(!v) v= int(sub.size()) + 1, sub.push_back(t);
    tab.push_back(v - 1);
   }
  }
  C->n= int(sub.size()), C->acc.resize(C->n);
  for(int i= 0; i < C->n; ++i) C->acc[i]= (sub[i] & am) != 0;
 } else {
  std::unordered_map<std::vector<int>, int, Hasher> h;
  std::deque<std::vector<int>> sub{start};
  h.emplace(start, 0);
  std::vector<int> buf, mk(N, -1);
  int tick= 0;
  for(std::size_t i= 0; i < sub.size(); ++i) {
   const std::vector<int>& cur= sub[i];
   for(int j= 0; j < m; ++j, ++tick) {
    buf.clear();
    for(int x: cur)
     for(int e= toff[x * m + j]; e < toff[x * m + j + 1]; ++e)
      if(const int z= tlst[e]; mk[z] != tick) mk[z]= tick, buf.push_back(z);
    if(buf.empty()) {
     tab.push_back(-1);
     continue;
    }
    std::sort(buf.begin(), buf.end());
    auto [it, ins]= h.try_emplace(buf, int(sub.size()));
    if(ins) sub.push_back(buf);
    tab.push_back(it->second);
   }
  }
  C->n= int(sub.size()), C->acc.resize(C->n);
  for(int i= 0; i < C->n; ++i) C->acc[i]= std::any_of(sub[i].begin(), sub[i].end(), [&](int x) { return nacc[x]; });
 }
 C->tab= std::move(tab);
 return C;
}
// 2 つの成分の積を、初期状態から来られる組だけの表にする。
inline CompPtr product(const Comp& A, const Comp& B, int m) {
 HashMap<int> h;
 std::vector<std::pair<int, int>> dec{{0, 0}};
 std::vector<int> tab;
 h[0]= 1;
 for(std::size_t i= 0; i < dec.size(); ++i) {
  const auto [a, b]= dec[i];
  tab.resize(tab.size() + m);
  int* row= tab.data() + i * m;
  for(int c= 0; c < m; ++c) {
   const int x= A.tab[std::size_t(a) * m + c], y= B.tab[std::size_t(b) * m + c];
   if(x < 0 || y < 0) {
    row[c]= -1;
    continue;
   }
   int& v= h[u64(x) << 32 | u64(y)];
   if(!v) v= int(dec.size()) + 1, dec.emplace_back(x, y);
   row[c]= v - 1;
  }
 }
 auto C= std::make_shared<Comp>();
 C->n= int(dec.size()), C->tab= std::move(tab), C->acc.resize(C->n);
 for(int i= 0; i < C->n; ++i) C->acc[i]= A.acc[dec[i].first] && B.acc[dec[i].second];
 return C;
}
struct Access;
}
template <class symbol_t> class DFA {
 template <class> friend class DFA;
 friend struct dfa_internal::Access;
 using Comp= dfa_internal::Comp;
 using CompPtr= dfa_internal::CompPtr;
 std::vector<symbol_t> alph;
 std::vector<CompPtr> comps;
 DFA()= default;
 static std::vector<symbol_t> sorted(std::vector<symbol_t> a) { return std::sort(a.begin(), a.end()), a; }
 // 成分が 1 つのとき。状態の番号をそのまま来た組の置き場所を引く添字にする。
 template <bool ALL, class T, class A, class F> std::vector<T> run1(const Comp& C, int n, const A& op, const T& e, const F& f, const T& init) const {
  const int m= int(alph.size());
  std::vector<T> out, val{init}, nval;
  std::vector<int> slot(C.n), key{0}, nkey;
  const int* tab= C.tab.data();
  for(int i= 0;; ++i) {
   if(ALL || i == n) {
    T r= e;
    for(std::size_t s= 0; s < key.size(); ++s)
     if(C.acc[key[s]]) r= op(r, val[s]);
    out.push_back(r);
   }
   if(i == n) break;
   nkey.clear(), nval.clear();
   for(std::size_t s= 0; s < key.size(); ++s) {
    const int* row= tab + std::size_t(key[s]) * m;
    const T& v= val[s];
    for(int c= 0; c < m; ++c) {
     const int t= row[c];
     if(t < 0) continue;
     if(int& p= slot[t]; !p) p= int(nkey.size()) + 1, nkey.push_back(t), nval.push_back(f(v, alph[c], i));
     else nval[p - 1]= op(nval[p - 1], f(v, alph[c], i));
    }
   }
   for(int t: nkey) slot[t]= 0;
   key.swap(nkey), val.swap(nval);
   if(key.empty()) {
    out.resize(ALL ? n + 1 : 1, e);
    break;
   }
  }
  return out;
 }
 // 成分が 2 つ以上のとき。添字の bit が足りなくなったら、混ぜる 2 つの成分を ma と mb に入れて false を返す。
 // 添字の bit は、集合が止まった (frozen) 成分を下位に、まだ変わる成分を上位に並べる。移る先の添字を求めるときは、
 // 下位の成分どうしの表の和を前もって 1 つの表にまとめ (止まっているので作り直しはまれ)、上位の成分どうしも毎歩小さい
 // 表にまとめて、遷移 1 回で引く表を 2 つ以下にする。
 template <bool ALL, class T, class A, class F> bool runk(const std::vector<CompPtr>& cs, int n, const A& op, const T& e, const F& f, const T& init, std::vector<T>& out, int& ma, int& mb) const {
  constexpr int LIMIT= 24, SBITS= 16, DBITS= 12;
  const int m= int(alph.size()), k= int(cs.size());
  // 成分ごとの、今の歩で来られる状態 R (番号の順に並べたもの) と、次の歩への表 ltab
  struct Run {
   const Comp* C;
   std::vector<int> R, Rn, mark, loc, ltab;  // mark[g]: g が最後に入った歩、loc[g]: 最後に振った番号
   std::vector<char> lacc;
   int bits= 0, nbits= 0, sh= 0, nsh= 0, ltab_sh= -1, ver= 0;  // ver: ltab を作り直した回数
   bool frozen= false, was= false, ltab_ok= false;  // was: この歩の初めに止まっていた (今の並びで下位にある)
  };
  std::vector<Run> rs(k);
  for(int j= 0; j < k; ++j) {
   Run& r= rs[j];
   r.C= cs[j].get(), r.R= {0}, r.mark.assign(r.C->n, -1), r.loc.assign(r.C->n, 0), r.mark[0]= 0, r.lacc= {r.C->acc[0]};
  }
  struct Row {
   const int* t;
   int sh, mk;
  };
  std::vector<Row> rows;
  std::vector<int> gs, gd;
  std::vector<std::array<int, 4>> gs_sig, sig;  // gs にまとめた成分の (番号, bit の位置, 幅, ltab の版)
  // 添字の bit の [base, base + b) にある成分 (どれも in(r) が真) の表の和を、その bit の値ごとの表 g にまとめる
  auto combine= [&](std::vector<int>& g, auto in, int base, int b) {
   g.assign((std::size_t(1) << b) * m, 0);
   for(std::size_t u= 0; u < (std::size_t(1) << b); ++u) {
    int* gr= g.data() + u * m;
    bool valid= true;
    for(const Run& r: rs)
     if(in(r) && ((u >> (r.sh - base)) & ((std::size_t(1) << r.bits) - 1)) >= r.R.size()) valid= false;
    if(!valid) {
     std::fill(gr, gr + m, -1);
     continue;
    }
    for(const Run& r: rs)
     if(in(r)) {
      const int* lt= r.ltab.data() + ((u >> (r.sh - base)) & ((std::size_t(1) << r.bits) - 1)) * m;
      for(int c= 0; c < m; ++c) {
       const int a= gr[c], x= lt[c];
       gr[c]= (a | x) < 0 ? -1 : a + x;
      }
     }
   }
  };
  out.clear();
  std::vector<T> val{init}, nval;
  std::vector<int> slot, key{0}, nkey, tg(m);
  for(int i= 0;; ++i) {
   if(ALL || i == n) {
    T res= e;
    for(std::size_t s= 0; s < key.size(); ++s) {
     bool ok= true;
     for(int j= 0; j < k && ok; ++j) ok= rs[j].lacc[(key[s] >> rs[j].sh) & ((1 << rs[j].bits) - 1)];
     if(ok) res= op(res, val[s]);
    }
    out.push_back(res);
   }
   if(i == n) break;
   int ns= 0, bs= 0, bd= 0;
   for(Run& r: rs) r.was= r.frozen, ns+= r.frozen, (r.frozen ? bs : bd)+= r.bits;
   // 成分ごとに、次の歩で来られる状態を求めて番号を振る
   for(Run& r: rs) {
    if(!r.frozen) {
     r.Rn.clear();
     std::size_t old= 0;
     for(int g: r.R) {
      const int* rw= r.C->tab.data() + std::size_t(g) * m;
      for(int c= 0; c < m; ++c) {
       const int t= rw[c];
       if(t < 0 || r.mark[t] == i + 1) continue;
       old+= r.mark[t] == i, r.mark[t]= i + 1, r.loc[t]= int(r.Rn.size()), r.Rn.push_back(t);
      }
     }
     if(r.Rn.size() == r.R.size() && old == r.R.size()) {  // 前の歩と同じ集合なので、以後は変わらない
      r.frozen= true;
      for(std::size_t u= 0; u < r.R.size(); ++u) r.loc[r.R[u]]= int(u);
      r.Rn= r.R;
     }
     r.ltab_ok= false;
     int b= 0;
     while((std::size_t(1) << b) < r.Rn.size()) ++b;
     r.nbits= std::max(r.bits, b);
    } else r.nbits= r.bits;
   }
   // 次の歩の並び: 止まった成分を下位に
   int total= 0;
   for(Run& r: rs)
    if(r.frozen) r.nsh= total, total+= r.nbits;
   for(Run& r: rs)
    if(!r.frozen) r.nsh= total, total+= r.nbits;
   if(total > LIMIT) {
    ma= 0, mb= 1;
    for(int j= 0; j < k; ++j)
     if(rs[j].nbits > rs[ma].nbits) mb= ma, ma= j;
     else if(j != ma && (mb == ma || rs[j].nbits > rs[mb].nbits)) mb= j;
    if(ma > mb) std::swap(ma, mb);
    return false;
   }
   for(Run& r: rs)
    if(!r.ltab_ok || r.ltab_sh != r.nsh) {
     r.ltab.resize(r.R.size() * m);
     for(std::size_t u= 0; u < r.R.size(); ++u) {
      const int* rw= r.C->tab.data() + std::size_t(r.R[u]) * m;
      int* lt= r.ltab.data() + u * m;
      for(int c= 0; c < m; ++c) lt[c]= rw[c] < 0 ? -1 : r.loc[rw[c]] << r.nsh;
     }
     r.ltab_ok= true, r.ltab_sh= r.nsh, ++r.ver;
    }
   // 移る先の添字を求めるのに引く表。下位の成分は、bit の下から SBITS に収まるところまでを 1 つの表にまとめる
   rows.clear();
   int cs_bits= 0;
   sig.clear();
   for(int j= 0; j < k; ++j)
    if(const Run& r= rs[j]; r.was && r.sh + r.bits <= SBITS) cs_bits= std::max(cs_bits, r.sh + r.bits), sig.push_back({j, r.sh, r.bits, r.ver});
   const bool use_gs= sig.size() >= 2;
   if(use_gs) {
    if(sig != gs_sig) combine(gs, [&](const Run& r) { return r.was && r.sh + r.bits <= SBITS; }, 0, cs_bits), gs_sig= sig;
    rows.push_back({gs.data(), 0, (1 << cs_bits) - 1});
   }
   for(const Run& r: rs)
    if(r.was && (!use_gs || r.sh + r.bits > SBITS)) rows.push_back({r.ltab.data(), r.sh, (1 << r.bits) - 1});
   if(k - ns >= 2 && bd <= DBITS && (std::size_t(1) << bd) <= std::max<std::size_t>(16, key.size())) {
    combine(gd, [](const Run& r) { return !r.was; }, bs, bd);
    rows.push_back({gd.data(), bs, (1 << bd) - 1});
   } else
    for(const Run& r: rs)
     if(!r.was) rows.push_back({r.ltab.data(), r.sh, (1 << r.bits) - 1});
   if(slot.size() < (std::size_t(1) << total)) slot.assign(std::size_t(1) << total, 0);
   // 1 歩進める。移る先の添字は表の値の和で、どれかが捨てる (-1) なら捨てる
   nkey.clear(), nval.clear();
   auto upd= [&](int t, const T& v, int c) __attribute__((always_inline)) {
    if(int& p= slot[t]; !p) p= int(nkey.size()) + 1, nkey.push_back(t), nval.push_back(f(v, alph[c], i));
    else nval[p - 1]= op(nval[p - 1], f(v, alph[c], i));
   };
   if(rows.size() == 1) {
    const Row r0= rows[0];
    for(std::size_t s= 0; s < key.size(); ++s) {
     const int* a= r0.t + std::size_t((key[s] >> r0.sh) & r0.mk) * m;
     const T& v= val[s];
     for(int c= 0; c < m; ++c)
      if(const int t= a[c]; t >= 0) upd(t, v, c);
    }
   } else if(rows.size() == 2) {
    const Row r0= rows[0], r1= rows[1];
    for(std::size_t s= 0; s < key.size(); ++s) {
     const int kk= key[s];
     const int *a= r0.t + std::size_t((kk >> r0.sh) & r0.mk) * m, *b= r1.t + std::size_t((kk >> r1.sh) & r1.mk) * m;
     const T& v= val[s];
     for(int c= 0; c < m; ++c)
      if(const int x= a[c], y= b[c]; (x | y) >= 0) upd(x + y, v, c);
    }
   } else {
    for(std::size_t s= 0; s < key.size(); ++s) {
     const int kk= key[s];
     const int* a= rows[0].t + std::size_t((kk >> rows[0].sh) & rows[0].mk) * m;
     for(int c= 0; c < m; ++c) tg[c]= a[c];
     for(std::size_t j= 1; j < rows.size(); ++j) {
      const int* b= rows[j].t + std::size_t((kk >> rows[j].sh) & rows[j].mk) * m;
      for(int c= 0; c < m; ++c) {
       const int x= tg[c], y= b[c];
       tg[c]= (x | y) < 0 ? -1 : x + y;
      }
     }
     const T& v= val[s];
     for(int c= 0; c < m; ++c)
      if(const int t= tg[c]; t >= 0) upd(t, v, c);
    }
   }
   for(int t: nkey) slot[t]= 0;
   key.swap(nkey), val.swap(nval);
   for(Run& r: rs) {
    if(r.Rn.size() && !r.frozen) {
     r.R.swap(r.Rn), r.lacc.resize(r.R.size());
     for(std::size_t u= 0; u < r.R.size(); ++u) r.lacc[u]= r.C->acc[r.R[u]];
    }
    r.bits= r.nbits, r.sh= r.nsh;
   }
   if(key.empty()) {
    out.resize(ALL ? n + 1 : 1, e);
    break;
   }
  }
  return true;
 }
 template <bool ALL, class T, class A, class F> std::vector<T> run(int n, const A& op, const T& e, const F& f, const T& init) const {
  if(n < 0) return std::vector<T>(ALL ? 0 : 1, e);
  std::vector<CompPtr> cs= comps;
  std::vector<T> out;
  for(int a, b;;) {
   if(cs.size() == 1) return run1<ALL>(*cs[0], n, op, e, f, init);
   if(runk<ALL>(cs, n, op, e, f, init, out, a, b)) return out;
   cs[a]= dfa_internal::product(*cs[a], *cs[b], int(alph.size())), cs.erase(cs.begin() + b);
  }
 }
public:
 // DFA。tr(状態, 記号) が次の状態を返す。
 template <class state_t, class F, class G>
  requires dfa_internal::dfa_transition<F, state_t, symbol_t>
 DFA(const std::vector<symbol_t>& alphabet, const state_t& init, const F& tr, const G& ac): alph(sorted(alphabet)) {
  comps.push_back(dfa_internal::build(int(alph.size()), init, [&](const state_t& s, int j) -> state_t { return tr(s, alph[j]); }, ac, [](const state_t&) { return false; }));
 }
 // DFA。reject に移る遷移は捨てる (reject からはどの列も受理されないときに渡すと速い)。
 template <class state_t, class F, class G>
  requires dfa_internal::dfa_transition<F, state_t, symbol_t>
 DFA(const std::vector<symbol_t>& alphabet, const state_t& init, const F& tr, const G& ac, const state_t& reject): alph(sorted(alphabet)) {
  comps.push_back(dfa_internal::build(int(alph.size()), init, [&](const state_t& s, int j) -> state_t { return tr(s, alph[j]); }, ac, [&](const state_t& t) { return t == reject; }));
 }
 // NFA。tr(状態, 記号) が次の状態の集まり (std::set や std::vector) を返す。部分集合構成で DFA に直す。
 template <class state_t, class F, class G>
  requires dfa_internal::nfa_transition<F, state_t, symbol_t>
 DFA(const std::vector<symbol_t>& alphabet, const state_t& init, const F& tr, const G& ac): alph(sorted(alphabet)) {
  comps.push_back(dfa_internal::build_nfa<false>(int(alph.size()), init, [&](const state_t& s, int j) { return tr(s, alph[j]); }, ac, 0));
 }
 // ε 遷移のある NFA。eps(状態) が ε で移れる状態の集まりを返す。
 template <class state_t, class F, class G, class H>
  requires dfa_internal::nfa_transition<F, state_t, symbol_t> && dfa_internal::eps_transition<H, state_t>
 DFA(const std::vector<symbol_t>& alphabet, const state_t& init, const F& tr, const G& ac, const H& eps): alph(sorted(alphabet)) {
  comps.push_back(dfa_internal::build_nfa<true>(int(alph.size()), init, [&](const state_t& s, int j) { return tr(s, alph[j]); }, ac, eps));
 }
 const std::vector<symbol_t>& alphabet() const { return alph; }
 // 両方に受理される列を受理する DFA。表は作らず、成分を並べるだけ。
 DFA operator&(const DFA& r) const {
  assert(alph == r.alph);
  DFA ret;
  ret.alph= alph, ret.comps= comps, ret.comps.insert(ret.comps.end(), r.comps.begin(), r.comps.end());
  return ret;
 }
 // 記号 c を g(c) と読む DFA (記号の集まりは alphabet)。g(c) がもとの記号に無ければ、c の遷移は捨てる。
 template <class S2, class G> DFA<S2> map_symbols(const std::vector<S2>& alphabet, const G& g) const {
  DFA<S2> ret;
  ret.alph= DFA<S2>::sorted(alphabet);
  const int m= int(alph.size()), m2= int(ret.alph.size());
  std::vector<int> to(m2);
  for(int j= 0; j < m2; ++j) {
   const symbol_t a= g(ret.alph[j]);
   auto it= std::lower_bound(alph.begin(), alph.end(), a);
   to[j]= it != alph.end() && *it == a ? int(it - alph.begin()) : -1;
  }
  for(const auto& C: comps) {
   auto D= std::make_shared<Comp>();
   D->n= C->n, D->acc= C->acc, D->tab.resize(std::size_t(C->n) * m2);
   for(int s= 0; s < C->n; ++s)
    for(int j= 0; j < m2; ++j) D->tab[std::size_t(s) * m2 + j]= to[j] < 0 ? -1 : C->tab[std::size_t(s) * m + to[j]];
   ret.comps.push_back(std::move(D));
  }
  return ret;
 }
 // 受理される長さ n の列についての DP。値は f(値, 記号, 位置) で進め、同じ状態に来た値と受理される状態の値を op
 // (可換なモノイドで、単位元は e) でまとめる。長さ 0 の列の値は init。
 template <class T, class A, class F> T dp(int n, const A& op, const T& e, const F& f, const T& init) const { return run<false>(n, op, e, f, init)[0]; }
 // dp の長さ 0 から n までの答え。
 template <class T, class A, class F> std::vector<T> dp_all(int n, const A& op, const T& e, const F& f, const T& init) const { return run<true>(n, op, e, f, init); }
 // 受理される長さ n の列の数。
 template <class T> T count(int n) const { return dp(n, std::plus<>(), T(), [](const T& x, const symbol_t&, int) -> const T& { return x; }, T(1)); }
 // 受理される長さ 0 から n までの列の数。
 template <class T> std::vector<T> count_all(int n) const { return dp_all(n, std::plus<>(), T(), [](const T& x, const symbol_t&, int) -> const T& { return x; }, T(1)); }
};
namespace dfa_internal {
struct Access {
 static DFA<int> make(int base, CompPtr C) {
  DFA<int> r;
  r.alph.resize(base);
  for(int c= 0; c < base; ++c) r.alph[c]= c;
  r.comps= {std::move(C)};
  return r;
 }
};
// 上の桁から読んだ長さ |d| の数字の列を、d を上の桁から並べた数と比べる。leq なら d 以下、でなければ d 以上を受理する。
// 状態 i (< L) は上から i 桁が等しい、L は全部等しい、L + 1 は受理する向きに大小が決まった。
inline DFA<int> cmp_msb(const std::vector<int>& d, int base, bool leq) {
 const int L= int(d.size());
 auto C= std::make_shared<Comp>();
 C->n= L + 2, C->tab.assign(std::size_t(L + 2) * base, -1), C->acc.assign(L + 2, 1);
 for(int i= 0; i < L; ++i)
  for(int c= 0; c < base; ++c) {
   int& t= C->tab[std::size_t(i) * base + c];
   if(c == d[i]) t= i + 1;
   else if((c < d[i]) == leq) t= L + 1;
  }
 for(int c= 0; c < base; ++c) C->tab[std::size_t(L) * base + c]= L, C->tab[std::size_t(L + 1) * base + c]= L + 1;
 return Access::make(base, std::move(C));
}
inline DFA<int> cmp_msb(unsigned long long N, int base, int len, bool leq) {
 std::vector<int> d(len);
 for(int i= len; i--;) d[i]= int(N % base), N/= base;
 if(N) {  // N が base^len 以上なので、leq ならすべて、でなければどれも受理しない
  if(leq) return cmp_msb(std::vector<int>(len, base - 1), base, true);
  auto C= std::make_shared<Comp>();
  C->n= 1, C->tab.assign(base, -1), C->acc= {0};
  return Access::make(base, std::move(C));
 }
 return cmp_msb(d, base, leq);
}
// 下の桁から読んだ列 (長さは問わない) の値を N と比べる。N は base 進で L 桁。状態 2i + b は、下から i 桁 (L 桁で止める)
// を読み、その値が N の下 i 桁と比べて受理する向き (leq なら以下) なら b = 0、逆なら b = 1。i < L なら値は N より小さい
// ので、leq はすべて受理し、geq はどれも受理しない。
inline DFA<int> cmp_lsb(unsigned long long N, int base, bool leq) {
 std::vector<int> d;
 for(; N; N/= base) d.push_back(int(N % base));
 const int L= int(d.size());
 auto C= std::make_shared<Comp>();
 C->n= 2 * L + 2, C->tab.assign(std::size_t(2 * L + 2) * base, -1), C->acc.assign(2 * L + 2, 0);
 for(int i= 0; i <= L; ++i)
  for(int b= 0; b < 2; ++b) {
   int* t= C->tab.data() + std::size_t(2 * i + b) * base;
   C->acc[2 * i + b]= i < L ? leq : b == 0;
   for(int c= 0; c < base; ++c) {
    if(i < L) t[c]= 2 * (i + 1) + (c == d[i] ? b : (c < d[i]) != leq);
    else if(c == 0) t[c]= 2 * L + b;
    else t[c]= leq ? -1 : 2 * L;
   }
  }
 if(leq)
  for(int& t: C->tab)
   if(t == 2 * L + 1) t= -1;  // N を超えた
 return Access::make(base, std::move(C));
}
}
// 上の桁から読んだ |N| 桁の数字の列 (記号は int の 0 から 9) のうち、10 進の文字列 N の表す数以下のものを受理する。
inline DFA<int> leq_digits(const std::string& N) {
 std::vector<int> d(N.size());
 for(std::size_t i= 0; i < N.size(); ++i) assert('0' <= N[i] && N[i] <= '9'), d[i]= N[i] - '0';
 return dfa_internal::cmp_msb(d, 10, true);
}
// 上の桁から読んだ |digits| 桁の列 (記号は 0 以上 base 未満) のうち、digits を上の桁から並べた数以下のものを受理する。
inline DFA<int> leq_digits(const std::vector<int>& digits, int base) { return dfa_internal::cmp_msb(digits, base, true); }
// 上の桁から読んだ len 桁の base 進の列のうち、N 以下のものを受理する。
inline DFA<int> leq_digits(unsigned long long N, int base, int len) { return dfa_internal::cmp_msb(N, base, len, true); }
// 下の桁から読んだ base 進の列 (長さは問わない) のうち、値が N 以下のものを受理する。
inline DFA<int> leq_digits_lsb(unsigned long long N, int base) { return dfa_internal::cmp_lsb(N, base, true); }
// 以下を以上にしたもの。
inline DFA<int> geq_digits(const std::string& N) {
 std::vector<int> d(N.size());
 for(std::size_t i= 0; i < N.size(); ++i) assert('0' <= N[i] && N[i] <= '9'), d[i]= N[i] - '0';
 return dfa_internal::cmp_msb(d, 10, false);
}
inline DFA<int> geq_digits(const std::vector<int>& digits, int base) { return dfa_internal::cmp_msb(digits, base, false); }
inline DFA<int> geq_digits(unsigned long long N, int base, int len) { return dfa_internal::cmp_msb(N, base, len, false); }
inline DFA<int> geq_digits_lsb(unsigned long long N, int base) { return dfa_internal::cmp_lsb(N, base, false); }
