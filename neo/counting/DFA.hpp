#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <deque>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <ranges>
#include <string>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <unordered_map>
#include <utility>
#include <vector>
#include "neo/data_structure/HashMap.hpp"
// 有限オートマトン (DFA) と、受理される長さ n の記号の列の上の DP。今の Library の Automaton (mylib/counting/Automaton.hpp)
// の作り直しで、口 (遷移の関数から作る、NFA と ε 遷移、&、DP の op と f(値, 記号, 位置)) は同じにしてある。
// 構築は、初期状態から幅優先にたどって状態に番号を振り、状態の数 × 記号の数の表にする。番号を引く表は、状態の型が
// 64 bit に詰められる (整数、bool、それらの pair や tuple や array) なら、整数の葉ごとの値の範囲を広げる密な配列
// (まばらなら NeoLibrary の HashMap)、ほかは要素のハッシュを混ぜる unordered_map、ハッシュを作れない型だけ std::map に
// する。NFA は、NFA の状態に先に番号を振ってから部分集合構成をする (状態が 64 個以下なら部分集合を u64 の bit で持つ)。
// 積 (&) は表を作らず、成分の表を共有して並べるだけにする。DP は 1 歩ずつ進め、成分ごとにその歩で来られる状態の集合を
// その成分だけで求めて 0 から番号を振り、成分ごとの番号を 2 の冪に切り上げた幅の bit に並べた値を、密な配列 (来た組の
// 置き場所を引く表) の添字にする。ハッシュも待ち行列も使わない。成分の集合が前の歩と同じになったら、以後は求め直さない。
// 集合が止まった成分を bit の下位に並べ、その表の和を前もって 1 つの表にまとめ、変わる成分どうしも毎歩まとめるので、
// 移る先を 1 つ求めるのに引く表はたいてい 2 つになる。
// 「N 以下」の DFA (位置を状態に持つ) との積では、その成分の集合が毎歩 2 つになるので、手で書く桁 DP の
// dp[i][tight][...] と同じ大きさの配列で進む。添字の bit が 24 を超えたら (成分どうしが連動していて組の大半には来られない
// とき)、bit の多い 2 つの成分を、来られる組だけの表にした積に替えて最初からやり直す。
// ただし、長さが 32 以上で、積の来られる組が少なく DP で同じ組を何歩も使い回すとき (「N 以下」の DFA と小さい成分の積で、
// 小さいと決まった側の組が毎歩同じもの) は、はじめに積を来られる組だけの表にして、成分が 1 つのときと同じ道で進める。
// 表にするかは、積を幅優先で作りながら組の増え方から見積もる (product_table)。
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
// 詰められる型の状態に番号を振る表。詰めた値を整数の葉 (bool、enum、整数) ごとの値に分け、葉ごとの値の範囲
// [lo, lo + 2^lb) を広げながら、葉の値を並べた添字の密な配列で番号を引く。位置や余りを状態に持つ DFA は値が詰まって
// いて、HashMap より速い (ハッシュで散らさないので、順に増える状態がメモリを順に触る)。範囲を広げるときは配列を作り
// 直し、配列が番号の数の 8 倍 (と 2^10) を超えるほどまばらなら HashMap で引く。HashMap のあいだも範囲は追い、番号の
// 数が 2 倍になるたび、密な配列へ戻せるか見直す。
struct Leaf {
 int off, w;
 bool sgn;
};
template <class T> constexpr int leaf_count() {
 if constexpr(std::is_same_v<T, bool> || std::is_enum_v<T> || std::is_integral_v<T>) return 1;
 else if constexpr(is_pair<T>::value) return leaf_count<typename T::first_type>() + leaf_count<typename T::second_type>();
 else if constexpr(is_array<T>::value) return leaf_count<typename T::value_type>() * int(std::tuple_size_v<T>);
 else return []<std::size_t... I>(std::index_sequence<I...>) { return (0 + ... + leaf_count<std::tuple_element_t<I, T>>()); }(std::make_index_sequence<std::tuple_size_v<T>>());
}
// pack と同じ並びで、葉の bit の位置、幅、符号の有無を out に書く
template <class T> constexpr void leaves(Leaf* out, int& i, int off) {
 if constexpr(std::is_enum_v<T>) leaves<std::underlying_type_t<T>>(out, i, off);
 else if constexpr(std::is_same_v<T, bool>) out[i++]= {off, 1, false};
 else if constexpr(std::is_integral_v<T>) out[i++]= {off, int(sizeof(T)) * 8, std::is_signed_v<T>};
 else if constexpr(is_pair<T>::value) leaves<typename T::first_type>(out, i, off), leaves<typename T::second_type>(out, i, off + pack_bits<typename T::first_type>());
 else if constexpr(is_array<T>::value) {
  for(std::size_t j= 0; j < std::tuple_size_v<T>; ++j) leaves<typename T::value_type>(out, i, off + pack_bits<typename T::value_type>() * int(j));
 } else
  [&]<std::size_t... I>(std::index_sequence<I...>) {
   int o= off;
   ((leaves<std::tuple_element_t<I, T>>(out, i, o), o+= pack_bits<std::tuple_element_t<I, T>>()), ...);
  }(std::make_index_sequence<std::tuple_size_v<T>>());
}
template <class T> constexpr auto leaf_list() {
 std::array<Leaf, leaf_count<T>()> a{};
 int i= 0;
 leaves<T>(a.data(), i, 0);
 return a;
}
template <class K> struct NumPack {
 static constexpr auto LF= leaf_list<K>();
 static constexpr int F= int(LF.size()), F1= F > 0 ? F : 1;
 using Vals= std::array<u64, F1>;
 u64 lo[F1]{};
 int lb[F1]{}, sh[F1]{}, tot= 0;
 std::vector<int> id;  // 密なとき、葉の値を並べた添字から番号 + 1
 std::vector<u64> keys;  // 番号の順の詰めた値 (表を作り直すときに使う)
 HashMap<int> h;  // まばらなときの表
 int n= 0, next_check= 1024;
 char mode= 0;  // 0: 空、1: 密な配列、2: HashMap
 bool never= false;  // 範囲が広すぎて密にできない
 // 葉 f の値を、大小の順を保つ符号なしの値にしたもの (位置と幅は定数)
 template <int f> static u64 val(u64 p) {
  constexpr Leaf L= LF[f];
  u64 x= p >> L.off;
  if constexpr(L.w < 64) x&= (u64(1) << L.w) - 1;
  if constexpr(L.sgn) x^= u64(1) << (L.w - 1);
  return x;
 }
 static Vals vals(u64 p) {
  Vals v{};
  [&]<std::size_t... I>(std::index_sequence<I...>) { ((v[I]= val<I>(p)), ...); }(std::make_index_sequence<F>());
  return v;
 }
 // p が今の範囲にあれば、その添字を r に入れて true
 [[gnu::always_inline]] bool find(u64 p, std::size_t& r) const {
  u64 bad= 0;
  r= 0;
  [&]<std::size_t... I>(std::index_sequence<I...>) { ((bad|= (val<I>(p) - lo[I]) >> lb[I], r|= std::size_t(val<I>(p) - lo[I]) << sh[I]), ...); }(std::make_index_sequence<F>());
  return bad == 0;
 }
 std::size_t index(u64 p) const {
  std::size_t r;
  find(p, r);
  return r;
 }
 // 範囲を p が入るように広げる (幅は 2 倍以上ずつ)。2^40 を超える幅が要るなら never にする。
 void widen(u64 p) {
  const Vals v= vals(p);
  tot= 0;
  for(int f= 0; f < F; ++f) {
   const u64 size= u64(1) << lb[f];
   if(v[f] - lo[f] >= size) {
    const u64 span= v[f] < lo[f] ? lo[f] + size - v[f] : v[f] - lo[f] + 1;
    int b= lb[f] + 1;
    while(b < 40 && (u64(1) << b) < span) ++b;
    if((u64(1) << b) < span) {
     never= true;
     return;
    }
    if(v[f] < lo[f]) lo[f]= lo[f] + size - (u64(1) << b);  // 下へ広げるときは上の端を保つ
    lb[f]= b;
   }
   sh[f]= tot, tot+= lb[f];
  }
 }
 bool fits() const { return !never && tot <= 34 && (std::size_t(1) << tot) <= std::max<std::size_t>(1024, std::size_t(8) * (n + 1)); }
 void to_dense() {
  mode= 1, id.assign(std::size_t(1) << tot, 0), h.clear();
  for(int i= 0; i < n; ++i) id[index(keys[i])]= i + 1;
 }
 void to_hash() {
  mode= 2, std::vector<int>().swap(id);
  for(int i= 0; i < n; ++i) h[keys[i]]= i + 1;
 }
 // 新しい状態を入れる、範囲を広げる、表を移すなどの遅い道
 [[gnu::noinline]] std::pair<int, bool> insert_slow(u64 p) {
  if(mode == 0) {
   const Vals v= vals(p);
   for(int f= 0; f < F; ++f) lo[f]= v[f];
   mode= 1, id.assign(1, 0);
  }
  std::size_t r;
  if(mode == 1) {
   if(!find(p, r)) {
    widen(p);
    if(fits()) to_dense(), r= index(p);
    else to_hash();
   }
   if(mode == 1) {
    int& v= id[r];
    if(v) return {v - 1, false};
    return v= ++n, keys.push_back(p), std::pair{n - 1, true};
   }
  }
  // まばらなとき。範囲は追い続け、番号の数が 2 倍になるごとに密にできるか見直す (はじめに離れた状態が出て、あとで
  // その間が埋まる DFA のため)
  int& v= h[p];
  if(v) return {v - 1, false};
  v= ++n, keys.push_back(p);
  if(!never && !find(p, r)) widen(p);
  if(n >= next_check) {
   next_check*= 2;
   if(fits()) to_dense();
  }
  return {n - 1, true};
 }
 // 速い道 (今の範囲の中か、HashMap にある状態) だけを展開させる。gcc は展開しないと、呼び出しと葉の位置の読み直しで
 // 番号付けが遅くなる。
 [[gnu::always_inline]] std::pair<int, bool> insert(const K& k) {
  const u64 p= pack(k);
  if(mode == 1) {
   std::size_t r;
   if(find(p, r)) [[likely]] {
    int& v= id[r];
    if(v) return {v - 1, false};
    return v= ++n, keys.push_back(p), std::pair{n - 1, true};
   }
  } else if(mode == 2)
   if(const int v= h.get(p)) return {v - 1, false};
  return insert_slow(p);
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
// st は利用者の状態の値 (std::deque<状態の型>、番号の順)。st_type はその型。構築で作った成分だけが持つ。
struct Comp {
 int n= 0;
 std::vector<int> tab;
 std::vector<char> acc;
 std::shared_ptr<const void> st;
 const std::type_info* st_type= nullptr;
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
 C->st= std::make_shared<const std::deque<S>>(std::move(dec)), C->st_type= &typeid(S);
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
// k 個の成分の積を、初期状態から n - 1 歩までに来られる組だけの表にする (n 歩目に初めて来る組の行は -1)。表にすると
// 得かを見積もり、得でなければ途中でやめて nullptr を返す。見積もりは、幅優先で段ごとに進め、L_MIN 段目とその 2 倍、
// 4 倍…の段と、組の数が allow を超えてから前の見積もりの 2 倍を超えた段で行う。後ろ半分の段の組の増え方 rate を毎段
// 増え続ける分とみなし、残り (P) を毎歩使われる組とみなして、DP で積の 1 状態を使い回す回数を (P + rate) n / (P + rate n)
// とする。これが REUSE より少なければやめる。ただし、増え方が段に WAIT 個以下で組が allow 以下のうちは続ける (はじめは
// 上限と等しい道しか無く、2 つの数の上の桁が分かれてから組が広がる形のため)。位置を状態に持つ DFA との積で、どの位置の
// 組も多くの組と並ぶもの (毎段同じだけ組が増える) は、使い回しが無いので表にしない。組が lim を超えてもやめる。allow と
// lim は、表の升目 (組の数 × 記号の数) が 2^16 と 2^20 になる組の数。
inline CompPtr product_table(const std::vector<CompPtr>& cs, int m, int n) {
 constexpr int L_MIN= 16;
 constexpr double REUSE= 32, WAIT= 4;
 const std::size_t lim= (std::size_t(1) << 20) / m, allow= (std::size_t(1) << 16) / m;
 const int k= int(cs.size());
 std::vector<int> sh(k);
 int tb= 0;
 for(int j= 0; j < k; ++j) {
  int b= 0;
  while((1ll << b) < cs[j]->n) ++b;
  sh[j]= tb, tb+= b;
 }
 if(tb > 64) return nullptr;
 HashMap<int> h;
 std::vector<int> dec(k, 0), tab, tg(k);
 std::vector<std::size_t> lv{1};  // lv[d]: 段 d までの組の数
 h[0]= 1;
 std::size_t head= 0, cnt= 1, last= 1;
 int next_check= L_MIN;
 for(int d= 0; d < n && head < cnt; ++d) {
  for(; head < lv[d]; ++head) {
   tab.resize(tab.size() + m);
   int* row= tab.data() + head * m;
   const int* tp= dec.data() + head * k;
   for(int c= 0; c < m; ++c) {
    u64 key= 0;
    bool ok= true;
    for(int j= 0; j < k; ++j) {
     const int t= cs[j]->tab[std::size_t(tp[j]) * m + c];
     if(t < 0) {
      ok= false;
      break;
     }
     tg[j]= t, key|= u64(t) << sh[j];
    }
    if(!ok) {
     row[c]= -1;
     continue;
    }
    int& v= h[key];
    if(!v) {
     v= int(++cnt);
     dec.insert(dec.end(), tg.begin(), tg.end());
     tp= dec.data() + head * k;
    }
    row[c]= v - 1;
   }
  }
  lv.push_back(cnt);
  if(cnt > lim) return nullptr;
  const int L= d + 1;
  if(head < cnt && L >= L_MIN && (L == next_check || (cnt > allow && cnt > 2 * last))) {
   const double rate= double(lv[L] - lv[L / 2]) / (L - L / 2);
   const double P= std::max(0.0, double(cnt) - rate * L), total= P + rate * (n + 1);
   const bool good= total <= double(lim) && (P + rate) * n >= REUSE * total, wait= cnt <= allow && rate <= WAIT;
   if(!good && !wait) return nullptr;
   last= cnt;
   while(next_check <= L) next_check*= 2;
  }
 }
 tab.resize(cnt * m, -1);
 auto C= std::make_shared<Comp>();
 C->n= int(cnt), C->tab= std::move(tab), C->acc.resize(cnt);
 for(std::size_t i= 0; i < cnt; ++i) {
  bool a= true;
  for(int j= 0; j < k && a; ++j) a= cs[j]->acc[dec[i * k + j]];
  C->acc[i]= a;
 }
 return C;
}
// 成分の積を、初期状態から来られる組だけの 1 つの表にしたもの。parts[s * k + j] は表の状態 s の成分 j の状態 (k = 1 なら
// parts は空で、表は成分そのもの)。cs は成分 (状態の値を引くために持つ)。
struct Tab {
 CompPtr C;
 int k= 1;
 std::vector<int> parts;
 std::vector<CompPtr> cs;
};
using TabPtr= std::shared_ptr<const Tab>;
inline TabPtr tabulate(const std::vector<CompPtr>& cs, int m) {
 auto T= std::make_shared<Tab>();
 T->cs= cs, T->k= int(cs.size());
 if(T->k == 1) return T->C= cs[0], T;
 const int k= T->k;
 std::vector<int> sh(k);
 int tb= 0;
 for(int j= 0; j < k; ++j) {
  int b= 0;
  while((1ll << b) < cs[j]->n) ++b;
  sh[j]= tb, tb+= b;
 }
 HashMap<int> h;
 std::unordered_map<std::vector<int>, int, Hasher> hv;
 std::vector<int> dec(k, 0), tab, tg(k);
 int cnt= 0;
 auto id= [&]() -> int {  // tg の組の番号 (無ければ足す)
  int* v;
  if(tb <= 64) {
   u64 key= 0;
   for(int j= 0; j < k; ++j) key|= u64(tg[j]) << sh[j];
   v= &h[key];
  } else v= &hv[tg];
  if(!*v) *v= ++cnt, dec.insert(dec.end(), tg.begin(), tg.end());
  return *v - 1;
 };
 for(int j= 0; j < k; ++j) tg[j]= 0;
 dec.clear(), id();
 for(int s= 0; s < cnt; ++s) {
  tab.resize(tab.size() + m);
  for(int c= 0; c < m; ++c) {
   bool ok= true;
   for(int j= 0; j < k && ok; ++j) {
    const int t= cs[j]->tab[std::size_t(dec[std::size_t(s) * k + j]) * m + c];
    ok= t >= 0, tg[j]= t;
   }
   tab[std::size_t(s) * m + c]= ok ? id() : -1;
  }
 }
 auto C= std::make_shared<Comp>();
 C->n= cnt, C->tab= std::move(tab), C->acc.resize(cnt);
 for(int s= 0; s < cnt; ++s) {
  bool a= true;
  for(int j= 0; j < k && a; ++j) a= cs[j]->acc[dec[std::size_t(s) * k + j]];
  C->acc[s]= a;
 }
 T->C= std::move(C), T->parts= std::move(dec);
 return T;
}
// 表の状態 1 つ。st[j] は成分 j の状態の番号 (成分の中での番号)、get<S>(j) は成分 j の利用者の状態の値。
struct StateRef {
 const Tab* T= nullptr;
 int id= -1;
 int operator[](int j) const { return T->k == 1 ? id : T->parts[std::size_t(id) * T->k + j]; }
 template <class S> const S& get(int j= 0) const {
  const Comp& c= *T->cs[j];
  assert(c.st && *c.st_type == typeid(S));  // 状態の値を持たない成分 (NFA、leq_digits など) か、型が違う
  return (*static_cast<const std::deque<S>*>(c.st.get()))[(*this)[j]];
 }
 bool accepted() const { return T->C->acc[id]; }
};
// DP の遷移 1 回の前後の状態。
struct Step {
 StateRef from, to;
};
// 利用者に返す状態。表を持ち続けるので、DFA が消えても使える。
struct State: StateRef {
 TabPtr hold;
};
// count の f (値をそのまま返す)。run1 はこの型を見て、同じ先へ移る記号をまとめて本数を掛ける。
struct CountF {
 template <class T, class X> const T& operator()(const T& x, const X&, int) const { return x; }
};
template <class T> T sat_add(const T& a, const T& b) {
 if constexpr(std::is_integral_v<T>) {
  T r;
  return __builtin_add_overflow(a, b, &r) ? std::numeric_limits<T>::max() : r;
 } else return a + b;
}
// Berlekamp–Massey。a_i = Σ_{j=1}^{L} c_j a_{i-j} を満たす最短の (c_1, ..., c_L) を返す。T は体。
template <class T> std::vector<T> berlekamp_massey(const std::vector<T>& a) {
 std::vector<T> C{T(1)}, B{T(1)};
 int L= 0, m= 1;
 T b= T(1);
 for(std::size_t i= 0; i < a.size(); ++i, ++m) {
  T d= a[i];
  for(int j= 1; j <= L; ++j) d+= C[j] * a[i - j];
  if(d == T(0)) continue;
  const std::vector<T> Tm= C;
  const T coef= d / b;
  if(C.size() < B.size() + m) C.resize(B.size() + m, T(0));
  for(std::size_t j= 0; j < B.size(); ++j) C[j + m]-= coef * B[j];
  if(2 * L <= int(i)) L= int(i) + 1 - L, B= Tm, b= d, m= 0;
 }
 C.resize(L + 1, T(0));
 std::vector<T> r(L);
 for(int j= 0; j < L; ++j) r[j]= -C[j + 1];
 return r;
}
// 漸化式 c (berlekamp_massey の形) と最初の項 a から、a_n を x^n mod (特性多項式) で求める (掛け算は素朴に O(L^2))。
template <class T> T linear_recurrence_nth(const std::vector<T>& a, const std::vector<T>& c, unsigned long long n) {
 const int L= int(c.size());
 if(n < a.size()) return a[n];
 if(L == 0) return T(0);
 // p, q は x^k mod P の係数 (次数 < L)。P(x) = x^L - Σ c_j x^{L-j}
 auto mul= [&](const std::vector<T>& p, const std::vector<T>& q) {
  std::vector<T> r(2 * L - 1, T(0));
  for(int i= 0; i < L; ++i)
   if(!(p[i] == T(0)))
    for(int j= 0; j < L; ++j) r[i + j]+= p[i] * q[j];
  for(int i= 2 * L - 2; i >= L; --i)
   if(!(r[i] == T(0)))
    for(int j= 1; j <= L; ++j) r[i - j]+= r[i] * c[j - 1];
  r.resize(L);
  return r;
 };
 std::vector<T> res(L, T(0)), base(L, T(0));
 res[0]= T(1);
 if(L == 1) base[0]= c[0];
 else base[1]= T(1);
 for(; n; n>>= 1, base= mul(base, base))
  if(n & 1) res= mul(res, base);
 T ans= T(0);
 for(int i= 0; i < L; ++i) ans+= res[i] * a[i];
 return ans;
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
 mutable dfa_internal::TabPtr tabc;  // 積を 1 つの表にしたもの (状態を見る口で作る)
 const dfa_internal::TabPtr& table() const {
  if(!tabc) tabc= dfa_internal::tabulate(comps, int(alph.size()));
  return tabc;
 }
 DFA()= default;
 static std::vector<symbol_t> sorted(std::vector<symbol_t> a) { return std::sort(a.begin(), a.end()), a; }
 // 成分が 1 つのとき。状態の番号をそのまま来た組の置き場所を引く添字にする。
 // gcc は run1 を run に展開すると、積を表にした道 (aoj-0570) の DP が 2 割ほど遅くなるので、展開させない。
 // 数え上げを、状態ごとの (移る先, 記号の本数) の組 (gofs、gto、gmul) で進める。run1 の速い道。
 template <bool ALL, class T, class A> [[gnu::noinline]] std::vector<T> run1g(const Comp& C, int n, const A& op, const T& e, const T& init, const std::vector<int>& gofs, const std::vector<int>& gto, const std::vector<int>& gmul) const {
  std::vector<T> out, val{init}, nval;
  std::vector<int> slot(C.n), key{0}, nkey;
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
    const int u= key[s];
    const T& v= val[s];
    for(int j= gofs[u]; j < gofs[u + 1]; ++j) {
     const int t= gto[j];
     const T x= gmul[j] == 1 ? v : v * T(gmul[j]);
     if(int& p= slot[t]; !p) p= int(nkey.size()) + 1, nkey.push_back(t), nval.push_back(x);
     else nval[p - 1]= op(nval[p - 1], x);
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
 // STEP なら f に遷移の前後の状態 (Step) も渡す (C は tp の表)。fk が null でなければ、長さ n で来た状態と値を fk と fv に入れる。
 template <bool ALL, bool STEP, class T, class A, class F> [[gnu::noinline]] std::vector<T> run1(const Comp& C, int n, const A& op, const T& e, const F& f, const T& init, const dfa_internal::Tab* tp= nullptr, std::vector<int>* fk= nullptr, std::vector<T>* fv= nullptr) const {
  const int m= int(alph.size());
  std::vector<T> out, val{init}, nval;
  std::vector<int> slot(C.n), key{0}, nkey;
  const int* tab= C.tab.data();
  // 数え上げで記号が多いときは、同じ先へ移る記号をまとめて、移る先と本数の組で進める (run1g)。始めに全部の状態をまとめてみて、
  // 組の数が遷移の数の 65% 以下のときだけ使う。同じ先へ多くの記号が集まると、その先の値への足し込みが 1 本の待ちの鎖になる
  // ので、まとめると速い (26 文字で 20 語の Aho–Corasick で 2 倍、3 語で 7 倍)。先が分かれている DFA では、まとめた表を引く
  // 分で 2〜4 割遅く、数字の 10 個では本数を掛ける手間で遅くなった (gcc で yuki-1417 が 1.6 倍) ので、記号が 16 個以上に限る。
  if constexpr(!STEP && std::is_same_v<F, dfa_internal::CountF> && requires(const T& a) { a * T(2); })
   if(m >= 16 && !fk) {
    std::vector<int> gpos(C.n, 0), gofs(C.n + 1, 0), gto, gmul;
    std::size_t valid= 0;
    for(int u= 0; u < C.n; ++u) {
     const int* row= tab + std::size_t(u) * m;
     for(int c= 0; c < m; ++c)
      if(const int t= row[c]; t >= 0) {
       ++valid;
       if(int& p= gpos[t]; p) ++gmul[p - 1];
       else p= int(gto.size()) + 1, gto.push_back(t), gmul.push_back(1);
      }
     for(std::size_t j= gofs[u]; j < gto.size(); ++j) gpos[gto[j]]= 0;
     gofs[u + 1]= int(gto.size());
    }
    if(gto.size() * 100 <= valid * 65) return run1g<ALL>(C, n, op, e, init, gofs, gto, gmul);
   }
  for(int i= 0;; ++i) {
   if(ALL || i == n) {
    T r= e;
    for(std::size_t s= 0; s < key.size(); ++s)
     if(C.acc[key[s]]) r= op(r, val[s]);
    out.push_back(r);
   }
   if(i == n) {
    if(fk) *fk= key, *fv= val;
    break;
   }
   nkey.clear(), nval.clear();
   for(std::size_t s= 0; s < key.size(); ++s) {
    const int* row= tab + std::size_t(key[s]) * m;
    const T& v= val[s];
    for(int c= 0; c < m; ++c) {
     const int t= row[c];
     if(t < 0) continue;
     auto g= [&]() -> T {
      if constexpr(STEP) return f(v, alph[c], i, Step{{tp, key[s]}, {tp, t}});
      else return f(v, alph[c], i);
     };
     if(int& p= slot[t]; !p) p= int(nkey.size()) + 1, nkey.push_back(t), nval.push_back(g());
     else nval[p - 1]= op(nval[p - 1], g());
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
  if constexpr(std::is_invocable_v<const F&, const T&, const symbol_t&, int, const Step&>) {  // f が状態を見る
   const dfa_internal::Tab& T0= *table();
   return run1<ALL, true>(*T0.C, n, op, e, f, init, &T0);
  } else {
   std::vector<CompPtr> cs= comps;
   if(cs.size() >= 2 && n >= 32)
    if(CompPtr P= dfa_internal::product_table(cs, int(alph.size()), n)) cs= {std::move(P)};
   std::vector<T> out;
   for(int a, b;;) {
    if(cs.size() == 1) return run1<ALL, false>(*cs[0], n, op, e, f, init);
    if(runk<ALL>(cs, n, op, e, f, init, out, a, b)) return out;
    cs[a]= dfa_internal::product(*cs[a], *cs[b], int(alph.size())), cs.erase(cs.begin() + b);
   }
  }
 }
public:
 // 遷移 1 回の前後の状態 (dp の f の 4 つ目の引数)。t.to[j] は & の左から j 番目の成分の状態の番号、t.to.get<S>(j) は
 // その成分の利用者の状態の値 (S は構築に渡した状態の型)。
 using Step= dfa_internal::Step;
 // 状態 (count_states と dp_states の答え)。使い方は Step の from と to と同じで、accepted() は受理されるか。
 using State= dfa_internal::State;
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
   D->n= C->n, D->acc= C->acc, D->st= C->st, D->st_type= C->st_type, D->tab.resize(std::size_t(C->n) * m2);
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
 template <class T> T count(int n) const { return dp(n, std::plus<>(), T(), dfa_internal::CountF(), T(1)); }
 // 受理される長さ 0 から n までの列の数。
 template <class T> std::vector<T> count_all(int n) const { return dp_all(n, std::plus<>(), T(), dfa_internal::CountF(), T(1)); }
 // 長さ n の列で来られる状態ごとの dp の値 (受理されない状態も含む)。積は来られる組だけの表にしてから進める。
 template <class T, class A, class F> std::vector<std::pair<State, T>> dp_states(int n, const A& op, const T& e, const F& f, const T& init) const {
  std::vector<std::pair<State, T>> ret;
  if(n < 0) return ret;
  const dfa_internal::TabPtr& T0= table();
  std::vector<int> fk;
  std::vector<T> fv;
  if constexpr(std::is_invocable_v<const F&, const T&, const symbol_t&, int, const Step&>) run1<false, true>(*T0->C, n, op, e, f, init, T0.get(), &fk, &fv);
  else run1<false, false>(*T0->C, n, op, e, f, init, T0.get(), &fk, &fv);
  for(std::size_t i= 0; i < fk.size(); ++i) {
   State st;
   st.T= T0.get(), st.id= fk[i], st.hold= T0;
   ret.emplace_back(std::move(st), std::move(fv[i]));
  }
  return ret;
 }
 // 長さ n の列で来られる状態ごとの、そこに来る列の数。
 template <class T> std::vector<std::pair<State, T>> count_states(int n) const { return dp_states(n, std::plus<>(), T(), dfa_internal::CountF(), T(1)); }
 // 受理される列を、記号を並べた順の辞書順で数える表 (ranker が返す)。状態と残りの長さごとの受理される列の数を後ろからの
 // DP で前もって作る (大きさは (長さの上限 + 1) × 積の来られる組の数)。T が整数型なら、和は型の最大で止める。
 template <class T> class Ranker {
  friend class DFA;
  dfa_internal::TabPtr T0;
  std::vector<symbol_t> alph;
  int L= 0;
  std::vector<T> cnt;  // cnt[r * 状態の数 + s]: 状態 s から長さ r の列で受理に至る数
  const T& at(int r, int s) const { return cnt[std::size_t(r) * T0->C->n + s]; }
  const int* row(int s) const { return T0->C->tab.data() + std::size_t(s) * alph.size(); }
 public:
  // 受理される長さ len (≤ 作るときの長さ) の列の数。
  T count(int len) const { return assert(0 <= len && len <= L), at(len, 0); }
  // 受理される長さ |w| の列のうち、w より辞書順で小さいもの。w に記号でないものがあれば、それより小さい記号で比べる。
  T count_less(const std::vector<symbol_t>& w) const {
   const int n= int(w.size());
   assert(n <= L);
   T r= T(0);
   for(int i= 0, s= 0; i < n; ++i) {
    const int p= int(std::lower_bound(alph.begin(), alph.end(), w[i]) - alph.begin());
    const int* rw= row(s);
    for(int c= 0; c < p; ++c)
     if(rw[c] >= 0) r= dfa_internal::sat_add(r, at(n - 1 - i, rw[c]));
    if(p == int(alph.size()) || !(alph[p] == w[i]) || rw[p] < 0) return r;
    s= rw[p];
   }
   return r;
  }
  // w が受理されるか。
  bool accepts(const std::vector<symbol_t>& w) const {
   int s= 0;
   for(const symbol_t& x: w) {
    const int p= int(std::lower_bound(alph.begin(), alph.end(), x) - alph.begin());
    if(p == int(alph.size()) || !(alph[p] == x) || (s= row(s)[p]) < 0) return false;
   }
   return T0->C->acc[s];
  }
  // 受理される長さ |w| の列のうち、w 以下のもの。
  T count_leq(const std::vector<symbol_t>& w) const { return dfa_internal::sat_add(count_less(w), T(accepts(w) ? 1 : 0)); }
  // 受理される長さ len の列のうち、辞書順で k 番目 (0 始まり) のもの。k が count(len) 以上なら空を返す。
  std::vector<symbol_t> kth(T k, int len) const {
   assert(0 <= len && len <= L);
   if(!(k < at(len, 0))) return {};
   std::vector<symbol_t> w;
   for(int i= 0, s= 0; i < len; ++i) {
    const int* rw= row(s);
    for(int c= 0;; ++c) {
     if(rw[c] < 0) continue;
     const T& x= at(len - 1 - i, rw[c]);
     if(k < x) {
      w.push_back(alph[c]), s= rw[c];
      break;
     }
     k-= x;
    }
   }
   return w;
  }
 };
 // 長さ L 以下の列についての Ranker を作る。
 template <class T> Ranker<T> ranker(int L) const {
  Ranker<T> R;
  R.T0= table(), R.alph= alph, R.L= L;
  const dfa_internal::Comp& C= *R.T0->C;
  const int S= C.n, m= int(alph.size());
  R.cnt.assign(std::size_t(L + 1) * S, T(0));
  for(int s= 0; s < S; ++s) R.cnt[s]= T(C.acc[s] ? 1 : 0);
  for(int r= 1; r <= L; ++r) {
   const T* prv= R.cnt.data() + std::size_t(r - 1) * S;
   T* cur= R.cnt.data() + std::size_t(r) * S;
   for(int s= 0; s < S; ++s) {
    const int* rw= C.tab.data() + std::size_t(s) * m;
    T x= T(0);
    for(int c= 0; c < m; ++c)
     if(rw[c] >= 0) x= dfa_internal::sat_add(x, prv[rw[c]]);
    cur[s]= x;
   }
  }
  return R;
 }
 // 受理される列が長さによらず有限個か (初期状態から来られて受理に至れる状態に閉路が無いか)。
 bool is_finite() const { return finite_order(nullptr); }
 // 受理される列の、長さによらない数。有限個でなければならない (is_finite)。部分列の自動機のような閉路の無い DFA で、長さごとに
 // 進めずに状態 × 記号の手間で数える。
 template <class T> T count_any() const {
  std::vector<int> ord;
  const bool fin= finite_order(&ord);
  assert(fin);
  if(!fin) return T(0);
  const dfa_internal::Comp& C= *table()->C;
  const int m= int(alph.size());
  std::vector<T> w(C.n, T(0));
  w[0]= T(1);
  T r= T(0);
  for(int s: ord) {
   if(C.acc[s]) r= dfa_internal::sat_add(r, w[s]);
   const int* rw= C.tab.data() + std::size_t(s) * m;
   for(int c= 0; c < m; ++c)
    if(rw[c] >= 0) w[rw[c]]= dfa_internal::sat_add(w[rw[c]], w[s]);
  }
  return r;
 }
 // 受理される長さ n の列の数を、n が大きくても求める。count_all で最初の 2 × (積の来られる組の数) 項を出し、Berlekamp–Massey で
 // 漸化式にしてから n 項目を求める。T は体 (ZMod の素数の法など)。手間は組の数 S について、S^2 × 記号の数と S^2 log n。
 template <class T> T count_large(unsigned long long n) const {
  const int S= table()->C->n;
  const int K= 2 * S + 2;
  if(n <= (unsigned long long)K) return count<T>(int(n));
  const std::vector<T> a= count_all<T>(K);
  return dfa_internal::linear_recurrence_nth(a, dfa_internal::berlekamp_massey(a), n);
 }
private:
 // 初期状態から来られて受理に至れる状態だけで閉路が無ければ、その位相順を ord に入れて true。
 bool finite_order(std::vector<int>* ord) const {
  const dfa_internal::Comp& C= *table()->C;
  const int S= C.n, m= int(alph.size());
  std::vector<std::vector<int>> rev(S);
  for(int s= 0; s < S; ++s)
   for(int c= 0; c < m; ++c)
    if(const int t= C.tab[std::size_t(s) * m + c]; t >= 0) rev[t].push_back(s);
  std::vector<char> co(S, 0);  // 受理に至れる
  std::vector<int> st;
  for(int s= 0; s < S; ++s)
   if(C.acc[s]) co[s]= 1, st.push_back(s);
  while(!st.empty()) {
   const int t= st.back();
   st.pop_back();
   for(int s: rev[t])
    if(!co[s]) co[s]= 1, st.push_back(s);
  }
  // 初期状態から受理に至れる状態だけをたどって来られる状態 q を集め、その中の辺で Kahn の方法の位相順を作る
  std::vector<int> indeg(S, 0), q;
  std::vector<char> reach(S, 0);
  if(co[0]) q.push_back(0), reach[0]= 1;
  for(std::size_t i= 0; i < q.size(); ++i)
   for(int c= 0; c < m; ++c)
    if(const int t= C.tab[std::size_t(q[i]) * m + c]; t >= 0 && co[t] && !reach[t]) reach[t]= 1, q.push_back(t);
  for(int s: q)
   for(int c= 0; c < m; ++c)
    if(const int t= C.tab[std::size_t(s) * m + c]; t >= 0 && reach[t]) ++indeg[t];
  std::vector<int> o;
  if(!q.empty() && indeg[0] == 0) o.push_back(0);
  for(std::size_t i= 0; i < o.size(); ++i)
   for(int c= 0; c < m; ++c)
    if(const int t= C.tab[std::size_t(o[i]) * m + c]; t >= 0 && reach[t] && --indeg[t] == 0) o.push_back(t);
  if(o.size() != q.size()) return false;
  if(ord) *ord= std::move(o);
  return true;
 }
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
