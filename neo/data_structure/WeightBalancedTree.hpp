#pragma once
#include <algorithm>
#include <cstddef>
#include <new>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#ifdef __linux__
#include <sys/mman.h>
#endif
// 列の重み平衡木 (今の Library の WeightBalancedTree の中身を写したもの)。値は葉にだけ持ち、左右の葉の数の比を 4 以下に保つ。
// split と併合は、重い方の背を下りて回転 1 回で直す join で書く。構造体 M の決まりは SegmentTree と同じで、モノイドは T、ti、op、
// 作用は E、mp、cp (mp は区間の長さを受けてもよい)。M の代わりに int などの型を渡すと値を並べるだけの列になる。REV で反転を
// 選び、可換でないモノイドを反転するときは逆順の積も持つ (M に using commute = void を書くと省く)。
// PERSIST が真なら永続で、どの操作も道の上の節を写して新しい根を作るので、木の値を写せば版の写しになる。永続しないときは節を
// その場で書き換え、木の値を写すと同じ木を指すだけになる (切ったり併合したりしたあとは元の値を使わない)。
// - 同じ型の木は節の置き場を共有するので、別の木どうしで切り貼りできる。置き場は伸びる配列で、1 MB 以上は 2 MB 境界に置いて
//   huge page を頼む (Linux のみ)。永続では操作のたびに節を新しく作るので、全体の大きさが分かれば reserve で先に取ると、配列を
//   作り直して写す分を省ける (procon-judge の aoj-3024 で、先に取らないと今の Library の 1.17 倍、取ると並んだ)。
// - 永続のときは、pool_empty() が置き場の伸び (前の gc で生きていた節の 2 倍、少なくとも 2^22 個) を知らせ、gc(版の列) が渡した
//   版から届く節に印を付けて置き場の前へ詰め、版の根を書き換える。版どうしで共有する節は共有したまま残り、渡さなかった木は
//   使えなくなる。永続しないときは、split で外れた節と消した葉を使い回すので gc は要らない (pool_empty() は偽を返す)。
// - 遅延の印 (作用と反転) は葉の数の上位 2 bit に置く。区間の積と 1 点の取得は、遅延を配らずに根から合わせて下りるので、永続の
//   ときにも節を写さない。
namespace weight_balanced_tree_internal {
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
// 1 MB 以上の確保は 2 MB 境界に寄せ、Linux では huge page を頼む。
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
}  // namespace weight_balanced_tree_internal
template <class M, bool REV= false, bool PERSIST= false> class WeightBalancedTree {
 using T= typename weight_balanced_tree_internal::value_of<M>::type;
 using E= typename weight_balanced_tree_internal::action_of<M>::type;
 static constexpr bool MO= weight_balanced_tree_internal::monoid<M>::value, AC= weight_balanced_tree_internal::action<M>::value;
 static constexpr bool RS= MO && REV && !weight_balanced_tree_internal::commute<M>::value;  // 逆順の積も持つ
 static constexpr bool FL= AC || REV;  // 葉の数の上位 bit に遅延の印を置く
 static constexpr unsigned LZ= 0x80000000u, TG= 0x40000000u, SZ= 0x3fffffffu;
 using S= std::conditional_t<MO, T, char>;
 using SR= std::conditional_t<RS, T, char>;
 using EA= std::conditional_t<AC, E, char>;
 struct Node {
  int ch[2];    // 正なら内側の節、負なら葉 (値は nl[-c])
  unsigned sz;  // 葉の数。FL なら上位 2 bit が作用 (LZ) と反転 (TG) の印 (子へ配る分で、自分の sum と rsum には当ててある)
  S sum;
  SR rsum;
  EA laz;
 };
 template <class U> using Vec= std::vector<U, weight_balanced_tree_internal::Alloc<U>>;
 static inline Vec<Node> nm= Vec<Node>(1);  // 0 は空の木
 static inline Vec<T> nl= Vec<T>(1);
 static inline std::vector<int> free_nm, free_nl;
 static inline std::size_t gc_at= std::size_t(1) << 22;
 static void map1(T& v, const E& e, int len) {
  if constexpr(std::is_invocable_v<decltype(M::mp), T&, const E&, int>) M::mp(v, e, len);
  else M::mp(v, e);
 }
 static unsigned msize(int i) {
  if constexpr(FL) return nm[i].sz & SZ;
  else return nm[i].sz;
 }
 static unsigned size(int i) { return i < 0 ? 1 : msize(i); }
 static T sum(int i) { return i < 0 ? nl[-i] : T(nm[i].sum); }
 static T rsum(int i) {
  if constexpr(RS) return i < 0 ? nl[-i] : T(nm[i].rsum);
  else return sum(i);
 }
 static int new_node() {
  if constexpr(!PERSIST)
   if(!free_nm.empty()) {
    const int i= free_nm.back();
    free_nm.pop_back();
    return i;
   }
  nm.emplace_back();
  return int(nm.size() - 1);
 }
 static int copy_node(int i) {
  const Node n= nm[i];
  const int j= new_node();
  nm[j]= n;
  return j;
 }
 static int new_leaf(T v) {
  if constexpr(!PERSIST)
   if(!free_nl.empty()) {
    const int i= free_nl.back();
    free_nl.pop_back();
    nl[i]= v;
    return -i;
   }
  nl.push_back(v);
  return -int(nl.size() - 1);
 }
 static void release(int c) {
  if constexpr(!PERSIST) {
   if(c > 0) free_nm.push_back(c);
   else if(c < 0) free_nl.push_back(-c);
  }
 }
 // 子から求め直す。SIZE が真なら葉の数も求め直す (印は消えるので、遅延を配ったあとの節だけ)。
 template <bool SIZE= true> static void update(int i) {
  const int l= nm[i].ch[0], r= nm[i].ch[1];
  if constexpr(SIZE) nm[i].sz= size(l) + size(r);
  if constexpr(MO) {
   nm[i].sum= M::op(sum(l), sum(r));
   if constexpr(RS) nm[i].rsum= M::op(rsum(r), rsum(l));
  }
 }
 static void propagate(int i, const E& x) {
  Node& t= nm[i];
  if(t.sz & LZ) M::cp(t.laz, x);
  else t.laz= x;
  t.sz|= LZ;
  if constexpr(MO) {
   map1(t.sum, x, t.sz & SZ);
   if constexpr(RS) map1(t.rsum, x, t.sz & SZ);
  }
 }
 static void toggle(int i) {
  Node& t= nm[i];
  if constexpr(RS) std::swap(t.sum, t.rsum);
  std::swap(t.ch[0], t.ch[1]), t.sz^= TG;
 }
 // 子 c へ印 f と作用 lz を配る。永続なら c を写してから当てる。
 static void push_child(int& c, unsigned f, const EA& lz) {
  if(c > 0) {
   if constexpr(PERSIST) c= copy_node(c);
   if constexpr(AC)
    if(f & LZ) propagate(c, lz);
   if constexpr(REV)
    if(f & TG) toggle(c);
  } else if constexpr(AC) {
   if(c < 0 && (f & LZ)) {
    if constexpr(PERSIST) c= new_leaf(nl[-c]);
    map1(nl[-c], lz, 1);
   }
  }
 }
 // 節 i の印を子へ配る。i の部分木の表す列は変わらないので、永続でも i 自身はその場で書き換えてよい。
 static void push(int i) {
  if constexpr(FL) {
   const unsigned f= nm[i].sz;
   if(!(f >> 30)) return;
   EA lz{};
   if constexpr(AC) lz= nm[i].laz;
   int l= nm[i].ch[0], r= nm[i].ch[1];
   push_child(l, f, lz), push_child(r, f, lz);
   Node& t= nm[i];
   t.ch[0]= l, t.ch[1]= r, t.sz= f & SZ;
  }
 }
 // a が重い側。R が偽なら a の右の背を下りて b をつなぎ、真なら a の左の背を下りて b をつなぐ。t はつなぎ目に使う節。
 template <bool R> static int join_(int t, int a, int b) {
  if constexpr(FL) push(a);
  if constexpr(R) b= join(b, t, nm[a].ch[0]);
  else b= join(nm[a].ch[1], t, b);
  if constexpr(PERSIST) a= copy_node(a);
  if(size(nm[a].ch[R]) * 4 >= msize(b)) return nm[a].ch[!R]= b, update(a), a;
  return nm[a].ch[!R]= nm[b].ch[R], update(a), nm[b].ch[R]= a, update(b), b;
 }
 static int join(int l, int t, int r) {
  const unsigned lsz= size(l), rsz= size(r);
  if(lsz > rsz * 4) return join_<0>(t, l, r);
  if(rsz > lsz * 4) return join_<1>(t, r, l);
  nm[t].ch[0]= l, nm[t].ch[1]= r, update(t);
  return t;
 }
 static int merge(int l, int r) { return !l ? r : !r ? l : join(l, new_node(), r); }
 // 0 < k < size(i) のとき、k の位置で割った左と右と、割れ目にあった節を返す。
 static std::tuple<int, int, int> split_(int i, unsigned k) {
  if constexpr(FL) push(i);
  const int l= nm[i].ch[0], r= nm[i].ch[1];
  const unsigned lsz= size(l);
  if(k == lsz) return {l, i, r};
  if constexpr(PERSIST) i= new_node();
  if(k < lsz) {
   const auto [a, b, c]= split_(l, k);
   return {a, b, join(c, i, r)};
  }
  const auto [a, b, c]= split_(r, k - lsz);
  return {join(l, i, a), b, c};
 }
 static std::pair<int, int> split(int i, unsigned k) {
  if(k == 0) return {0, i};
  if(k >= size(i)) return {i, 0};
  const auto [l, c, r]= split_(i, k);
  release(c);
  return {l, r};
 }
 // i の部分木の列の順の [l, r) の積を、遅延を配らずに求める。fl は祖先から受けた反転の印。
 static T prod_(int i, unsigned l, unsigned r, bool fl) {
  if(i < 0) return nl[-i];
  const unsigned f= nm[i].sz, s= f & SZ;
  if(l == 0 && r == s) return fl ? rsum(i) : sum(i);
  const int a= nm[i].ch[fl], b= nm[i].ch[!fl];
  const bool cf= fl ^ bool(f & TG);
  const unsigned asz= size(a);
  T res;
  if(r <= asz) res= prod_(a, l, r, cf);
  else if(asz <= l) res= prod_(b, l - asz, r - asz, cf);
  else res= M::op(prod_(a, l, asz, cf), prod_(b, 0, r - asz, cf));
  if constexpr(AC)
   if(f & LZ) map1(res, nm[i].laz, r - l);
  return res;
 }
 // 根から葉まで下り、葉の値に下の節の作用から順に当てる (下の節の作用ほど古い)。
 static T get_(int i, unsigned k) {
  static std::vector<int> ps;
  ps.clear();
  bool fl= false;
  while(i > 0) {
   ps.push_back(i);
   const unsigned f= nm[i].sz;
   const int a= nm[i].ch[fl], b= nm[i].ch[!fl];
   fl^= bool(f & TG);
   const unsigned asz= size(a);
   if(k < asz) i= a;
   else k-= asz, i= b;
  }
  T v= nl[-i];
  if constexpr(AC)
   for(int j= (int)ps.size(); j--;)
    if(nm[ps[j]].sz & LZ) map1(v, nm[ps[j]].laz, 1);
  return v;
 }
 static int apply_(int i, unsigned l, unsigned r, const E& x) {
  if(i < 0) {
   if constexpr(PERSIST) i= new_leaf(nl[-i]);
   map1(nl[-i], x, 1);
   return i;
  }
  if constexpr(PERSIST) i= copy_node(i);
  if(l == 0 && msize(i) <= r) return propagate(i, x), i;
  push(i);
  int a= nm[i].ch[0], b= nm[i].ch[1];
  const unsigned asz= size(a);
  if(r <= asz) a= apply_(a, l, r, x);
  else if(asz <= l) b= apply_(b, l - asz, r - asz, x);
  else a= apply_(a, l, asz, x), b= apply_(b, 0, r - asz, x);
  nm[i].ch[0]= a, nm[i].ch[1]= b;
  if constexpr(MO) update<false>(i);
  return i;
 }
 // k 番目の葉を f(今の値) にする。道の上の遅延は配ってから下りるので、葉の値は今の値になっている。
 template <class F> static int modify_(int i, unsigned k, const F& f) {
  if(i < 0) {
   const T v= f(nl[-i]);
   if constexpr(PERSIST) return new_leaf(v);
   nl[-i]= v;
   return i;
  }
  push(i);
  if constexpr(PERSIST) i= copy_node(i);
  int a= nm[i].ch[0], b= nm[i].ch[1];
  const unsigned asz= size(a);
  if(k < asz) a= modify_(a, k, f);
  else b= modify_(b, k - asz, f);
  nm[i].ch[0]= a, nm[i].ch[1]= b;
  if constexpr(MO) update<false>(i);
  return i;
 }
 template <class F> static int build(const F& at, int lo, int hi) {
  if(hi - lo == 1) return new_leaf(at(lo));
  const int m= (lo + hi) / 2;
  const int a= build(at, lo, m), b= build(at, m, hi);
  const int i= new_node();
  nm[i].ch[0]= a, nm[i].ch[1]= b, update(i);
  return i;
 }
 // 置き場に n 個の空きを作る。足りなければ倍に伸ばす (ちょうどの大きさで取ると、小さい木を何本も作るたびに作り直すため)。
 template <class V> static void make_room(V& v, std::size_t n) {
  if(v.capacity() < v.size() + n) v.reserve(std::max(v.size() + n, 2 * v.capacity()));
 }
 // n 個の葉と n - 1 個の節を作るので、先に空きを作ってから組む。
 template <class F> static int build(const F& at, int n) {
  if(!n) return 0;
  make_room(nm, n), make_room(nl, n);
  return build(at, 0, n);
 }
 static void dump_(int i, bool fl, std::vector<T>& out, std::vector<E>& acts) {
  if(i < 0) {
   T v= nl[-i];
   if constexpr(AC)
    for(int j= (int)acts.size(); j--;) map1(v, acts[j], 1);
   out.push_back(v);
   return;
  }
  const unsigned f= nm[i].sz;
  const int a= nm[i].ch[fl], b= nm[i].ch[!fl];
  const bool cf= fl ^ bool(f & TG);
  if constexpr(AC)
   if(f & LZ) acts.push_back(nm[i].laz);
  dump_(a, cf, out, acts), dump_(b, cf, out, acts);
  if constexpr(AC)
   if(f & LZ) acts.pop_back();
 }
 // 渡した根から届く節に印を付け、置き場の前へ詰める。根の番号は書き換える。
 static void gc_roots(const std::vector<int*>& roots) {
  std::vector<char> mi(nm.size()), ml(nl.size());
  std::vector<int> st;
  for(int* r: roots) st.push_back(*r);
  while(!st.empty()) {
   const int c= st.back();
   st.pop_back();
   if(c > 0) {
    if(mi[c]) continue;
    mi[c]= 1, st.push_back(nm[c].ch[0]), st.push_back(nm[c].ch[1]);
   } else if(c < 0) ml[-c]= 1;
  }
  std::vector<int> ni(nm.size()), nj(nl.size());
  int k= 1, kl= 1;
  for(std::size_t i= 1; i < nm.size(); ++i)
   if(mi[i]) ni[i]= k++;
  for(std::size_t i= 1; i < nl.size(); ++i)
   if(ml[i]) nj[i]= kl++;
  auto remap= [&](int c) { return c > 0 ? ni[c] : c < 0 ? -nj[-c] : 0; };
  for(std::size_t i= 1; i < nm.size(); ++i)
   if(mi[i]) {
    Node n= nm[i];
    n.ch[0]= remap(n.ch[0]), n.ch[1]= remap(n.ch[1]);
    nm[ni[i]]= n;
   }
  for(std::size_t i= 1; i < nl.size(); ++i)
   if(ml[i]) nl[nj[i]]= nl[i];
  nm.resize(k), nl.resize(kl);
  free_nm.clear(), free_nl.clear();
  for(int* r: roots) *r= remap(*r);
  gc_at= std::max<std::size_t>(std::size_t(1) << 22, 2 * (nm.size() + nl.size()));
 }
 int root= 0;
 static WeightBalancedTree of(int rt) {
  WeightBalancedTree s;
  s.root= rt;
  return s;
 }
public:
 WeightBalancedTree()= default;
 explicit WeightBalancedTree(int n, const T& v= T()) { root= build([&](int) { return v; }, n); }
 WeightBalancedTree(const std::vector<T>& v) { root= build([&](int i) { return v[i]; }, (int)v.size()); }
 template <class It, std::enable_if_t<!std::is_integral_v<It>, int> = 0> WeightBalancedTree(It bg, It ed) {
  const std::vector<T> v(bg, ed);
  root= build([&](int i) { return v[i]; }, (int)v.size());
 }
 // 置き場を節と葉それぞれ n 個ぶん先に取る (配列を作り直して写す分を省く)。
 static void reserve(std::size_t n) { nm.reserve(n), nl.reserve(n); }
 int size() const { return root ? (int)size(root) : 0; }
 bool empty() const { return !root; }
 void clear() { root= 0; }
 std::pair<WeightBalancedTree, WeightBalancedTree> split(int k) const {
  const auto [l, r]= split(root, k);
  return {of(l), of(r)};
 }
 std::tuple<WeightBalancedTree, WeightBalancedTree, WeightBalancedTree> split3(int a, int b) const {
  const auto [tmp, r]= split(root, b);
  const auto [l, c]= split(tmp, a);
  return {of(l), of(c), of(r)};
 }
 WeightBalancedTree& operator+=(const WeightBalancedTree& rhs) { return root= merge(root, rhs.root), *this; }
 WeightBalancedTree operator+(const WeightBalancedTree& rhs) const { return of(merge(root, rhs.root)); }
 void insert(int k, const T& v) {
  const auto [l, r]= split(root, k);
  root= merge(merge(l, new_leaf(v)), r);
 }
 void push_back(const T& v) { root= merge(root, new_leaf(v)); }
 void push_front(const T& v) { root= merge(new_leaf(v), root); }
 T erase(int k) {
  const auto [l, tmp]= split(root, k);
  const auto [t, r]= split(tmp, 1);
  const T v= nl[-t];
  release(t);
  root= merge(l, r);
  return v;
 }
 T pop_back() { return erase(size() - 1); }
 T pop_front() { return erase(0); }
 T get(int k) const { return get_(root, k); }
 T operator[](int k) const { return get_(root, k); }
 void set(int k, const T& v) {
  root= modify_(root, k, [&](const T&) { return v; });
 }
 // k 番目を op(今の値, v) にする。
 void mul(int k, const T& v) {
  static_assert(MO, "mul はモノイドがあるときだけ");
  root= modify_(root, k, [&](const T& x) { return M::op(x, v); });
 }
 T prod(int l, int r) const {
  static_assert(MO, "prod はモノイドがあるときだけ");
  if(l >= r) return M::ti();
  return prod_(root, l, r, false);
 }
 T all_prod() const {
  static_assert(MO, "all_prod はモノイドがあるときだけ");
  return root ? sum(root) : M::ti();
 }
 void apply(int l, int r, const E& x) {
  static_assert(AC, "apply は作用があるときだけ");
  if(l < r) root= apply_(root, l, r, x);
 }
 void reverse() {
  static_assert(REV, "reverse は REV のときだけ");
  if(root > 0) {
   if constexpr(PERSIST) root= copy_node(root);
   toggle(root);
  }
 }
 void reverse(int a, int b) {
  static_assert(REV, "reverse は REV のときだけ");
  if(b - a < 2) return;
  const auto [tmp, r]= split(root, b);
  auto [l, c]= split(tmp, a);
  if constexpr(PERSIST) c= copy_node(c);
  toggle(c);
  root= merge(merge(l, c), r);
 }
 std::vector<T> dump() const {
  std::vector<T> out;
  std::vector<E> acts;
  out.reserve(size());
  if(root) dump_(root, false, out, acts);
  return out;
 }
 // 永続のとき、置き場が前の gc のあとに生きていた節の 2 倍 (少なくとも 2^22 個) まで伸びたら真。
 static bool pool_empty() { return PERSIST && nm.size() + nl.size() >= gc_at; }
 // 渡した版から届く節だけを残して置き場を詰める。渡さなかった木は使えなくなる。
 static void gc(std::vector<WeightBalancedTree>& vs) {
  std::vector<int*> rs;
  rs.reserve(vs.size());
  for(auto& v: vs) rs.push_back(&v.root);
  gc_roots(rs);
 }
 template <class... Ps> static void gc(WeightBalancedTree& v, Ps&... vs) {
  std::vector<int*> rs{&v.root, &vs.root...};
  gc_roots(rs);
 }
};
