#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <cstring>
#include <optional>
#include <type_traits>
#include <vector>
// int の順序付き集合。B+ 木で、葉はキーを 63 個まで昇順に持って前後の葉とつなぎ、内部の節点は子を 63 個まで持って、子ごとに
// 部分木のキーの最大を持つ。キーで下りるときは、最大の列のうち x 未満 (または x 以下) のものの数を AVX2 で 8 個ずつ比べて
// 数え、子を選ぶ。入れて 64 個になった節点は半分に割り、除いて空になった節点は親から外す (隣と併せることはしない)。
// Counted = true なら k 番目と x 以下の個数も引ける。内部の節点は子ごとの部分木の個数と、子を 8 個ずつに分けた組の累積を
// 持ち、k 番目は組の累積を 1 回比べて組を選んでから組の中の累積をレジスタで作って選び、x 以下の個数は組の累積に組の中の
// 手前の和を足す。どちらも分岐を使わない。入れるときと除くときは、子の個数 1 つと組の累積 8 個を足し引きする。
// Counted = false はこれらを持たず、その分だけ入れる・除くが安い。procon-judge の yosupo-ordered-set で書き比べた
// bptree_c8b を写したもの。
namespace ordered_set_internal {
constexpr int B= 64, FILL= 48;
struct alignas(64) Leaf {
 int key[B];
 int n, prv, nxt;
};
struct alignas(64) InnerN {
 int mx[B], ch[B];
 int n;
};
struct alignas(64) InnerC {
 int mx[B], ch[B];
 int cnt[B];  // 子ごとの部分木の個数。n より後ろは 0
 int gc[9];   // gc[g] は cnt[0, 8g) の和
 int n;
};
// 昇順の a[0, n) のうち x 未満のものの個数。
inline int count_lt(const int* a, int n, int x) {
 const __m256i vx= _mm256_set1_epi32(x);
 int c= 0;
 for(int i= 0; i < n; i+= 8) {
  unsigned m= _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpgt_epi32(vx, _mm256_loadu_si256(reinterpret_cast<const __m256i*>(a + i)))));
  if(n - i < 8) m&= (1u << (n - i)) - 1;
  c+= __builtin_popcount(m);
 }
 return c;
}
// 昇順の a[0, n) のうち x 以下のものの個数。
inline int count_le(const int* a, int n, int x) {
 const __m256i vx= _mm256_set1_epi32(x);
 int c= 0;
 for(int i= 0; i < n; i+= 8) {
  unsigned m= ~_mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpgt_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(a + i)), vx))) & 0xff;
  if(n - i < 8) m&= (1u << (n - i)) - 1;
  c+= __builtin_popcount(m);
 }
 return c;
}
inline int sum8(const int* a) {
 int s= 0;
 for(int i= 0; i < 8; ++i) s+= a[i];
 return s;
}
// 子の個数を見直して組の累積を作る。n より後ろの個数は 0 にする。
inline void regroup(InnerC& t) {
 std::memset(t.cnt + t.n, 0, (B - t.n) * sizeof(int));
 t.gc[0]= 0;
 for(int g= 0, c= 0; g < 8; ++g) t.gc[g + 1]= c+= sum8(t.cnt + g * 8);
}
// 子 i を含む組から後ろの累積 gc[i / 8 + 1, 9) に d を足す。
inline void add_group(InnerC& t, int i, int d) {
 const __m256i idx= _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7), from= _mm256_set1_epi32(i / 8 - 1);
 __m256i* p= reinterpret_cast<__m256i*>(t.gc + 1);
 _mm256_storeu_si256(p, _mm256_add_epi32(_mm256_loadu_si256(p), _mm256_and_si256(_mm256_cmpgt_epi32(idx, from), _mm256_set1_epi32(d))));
}
// p[0, r) の和 (r は 0 以上 7 以下)。
inline int group_prefix(const int* p, int r) {
 const __m256i idx= _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
 const __m256i v= _mm256_and_si256(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(p)), _mm256_cmpgt_epi32(_mm256_set1_epi32(r), idx));
 __m128i s= _mm_add_epi32(_mm256_castsi256_si128(v), _mm256_extracti128_si256(v, 1));
 s= _mm_add_epi32(s, _mm_shuffle_epi32(s, 0x4e));
 s= _mm_add_epi32(s, _mm_shuffle_epi32(s, 0xb1));
 return _mm_cvtsi128_si32(s);
}
// p[0, 8) の先頭からの和が k を超える最初の位置 j を返し、k から p[0, j) の和を引く。p の和は k より大きい。
inline int group_find(const int* p, int& k) {
 const __m256i x= _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p));
 __m256i v= _mm256_add_epi32(x, _mm256_slli_si256(x, 4));
 v= _mm256_add_epi32(v, _mm256_slli_si256(v, 8));
 v= _mm256_add_epi32(v, _mm256_blend_epi32(_mm256_setzero_si256(), _mm256_permutevar8x32_epi32(v, _mm256_set1_epi32(3)), 0xf0));
 const int j= __builtin_ctz(_mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpgt_epi32(v, _mm256_set1_epi32(k)))));
 k-= _mm_cvtsi128_si32(_mm256_castsi256_si128(_mm256_permutevar8x32_epi32(_mm256_sub_epi32(v, x), _mm256_set1_epi32(j))));
 return j;
}
}
template <bool Counted= true> class OrderedSet {
 using Leaf= ordered_set_internal::Leaf;
 using Inner= std::conditional_t<Counted, ordered_set_internal::InnerC, ordered_set_internal::InnerN>;
 static constexpr int B= ordered_set_internal::B, FILL= ordered_set_internal::FILL;
 std::vector<Leaf> lf;
 std::vector<Inner> in;
 std::vector<int> flf, fin;  // 空いた葉と節点
 int root, h, total;         // h は根の高さ (0 なら根が葉)
 int path[16], pos[16];      // 段 d の通った節点と、その中で選んだ子
 int new_leaf() {
  if(!flf.empty()) {
   const int v= flf.back();
   flf.pop_back();
   return v;
  }
  lf.emplace_back();
  return int(lf.size()) - 1;
 }
 int new_inner() {
  if(!fin.empty()) {
   const int v= fin.back();
   fin.pop_back();
   return v;
  }
  in.emplace_back();
  return int(in.size()) - 1;
 }
 static void regroup(Inner& t) {
  if constexpr(Counted) ordered_set_internal::regroup(t);
 }
 // 段 d - 1 の節点 l を割って r を作ったので、段 d の親に r を足す。lc と rc は 2 つの部分木の個数 (Counted のときだけ使う)。
 void push_up(int d, int l, int lc, int lm, int r, int rc, int rm) {
  if(d > h) {
   const int v= new_inner();
   Inner& t= in[v];
   t.n= 2, t.ch[0]= l, t.ch[1]= r, t.mx[0]= lm, t.mx[1]= rm;
   if constexpr(Counted) t.cnt[0]= lc, t.cnt[1]= rc;
   regroup(t);
   root= v, ++h;
   return;
  }
  const int v= path[d], i= pos[d];
  {
   Inner& t= in[v];
   const int k= t.n - i - 1;
   std::memmove(t.ch + i + 2, t.ch + i + 1, k * sizeof(int));
   std::memmove(t.mx + i + 2, t.mx + i + 1, k * sizeof(int));
   t.mx[i]= lm, t.ch[i + 1]= r, t.mx[i + 1]= rm;
   if constexpr(Counted) std::memmove(t.cnt + i + 2, t.cnt + i + 1, k * sizeof(int)), t.cnt[i]= lc, t.cnt[i + 1]= rc;
   if(++t.n < B) return regroup(t);
  }
  const int u= new_inner();
  Inner &T= in[v], &U= in[u];
  T.n= U.n= B / 2;
  std::memcpy(U.ch, T.ch + B / 2, B / 2 * sizeof(int));
  std::memcpy(U.mx, T.mx + B / 2, B / 2 * sizeof(int));
  int tc= 0, uc= 0;
  if constexpr(Counted) {
   std::memcpy(U.cnt, T.cnt + B / 2, B / 2 * sizeof(int));
   regroup(T), regroup(U);
   tc= T.gc[8], uc= U.gc[8];
  }
  push_up(d + 1, v, tc, T.mx[B / 2 - 1], u, uc, U.mx[B / 2 - 1]);
 }
 void split_leaf(int v) {
  const int u= new_leaf();
  Leaf &L= lf[v], &R= lf[u];
  R.n= L.n= B / 2;
  std::memcpy(R.key, L.key + B / 2, B / 2 * sizeof(int));
  R.prv= v, R.nxt= L.nxt, L.nxt= u;
  if(R.nxt >= 0) lf[R.nxt].prv= u;
  push_up(1, v, B / 2, L.key[B / 2 - 1], u, B / 2, R.key[B / 2 - 1]);
 }
 // 段 d の pos[d] 番目の子の部分木の最大が m に変わったので、上へ直す。
 void fix_max(int d, int m) {
  for(; d <= h; ++d) {
   Inner& t= in[path[d]];
   t.mx[pos[d]]= m;
   if(pos[d] != t.n - 1) return;
  }
 }
 // 段 d の pos[d] 番目の子 (空になった) を外す。
 void remove_child(int d) {
  const int v= path[d], i= pos[d];
  Inner& t= in[v];
  const int k= t.n - i - 1;
  std::memmove(t.ch + i, t.ch + i + 1, k * sizeof(int));
  std::memmove(t.mx + i, t.mx + i + 1, k * sizeof(int));
  if constexpr(Counted) std::memmove(t.cnt + i, t.cnt + i + 1, k * sizeof(int));
  --t.n;
  regroup(t);
  if(t.n == 0) {
   fin.push_back(v);
   if(d < h) return remove_child(d + 1);
   root= new_leaf(), h= 0;
   lf[root].n= 0, lf[root].prv= lf[root].nxt= -1;
   return;
  }
  if(i == t.n) fix_max(d + 1, t.mx[t.n - 1]);
  while(h > 0 && in[root].n == 1) fin.push_back(root), root= in[root].ch[0], --h;
 }
 // x 未満 (LE なら x 以下) の要素の個数。
 template <bool LE> int rank(int x) const {
  static_assert(Counted, "rank needs OrderedSet<true>");
  int c= 0, v= root;
  for(int d= h; d > 0; --d) {
   const Inner& t= in[v];
   const int i= LE ? ordered_set_internal::count_le(t.mx, t.n, x) : ordered_set_internal::count_lt(t.mx, t.n, x);
   if(i == t.n) return c + t.gc[8];
   c+= t.gc[i / 8] + ordered_set_internal::group_prefix(t.cnt + i / 8 * 8, i % 8), v= t.ch[i];
  }
  const Leaf& L= lf[v];
  return c + (LE ? ordered_set_internal::count_le(L.key, L.n, x) : ordered_set_internal::count_lt(L.key, L.n, x));
 }
public:
 OrderedSet(): OrderedSet(std::vector<int>()) {}
 // a の要素からなる集合を作る。a は狭義単調増加。
 explicit OrderedSet(const std::vector<int>& a): h(0), total(a.size()) {
  const int n= a.size();
  lf.reserve(n / FILL + 2), in.reserve(n / FILL / FILL + 8);
  std::vector<int> ids, mxs, cnts;
  for(int i= 0; i < n || i == 0; i+= FILL) {
   const int v= new_leaf(), m= std::min(FILL, n - i);
   Leaf& L= lf[v];
   L.n= m, L.prv= ids.empty() ? -1 : ids.back(), L.nxt= -1;
   if(m > 0) std::memcpy(L.key, a.data() + i, m * sizeof(int));
   if(!ids.empty()) lf[ids.back()].nxt= v;
   ids.push_back(v), mxs.push_back(m ? L.key[m - 1] : 0), cnts.push_back(m);
  }
  while(ids.size() > 1) {
   std::vector<int> ids2, mxs2, cnts2;
   for(std::size_t i= 0; i < ids.size(); i+= FILL) {
    const int v= new_inner(), m= std::min<int>(FILL, ids.size() - i);
    Inner& t= in[v];
    t.n= m;
    int c= 0;
    for(int j= 0; j < m; ++j) {
     t.ch[j]= ids[i + j], t.mx[j]= mxs[i + j], c+= cnts[i + j];
     if constexpr(Counted) t.cnt[j]= cnts[i + j];
    }
    regroup(t);
    ids2.push_back(v), mxs2.push_back(t.mx[m - 1]), cnts2.push_back(c);
   }
   ids.swap(ids2), mxs.swap(mxs2), cnts.swap(cnts2), ++h;
  }
  root= ids[0];
 }
 int size() const { return total; }
 bool empty() const { return total == 0; }
 bool contains(int x) const {
  int v= root;
  for(int d= h; d > 0; --d) {
   const Inner& t= in[v];
   const int i= ordered_set_internal::count_lt(t.mx, t.n, x);
   if(i == t.n) return false;
   v= t.ch[i];
  }
  const Leaf& L= lf[v];
  const int j= ordered_set_internal::count_lt(L.key, L.n, x);
  return j < L.n && L.key[j] == x;
 }
 // x を入れる。入ったら true、既にあれば false。
 bool insert(int x) {
  int v= root;
  for(int d= h; d > 0; --d) {
   const Inner& t= in[v];
   int i= ordered_set_internal::count_lt(t.mx, t.n, x);
   if(i == t.n) i= t.n - 1;
   path[d]= v, pos[d]= i, v= t.ch[i];
  }
  Leaf& L= lf[v];
  const int j= ordered_set_internal::count_lt(L.key, L.n, x);
  if(j < L.n && L.key[j] == x) return false;
  std::memmove(L.key + j + 1, L.key + j, (L.n - j) * sizeof(int));
  L.key[j]= x, ++L.n, ++total;
  for(int d= 1; d <= h; ++d) {
   Inner& t= in[path[d]];
   if constexpr(Counted) ++t.cnt[pos[d]], ordered_set_internal::add_group(t, pos[d], 1);
   if(t.mx[pos[d]] < x) t.mx[pos[d]]= x;
  }
  if(L.n == B) split_leaf(v);
  return true;
 }
 // x を除く。除いたら true、無ければ false。
 bool erase(int x) {
  int v= root;
  for(int d= h; d > 0; --d) {
   const Inner& t= in[v];
   const int i= ordered_set_internal::count_lt(t.mx, t.n, x);
   if(i == t.n) return false;
   path[d]= v, pos[d]= i, v= t.ch[i];
  }
  Leaf& L= lf[v];
  const int j= ordered_set_internal::count_lt(L.key, L.n, x);
  if(j == L.n || L.key[j] != x) return false;
  std::memmove(L.key + j, L.key + j + 1, (L.n - j - 1) * sizeof(int));
  --L.n, --total;
  if constexpr(Counted)
   for(int d= 1; d <= h; ++d) {
    Inner& t= in[path[d]];
    --t.cnt[pos[d]], ordered_set_internal::add_group(t, pos[d], -1);
   }
  if(L.n > 0) {
   if(j == L.n) fix_max(1, L.key[L.n - 1]);
   return true;
  }
  if(h == 0) return true;
  if(L.prv >= 0) lf[L.prv].nxt= L.nxt;
  if(L.nxt >= 0) lf[L.nxt].prv= L.prv;
  flf.push_back(v);
  remove_child(1);
  return true;
 }
 // x 以下で最大の要素。
 std::optional<int> prev(int x) const {
  int v= root;
  for(int d= h; d > 0; --d) {
   const Inner& t= in[v];
   const int i= ordered_set_internal::count_le(t.mx, t.n, x);
   if(i == t.n) return t.mx[t.n - 1];
   v= t.ch[i];
  }
  const Leaf& L= lf[v];
  if(const int j= ordered_set_internal::count_le(L.key, L.n, x)) return L.key[j - 1];
  if(L.prv >= 0) return lf[L.prv].key[lf[L.prv].n - 1];
  return std::nullopt;
 }
 // x 以上で最小の要素。
 std::optional<int> next(int x) const {
  int v= root;
  for(int d= h; d > 0; --d) {
   const Inner& t= in[v];
   const int i= ordered_set_internal::count_lt(t.mx, t.n, x);
   if(i == t.n) return std::nullopt;
   v= t.ch[i];
  }
  const Leaf& L= lf[v];
  const int j= ordered_set_internal::count_lt(L.key, L.n, x);
  if(j == L.n) return std::nullopt;
  return L.key[j];
 }
 // 小さい方から k 番目 (0 始まり) の要素。Counted のときだけ使える。
 std::optional<int> kth(int k) const {
  static_assert(Counted, "kth needs OrderedSet<true>");
  if(k < 0 || k >= total) return std::nullopt;
  int v= root;
  for(int d= h; d > 0; --d) {
   const Inner& t= in[v];
   const int g= ordered_set_internal::count_le(t.gc + 1, 8, k);
   k-= t.gc[g];
   v= t.ch[g * 8 + ordered_set_internal::group_find(t.cnt + g * 8, k)];
  }
  return lf[v].key[k];
 }
 // x 未満の要素の個数。Counted のときだけ使える。
 int count_lt(int x) const { return rank<false>(x); }
 // x 以下の要素の個数。Counted のときだけ使える。
 int count_le(int x) const { return rank<true>(x); }
};
