#pragma once
#include <algorithm>
#include <cstdint>
#include <memory>
#include <numeric>
#include <string>
#include <vector>
// SA-IS による接尾辞配列。1 回目の induced sorting は libsais (Ilya Grebnov) の組み立て方を、文字の種類が少ない段で使い、
// 最後の induced sorting は sais-lite (Yuta Mori) の符号の手で回す。procon-judge の yosupo-suffixarray で書き比べた sais_lr13 と同じ。
namespace suffix_array_internal {
constexpr int MARK= int(0x80000000u), POS= 0x7fffffff;
// induced sorting の 2 回の走査 (L 型を前から、S 型を後ろから)。cnt[c] はバケット c の先頭、cnt[c + 1] は末尾の次。置く値の符号に
// 「この接尾辞の 1 つ前を、今の走査で induce するか」を持たせ、L 型と S 型は隣の 2 文字から決める (sais-lite の手)。First は 1 回目
// (LMS の部分文字列を並べる) で、読んだ値を 0 に戻し、LMS だけを負の値で残す。Table は最後の induced sorting で、バケットの位置を
// レジスタに置かず、置くたびに表を読み書きする (置くたびにバケットが替わる入力で速い)。
template <bool First, bool Table, class C> inline void induce_scan(const C* s, int n, int K, const int* cnt, int* sa, int* bkt) {
 for(int c= 0; c < K; ++c) bkt[c]= cnt[c];
 if constexpr(Table) {
  {
   const int j= n - 1, c= int(s[j]);  // 番兵から最後の文字の接尾辞を induce する
   sa[bkt[c]++]= j ^ -int((j > 0) & (int(s[j - (j > 0)]) < c));
  }
  for(int i= 0; i < n; ++i) {
   int j= sa[i];
   sa[i]= ~j;
   if(j > 0) {
    --j;
    const int c= int(s[j]);
    sa[bkt[c]++]= j ^ -int((j > 0) & (int(s[j - (j > 0)]) < c));
   }
  }
  for(int c= 0; c < K; ++c) bkt[c]= cnt[c + 1];
  for(int i= n - 1; i >= 0; --i) {
   int j= sa[i];
   if(j > 0) {
    --j;
    const int c= int(s[j]);
    sa[--bkt[c]]= j ^ -int((j == 0) | (int(s[j - (j > 0)]) > c));
   } else {
    sa[i]= ~j;
   }
  }
  return;
 }
 int c1= int(s[n - 1]);
 int* b= sa + bkt[c1];
 {
  const int j= n - 1;  // 番兵から最後の文字の接尾辞を induce する
  *b++= j ^ -int((j > 0) & (int(s[j - (j > 0)]) < c1));
 }
 for(int i= 0; i < n; ++i) {
  int j= sa[i];
  sa[i]= First ? ~j & (j >> 31) : ~j;
  if(j > 0) {
   --j;
   const int c0= int(s[j]);
   if(c0 != c1) bkt[c1]= int(b - sa), b= sa + bkt[c1= c0];
   *b++= j ^ -int((j > 0) & (int(s[j - (j > 0)]) < c1));  // ~j は j ^ -1。j = 0 のときは s[0] を読んで捨てる
  }
 }
 for(int c= 0; c < K; ++c) bkt[c]= cnt[c + 1];
 c1= 0;
 b= sa + bkt[0];
 for(int i= n - 1; i >= 0; --i) {
  int j= sa[i];
  if(j > 0) {
   if constexpr(First) sa[i]= 0;
   --j;
   const int c0= int(s[j]);
   if(c0 != c1) bkt[c1]= int(b - sa), b= sa + bkt[c1= c0];
   const int lt= int(s[j - (j > 0)]) > c1;  // j - 1 が L 型
   *--b= j ^ -(First ? int(j > 0) & lt : int(j == 0) | lt);
  } else if constexpr(!First) {
   sa[i]= ~j;
  }
 }
}
// sa[0, m) に並んだ LMS を、その場でバケットの末尾へ移して最後の induced sorting をする (k 番目の LMS の行き先は k 以上)。
template <bool Table, class C> inline void induce_inplace(const C* s, int n, int K, const int* cnt, int m, int* sa, int* bkt) {
 std::fill(sa + m, sa + n, 0);
 for(int c= 0; c < K; ++c) bkt[c]= cnt[c + 1];
 for(int k= m - 1; k >= 0; --k) {
  const int p= sa[k];
  sa[k]= 0;
  sa[--bkt[s[p]]]= p;
 }
 induce_scan<false, Table>(s, n, K, cnt, sa, bkt);
}
// LMS の位置を後ろからの走査で見つけて bits に印を付け、その数を返す。Four なら位置 p ごとに
// q[4 * s[p] + 2 * (p が L 型) + (p - 1 が L 型)] を、そうでなければ q[s[p] + 1] を数える。
template <bool Four, class C> inline int gather_lms(const C* s, int n, uint64_t* bits, int* q) {
 int m= 0, c0= int(s[n - 1]), f0= 1;  // 末尾の番兵は最も小さいので、最後の文字は L 型
 uint64_t w= 0;
 for(int p= n - 1; p > 0; --p) {  // 位置 p の種類 f0 と、p - 1 の種類 f1
  const int c1= int(s[p - 1]), f1= c1 > c0 - f0, is= f1 & (f0 ^ 1);
  ++q[Four ? 4 * c0 + 2 * f0 + f1 : c0 + 1];
  w|= uint64_t(is) << (p & 63), m+= is;
  if((p & 63) == 0) bits[p >> 6]= w, w= 0;
  c0= c1, f0= f1;
 }
 ++q[Four ? 4 * c0 + 2 * f0 : c0 + 1];
 bits[0]= w;  // 位置 0 は LMS にならない
 return m;
}
// 印を付けた LMS の位置を、位置の順に dst に書く。
inline void expand_lms(const uint64_t* bits, int n, int* dst) {
 for(int b= 0; b <= (n - 1) >> 6; ++b)
  for(uint64_t w= bits[b]; w; w&= w - 1) *dst++= 64 * b + __builtin_ctzll(w);
}
// sort_lms_lr の区画の 1 回の走査 (Up なら前から読んで前から置く)。読んだ要素の 1 つ前を置き先 ptr[2c + t] (t は置く要素の 1 つ前が
// L 型なら 1) に置き、元の要素の組 (d、印を読むたびに増える) が同じ置き先の前の要素と違えば印を付ける。直前の置き先の位置と組は
// レジスタに置く。
template <bool Up, class C> inline int scan_lr(const C* s, int* sa, int from, int to, int* ptr, int* D, int d, int lv) {
 int *lb= sa + ptr[lv], ld= D[lv];
 for(int i= from; Up ? i < to : i > to; Up ? ++i : --i) {
  const int x= sa[i];
  d+= x < 0;
  const int p= (x & POS) - 1, c= int(s[p]), v= 2 * c + (Up ? int(s[p - 1]) >= c : int(s[p - 1]) > c);
  if(v != lv) ptr[lv]= int(lb - sa), D[lv]= ld, lv= v, lb= sa + ptr[v], ld= D[v];
  const int y= p | int(unsigned(ld != d) << 31);
  if constexpr(Up) *lb++= y;
  else *--lb= y;
  ld= d;
 }
 ptr[lv]= int(lb - sa), D[lv]= ld;
 return d;
}
// 文字の種類が少ない段で、LMS の部分文字列を並べて名前を付ける (libsais の組み立て方)。位置を (L 型か、1 つ前が L 型か) の 4 種類に
// 分けて文字ごとに数え、前から読む走査で読む要素を左の区画に、後ろから読む走査で読む要素を右の区画に隙間なく置くので、走査に分岐が
// 要らない。いちばん左の LMS とそれより左の位置は数から外す。種はバケットごとにいちばん前のものに印を付ける (名前を正確にする)。
// 並んだ LMS を sa[0, m) に、位置 p の名前を印を付けて sa[m + p / 2] に置き、名前の数を返す。
template <class C> inline int sort_lms_lr(const C* s, int n, int K, int* q, const int* lms, int m, int* sa) {
 const int first= lms[0];
 {
  int c0= int(s[first]), f0= 0;
  for(int p= first; p > 0; --p) {
   const int c1= int(s[p - 1]), f1= c1 > c0 - f0;
   --q[4 * c0 + 2 * f0 + f1];
   c0= c1, f0= f1;
  }
  --q[4 * c0 + 2 * f0];
 }
 std::vector<int> ptr(4 * K);
 int *lp= ptr.data(), *rp= lp + 2 * K, left= 0;
 for(int c= 0; c < K; ++c) lp[2 * c + 1]= left, left+= q[4 * c + 3] + q[4 * c + 1], rp[2 * c + 1]= left;
 for(int k= m - 1; k > 0; --k) sa[--rp[2 * int(s[lms[k]]) + 1]]= lms[k];
 for(int c= 0; c < K; ++c)
  if(q[4 * c + 1] > 0) sa[rp[2 * c + 1]]|= MARK;
 for(int c= 0, r= left + 1, l= 0, cf= int(s[first]); c < K; ++c) {
  lp[2 * c]= r, r+= q[4 * c] + q[4 * c + 2], rp[2 * c]= r;
  l+= q[4 * c + 1] + (c == cf), rp[2 * c + 1]= l;  // first も LMS として置く
 }
 int* D= q;  // 数はもう使わないので、組の番号の表に q の領域を使う
 std::fill(D, D + 2 * K, 0);
 const int c= int(s[n - 1]), lv= 2 * c + (int(s[n - 2]) >= c);  // 番兵から最後の文字の接尾辞を induce する
 sa[lp[lv]++]= (n - 1) | MARK, D[lv]= 1;
 const int d= scan_lr<true>(s, sa, 0, left, lp, D, 1, lv);
 for(int c= 0; c < K; ++c) {
  int f= MARK;
  for(int i= lp[2 * c] - 1, lo= c ? rp[2 * c - 2] : left + 1; i >= lo; --i) {
   const int x= sa[i], y= (x & MARK) ^ f;
   f^= y, sa[i]= x ^ y;
  }
 }
 scan_lr<false>(s, sa, n - first - 1, left, rp, D, d, 0);
 std::fill(sa + m, sa + m + (n >> 1), 0);
 int names= 0;
 for(int k= 0; k < m; ++k) {
  const int x= sa[k];
  sa[m + ((x & POS) >> 1)]= names | MARK, names+= x < 0;
 }
 return names;
}
// 文字の種類が多い段 (K > n / 32) で、LMS の部分文字列を並べて名前を付ける。4 種類ずつの表は大きすぎるので、表 1 本の induced
// sorting で並べ、長さと文字の並びで比べる。出し方は sort_lms_lr と同じ。
template <class C> inline int sort_lms_wide(const C* s, int n, int K, const int* cnt, int* bkt, const int* lms, int m, int* sa) {
 std::fill(sa, sa + n, 0);
 for(int c= 0; c < K; ++c) bkt[c]= cnt[c + 1];
 for(int k= m - 1; k >= 0; --k) sa[--bkt[s[lms[k]]]]= lms[k];  // LMS の 1 つ前は L 型なので、正のまま置く
 induce_scan<true, false>(s, n, K, cnt, sa, bkt);
 for(int i= 0, k= 0; i < n; ++i) {  // 並んだ順に LMS の位置 (負の値で残っている) を sa[0, m) に詰める
  const int v= sa[i];
  sa[k]= ~v, k+= v < 0;
 }
 std::fill(sa + m, sa + m + (n >> 1), 0);
 for(int j= 0; j + 1 < m; ++j) sa[m + (lms[j] >> 1)]= lms[j + 1] - lms[j] + 1;
 int names= 0;
 for(int k= 0, prev= 0, prev_len= 0; k < m; ++k) {
  const int p= sa[k], l= sa[m + (p >> 1)];
  int r= 0;
  if(l != 0 && l == prev_len)
   while(r < l && s[p + r] == s[prev + r]) ++r;
  names+= r < l || l == 0;
  sa[m + (p >> 1)]= (names - 1) | MARK, prev= p, prev_len= l;
 }
 return names;
}
// s[0, n) (値は 0 以上 K 未満) の接尾辞配列を sa[0, n) に書く。末尾には、どの文字より小さい番兵があるものとして扱う。
template <class C> void sa_is(const C* s, int n, int K, int* sa) {
 if(n < 2) {  // n が 2 以上だと分かる形にしておくと、gcc が std::fill の長さを負と疑わない
  if(n == 1) sa[0]= 0;
  return;
 }
 const bool wide= K > n / 32;
 std::vector<int> q(wide ? 0 : 4 * K, 0), cnt(K + 1, 0), bkt(K);
 std::unique_ptr<uint64_t[]> bits(new uint64_t[(n + 63) / 64]);
 const int m= wide ? gather_lms<false>(s, n, bits.get(), cnt.data()) : gather_lms<true>(s, n, bits.get(), q.data());
 for(int c= 0; c < K; ++c) cnt[c + 1]= wide ? cnt[c + 1] + cnt[c] : cnt[c] + q[4 * c] + q[4 * c + 1] + q[4 * c + 2] + q[4 * c + 3];
 if(m == 0) {  // すべて L 型なら、番兵からの induce だけで並ぶ
  induce_inplace<false>(s, n, K, cnt.data(), 0, sa, bkt.data());
  return;
 }
 int names;
 if(wide) {
  std::unique_ptr<int[]> lms(new int[m]);
  expand_lms(bits.get(), n, lms.get());
  names= sort_lms_wide(s, n, K, cnt.data(), bkt.data(), lms.get(), m, sa);
 } else {
  expand_lms(bits.get(), n, sa + n - m);
  names= sort_lms_lr(s, n, K, q.data(), sa + n - m, m, sa);
 }
 if(names < m) {
  for(int i= m + (n >> 1) - 1, l= n; i >= m; --i) {
   const int x= sa[i];
   sa[l - 1]= x & POS, l-= x < 0;
  }
  sa_is<int>(sa + n - m, m, names, sa);
  expand_lms(bits.get(), n, sa + n - m);
  for(int k= 0; k < m; ++k) sa[k]= sa[n - m + sa[k]];
 } else {
  for(int k= 0; k < m; ++k) sa[k]&= POS;
 }
 // 最後の induced sorting の形を、並んだ LMS の隣どうしで 1 つ前の文字が違う割合 (置くたびにバケットが替わる割合の見積もり) で選ぶ。
 int switched= 0;
 for(int k= 0; k < std::min(m - 1, 4096); ++k) switched+= s[sa[k] - 1] != s[sa[k + 1] - 1];
 if(2 * switched > std::min(m - 1, 4096)) induce_inplace<true>(s, n, K, cnt.data(), m, sa, bkt.data());
 else induce_inplace<false>(s, n, K, cnt.data(), m, sa, bkt.data());
}
}
// 文字列 s の接尾辞配列 (辞書順で i 番目の接尾辞の開始位置) を返す。文字は unsigned char として比べる。
inline std::vector<int> suffix_array(const std::string& s) {
 std::vector<int> sa(s.size());
 suffix_array_internal::sa_is(reinterpret_cast<const unsigned char*>(s.data()), int(s.size()), 256, sa.data());
 return sa;
}
// 値が 0 以上 K 未満の列 s の接尾辞配列を返す。
inline std::vector<int> suffix_array(const std::vector<int>& s, int K) {
 std::vector<int> sa(s.size());
 suffix_array_internal::sa_is(s.data(), int(s.size()), K, sa.data());
 return sa;
}
// 比べられる値の列 s の接尾辞配列を返す。値を座標圧縮してから解く。
template <class T> std::vector<int> suffix_array(const std::vector<T>& s) {
 const int n= s.size();
 std::vector<int> idx(n), t(n);
 std::iota(idx.begin(), idx.end(), 0);
 std::sort(idx.begin(), idx.end(), [&](int a, int b) { return s[a] < s[b]; });
 int K= 0;
 for(int i= 0; i < n; ++i) t[idx[i]]= K+= i && s[idx[i - 1]] < s[idx[i]];
 return suffix_array(t, K + 1);
}
