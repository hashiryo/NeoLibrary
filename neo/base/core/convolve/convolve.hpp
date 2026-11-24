#pragma once
#include <algorithm>
#include <cassert>
#include <vector>
// 畳み込みの共通の口。積を求める核 (static 関数 convolve_range を持つ構造体) をテンプレート引数 Conv で受け、全体の積、範囲の積、
// middle product を出す。このヘッダは核を include しない。核は使う側が include する (neo/base/core/convolve/naive.hpp の Naive など)。
//   convolve<Conv>(a, b): 積の全体。
//   convolve_range<Conv>(a, b, lo, hi): 積の lo 次から hi - 1 次の係数。先頭 len 項は (0, len)。
//   middle_product<Conv>(a, b): c[i] = Σ_j a[i + j] b[j] (0 <= i <= |a| - |b|)。
// Conv を省くと、係数の型 K の既定の核になる。既定は、核のヘッダが宣言する default_kernel(const K*) の戻り値の型で決まる。
// include した核のうち K を受け持つものが 1 つならそれになる。無いか 2 つ以上あるとコンパイルエラーになるので、Conv を書いて選ぶ。
// 核が持つ関数は
//   static void convolve_range(const K* a, int n, const K* b, int m, int lo, int hi, K* out);
// で、out[k - lo] = Σ_{i+j=k} a[i] b[j] (lo <= k < hi) を書く。n >= 1、m >= 1、0 <= lo < hi <= n + m - 1 で呼ばれ、out は a、b と
// 重ならない。a と b が同じ配列 (a == b かつ n == m) のこともある。
namespace convolve_internal {
struct Auto {};
template <class Conv, class K> struct pick {
 using type= Conv;
};
template <class K> struct pick<Auto, K> {
 // ここでエラーになるときは、K を受け持つ核を include していないか、既定の核が 1 つに決まらない。
 using type= decltype(default_kernel(static_cast<const K*>(nullptr)));
};
}
// Conv が省かれていれば K の既定の核、書かれていれば Conv。
template <class Conv, class K> using kernel_t= typename convolve_internal::pick<Conv, K>::type;
// out[k - lo] = Σ_{i+j=k} a[i] b[j] (lo <= k < hi)。積の次数を超える k には 0 を書く。0 <= lo で、n や m は 0 でもよい。
// out は a、b と重ならないこと。範囲に届かない項を落としてから核を呼ぶ。
template <class Conv= convolve_internal::Auto, class K> void convolve_range(const K* a, int n, const K* b, int m, int lo, int hi, K* out) {
 if(lo >= hi) return;
 const int top= n > 0 && m > 0 ? n + m - 1 : 0;  // 積の項の数
 if(hi > top) {
  std::fill(out + (std::max(lo, top) - lo), out + (hi - lo), K());
  hi= top;
  if(lo >= hi) return;
 }
 n= std::min(n, hi), m= std::min(m, hi);  // i >= hi の a[i] は k >= hi にしか入らない
 const int sa= std::max(0, lo - (m - 1));  // i < sa の a[i] は k < lo にしか入らない
 a+= sa, n-= sa, lo-= sa, hi-= sa;
 const int sb= std::max(0, lo - (n - 1));
 b+= sb, m-= sb, lo-= sb, hi-= sb;
 n= std::min(n, hi), m= std::min(m, hi);
 kernel_t<Conv, K>::convolve_range(a, n, b, m, lo, hi, out);
}
// 積の lo 次から hi - 1 次の係数 (0 <= lo <= hi)。積の次数を超える部分は 0。
template <class Conv= convolve_internal::Auto, class K> std::vector<K> convolve_range(const std::vector<K>& a, const std::vector<K>& b, int lo, int hi) {
 assert(0 <= lo && lo <= hi);
 std::vector<K> c(hi - lo);
 convolve_range<Conv>(a.data(), int(a.size()), b.data(), int(b.size()), lo, hi, c.data());
 return c;
}
// 積の全体。どちらかが空なら空を返す。
template <class Conv= convolve_internal::Auto, class K> std::vector<K> convolve(const std::vector<K>& a, const std::vector<K>& b) {
 if(a.empty() || b.empty()) return {};
 return convolve_range<Conv>(a, b, 0, int(a.size() + b.size()) - 1);
}
// out[i] = Σ_j a[i + j] b[j] (0 <= i <= n - m)。n >= m >= 1 で、out は a、b と重ならないこと。
// 積の転置で、b を反転した積の m - 1 次から n - 1 次と同じ。
template <class Conv= convolve_internal::Auto, class K> void middle_product(const K* a, int n, const K* b, int m, K* out) {
 std::vector<K> r(b, b + m);
 std::reverse(r.begin(), r.end());
 convolve_range<Conv>(a, n, static_cast<const K*>(r.data()), m, m - 1, n, out);
}
// c[i] = Σ_j a[i + j] b[j] (0 <= i <= |a| - |b|)。b が空か、a が b より短ければ空を返す。
template <class Conv= convolve_internal::Auto, class K> std::vector<K> middle_product(const std::vector<K>& a, const std::vector<K>& b) {
 const int n= int(a.size()), m= int(b.size());
 if(m == 0 || n < m) return {};
 std::vector<K> c(n - m + 1);
 middle_product<Conv>(a.data(), n, b.data(), m, c.data());
 return c;
}
