#pragma once
// relaxed (オンライン) の畳み込み。f と g の係数を 1 つずつ受け取り、積 f g の係数を次数の順に返す。計算量は O(n log^2 n)
// (n は next を呼んだ回数)。核はテンプレート引数 Conv で受け、convolve.hpp の kernel_t で既定が決まる
// (core/convolve/ntt.hpp などを include すると、その核が既定になる)。内部では核の convolve_range だけを呼ぶので、
// dft や idft を持たない核 (Naive など) でも動く。
// クラステンプレートは既定の引数のあとに既定の無い引数を置けないので、係数の型 K を先に書く。核は 2 つ目の引数。
//   relaxed::SemiRelaxedConvolution<K, Conv= Auto>: g を構築で全部受け取り (オフライン)、f だけを next で 1 つずつ受ける。
//   relaxed::Convolution<K, Conv= Auto>: f と g の両方を next で 1 つずつ受ける。
// どちらも、まだ渡していない次数の係数は 0 として扱う。next は何回でも呼べて、i 回目の呼び出し (0 始まり) が i 次の
// 係数を返す。長さ 0 は next を 1 回も呼ばないこと (SemiRelaxedConvolution の g は空でもよい)。長さ 1、2 冪でない長さ、
// f と g の長さが違う場合も、短い側を 0 で埋めた積の係数になる。返り値は「ここまでに渡した係数だけを持つ多項式の積」の
// 係数で、あとから渡す係数は前の次数には影響しない。
// 使い方:
//   ZMod<998244353> の g がオフラインで、f を先頭から与える:
//     relaxed::SemiRelaxedConvolution<Z> s(g);
//     for (int i= 0; i < n; ++i) c[i]= s.next(f[i]);
//   f と g の両方がオンライン:
//     relaxed::Convolution<Z> s;
//     for (int i= 0; i < n; ++i) c[i]= s.next(f[i], g[i]);
// 分割は hly1204 さんの記事 (algo-notes の sources/blogs/一种简单的在线卷积理解.md) の平方分割を、再帰ではなく next の
// 反復に落としたもの。時間 t (next を呼んだ回数) で 2^k | t となる大きさ 2^k の f のブロック [t - 2^k, t) と、g の
// ブロック [2^k, 2^{k+1}) の積を 1 回の畳み込みで求めて、まだ返していない次数に足す。記事の P_{0,0} が再帰、P_{0,1} と
// P_{1,0} が半オンライン、P_{1,1} が最後のオフラインの畳み込みに当たり、この反復はその 3 つを t ごとのブロックの対として
// まとめて回す形になる。g のブロックがオフラインなので、t の約数となるすべての 2^k で対を取れる。
// 両方がオンラインの Convolution は、f のブロック [t - 2^k, t) と g のブロック [2^k, 2^{k+1})、およびその転置 (g の
// ブロックを [t - 2^k, t)、f のブロックを [2^k, 2^{k+1}) にしたもの) を、t >= 2^{k+1} のときに取る。t = 2^{k+1} の
// ときは 2 つが同じ正方形になるので 1 回だけ取る。0 次の行と列 (f_0 g_j と f_i g_0) は、どの正方形にも入らないので、
// next の中で直接足す。
// どちらの構造体も作業の配列を持ち回るので、コピーすると状態ごと複製される (続きを別々に進められる)。ふつうは
// コピーしない。
#include <cassert>
#include <utility>
#include <vector>
#include "neo/base/core/convolve/convolve.hpp"

namespace relaxed {
namespace relaxed_internal {
// c を need 項以上に伸ばす (伸ばした分は K() で埋まる)。
template <class K> void ensure(std::vector<K>& c, int need) {
 if (int(c.size()) < need) c.resize(need);
}
// c[pos], c[pos + 1], ... に out の先頭 len 項を足す。
template <class K> void add(std::vector<K>& c, int pos, const K* out, int len) {
 ensure(c, pos + len);
 for (int i= 0; i < len; ++i) c[pos + i]+= out[i];
}
// 幅 len の畳み込み 1 回の結果を置く作業の配列を、足りなければ伸ばす (縮めない)。
template <class K> K* work(std::vector<K>& buf, int len) {
 if (int(buf.size()) < len) buf.resize(len);
 return buf.data();
}
}  // namespace relaxed_internal

// g がオフラインの半オンライン畳み込み。next(f_i) が (f g) の i 次の係数を返す。
template <class K, class Conv= convolve_internal::Auto> struct SemiRelaxedConvolution {
 // g はオフラインの係数 (長さは自由)。空でもよい。
 explicit SemiRelaxedConvolution(std::vector<K> g_): g(std::move(g_)) {}
 // f の次の係数を受け取り、(f g) のその次数の係数を返す。
 K next(K f_i) {
  const int i= int(f.size()), t= i + 1;
  f.push_back(f_i);
  relaxed_internal::ensure(c, t);
  if (!g.empty()) c[i]+= f[i] * g[0];
  // 2^k | t となるすべての大きさ 2^k で、f の [t - 2^k, t) と g の [2^k, 2^{k+1}) の積を足す。
  const int m= int(g.size());
  for (int L= 1; L <= t && (t & (L - 1)) == 0; L<<= 1) {
   if (m <= L) continue;                              // g の [L, 2L) が空
   const int len= L < m - L ? L : m - L;              // g の [L, L + len) (2L を超えない)
   const int outlen= L + len - 1;
   K* out= relaxed_internal::work(buf, outlen);
   convolve_range<Conv>(f.data() + (t - L), L, g.data() + L, len, 0, outlen, out);
   relaxed_internal::add(c, t, out, outlen);
  }
  return c[i];
 }
 private:
 std::vector<K> g;    // オフラインの係数
 std::vector<K> f;    // 受け取った f の係数
 std::vector<K> c;    // まだ返していない次数への寄与
 std::vector<K> buf;  // 畳み込み 1 回分の作業
};

// f と g の両方がオンラインの畳み込み。next(f_i, g_i) が (f g) の i 次の係数を返す。
template <class K, class Conv= convolve_internal::Auto> struct Convolution {
 // f と g の次の係数を受け取り、(f g) のその次数の係数を返す。
 K next(K f_i, K g_i) {
  const int i= int(f.size()), t= i + 1;
  f.push_back(f_i), g.push_back(g_i);
  relaxed_internal::ensure(c, t);
  if (t == 1) c[0]+= f[0] * g[0];
  else c[i]+= f[i] * g[0] + f[0] * g[i];  // 0 次の行と列
  // 2^k | t かつ t >= 2^{k+1} の大きさ 2^k で、正方形を足す。
  for (int L= 1; 2 * L <= t && (t & (L - 1)) == 0; L<<= 1) {
   const int outlen= 2 * L - 1;
   K* out= relaxed_internal::work(buf, outlen);
   if (t == 2 * L) {  // 2 つが同じ正方形になる対角の場合
    convolve_range<Conv>(f.data() + L, L, g.data() + L, L, 0, outlen, out);
    relaxed_internal::add(c, t, out, outlen);
   } else {
    convolve_range<Conv>(f.data() + (t - L), L, g.data() + L, L, 0, outlen, out);
    relaxed_internal::add(c, t, out, outlen);
    convolve_range<Conv>(g.data() + (t - L), L, f.data() + L, L, 0, outlen, out);
    relaxed_internal::add(c, t, out, outlen);
   }
  }
  return c[i];
 }
 private:
 std::vector<K> f, g;  // 受け取った係数
 std::vector<K> c;     // まだ返していない次数への寄与
 std::vector<K> buf;   // 畳み込み 1 回分の作業
};

// f の先頭 n 項の exp (f[0] = 0)。g = exp f は g' = f' g を満たすので、i 次までの係数から (i + 1) g_{i+1} が決まる。
// (k f_k) と g のオンライン畳み込み 1 本で、係数を 1 つずつ決めていく。O(n log^2 n)。値は fps::exp と同じ。
// 1 から n - 1 が K で逆元を持つこと (法が素数で n <= 法なら満たされる)。f は n 項より短くてもよい (0 で埋める)。
template <class Conv= convolve_internal::Auto, class K> std::vector<K> exp(const std::vector<K>& f, int n) {
 assert(n >= 1);
 assert(f.empty() || f[0] == K());
 std::vector<K> g(n);
 g[0]= K(1);
 if (n == 1) return g;
 Convolution<K, Conv> cv;
 for (int i= 0; i + 1 < n; ++i) {
  const K fi1= i + 1 < int(f.size()) ? f[i + 1] : K();
  g[i + 1]= cv.next(K(i + 1) * fi1, g[i]) / K(i + 1);
 }
 return g;
}

// f の先頭 n 項の逆数 g (f g = 1 mod x^n)。f[0] は逆元を持つこと。
// g_{i+1} = -(1/f_0) Σ_{k=1..i+1} f_k g_{i+1-k} なので、(f_{k+1}) と g のオンライン畳み込み 1 本で係数を決めていける。
// O(n log^2 n)。値は fps::inv と同じ。
template <class Conv= convolve_internal::Auto, class K> std::vector<K> inv(const std::vector<K>& f, int n) {
 assert(n >= 1);
 assert(!f.empty() && f[0] != K());
 std::vector<K> g(n);
 const K i0= f[0].inv();
 g[0]= i0;
 if (n == 1) return g;
 Convolution<K, Conv> cv;
 for (int i= 0; i + 1 < n; ++i) {
  const K fi1= i + 1 < int(f.size()) ? f[i + 1] : K();
  g[i + 1]= -cv.next(fi1, g[i]) * i0;
 }
 return g;
}
}  // namespace relaxed
