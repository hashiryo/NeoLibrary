#pragma once
// 多項式の商と余り (fast)。a = q b + r、deg r < deg b。K は体で、b の先頭の係数は 0 でないこと。
// 商は q = rev(a) / rev(b) の先頭 deg a - deg b + 1 項で、fps::inv の上で作る。余りは a - q b の deg b 未満の部分。
// 多項式は係数の vector (k 次の係数が [k]) で、末尾の 0 は取り除いて返す (0 は空の vector)。素朴な版は poly_divmod.hpp。
#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>
#include "neo/base/core/convolve/convolve.hpp"
#include "neo/base/core/poly/fps.hpp"

namespace poly_divmod_fast_internal {
template <class K> void trim(std::vector<K>& a) {
 while (!a.empty() && a.back() == K()) a.pop_back();
}
}  // namespace poly_divmod_fast_internal

template <class Conv= convolve_internal::Auto, class K>
std::pair<std::vector<K>, std::vector<K>> poly_divmod_fast(const std::vector<K>& a, const std::vector<K>& b) {
 using poly_divmod_fast_internal::trim;
 std::vector<K> aa= a, bb= b;
 trim(aa), trim(bb);
 assert(!bb.empty());
 const int n= int(aa.size()), m= int(bb.size());
 if (n < m) return {std::vector<K>(), aa};
 const int k= n - m + 1;  // 商の項の数 (先頭の係数は a の先頭 / b の先頭で 0 でない)
 // rev(a) の先頭 k 項と rev(b) (m 項)。q = rev(a) / rev(b) の先頭 k 項。
 std::vector<K> ar(k), br(m);
 for (int i= 0; i < k; ++i) ar[i]= aa[n - 1 - i];
 for (int i= 0; i < m; ++i) br[i]= bb[m - 1 - i];
 const std::vector<K> inv_br= fps::inv<Conv>(br, k);
 std::vector<K> qr(k);
 convolve_range<Conv>(ar.data(), k, inv_br.data(), k, 0, k, qr.data());
 std::vector<K> q(k);
 for (int i= 0; i < k; ++i) q[k - 1 - i]= qr[i];
 // r = a - q b の deg b 未満の部分 (deg b = m - 1)。
 std::vector<K> r;
 const int rlen= std::min(n, m - 1);
 if (rlen > 0) {
  std::vector<K> bq(rlen);
  convolve_range<Conv>(bb.data(), m, q.data(), k, 0, rlen, bq.data());
  r.resize(rlen);
  for (int i= 0; i < rlen; ++i) r[i]= aa[i] - bq[i];
  trim(r);
 }
 return {std::move(q), std::move(r)};
}
