#pragma once
#include <algorithm>
#include <string>
#include <vector>
#include "neo/data_structure/RangeMin.hpp"
#include "neo/string/suffix_array.hpp"
// 接尾辞どうしの最長共通接頭辞 (LCP) と、部分文字列の大小を答える。接尾辞配列と LCP 配列 (neo/string/suffix_array.hpp) を作り、
// LCP 配列に区間最小 (neo/data_structure/RangeMin.hpp) を載せる。lcp(i, j) は 2 つの接尾辞の順位の間の LCP 配列の最小で、順位の表を
// 2 回と区間最小の表を引く。procon-judge の self-string-suffix-lcp で書き比べた sa_rmq を写したもの。
class SuffixLCP {
 int n= 0;
 std::string s;
 std::vector<int> rnk;
 RangeMin rm;
public:
 SuffixLCP()= default;
 explicit SuffixLCP(const std::string& str): n(str.size()), s(str), rnk(str.size()) {
  const std::vector<int> sa= suffix_array(s), lcp= lcp_array(s, sa);
  for(int i= 0; i < n; ++i) rnk[sa[i]]= i;
  rm= RangeMin(std::vector<unsigned>(lcp.begin(), lcp.end()));
 }
 int size() const { return n; }
 // s[i, n) と s[j, n) の最長共通接頭辞の長さ。0 <= i, j <= n。
 [[gnu::always_inline]] int lcp(int i, int j) const {
  if(i == j) return n - i;
  if(i == n || j == n) return 0;
  int l= rnk[i], r= rnk[j];
  if(l > r) std::swap(l, r);
  return int(rm.min(l, r));
 }
 // s[l1, r1) と s[l2, r2) の辞書順の大小。前が小さければ負、等しければ 0、大きければ正。
 int compare(int l1, int r1, int l2, int r2) const {
  const int a= r1 - l1, b= r2 - l2, m= std::min(a, b), h= std::min(lcp(l1, l2), m);
  if(h == m) return (a > b) - (a < b);
  return (unsigned char)s[l1 + h] < (unsigned char)s[l2 + h] ? -1 : 1;
 }
};
