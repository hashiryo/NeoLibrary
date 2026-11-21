#pragma once
#include <cstring>
#include <string>
#include <vector>
// 接頭辞関数 (KMP の失敗関数)。pi[i] は s[0, i] の、自分より短い最長の border の長さ。KMP の一致の判定は分岐で書く
// (k += (s[i] == s[k]) と書くと、GCC が sete と add にして、次に読む s[k] の位置が前の読み込みを待つ鎖になる)。文字列の版は、
// 使われている文字が 4 種類以下なら、状態 k が 64 未満の間は KMP のオートマトンの遷移の表 (64 行 4 列) を引いて次の状態を求める。
// 2 文字の一様な列では border をたどる分岐が半分ほど外れるので、表の 1 回の読み込みの方が速い。状態 64 以上と、5 種類以上の文字の
// 列は分岐の KMP にする。procon-judge の self-string-prefix-function で書き比べた kmp_tab2 を写したもの。
namespace prefix_function_internal {
template <class S> inline void kmp(const S& s, int n, int* p) {
 p[0]= 0;
 for(int i= 1, k= 0; i < n; ++i) {
  if(s[i] == s[k]) {
   ++k;
  } else if(k > 0) {
   do k= p[k - 1];
   while(k > 0 && !(s[i] == s[k]));
   if(s[i] == s[k]) ++k;
  }
  p[i]= k;
 }
}
}
inline std::vector<int> prefix_function(const std::string& str) {
 const int n= str.size();
 std::vector<int> pi(n);
 if(n == 0) return pi;
 const unsigned char* t= reinterpret_cast<const unsigned char*>(str.data());
 int* p= pi.data();
 unsigned char code[256]= {};
 int sigma= 0;
 {
  bool seen[256]= {};
  for(int i= 0; i < n; ++i) seen[t[i]]= true;
  for(int c= 0; c < 256; ++c)
   if(seen[c]) code[c]= (unsigned char)(sigma++ & 3);
 }
 if(sigma > 4) {
  prefix_function_internal::kmp(t, n, p);
  return pi;
 }
 constexpr int R= 64, W= 4;
 alignas(64) int d[R * W]= {};
 d[code[t[0]]]= 1;
 p[0]= 0;
 int k= 0;
 for(int i= 1; i < n; ++i) {
  const int c= code[t[i]];
  if(i < R) std::memcpy(d + i * W, d + k * W, sizeof(int) * W), d[i * W + c]= i + 1;
  if(k < R) {
   k= d[k * W + c];
  } else if(t[i] == t[k]) {
   ++k;
  } else {
   do k= p[k - 1];
   while(k >= R && t[i] != t[k]);
   k= k >= R ? k + 1 : d[k * W + c];
  }
  p[i]= k;
 }
 return pi;
}
// 値を == で比べる版。分岐の KMP。
template <class T> std::vector<int> prefix_function(const std::vector<T>& s) {
 const int n= s.size();
 std::vector<int> pi(n);
 if(n) prefix_function_internal::kmp(s, n, pi.data());
 return pi;
}
