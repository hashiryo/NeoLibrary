---
title: ローリングハッシュ
documentation_of: ../../neo/string/RollingHash.hpp
---

列の部分列のハッシュを、法 $2^{61} - 1$ の多項式で求めます。

```cpp
explicit RollingHash(const std::string& s);
template <class T> explicit RollingHash(const std::vector<T>& s);
int size() const;
std::uint64_t get(int l, int r) const;    // s[l, r) のハッシュ
bool equal(int a, int b, int len) const;  // s[a, a + len) と s[b, b + len) が等しいか
int lcp(int a, int b) const;              // s[a, n) と s[b, n) の最長共通接頭辞の長さ
static std::uint64_t hash(const std::string& t);
template <class T> static std::uint64_t hash(const std::vector<T>& t);
static std::uint64_t concat(std::uint64_t x, std::uint64_t y, std::size_t len_y);
```

基数はプロセスごとに乱択し、すべてのインスタンスで同じものを使います。そのため、別の列の部分列のハッシュや、`hash(t)` で求めた列のハッシュと比べられます。文字列の版は文字を unsigned char の値として、`std::vector<T>` の版は値を $2^{61} - 1$ で割った余りとして扱います。値が 0 の項が先頭にあっても列のハッシュは変わらないので、長さの違う列を比べるときは長さも比べてください。長さの同じ異なる 2 つの列のハッシュが一致する確率は、長さを $\ell$ として、基数の選び方について $\ell / (2^{61} - 1)$ 以下です。

- 構築は $O(n)$ です。
- `get` と `equal` は $O(1)$ で、`equal` の掛け算は 1 回です。
- `lcp` は、長さを 1、2、4、... と倍にして範囲を絞ってから二分探索します。LCP が $L$ なら $O(\log L)$ です。
- `concat(x, y, len_y)` は、ハッシュ `x` の列の後ろに、ハッシュ `y` で長さ `len_y` の列をつないだもののハッシュを返します。$O(\log \mathit{len}_y)$ です。

$h_i$ を `s[0, i)` のハッシュ、$p_i$ を基数の $i$ 乗として持ちます。`s[l, r)` のハッシュは $h_r - h_l p_{r - l}$ です。2 つの部分列が等しいかは、$h_{a + \mathit{len}} - h_{b + \mathit{len}} = (h_a - h_b) p_{\mathit{len}}$ で見ます。表は 2 MB 以上なら、Linux では 2 MB 境界の mmap に置いて huge page を頼みます。$h$ は 512 個ごとに 1 つずらして置きます。2 MB 境界から置いた $h$ と $p$ の同じ添字の項や、2 の冪おきの項が、キャッシュの同じ組へ集まるのを防ぐためです。

procon-judge の self-string-substring-equal (長さ $10^6$ の文字列に $10^6$ 個の一致判定) で書き比べました。EPYC 7763 の x64-gcc では、構築とクエリを合わせた時間がケースごとの最大で 16.1 ms でした。同じ回に、今の Library の RollingHash (法 998244353 の 2 つの組) は 53.8 ms、表を 4 KB のページに置いた形は 28.5 ms でした。接尾辞配列と区間最小で LCP を引く形は 69.2 ms です。法を $2^{64}$ にした多項式のハッシュは、Thue–Morse 列の長さ $2^{12}$ の塊とその a と b を入れ替えた塊を、どの基数でも区別できません。記録は algo-notes の notes/string-basics.md にあります。
