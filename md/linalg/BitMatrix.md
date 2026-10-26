---
title: GF(2) の行列
documentation_of: ../../neo/linalg/BitMatrix.hpp
---

成分が 0 と 1 で、和を XOR、積を AND とする行列 (GF(2) の上の行列) です。行を 64 bit の語の列に詰めて持ち、積、累乗、行列式、階数、連立方程式、核、逆行列を求めます。

```cpp
class BitMatrix;
BitMatrix(int h, int w);                                 // h × w の零行列
explicit BitMatrix(const std::vector<std::string>& s);  // 各行を '0' と '1' の文字列で
static BitMatrix identity(int n);
int height() const;
int width() const;
int words() const;                                       // 行の語の数 (⌈w / 64⌉)
unsigned long long* operator[](int i);                   // 行 i の語の列
bool get(int i, int j) const;
void set(int i, int j, bool v = true);
void flip(int i, int j);
std::vector<std::string> to_strings() const;
BitMatrix transpose() const;
BitMatrix& operator+=(const BitMatrix& r);               // 和 (XOR)。+ もある
BitMatrix operator*(const BitMatrix& l, const BitMatrix& r);
std::vector<int> operator*(const BitMatrix& l, const std::vector<int>& v);
BitMatrix pow(uint64_t k) const;
bool det() const;
int rank() const;
std::optional<std::vector<int>> solve(const std::vector<int>& b) const;
std::optional<std::pair<std::vector<int>, BitMatrix>> solve_with_kernel(const std::vector<int>& b) const;
BitMatrix kernel() const;
std::optional<BitMatrix> inverse() const;
```

行 i の j 列目は、`A[i]` の語 j / 64 の bit j % 64 です。行を語のまま読み書きしてもかまいませんが、w 列目より後ろの bit は 0 のままにしてください。文字列から作るときは文字の下位 bit を見て、行の長さはそろえます。行が 0 本の行列は `BitMatrix(0, w)` で作ります。

`solve(b)` は A x = b の解を 1 つ返し (自由な変数は 0)、解がなければ `std::nullopt` を返します。`solve_with_kernel(b)` は、解の 1 つと、A v = 0 の解空間の基底を行に並べた行列 (w - rank 行) の組を返します。`kernel()` は基底の行列だけを返します。ベクトルの要素は 0 か 1 の `int` です。`det()` と `inverse()` は正方行列で使い、`inverse()` は正則でなければ `std::nullopt` を返します。

n 次の正方行列の積、行列式、逆行列は、どれも $O(n^3 / 64)$ の語の演算で求まります。

## 中身

積 A B は Four Russians の方法で求めます。C の列を 1024 bit の帯に分け、帯ごとに、B の 64 行分 (A の 1 語分) を 8 行ずつ 256 通りに XOR した表を 8 枚作ります。そのあと、A の各行の C の帯 (ymm 4 本) に表の行を 8 本 XOR します。C の帯の読み書きは、B の 64 行につき 1 回で済みます。A の行が 64 本以下なら、A の 1 の bit ごとに B の行を足すだけにします。

行列式、階数、連立方程式、核、逆行列は、要 (ピボット) を 16 本ずつまとめる掃き出し (M4RI の形) で求めます。組の中では、残りの行ごとに、どの要の行を足すかの印と、要を探す列の語に貯めた更新を足した今の値だけを持ちます。列の bit はその値の bit なので、要を探すのも要を足すのも、4 行ずつ AVX2 で調べて書けます。要のない列では、残りの行の値の OR を見て、次に 1 のある列まで飛びます。16 本たまったら、要の行 8 本ずつの 256 通りの XOR の表を 2 枚作り、残りの各行に表の行を 1 本ずつ XOR します。行列の読み書きは要 16 本につき 1 回になります。

連立方程式の解の 1 つは、前進消去のあとの後退代入で求めます。核と逆行列は、要の行を下から 16 本ずつ組にし、組の表で上の要の行を消して RREF に直してから読みます。核の基底は、要の行を要の列の位置に、自由な列 g には単位ベクトル $e_g$ を置いた正方行列の列 g です。その行列の 0 でない 64 × 64 の小行列だけを転置して取り出します。逆行列は [A | I] を消しますが、I の列は初めには置きません。t 番目の要として選んだ行の 1 を、右の部分の列 t へ置きます。前進消去で読み書きする右の部分が、それまでの要の数の分だけになります。最後に、右の部分の列 t を元の行 prow[t] の列へ並べ替えます (転置して行を並べ替え、転置し直します)。

64 列以下の階数と 64 次以下の行列式は、1 語の行のまま掃き出します。行が 64 本以下の階数は行を 1 本ずつ入れ、行が列の 2 倍より多い行列の階数は転置してから消します。

procon-judge の yosupo の mod 2 の 5 問 (N ≤ 4096) で書き比べた gf2_mul_m4r8 と gf2_*_blk16 を写し、掃き出しの要を探す手順を上の形に書き直しました。写し元は、run 38110682413 と 38111772190 の x64-gcc で、5 問とも前の首位の 3.4 倍から 4.9 倍の速さでした。EPYC 7763 では積が 45.4 ms、階数が 39.3 ms、連立方程式が 52.0 ms、逆行列が 81.6 ms で、行列式は EPYC 9V74 で 23.9 ms です。書き直した掃き出しで、手元の M2 では 4096 × 4096 の階数が 28 ms から 19 ms になりました。階数の低い行列では列ごとに全部の行をなめなくなり、連立方程式の max_lowrank は 28 ms から 3 ms になっています。手元の M2 での 5 問の最大は、写し元の 37、30、42、54、88 ms に対して、積が 35 ms、行列式が 21 ms、階数が 40 ms、連立方程式が 35 ms、逆行列が 71 ms です。CI ではまだ測っていません。写し元の記録は algo-notes の notes/linear-algebra.md にあります。
