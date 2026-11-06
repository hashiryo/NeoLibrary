---
title: 行列
documentation_of: ../../neo/linalg/Matrix.hpp
---

大きさを実行時に決める行列です。行優先の 1 本の配列 `a` に持ち、`a[i * w + j]` が i 行 j 列です。`a` は $2^{19}$ byte 以上なら、Linux では 2 MB 境界の mmap に置いて huge page を頼みます (allocator の `linalg_internal::HugeAlloc`)。

```cpp
template <class T> struct Matrix {
 int h, w;
 std::vector<T, linalg_internal::HugeAlloc<T>> a;
 Matrix(int h, int w);               // 0 で埋める
 Matrix(int h, int w, const T& v);
 Matrix(const std::vector<std::vector<T>>& v);
 static Matrix identity(int n);
 int height() const;
 int width() const;
 T* operator[](int i);               // i 行目の先頭
 Matrix pow(uint64_t k) const;
 Matrix transpose() const;
};
Matrix operator+(Matrix l, const Matrix& r);   // -、+=、-= も
Matrix operator*(const Matrix& l, const Matrix& r);
Matrix operator*(Matrix l, const T& c);        // c * A、A *= c も
std::vector<T> operator*(const Matrix& l, const std::vector<T>& v);
bool operator==(const Matrix& l, const Matrix& r);
```

要素の型 `T` は、和 (`+`)、積 (`*`)、0、1 を持つ半環であればよく、積は可換でなくてかまいません。0 と 1 は、`T` に静的メンバ関数 `zero()` と `one()` があればそれを、なければ `T(0)` と `T(1)` を使います。`neo/algebra/semiring.hpp` の `MinPlus` などは前者です。引き算とスカラー倍は、`T` がそれを持つときだけ使えます。

`A * B` は $h \times m$ と $m \times w$ の積、`A * v` は列ベクトルとの積、`A.pow(k)` は正方行列の $k$ 乗です。一般の `T` では、積は i-k-j の順に足す $O(hmw)$、累乗は二分累乗です。

## ZMod の行列

`T` が `ZMod` のときは、積、ベクトルとの積、累乗を専用の核で求めます。法は `ZMod` と同じく $1 \le M < 2^{30}$ で、偶奇も素数かどうかも問わず、`ZMod<0>` も使えます。

積は、値を $M$ 未満に直して u64 の和にため、最後に 1 回だけ還元します。小さい行列や細長い行列 (タイルに広げた分が 2.5 倍を超えるもの) はスカラーで掛けます。それより大きければ AVX2 の 4 行 × 8 列のタイルで掛けます。$M$ が奇数で、1 辺を 2 の冪に広げても積和の数がタイルより 2 割以上少ない大きさ (256、500 から 512、1024 など) では、Strassen の再帰にします。

累乗は、$M$ が素数で、辺が 16 以上、指数も大きいとき (辺 200 なら $K \ge 2^{15}$ 程度) は Frobenius 標準形を経由して求めます。期待される手間は $O(n^3 + n^2 \log K)$ です。それ以外は二分累乗で、$M$ が奇数なら値を Montgomery の表現のまま掛け続けます。

## 中身

タイルの核は、procon-judge の yosupo-pow-of-matrix の zmod_tile の核を、実行時の法と長方形の行列に広げたものです。B を 8 列ずつの帯に詰め直し、B の行の 8 個を 1 回読んで、偶数番目はそのまま vpmuludq に渡します。奇数番目は、x86 では 4 byte ずらした番地から読み直して渡します。vpmuludq は 64 bit の各々の下の 32 bit しか見ないので、ずらす命令が要りません。yosupo-matrix-product の zmod_strassen_hpo で、Intel の 3 機種では 3% から 5% 速く、AMD では変わりませんでした。SIMDe の道では読み直すほうが遅いので、32 bit ずらします。和は u64 の 8 本 (4 行 × 偶奇) に足し、8 段ごとに、最上位の bit が立っていれば $2M \cdot 2^{32}$ を引いて $2^{63}$ 未満に保ちます。最後は、$M$ が奇数なら和を $2^{32}$ の位で 1 回畳んでから Montgomery の還元で $2^{32}$ で割り、$2^{64} \bmod M$ を掛けて戻します。$M$ が偶数なら和を 1 つずつ Barrett で還元します。

Strassen の再帰は yosupo-matrix-product の zmod_winograd_hpo を、実行時の奇数の法に広げたものです。1 辺 64 のブロックを Morton 順に並べ、ブロックの積を同じタイルの形で求めます。再帰の中は足し算と引き算だけなので、どのブロックも同じ $2^{-32}$ 倍を持ち、最後に戻します。7 回の積は Winograd の形で、1 段の足し算と引き算は 15 回です (Strassen の形は 18 回)。並べ方は Boyer、Dumas、Pernet、Zhou (2009) の表 1 の、作業の領域 2 つで $C = AB$ を求める順です。作業の配列 (1024 次で 16 MB) は 1 つの領域にまとめ、Linux では huge page を頼んで先にページを用意させます。頼まない領域を 4 KB のページに置く判定機では、このページを初めて触る分が 7 ms ほどかかっていました。

行列の値を持つ配列を huge page に置くのも同じ理由です。matrix-product の lib-neo は、4 MB の行列を 3 つ (入力の 2 つと答え) 作るので、`std::vector<T>` に持つ形では写し元より 9 ms 遅くなっていました。

Frobenius 標準形の道は yosupo-pow-of-matrix の zmod_frob を、静的な素数の法と `ZMod<0>` の素数の法に広げたものです。乱数のベクトルから Krylov の鎖を作り、それまでの鎖の張る空間に入ったところで止めて、関係から鎖の多項式を読みます。前の鎖の成分が残れば、それを鎖の多項式で割って始めのベクトルから引き、鎖を作り直します。A を同伴行列を並べたブロック対角の形 $F$ に直し、$A^K = T F^K T^{-1}$ を 2 回のタイルの積で求めます。写し元は鎖ごとに乱数を選び直していましたが、法が小さいと前の鎖が最小多項式の次数を外して割り切れなくなり、法 2 の乱択の行列で 3 割ほどが二分累乗に落ちていました。そこで、法が 1024 未満なら乱数のベクトルを 4 本試して鎖の最も長いものを選び、割り切れなければ最初からやり直します (法が小さければ 8 回、大きければ 2 回まで)。それでも駄目なら二分累乗にします。法が素数かどうかは、底 2、7、61 の Miller-Rabin で調べます。静的な法ならコンパイル時に、`ZMod<0>` なら最後に調べた法の結果を覚えておきます。

写し元の CI の記録は、x64-gcc で次のとおりです。yosupo-matrix-product の zmod_winograd_hpo は EPYC 7763 で 49.80 ms でした (run 41412829910)。同じ回の zmod_strassen は 62.77 ms、Library Checker の首位の写しは 72.59 ms、今の Library は 887.54 ms です。yosupo-pow-of-matrix の zmod_frob は EPYC 9V45 で 6.21 ms (106 倍) でした。この行列を呼ぶ lib-neo の提出は、手元の M2 (SIMDe) で写し元と並びました。matrix-product の最大のケースは 82.1 ms と 79.8 ms、pow-of-matrix は 11.0 ms と 10.8 ms です。matrix-product の 22 ケースの合計は 370.6 ms と 477.9 ms で、辺 600 から 800 の長方形を、1024 まで広げずタイルで掛ける分だけ短くなっています。手元で決めた切り替えの大きさは、x64 では測り直していません。

C++20 以上が要ります。
