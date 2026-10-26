---
title: xor 基底
documentation_of: ../../neo/linalg/XorBasis.hpp
---

符号なし整数を GF(2) 上のベクトルとみて、xor で作れる値の空間 (線形空間) を扱います。1 つの基底を育てる `XorBasis`、2 つの空間の共通部分を求める `intersection`、列の区間の元が張る空間を問う `RangeXorBasis` があります。

```cpp
template <class U = std::uint64_t> class XorBasis;
bool insert(U x);                 // 張る空間が大きくなったら true
void merge(const XorBasis& o);    // o の張る空間を足す
int size() const;                 // 基底の元の数 (次元)
auto count() const;               // 張る空間の元の数 2^size()
bool contains(U x) const;
U reduce(U x) const;              // x ^ (張る空間の元) の最小
U max_xor(U x = 0) const;         // x ^ (張る空間の元) の最大
U kth(U k, U x = 0) const;        // x ^ (張る空間の元) の、小さいほうから k 番目 (0 始まり)
U operator[](int i) const;        // 基底の元 (begin() と end() でも回せる)

template <class U> XorBasis<U> intersection(const XorBasis<U>& a, const XorBasis<U>& b);

template <class U = std::uint64_t> class RangeXorBasis;
explicit RangeXorBasis(const std::vector<U>& a);
bool contains(int l, int r, U x) const;      // x が a[l..r) の元の xor で作れるか
U max_xor(int l, int r, U x = 0) const;      // x ^ (a[l..r) の元の xor) の最大
int rank(int l, int r) const;                // a[l..r) の元が張る空間の次元
XorBasis<U> basis(int l, int r) const;       // a[l..r) の元が張る空間の基底
```

`U` は 128 bit までの符号なし整数です (`unsigned __int128` も使えます)。張る空間はいつも 0 を含むので、`contains(0)` は真で、`kth(0)` は 0 です。`kth` は $0 \le k < 2^{\mathrm{size}()}$ で使い、`x` を渡すと $x$ を足した値の並びの k 番目を返します (グラフの閉路の空間で、道の xor を最小や k 番目にするときの形です)。`count` は U が 64 bit 未満なら `std::uint64_t`、64 bit 以上なら `unsigned __int128` で返し、U が 128 bit で次元が 128 のときだけ 0 になります。

`XorBasis` の 1 回の操作は、`insert`、`reduce`、`contains`、`max_xor`、`kth` がどれも基底の元の数に比例する時間で、`merge` は o の元を 1 つずつ入れます。`intersection` は 2 つの次元の積に比例する時間です。基底の元は、入れた順に、前の元で消したものが並びます。

`RangeXorBasis` の区間は半開区間で、$0 \le l \le r \le N$ です。$l = r$ の空の区間は {0} を張るとみます。構築は $O(N L)$、1 回の問い合わせは $O(16 L)$ ほどです (L は列の値の最上位の bit の位置の最大に 1 を足したもの)。メモリは $(N / 16 + 1) L$ 組の (U, int) と列の写しで、U = u64、$N = 4 \times 10^5$、$L = 60$ で 21 MB ほどです。

## 中身

基底の元は入れた順に持ち、あとから入った元は前の元の最上位の bit (要) をどれも持たないようにしておきます。値 z を基底で消すときは、入れた順に $z \leftarrow \min(z, z \oplus b)$ を当てます。b の要が z に立っているときだけ $z \oplus b < z$ なのでその bit が消え、あとの元はその bit を持たないので、消した bit がまた立つことはありません。残るのは z に張る空間の元を足した値のうち要がどれも 0 の唯一の値で、それがそれらの最小です。$\max$ で同じことをすると要がどれも 1 の値、つまり最大が残ります。分岐がなく、要の位置を clz で探して消す形の 3.4 倍の速さでした (procon-judge の yosupo-intersection-of-f2-vector-spaces の zassenhaus_min)。

値の大小は要の bit の並びで決まるので、`kth` は k の bit を要の位置に配った K (x86 では BMI2 の pdep) を作り、要の bit が K になる値を `reduce` を 2 回使って求めます。

`intersection` は Zassenhaus の方法です。元の多いほうを A、少ないほうを B とし、2W bit の行 (x, t) を A の元 a から (a, 0)、B の元 b から (b, b) として並べ、x の側で消していきます。x が 0 まで消えた B の行の t は共通部分に入ります。B の元を要の小さい順に流すと、t の要は流した元の要と同じになるので、得られた t を逆に並べればそのまま基底の形になります。2W が 64 以下なら行を 1 語に詰めます。A が $2^{\dim A}$ 未満の値をすべて張り、B の元がどれも $2^{\dim A}$ 未満のときは、B をそのまま返します。

`RangeXorBasis` は、接頭辞ごとに、要 h ごとに最も右の元を残す基底を使います。新しい元を入れるときは、上の bit から見て、その欄の元より新しいほうを欄に残し、2 つの xor を古いほうの位置で下へ流します。こうすると、位置が l 以上の欄の元が、位置 l 以上 r 未満の元の張る空間の基底になります。写しは 16 個おきにだけ持ちます。問い合わせでは、r を 16 の倍数 kB に切り下げた写しのうち位置が l 以上の欄の元で、x と残りの元 (位置 max(kB, l) から r - 1 まで、15 個以下) を上の bit から消します。消したあとの値は写しの要をどれも持たないので、残りの元を入れた順に互いに消してから x に当てれば判定できます。U が 64 bit なら、残りの元は AVX2 の 4 本に持って欄ごとにまとめて消します。構築の挿入は、x が 0 になるか空の欄に入ったら止め、欄の入れ替えはマスクの算術で書きます (三項演算子では gcc が分岐にしました)。

procon-judge の 2 問を、手元の Apple M2 Max (Apple clang 21) で測りました (algo_time の最大)。CI の x64 ではまだ測っていません。

| 問題 | NeoLibrary | 比べた提出 |
| --- | ---: | ---: |
| atcoder-abc223-h ($N = 4 \times 10^5$、$Q = 2 \times 10^5$) | 133 ms | snap16_red 172 ms |
| yosupo-intersection-of-f2-vector-spaces ($T = 10^5$) | 116 ms | zassenhaus_min 88 ms |

abc223-h は、写しの欄を値の最上位の bit までに絞る分、挿入の打ち切り、残りの元の AVX2 で、首位の snap16_red より速くなりました。intersection-of-f2-vector-spaces は、2 つの空間をそれぞれ `XorBasis` に入れる手間の分だけ zassenhaus_min より遅くなりました。最大のケース (29 次元の部分空間の組) で 1.3 倍かかります。zassenhaus_min は入力の元が線形独立と分かっている問題に合わせた形で、片方の基底を作らずに消し、30 本の空間なら何も計算せずに返します。記録は algo-notes の notes/linear-algebra.md にあります。
