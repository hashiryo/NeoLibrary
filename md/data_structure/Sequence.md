---
title: 列の木 (塊の形)
documentation_of: ../../neo/data_structure/Sequence.hpp
---

列を切り貼りしながら、区間の積、区間への作用、区間の反転を扱う木です。1 つの節に連続する値の配列 (塊) を持たせ、木が扱う節の数を値の数の塊の大きさ分の 1 ほどにします。構造体 `M` は [SegmentTree](SegmentTree.md) と同じ決まりで渡します。

- モノイド: `using T`、`static T ti()`、`static T op(T a, T b)`
- 作用: `using E`、`static void mp(T& v, E x)` か `mp(T& v, E x, int len)` (len は区間の長さ)、`static void cp(E& pre, E suf)`

`M` の代わりに `int` などの型を渡すと、値を並べるだけの列になります。可換でないモノイドを反転するときは逆順の積も持ちます。`M` に `using commute = void;` を書くと、それを省きます。永続にはできません。永続にするなら [WeightBalancedTree](WeightBalancedTree.md) を使ってください。

```cpp
template <class M, bool REV = false, int B = 0> class Sequence;  // B は塊の大きさ (0 なら型から決める)
Sequence();
explicit Sequence(int n, const T& v = T());
Sequence(const std::vector<T>& v);
template <class It> Sequence(It bg, It ed);
int size() const;
bool empty() const;
void clear();
std::pair<Sequence, Sequence> split(int k);              // 先頭の k 個と残り
std::tuple<Sequence, Sequence, Sequence> split3(int a, int b);  // [0, a)、[a, b)、[b, n)
Sequence operator+(Sequence rhs) const;                  // つなげた列
Sequence& operator+=(Sequence rhs);
void insert(int k, const T& v);
T erase(int k);                      // k 番目を消して返す
void push_back(const T& v);
void push_front(const T& v);
T pop_back();
T pop_front();
T get(int k);
T operator[](int k);
void set(int k, const T& v);
void mul(int k, const T& v);         // k 番目を op(今の値, v) にする (モノイドがあるとき)
T prod(int l, int r);                // op(a[l], ..., a[r-1])。空なら ti() (モノイドがあるとき)
T all_prod() const;
void apply(int l, int r, const E& x);  // (作用があるとき)
void reverse();                      // 列の全体を反転 (REV のとき)
void reverse(int l, int r);          // [l, r) を反転 (REV のとき)
std::vector<T> dump();
```

どの操作も、木の深さと塊の大きさの和に比例する手間です。同じ型の木は節と塊の置き場を共有するので、別の木どうしで切り貼りできます。木の値を写すと同じ木を指すだけになるので、切ったりつなげたりしたあとは元の値を使わないでください。

`monoid::AffineSum<ZMod<P>>` (法が定数の 1 次関数と区間和、[segtree/Affine.hpp](segtree/Affine.md)) を渡すと、塊の大きさが 512 になり、塊の中の和と作用を AVX2 で回します。

## 使い分け

- op が軽く永続しないなら、この型が速いです。
- 永続にするときは [WeightBalancedTree](WeightBalancedTree.md) を使います。
- op が重く可換でない列 (行列の積など) を何度も反転するときも [WeightBalancedTree](WeightBalancedTree.md) の方が速くなります。反転で入れ替えた塊の、順と逆順の積を塊の長さだけ求め直すためです。

## 中身

木の部分は treap で、優先度の乱数の種は起動の時刻から取ります。塊の大きさ (CAP) は値の大きさから決めます (16 から 128)。塊の長さは、列の両端を除いて下限 (MIN、CAP の 4 分の 1) 以上に保ちます。併合では、境目の塊のどちらかが MIN より短いあいだ、2 つが 1 つの塊に収まればまとめ、収まらなければ 2 つの塊で長さを半分ずつに分け直します (各節に、部分木の左端と右端の塊の長さを持つ)。満杯の塊への挿入は塊を半分ずつに割り、空になった節は外します。

区間の積と作用は形を変えずに根から下り、端の塊だけを触ります。節には、子へ配る遅延と、自分の塊へ当てる遅延 (作用と反転の印) を別に持ちます。塊の反転の印は、配列を並べ替えずに添字を読み替えます。

区間の反転は、両端の塊のうち区間の外にある部分を逆順にして互いの塊へ移してから、塊の境目で割って印を付けます (Library Checker の nandhagk の 402977 の方法)。こうすると塊を割らずに済みます。移したあとの塊が CAP を超えるときだけ区間の端で塊を割り、MIN より短くなった塊は、つなぎ直すときに隣の塊とまとめます。移す 2 つの塊は、反転の印を残したまま外の部分だけを書き換えます。下限を保たないと、反転を続けるうちに短い塊が増えて節が増えます (dynamic-sequence の反転だけのケースで、平均の長さが 37 まで縮んだ)。

`monoid::AffineSum<ZMod<P>>` のときは、塊の大きさを 512、下限を 64、作るときの塊を 128 にします (ほかの型は下限が CAP の 4 分の 1、作るときは CAP の 4 分の 3)。核は、値をそのままの値で 2P 未満に持ち、掛け算は Shoup の方法で 8 個ずつ行います。区間の和は、塊を書き換えずに残っている塊の遅延を部分の和に当て、短い側の和から求めます。塊の一部への作用は、残っている塊の遅延と合わせて 1 回の走査で書きます (区間の内と外で掛ける数をレーンごとに選ぶ)。塊を割る、反転で外の部分を移す、削除、代入のときも、塊の和を短い側の和の足し引きで求めます。作用は恒等写像で初めて、印で分けずに合成します。

procon-judge で書き比べました (記録は algo-notes の notes/sequence-tree.md)。値は run 44442381691 の x64-gcc の algo の最大です。

- yosupo-dynamic-sequence-range-affine-range-sum (Xeon 6975P-C): 626.96 ms で、下限を保つ前のこの型 (794.71 ms) の 0.79 倍、同じ作りの書き比べ (chunk2_256、658.47 ms) の 0.95 倍、今の Library の WeightBalancedTree (884.12 ms) の 0.71 倍でした。
- aoj-1508 (Xeon 6975P-C): 汎用の道で 40.04 ms で、今の Library の WeightBalancedTree (141.82 ms) の 3.5 倍の速さでした。
