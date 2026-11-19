---
title: 列の重み平衡木 (永続も選べる)
documentation_of: ../../neo/data_structure/WeightBalancedTree.hpp
---

列を切り貼りしながら、区間の積、区間への作用、区間の反転を扱う重み平衡木です。中身は今の Library の WeightBalancedTree を写したもので、永続にもできます。構造体 `M` は [SegmentTree](SegmentTree.md) と同じ決まりで渡します。

- モノイド: `using T`、`static T ti()`、`static T op(T a, T b)`
- 作用: `using E`、`static void mp(T& v, E x)` か `mp(T& v, E x, int len)` (len は区間の長さ)、`static void cp(E& pre, E suf)`

`M` の代わりに `int` などの型を渡すと、値を並べるだけの列になります。可換でないモノイドを反転するときは逆順の積も持ちます。`M` に `using commute = void;` を書くと、それを省きます。

```cpp
template <class M, bool REV = false, bool PERSIST = false> class WeightBalancedTree;
WeightBalancedTree();
explicit WeightBalancedTree(int n, const T& v = T());
WeightBalancedTree(const std::vector<T>& v);
template <class It> WeightBalancedTree(It bg, It ed);
static void reserve(std::size_t n);  // 置き場を節と葉それぞれ n 個ぶん先に取る
int size() const;
bool empty() const;
void clear();
std::pair<WeightBalancedTree, WeightBalancedTree> split(int k) const;                         // 先頭の k 個と残り
std::tuple<WeightBalancedTree, WeightBalancedTree, WeightBalancedTree> split3(int a, int b) const;  // [0, a)、[a, b)、[b, n)
WeightBalancedTree operator+(const WeightBalancedTree& rhs) const;  // つなげた列
WeightBalancedTree& operator+=(const WeightBalancedTree& rhs);
void insert(int k, const T& v);
T erase(int k);                      // k 番目を消して返す
void push_back(const T& v);
void push_front(const T& v);
T pop_back();
T pop_front();
T get(int k) const;
T operator[](int k) const;
void set(int k, const T& v);
void mul(int k, const T& v);         // k 番目を op(今の値, v) にする (モノイドがあるとき)
T prod(int l, int r) const;          // op(a[l], ..., a[r-1])。空なら ti() (モノイドがあるとき)
T all_prod() const;
void apply(int l, int r, const E& x);  // (作用があるとき)
void reverse();                      // 列の全体を反転 (REV のとき)
void reverse(int l, int r);          // [l, r) を反転 (REV のとき)
std::vector<T> dump() const;
static bool pool_empty();            // 永続のとき、置き場を詰めるとよい頃なら真
static void gc(std::vector<WeightBalancedTree>& vs);  // 渡した版から届く節だけを残して置き場を詰める
static void gc(WeightBalancedTree& v, ...);           // 版を並べて渡す形
```

どの操作も $O(\log n)$ です。同じ型の木は節の置き場を共有するので、別の木どうしで切り貼りできます。

## 永続にする

`PERSIST` を真にすると、どの操作も道の上の節を写して新しい根を作るので、木の値を写せば版の写しになります ($O(1)$)。永続では操作のたびに節が増えるので、置き場が伸びたら、生きている版を渡して `gc` を呼びます。

```cpp
using P = WeightBalancedTree<M, true, true>;
vector<P> ver(q + 1);
ver[i] = ver[j];          // 版の写し
ver[i].insert(p, v);
if (P::pool_empty()) P::gc(ver);
```

`gc` は渡した版から届く節に印を付けて置き場の前へ詰め、版の根を書き換えます。版どうしで共有している節は共有したまま残り、手間は生きている節の数で済みます。渡さなかった木は使えなくなります。`pool_empty()` は、前の `gc` のあとに生きていた節の 2 倍 (少なくとも $2^{22}$ 個) まで置き場が伸びたら真を返します。`gc` で半分より多くが生き残ったときは、次までの伸びを 4 倍に広げます。版を全部残す使い方 (永続キューなど) では詰めても減らないので、`gc` を呼ばなくてかまいません。

永続では操作のたびに節を新しく作るので、全体の大きさの見当がつくなら、先に `reserve` で置き場を取ると速くなります。配列が伸びるたびに作り直して写す分を省けるためです。

永続しないとき (`PERSIST` が偽) は、節をその場で書き換え、split で外れた節と消した葉を使い回すので、`gc` は要りません (`pool_empty()` は偽を返す)。木の値を写すと同じ木を指すだけになるので、切ったりつなげたりしたあとは元の値を使わないでください。

## 中身

値は葉にだけ持ち (葉は負の番号で指す)、左右の葉の数の比を 4 以下に保ちます。split と併合は、重い方の背を下りて回転 1 回で直す join で書きます。split は、割れ目にあった節をつなぎ目に使い回します。遅延の印 (作用と反転) は、葉の数を入れる欄の上位 2 bit に置きます。区間の積と 1 点の取得は、遅延を配らずに根から合わせて下りるので、永続のときにも節を写しません。置き場は伸びる配列で、1 MB 以上は 2 MB 境界に置いて huge page を頼みます (Linux のみ)。

procon-judge で今の Library の WeightBalancedTree と比べました (記録は algo-notes の notes/sequence-tree.md)。

- 永続の aoj-3024 (EPYC 7763): 置き場を先に取ると 105.69 ms で今の Library (104.44 ms) に並びました。先に取らないと 122.27 ms でした。
- 永続しない形の、op が重く可換でない self-seqtree-matrix-composite (2×2 行列の積) では 465.50 ms で、今の Library (466.43 ms) と並びました。

op が軽く永続しない問題 (dynamic-sequence-range-affine-range-sum など) では、1 つの節に値の配列を持たせる塊の形の方が速くなります。
