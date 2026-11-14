---
title: k 値の燃やす埋める
documentation_of: ../../neo/combinatorial_optimization/KaryOptimization.hpp
---

変数 $x_i$ が 0 から $k_i - 1$ の値を取るとき、1 変数の項と 2 変数の項の和の最小を求めます。1 変数の項は何でもよく、2 変数の項は Monge である必要があります。項は値の表か、値を返す関数で渡します。

```cpp
template <class Cost> class KaryOptimization {
  explicit KaryOptimization(const std::vector<int>& k);  // 変数 i は 0 から k[i] - 1 の値を取る
  KaryOptimization(int n, int k);
  int add_variable(int k);
  void add(Cost c);                                         // 定数
  void add(int i, const std::vector<Cost>& c);              // x_i が a なら c[a]
  void add(int i, F f);                                     // x_i が a なら f(a)
  void add(int i, int j, const std::vector<Cost>& c);       // (x_i, x_j) が (a, b) なら c[a k_j + b]
  void add(int i, int j, F f);                              // (x_i, x_j) が (a, b) なら f(a, b)
  std::optional<Cost> solve();
  std::vector<int> argmin();
};
```

同じ変数や同じ組に何度 `add` してもよく、項はすべて足し合わされます。2 変数の項は、足した組の分だけ網の辺になります。今の Library の monge_mincut のように全組の費用を呼ばないので、辺のある組だけを足せば済みます。`solve()` は和の最小を返し、`argmin()` は直前の `solve()` の最小を与える割り当てを返します。最大化は、費用の符号を反転して渡し、`solve()` の符号を戻します。

<!-- textlint-disable ja-technical-writing/max-comma -->
2 変数の項 $f(a, b)$ は、$f(a, b) + f(a + 1, b + 1) \le f(a, b + 1) + f(a + 1, b)$ (Monge) を満たす必要があります。$w |A_a - B_b|$ ($A$ と $B$ は昇順) や、$a - b$ の凸関数がこの形です。逆向きの不等号の項 (逆 Monge) があっても、いくつかの変数の値の並びを逆にして全部の項を Monge にできるなら、それを探して解きます。どちらでもない項があるか、どう逆にしても直らないときは、`solve()` が `nullopt` を返します。

禁止したい割り当てには、大きな費用を置きます。ただし表は有限の値のまま Monge である必要があります。$(x_i, x_j) = (2, 0)$ のように角の 1 つを禁じる表や、$(0, 2)$ と $(2, 0)$ の 2 つを禁じる表は Monge です。$x_i - x_j \ge 2$ を禁じるなら、$\mathrm{INF} \cdot [a - b \ge 2]$ は Monge でないので、$\mathrm{INF} \cdot \max(0, a - b - 1)$ と書きます。`Cost` は整数の型で、和は中で 128 bit で数えます。最小が `Cost` に収まらないときは、`Cost` の最大か最小に丸めた値を返します。
<!-- textlint-enable ja-technical-writing/max-comma -->

## 中身

変数 i に頂点を $k_i - 1$ 個置き、$[x_i \ge a]$ が 1 なら t の側とします。頂点を鎖に並べ、1 変数の表の値を鎖の辺の容量に、並びを保つ大きな容量を逆向きに載せます (Ishikawa の形)。2 変数の項は、表の 2 階の差分が 0 でない升目ごとに両向きの辺 1 本にします。2 つの向きへの分け方は、1 次の部分の増分ができるだけ 0 に近くなるものを選びます。そうすると、禁止のために置いた大きな値は辺の容量に入り、定数の項へ漏れません。並びを逆にする変数は、Monge の項は両端を同じに、逆 Monge の項は逆にする、という制約を、重み付きの union-find で 2 色に塗って決めます。

最小カットは [最大流](MaxFlow.md) で求めます。容量の扱いは [燃やす埋める](BinaryOptimization.md) と同じで、自明なカットの値の最小を上界にして 64 bit で流し、上界が $2^{61}$ 以上のときだけ 128 bit の容量で流します。

小さい乱択の入力 (逆にすれば直るもの、直らないもの、禁止の角や凸な差の関数を混ぜたもの) で、総当たりと 90 万回ほど突き合わせました。procon-judge では、ARC 107 F (3 値、2 つの角の禁止)、ARC 129 E (M 値、$W |A_a - B_b|$)、yuki-119 (3 値、$2^{62}$ での角の禁止) で確かめています。記録は algo-notes の notes/network-flow.md にあります。
