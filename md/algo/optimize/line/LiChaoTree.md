---
title: Li Chao 木
documentation_of: ../../neo/algo/optimize/line/LiChaoTree.hpp
---

関数の集合に関数を足しながら、点 x での最小 (または最大) を引きます。直線 y = a x + b と線分を既定で受け、どの 2 つも高々 1 回しか交わらない関数の族 (擬直線) なら、直線でなくても受けます。追加も問いも O(log) です。

型は 2 つあります。

```cpp
template <class T= long long, bool MAX= false, class F= LiChaoLine<T, MAX>> class LiChaoTree;
template <class T= long long, bool MAX= false, class F= LiChaoLine<T, MAX>> class LiChaoTreePoints;
```

- `LiChaoTree` は、x の範囲 `[lo, hi)` を渡します。問いの x を前もって知らなくても使えます。
- `LiChaoTreePoints` は、問いの x を全部先に渡します。線分が多いときは、こちらが数倍速いです。

`T` は値の型、`MAX` が true なら最大を引きます。関数の値 (直線なら a x + b) は `T` に収まる必要があります。関数が 1 つも載っていない点の答えは、最小なら `T` の最大値、最大なら `T` の最小値です。

傾きが単調な順に直線を足すだけなら、MonotoneCHT のほうが 3〜6 倍速いです。

## 使い方

```cpp
LiChaoTree<long long> t(lo, hi);        // x の範囲 [lo, hi)、最小
t.add_line(a, b);                       // 直線 y = a x + b
t.add_segment(l, r, a, b);              // 線分 y = a x + b (x in [l, r))
long long y= t.query(x);

LiChaoTree<long long, true> u(lo, hi);  // 最大

LiChaoTreePoints<long long> p(xs);      // 問いの x を全部渡す (順は問わず、重複があってよい)
p.add_segment(l, r, a, b);
p.query(x);                             // 先に渡した点 x を二分探索で引く
p.query_index(i);                       // 並べて重複を除いた i 番目の点
```

| 口 | 中身 |
| --- | --- |
| `LiChaoTree(lo, hi, reserve_nodes= 0)` | x の範囲。`reserve_nodes` は節の数の見込み (全体に足す関数 1 つにつき高々 1 つ) |
| `LiChaoTreePoints(xs)` | 問いの x の列。並べて重複を除く |
| `LiChaoTreePoints(xs, SortedUnique{})` | 並べて重複を除いた x の列をそのまま使う |
| `add_line(a, b)` | 直線を足す |
| `add_segment(l, r, a, b)` | x in `[l, r)` の線分を足す |
| `add(f)`、`add(f, l, r)` | 関数 `f` を全体に、または x in `[l, r)` に足す |
| `query(x)` | x での最小 (最大) |
| `index(x)`、`add_index(f, l, r)`、`query_index(i)` | `LiChaoTreePoints` だけ。点の添字で足し引きする |

## 直線でない関数

`F` に、x (`long long`) を受けて `T` を返す `operator()` を持つ型を渡します。どの 2 つの関数も、x の範囲の中で高々 1 回しか交わらないことが条件です。たとえば、2 次の係数が共通の放物線は、差が 1 次式なのでこの条件を満たします。

```cpp
struct Curve {
 long long j, dp;
 long long operator()(long long x) const { return dp + (x - j) * (x - j - 1) / 2; }
};
LiChaoTree<long long, false, Curve> t(lo, hi);
t.add(Curve{j, dp});
```

`F` が static の `none()` で「どこでも答えにならない値 (最小なら `T` の最大値) を返す関数」を出せると、関数の載っていない節にそれを置くので、節ごとの印が要らず、少し速くなります。`LiChaoTreePoints` に渡す `F` は、既定のコンストラクタも持つ必要があります。

## 中身

追加は、節の両端で新しい関数と節の関数を比べます。両端とも新しい関数が勝てば置き換え、両端とも負ければ止めます。片方だけ勝つなら、中点で勝つほうを節に残し、負けたほうを、それが勝つ端のある側の子へ下ろします。交わりが高々 1 回なので、両端で勝つ関数は節の区間の全体で勝ちます。

`LiChaoTree` は、範囲を 2 の冪に切り上げ、節を使う所だけ作ります。線分は、範囲を覆う節に分けてから、それぞれの節の下へ下ろします。`LiChaoTreePoints` は、並べた点の上の完全二分木を配列に置き、問いは葉から根へ上ります。通る節の添字が前の節の結果を待たないので、問いが速く、木も浅く済みます。

procon-judge で書き比べた核の、x64-gcc (Xeon 6975P-C) の計測区間の最大です (run 45323501621)。

| 問題 | LiChaoTree | LiChaoTreePoints | 今の Library |
| --- | --- | --- | --- |
| yosupo-line-add-get-min (直線、任意の順) | 61.6 ms | 86.3 ms (添字で問うと 70.9 ms) | ConvexHullTrick 108.7 ms、LiChaoTree 123.7 ms |
| yosupo-segment-add-get-min (線分) | 651.7 ms | 119.6 ms | LiChaoTree 770.3 ms |
| yuki-2458 (直線、最大) | 26.3 ms | 72.7 ms | ConvexHullTrick 25.1 ms、LiChaoTree 48.9 ms |
| yuki-409 (曲線のまま) | 76.0 ms | 66.4 ms | LiChaoTree 72.1 ms |

記録は algo-notes の notes/dp-optimization.md にあります。
