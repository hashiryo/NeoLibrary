---
title: セグメント木
documentation_of: ../../neo/data_structure/SegmentTree.hpp
---

構造体 `M` でモノイドと作用を渡すセグメント木です。モノイドだけなら点の更新と区間の積、作用も渡せば遅延セグ木、作用だけなら双対セグ木になります。

```cpp
template <class M> class SegmentTree;
SegmentTree(int n);                                    // ti() で埋める (作用だけなら T())
SegmentTree(int n, T x);
SegmentTree(const std::vector<T>& a);
template <class F> SegmentTree(int n, const F& init);  // a[i] = init(i)
int size() const;
T get(int i);                   // operator[] も同じ
void set(int i, T x);
void mul(int i, T x);           // a[i] = op(a[i], x)
T prod(int l, int r);           // op(a[l], ..., a[r - 1])
T all_prod();
void apply(int l, int r, E x);  // 作用があるときだけ
template <class C> int max_right(int l, const C& check);
template <class C> int min_left(int r, const C& check);
```

`M` の決まりは次のとおりです。

| 持つもの | 意味 |
| --- | --- |
| `using T`、`static T ti()`、`static T op(T a, T b)` | モノイド (単位元と積) |
| `using E`、`static void mp(T& v, E x)` | 作用。区間の長さが要るなら `mp(T& v, E x, int len)` |
| `static void cp(E& pre, E suf)` | 作用の合成。pre を先に、suf をあとに当てる合成を pre に書く |

`mp` が bool を返す形は Segment Tree Beats です。偽を返すと子へ下りて作り直します (葉では必ず真を返すこと)。`max_right` と `min_left` は ACL と同じで、`check(ti())` が真、条件が区間の伸びに対して単調であることを前提にします。

[segtree/](segtree/) の標準の構造体 (`monoid::Sum`、`monoid::AddMin`、`monoid::AffineSum` など) を渡すと、そのヘッダが中身を子を 8〜16 個持つ AVX2 の木に替えます。口は同じです。それ以外の構造体は、このヘッダの二分木で受けます。

二分木は大きさを 2 の冪に揃えた ACL と同じ非再帰の形です。区間の積は、値の型が算術型なら、節が決まった段ですぐ畳みます。それ以外の型 (1 次関数の合成や行列) では、覆う節の番号を分岐なしで集めてから、左と右の 2 本の鎖で畳みます。重い演算では、分岐の予測の外れが無くなるぶん速くなります。作用が残っているかは別の配列の印で持つので、作用の型に恒等元は要りません。1 MB 以上の配列は 2 MB 境界に置いて huge page を頼みます (Linux のみ)。速さは今の Library の SegmentTree と同じくらいです。
