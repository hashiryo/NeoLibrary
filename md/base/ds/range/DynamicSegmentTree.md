---
title: 添字の大きいセグメント木
documentation_of: ../../neo/base/ds/range/DynamicSegmentTree.hpp
---

添字の範囲 `[0, n)` が大きい (n は long long) セグメント木です。触った所だけ節を作ります。構造体 M の決まりと口は [SegmentTree](SegmentTree.md) と同じで、添字と区間の端が long long になります。はじめはどの要素も `ti()` (作用だけなら `T()`) です。

```cpp
template <class M> class DynamicSegmentTree;
explicit DynamicSegmentTree(long long n);   // どの要素も ti() (作用だけなら T())
long long size() const;
T get(long long i);                         // operator[] も同じ
void set(long long i, T x);
void mul(long long i, T x);                 // a[i] = op(a[i], x)
T prod(long long l, long long r);           // op(a[l], ..., a[r - 1])
T all_prod();
void apply(long long l, long long r, E x);  // 作用を渡したときだけ
template <class C> long long max_right(long long l, const C& check);
template <class C> long long min_left(long long r, const C& check);
```

渡す構造体で、中身が次の 2 つに替わります。

| 構造体 | 中身 | 節の数 | 1 回の操作 |
| --- | --- | --- | --- |
| モノイドだけ | パトリシア木 (枝分かれの無い道を 1 つの節に畳んだ二分木) | set した添字の数の 2 倍未満 | 木の深さ (添字のばらつきで決まり、64 以下) |
| 作用がある | 深さ $\lceil \log_2 n \rceil$ の二分木で、子の無い節は「その下の要素は `ti()` のまま」を表す | 1 回の `apply` で $4 \lceil \log_2 n \rceil$ 個まで増える | $O(\log n)$ |

`mp` が bool を返す形 (Segment Tree Beats) は受けません。区間の長さを受ける `mp` には長さを long long のまま渡すので、n が $2^{31}$ を越えるときは `mp` も長さを long long で受けてください。はじめの値を `ti()` 以外にする形はありません。作用は、まだ触っていない要素 (値は `ti()`) にもそのまま当たります。そのため、作用を当てたあとの値に対して `ti()` が単位元でなくなる構造体は使えません。たとえば区間加算と区間最小で `ti()` を $10^{18}$ にすると、触っていない区間に正の数を足したところで最小が崩れます。添字が int に収まって配列を取れるなら、[SegmentTree](SegmentTree.md) の方が速く動きます。

## 中身

パトリシア木の節は `[lo, lo + 2^lv)` を受け持ち、葉は `lv = 0` です。内部の節は必ず子を 2 つ持ちます。`set` で範囲の外の添字が来たら、その節と新しい葉を子に持つ節を、その節の場所に置きます。区間の積は、l と r - 1 が同じ子にある間は根から下り、分かれた節から先は、左の道で l 以上の側、右の道で r - 1 以下の側を集めます。`max_right` と `min_left` は、l から先 (r より前) の節を順にたどり、全部を足しても判定が真の節は丸ごと足して、偽になる節だけを下ります。

作用がある形では、子を作るのは、下りる道で子が無いときと、節に残った作用を子へ配るときだけです。作用は値の等しい要素の並びを値の等しい並びへ移すので、区間に作用を当てても、子ができるのは区間の端の道だけになります。区間の積は書き込みをせず、配っていない作用は、子の積を畳んだ戻りに、区間との重なりの長さで当てます。

どちらも節は配列に置いて子を添字で持ち、1 MB 以上の配列は 2 MB 境界に置いて huge page を頼みます (Linux のみ)。

標準の構造体 `monoid::AffineSum<ZMod<P>>` (P は奇数) を渡すと、区間の 1 次関数と区間の和だけを扱う形へ中身を替えます。使うときは neo/base/core/algebra/ZMod.hpp と neo/base/ds/range/segtree/Affine.hpp も読みます。子の無い節は「区間の値がすべて等しい」ことを表してその値を持ち、作用は値をそのまま書き換えます。値は Montgomery で持ち、長さ $2^h$ の節の和は $2^h$ の表から引きます。

方式は procon-judge の次の 2 問で比べて選びました。どちらも N は $10^9$ までで、クエリを 1 つずつ渡します。

| 問題 | 操作 | Q |
| --- | --- | --- |
| yosupo-point-set-range-composite-large-array | 1 点の更新と区間の 1 次関数の合成 | $5 \times 10^5$ まで |
| yosupo-range-affine-range-sum-large-array | 区間の 1 次関数と区間の和 | $10^5$ まで |

Xeon 6975P-C の x64-gcc で、この型を使う提出は 211.2 ms と 142.7 ms でした。今の Library の SegmentTree_Patricia と SegmentTree_Dynamic は同じ回で 374.6 ms と 322.1 ms です。写し元の提出は 203.5 ms と 151.9 ms です (x64-clang では、この型が 224.1 ms と 141.9 ms、写し元が 198.3 ms と 151.4 ms)。区間の和の問題では上の `monoid::AffineSum<ZMod<P>>` の形が写し元より速く、メモリも写し元の 167 MB に対して 127 MB です。1 点の更新の問題では、写し元が値を Montgomery で持って区間の 1 次関数を x に 1 回ずつ当てるのに対し、この型は ZMod の合成で畳むので、掛け算が多くなります。CPU によって写し元の 1.0 倍から 1.4 倍かかります (EPYC 9V45 で 0.97 倍、9V74 で 1.38 倍)。記録は algo-notes の notes/segment-tree-design.md にあります。
