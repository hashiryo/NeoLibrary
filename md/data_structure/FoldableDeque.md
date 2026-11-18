---
title: 積を取れる両端キュー
documentation_of: ../../neo/data_structure/FoldableDeque.hpp
---

両端に足し、両端から消し、列の全体の積を返す両端キューです (両端のスライド窓の集約)。モノイドは [SegmentTree](SegmentTree.md) と同じ構造体 `M` (`using T`、`static T ti()`、`static T op(T a, T b)`) で渡し、可換でなくてかまいません。

```cpp
template <class M> class FoldableDeque;
FoldableDeque();
void reserve(int n);                 // 列が n 個までなら配列を伸ばさないよう、先に取る (中は 2n 個ぶん)
int size() const;
bool empty() const;
void push_front(const T& x);
void push_back(const T& x);
void pop_front();                    // 空なら使えない
void pop_back();                     // 空なら使えない
const T& front() const;
const T& back() const;
const T& operator[](int i) const;    // 先頭から i 番目
T all_prod() const;                  // op(先頭, ..., 末尾)。空なら ti()
void clear();
```

1 回の `push_front`、`push_back`、`all_prod` は $O(1)$ で、`pop_front` と `pop_back` もならして $O(1)$ です。列の長さの上限が分かっているなら、先に `reserve` を呼ぶと、列が伸びるたびに配列を取り直す分を省けます。足すのが末尾だけで消すのが先頭だけなら、[FoldableQueue](FoldableQueue.md) の方が速く動きます。

## 中身

向かい合わせの 2 本のスタックで持つ定番の形を、1 本の配列の上で持ちます。列を配列の区間 `[hd, tl)` に置き、境 `md` の前の `[hd, md)` には各位置から `md` の手前までの積を、後ろの `[md, tl)` には `md` から各位置までの積を、値と並べて置きます。片方が空のときにそちらから消すと、境を列の真ん中へ移して、両方の積をその場で数え直します。2 本のスタックで半分を移す形と同じ手間で、値を別の配列へ写しません。配列の端に空きが無くなったら、生きている要素が半分以下なら同じ配列の真ん中へ、そうでなければ 2 倍の配列の真ん中へ移します。1 MB 以上の配列は 2 MB 境界に置いて huge page を頼みます (Linux のみ)。

procon-judge の yosupo-deque-operate-all-composite で書き比べました。Intel の Xeon 6975P-C では、列が大きく伸びるケースで、2 本のスタックを別々の vector に持つ形 (swag) の 1.2 倍かかりました。伸ばすたびに倍の配列を初めて触って要素を写す分です。構築のときに `reserve` で容量を取っておく形は、どのケースでも swag より速くなりました (max_random で 4.90 ms と 6.48 ms)。記録は algo-notes の notes/segment-tree-design.md にあります。
