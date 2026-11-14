---
title: 積を取れるキュー (SWAG)
documentation_of: ../../neo/data_structure/FoldableQueue.hpp
---

末尾に足し、先頭から消し、列の全体の積を返すキューです (スライド窓の集約、SWAG)。モノイドは [SegmentTree](SegmentTree.md) と同じ構造体 `M` (`using T`、`static T ti()`、`static T op(T a, T b)`) で渡し、可換でなくてかまいません。

```cpp
template <class M> class FoldableQueue;
FoldableQueue();
int size() const;
bool empty() const;
void push(const T& x);  // 末尾に足す
void pop();             // 先頭を消す (空なら使えない)
T all_prod() const;     // op(先頭, ..., 末尾)。空なら ti()
void clear();
```

1 回の `push` と `all_prod` は $O(1)$ で、`pop` もならして $O(1)$ です。両端に足したり消したりするなら [FoldableDeque](FoldableDeque.md) を使います。

## 中身

2 本のスタックで持つ定番の形を、1 本の配列の上で持ちます。列を配列の区間 `[hd, tl)` に置き、境 `md` の前の `[hd, md)` には各位置から `md` の手前までの積を、後ろの `[md, tl)` には足した値をそのまま置きます。後ろの全体の積は別に 1 つ持ちます。先頭を消すときに前が空なら、後ろの値をその場で末尾の側から積に書き換えて、境を `tl` に移します。2 本のスタックで持つ形と違い、後ろから前へ値を写しません。配列の末尾に空きが無くなったら、生きている要素が半分以下なら配列の先頭へ詰め、そうでなければ 2 倍の配列へ移します。1 MB 以上の配列は 2 MB 境界に置いて huge page を頼みます (Linux のみ)。
