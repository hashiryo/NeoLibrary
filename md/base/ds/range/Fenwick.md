---
title: Fenwick 木
documentation_of: ../../neo/base/ds/range/Fenwick.hpp
---

1 点の更新と接頭辞の積を $O(\log n)$ で行う木です。

```cpp
template <class M> struct Fenwick {
  using T = typename M::T;
  explicit Fenwick(int n = 0);
  void add(int i, const T& x);  // i 番目 (0 始まり) の値を M::op で x と合わせる
  T prod(int r) const;          // [0, r) の積
};
```

`M` はセグメント木と同じモノイド (`using T`、`static T ti()`、`static T op(T a, T b)`) です。和にも最大にも使えます。逆元の要らないモノイド (最大、最小、gcd など) をそのまま渡せます。

`add` は、和のモノイドなら値を足し、最大のモノイドなら値を大きくします。`prod(r)` は接頭辞 `[0, r)` の積です。区間 `[l, r)` の積は、逆元のあるモノイドなら `prod(r)` と `prod(l)` から作れますが、この木はそこまでを持ちません。
