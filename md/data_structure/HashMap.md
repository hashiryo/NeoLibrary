---
title: ハッシュ表 (u64 のキー)
documentation_of: ../../neo/data_structure/HashMap.hpp
---

64 bit の符号なし整数のキーから値への連想配列です。

```cpp
template <class V> class HashMap;
HashMap();
int size() const;
bool empty() const;
bool contains(unsigned long long key) const;
V get(unsigned long long key, V d = V()) const;  // 無ければ d
V& operator[](unsigned long long key);          // 無ければ V() で入れる
void set(unsigned long long key, V v);
bool erase(unsigned long long key);             // あったら true
void clear();
template <class F> void for_each(F&& f) const;  // f(key, value) を入っている数だけ呼ぶ。順番は決まらない
void swap(HashMap& o);
```

キーは u64 の全域を使えます。int や long long のキーは u64 に直して渡します (負の値も別のキーになります)。`V` は、0 の bit 列が `V()` と等しい、trivially copyable な型に限ります (整数、浮動小数点数、それらの pair や配列など)。コピーはできず、move と swap はできます。1 回の操作は期待値で $O(1)$ です。

キーの順序 (前後の要素や k 番目) が要るなら [OrderedSet](OrderedSet.md) を使います。キーが 0 以上 n 未満の整数 (n は $10^7$ ほどまで) で、入っているかと前後の要素だけを見るなら、[FastSet](FastSet.md) の方が速く動きます。

## 中身

開番地法の線形探索で、{キー, 値} の組を 2 の冪の個数だけ並べ、埋まりが半分を超えたら 2 倍にします。キーは bit を反転して持ち、0 を空きの印にするので、calloc や mmap で取った 0 の領域がそのまま空の表になります。反転すると 0 になるキー ($2^{64} - 1$) だけは表の外に置きます。除くときは、後ろに続く組のうち動かせるものを詰め直すので、墓標を置かず、消したあとも探す長さが延びません。

位置は、表を作るときに steady_clock から取った乱数の種をキーに足してから splitmix64 で混ぜ、その上位 bit で決めます。乱数を掛けるだけの軽いハッシュは、連続する整数や等差数列のキーで、まれに 1 個あたり数百組を探す種を引くので使いません。2 MB 以上の表は、Linux では 2 MB 境界に取って huge page を頼みます。ばらばらの位置を引く表なので、TLB の外れが減ります。

方式は procon-judge の yosupo-associative-array で比べて選びました ($Q = 10^6$、最大は 99 万個の新しいキーを入れ続けるケース)。同じ方式の lp_hp は、EPYC 9V74 の x64-gcc で最大 32 ms ほどで、乱数を混ぜるハッシュを渡した std::unordered_map の 8.6 分の 1、std::map の 22 分の 1 でした。Swiss table の形は、全ケースの合計では 2 割近く速かったものの、最大のケースの順は CPU で入れ替わり、消去に墓標が要るので、短く書ける線形探索にしました。記録は algo-notes の notes/associative-containers.md にあります。
