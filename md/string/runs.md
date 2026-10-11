---
title: run の列挙
documentation_of: ../../neo/string/runs.hpp
---

列のすべての run を求めます。run は、最小の周期 $t$ が長さの半分以下で、同じ周期のまま左右へ伸ばせない区間です。

```cpp
std::vector<std::array<int, 3>> runs(const std::string& s);
template <class T> std::vector<std::array<int, 3>> runs(const std::vector<T>& s);
```

返す列の要素 $(t, l, r)$ は、`s[l, r)` の最小の周期が $t$ で、$r - l \ge 2t$ を満たし、`s[l - 1, r)` も `s[l, r + 1)` も周期 $t$ を持たないことを表します。列は $(t, l)$ の昇順です。長さを $n$ とすると、run の数は $n$ 未満です。

- `runs(s)` は文字列を受け、文字を unsigned char として比べます。$O(n)$ です。
- `std::vector<T>` の版は値を座標圧縮してから解きます。`T` は `<` で比べられれば何でもよく、$O(n \log n)$ です。

Bannai、I、Inenaga、Nakashima、Takeda、Tsuruta (2017) の Runs 定理によります。どの run も、文字の順とその逆のどちらかで、最長の Lyndon 接頭辞として現れる長さ $t$ の根を持ちます。位置 $i$ の最長の Lyndon 接頭辞の終わりは、$i$ より後ろで、接尾辞の順位が $i$ より小さい最初の位置です。接尾辞配列を順位の小さい順に読み、順位が今より小さい位置の集合と大きい位置の集合 (FastSet) から、2 つの順序の終わりを同時に引きます。

根は前後へ伸ばし、伸ばした長さが周期の 2 倍以上なら run です。伸ばすときは前後の 8 byte を読み、排他的論理和で一致の長さを数えます。一様な文字列の根は、ほとんどが初めの 8 byte で決まります。同じ周期の前に見つけた run の中に収まる根は飛ばすので、周期ごとには始まりの順に見つかり、周期で数え分けるだけで並びます。文字を直接比べる長さは、$a b a^2 b a^3 b \cdots$ のような列で $O(n^{1.5})$ に膨らみます。そこで、比べた長さが全体で $32n$ を超えたら、LCE を LCP 配列と RangeMin で引く形に切り替えます。

接尾辞配列、2 つの Lyndon 配列、写した列、run の列は 1 つの作業領域に置きます。OS は、確保した領域の 4 KB ごとのページへ初めて書き込まれたところで、メモリを割り当てて 0 で埋めます。GitHub の x64 のランナー (EPYC 7763) では、これにページ 1 枚あたり 1.7 µs かかりました。そこで作業領域が 1 MB 以上なら、2 MB 境界に寄せて `madvise(MADV_HUGEPAGE)` を呼び、2 MB のページを使わせます (Linux のみ)。この手間は 8 分の 1 ほどになります。

procon-judge の yosupo-runenumerate ($n \le 2 \times 10^5$) で書き比べました。EPYC 7763 の x64-gcc では、ケースの組ごとの中央値で最も遅い組が 7.70 ms でした。同じ回に、今の Library の run_enumerate (Main–Lorentz の分割統治と Z 配列) は 78.77 ms でした。Neoverse-N2 の arm-gcc では 5.53 ms と 67.73 ms です。記録は algo-notes の notes/runs.md にあります。
