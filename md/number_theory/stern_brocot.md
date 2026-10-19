---
title: Stern-Brocot 木
documentation_of: ../../neo/number_theory/stern_brocot.hpp
---

Stern-Brocot 木の節点 (正の既約分数) について、根 $1/1$ からの道、道の行き先、LCA、祖先、子孫の区間を求めます。

```cpp
unsigned* sb_encode(unsigned a, unsigned b, unsigned* out);
std::pair<unsigned, unsigned> sb_decode(const unsigned* first, const unsigned* last);
std::pair<unsigned, unsigned> sb_lca(unsigned a, unsigned b, unsigned c, unsigned d);
std::pair<unsigned, unsigned> sb_ancestor(unsigned long long k, unsigned a, unsigned b);
std::array<unsigned, 4> sb_range(unsigned a, unsigned b);
```

分数 $a/b$ は $1 \le a, b < 2^{32}$ で既約であることが前提です。道は連の長さを R、L、R、... の順に並べて表し、最初の R の連だけ 0 になりえます ($a/b < 1$ のとき)。たとえば $2/3$ への道 LR は $(0, 1, 1)$ で、$1/1$ への道は $(0)$ です。

- `sb_encode(a, b, out)` は、$a/b$ への道を out に書き、書き終わりを返します。連の数は 46 以下で、out は 46 個ぶんあれば足ります。
- `sb_decode(first, last)` は、first から last の手前までの道の行き先を $(a, b)$ で返します。答えは $2^{32}$ 未満であることが前提です。
- `sb_lca(a, b, c, d)` は、$a/b$ と $c/d$ の最も深い共通の祖先を返します。
- `sb_ancestor(k, a, b)` は、$a/b$ の祖先のうち深さ $k$ のものを返します ($1/1$ の深さは 0、$a/b$ 自身も祖先に含みます)。無ければ $(0, 0)$ を返します。
- `sb_range(a, b)` は、$a/b$ の子孫の集まりである開区間の下端 $f/g$ と上端 $h/k$ を $\{f, g, h, k\}$ で返します。$0$ を $0/1$、$+\infty$ を $1/0$ で表します。

どれも $a/b$ の連分数展開 (互除法の商の列) で求めます。$a/b$ への道は連分数の項を並べ、最後の項から 1 を引いたものです。道の途中の節点は近似分数から作ります。互除法は 32 bit の割り算で回し、1 段の待ちは割り算 1 回ぶんです。

- sb_encode は、互除法の 1 段 $a = qb + r$ で $r > b/2$ なら次の段の商が 1 なので、次の段を $(b, b - r)$ にして 2 段を 1 回の割り算で済ませます。一様な分数では割り算が 7 割ほどに減ります。
- sb_lca は 2 つの互除法を同じ段で進め、連が初めて食い違うところで止めます。
- sb_ancestor は、残りの深さが互除法の組 $(x, y)$ の $x$ 以上になったら、祖先が無いと分かるので止めます。
- sb_range は、$a g - b f = 1$ から $g = a^{-1} \bmod b$ ($b$ が偶数なら $a$ を法にした逆元) を割り算を使わない Stein の方法で求め、もう一方は割り切れる割り算で出します。
- sb_decode は、x86 では vpmuludq で掛けます。

手間は互除法の段の数で、$O(\log \max(a, b))$ です。procon-judge の yosupo-stern-brocot-tree ($10^5$ 問) では、x64 の CI のランナーで最大 5 ms から 7 ms ほどでした。yosupo の最速提出の写しの 2.8 倍ほどの速さです。
