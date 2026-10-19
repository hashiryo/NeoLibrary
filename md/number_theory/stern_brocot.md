---
title: Stern-Brocot 木
documentation_of: ../../neo/number_theory/stern_brocot.hpp
---

Stern-Brocot 木の節点 (正の既約分数) について、根 $1/1$ からの道、道の行き先、LCA、祖先、子孫の区間を求めます。分子と分母に上限のある分数での最良の近似も求めます。

```cpp
unsigned* sb_encode(unsigned a, unsigned b, unsigned* out);
std::pair<unsigned, unsigned> sb_decode(const unsigned* first, const unsigned* last);
std::pair<unsigned, unsigned> sb_lca(unsigned a, unsigned b, unsigned c, unsigned d);
std::pair<unsigned, unsigned> sb_ancestor(unsigned long long k, unsigned a, unsigned b);
std::array<unsigned, 4> sb_range(unsigned a, unsigned b);
std::array<unsigned, 4> sb_approx(unsigned n, unsigned x, unsigned y);
```

分数 $a/b$ は $1 \le a, b < 2^{32}$ で既約であることが前提です。道は連の長さを R、L、R、... の順に並べて表し、最初の R の連だけ 0 になりえます ($a/b < 1$ のとき)。たとえば $2/3$ への道 LR は $(0, 1, 1)$ で、$1/1$ への道は $(0)$ です。

- `sb_encode(a, b, out)` は、$a/b$ への道を out に書き、書き終わりを返します。連の数は 46 以下で、out は 46 個ぶんあれば足ります。
- `sb_decode(first, last)` は、first から last の手前までの道の行き先を $(a, b)$ で返します。答えは $2^{32}$ 未満であることが前提です。
- `sb_lca(a, b, c, d)` は、$a/b$ と $c/d$ の最も深い共通の祖先を返します。
- `sb_ancestor(k, a, b)` は、$a/b$ の祖先のうち深さ $k$ のものを返します ($1/1$ の深さは 0、$a/b$ 自身も祖先に含みます)。無ければ $(0, 0)$ を返します。
- `sb_range(a, b)` は、$a/b$ の子孫の集まりである開区間の下端 $f/g$ と上端 $h/k$ を $\{f, g, h, k\}$ で返します。$0$ を $0/1$、$+\infty$ を $1/0$ で表します。
- `sb_approx(n, x, y)` は、分子と分母がどちらも $n$ 以下の既約分数で、$x/y$ を下と上から挟む最良の近似を $\{a, b, c, d\}$ で返します。$a/b$ は $x/y$ 以下で最大のもの (無ければ $0/1$)、$c/d$ は $x/y$ 以上で最小のもの (無ければ $1/0$) です。$1 \le n, x, y < 2^{32}$ を受け、$x/y$ は既約でなくてかまいません。

どれも $a/b$ の連分数展開 (互除法の商の列) で求めます。$a/b$ への道は連分数の項を並べ、最後の項から 1 を引いたものです。道の途中の節点は近似分数から作ります。互除法は 32 bit の割り算で回し、1 段の待ちは割り算 1 回ぶんです。

- sb_encode は、互除法の 1 段 $a = qb + r$ で $r > b/2$ なら次の段の商が 1 なので、次の段を $(b, b - r)$ にして 2 段を 1 回の割り算で済ませます。一様な分数では割り算が 7 割ほどに減ります。
- sb_lca は 2 つの互除法を同じ段で進め、連が初めて食い違うところで止めます。
- sb_ancestor は、残りの深さが互除法の組 $(x, y)$ の $x$ 以上になったら、祖先が無いと分かるので止めます。
- sb_range は、$a g - b f = 1$ から $g = a^{-1} \bmod b$ ($b$ が偶数なら $a$ を法にした逆元) を割り算を使わない Stein の方法で求め、もう一方は割り切れる割り算で出します。
- sb_decode は、x86 では vpmuludq で掛けます。
- sb_approx は、$x/y$ への道を連ごとに下り、今いる節点が $n$ を超えたところの端を返します。連を全部進めても $n$ 以下かは掛け算と比べで分かるので、超える連でだけ割り算で歩数を出します。整数部分の次からは、商 2 以下の段を引き算と cmov で進め、商 3 以上が出たら残りを割り算で進めます。隣り合う Fibonacci 数の比のように商の小さい分数を、割り算なしで速く下るためです。

手間は互除法の段の数で、$O(\log \max(a, b))$ です。procon-judge の yosupo-stern-brocot-tree ($10^5$ 問) では、x64 の CI のランナーで最大 5 ms から 7 ms ほどでした。yosupo の最速提出の写しの 2.8 倍ほどの速さです。sb_approx は yosupo-rational-approximation ($10^5$ 問) で最大 9 ms から 13 ms ほどで、yosupo の最速提出の写しの 3 倍ほどの速さです。
