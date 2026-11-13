---
title: 有限オートマトンの上の DP
documentation_of: ../../neo/counting/DFA.hpp
---

遷移の関数で表した有限オートマトン (DFA や NFA) について、受理される長さ $n$ の記号の列を数えたり、列の上で DP をしたりします。桁 DP の多くは、条件ごとの DFA の積の上の DP として書けます。

```cpp
template <class symbol_t> class DFA;
DFA(const std::vector<symbol_t>& alphabet, const state_t& init, F tr, G ac);
DFA(const std::vector<symbol_t>& alphabet, const state_t& init, F tr, G ac, const state_t& reject);
DFA(const std::vector<symbol_t>& alphabet, const state_t& init, F tr, G ac);         // tr が状態の集まりを返す NFA
DFA(const std::vector<symbol_t>& alphabet, const state_t& init, F tr, G ac, H eps);  // ε 遷移のある NFA
DFA operator&(const DFA& r) const;
template <class S2, class G> DFA<S2> map_symbols(const std::vector<S2>& alphabet, G g) const;
const std::vector<symbol_t>& alphabet() const;
template <class T> T count(int n) const;
template <class T> std::vector<T> count_all(int n) const;
template <class T, class A, class F> T dp(int n, A op, T e, F f, T init) const;
template <class T, class A, class F> std::vector<T> dp_all(int n, A op, T e, F f, T init) const;

DFA<int> leq_digits(const std::string& N);
DFA<int> leq_digits(const std::vector<int>& digits, int base);
DFA<int> leq_digits(unsigned long long N, int base, int len);
DFA<int> leq_digits_lsb(unsigned long long N, int base);
// geq_digits と geq_digits_lsb も同じ引数
```

コンストラクタは、記号の集まり alphabet、初期状態 init、遷移の関数 `tr(状態, 記号)`、受理の判定 `ac(状態)` を受けます。初期状態から幅優先にたどり、来られる状態に番号を振って表にします。来られる状態が有限でなければ終わりません。状態の型は `==` で比べられるもので、ハッシュを作れる型 (整数、bool、それらの pair、tuple、array、vector、set、string など) か、`<` で比べられる型です。`reject` を渡すと、reject に移る遷移を捨てます。reject からはどの列も受理されないときに渡すと、DP で回る状態が減ります。

`tr` が状態の集まり (std::set や std::vector) を返すときは NFA として受け、部分集合構成で DFA に直します。`eps(状態)` を渡すと、ε 遷移で移れる状態の集まりとして扱います。部分集合の数は、NFA の状態の数の指数で増えることがあります。

`a & b` は、a と b の両方に受理される列を受理する DFA です。記号の集まりが同じものどうしでだけ取れます。表は作らずに成分の DFA を並べて持つだけなので、何度取っても軽く済みます。

`map_symbols(alphabet, g)` は、記号 c を g(c) と読む DFA を返します。記号の集まりは alphabet になり、g(c) がもとの記号の集まりに無ければ、c の遷移を捨てます。2 つの数の桁の組を記号にする問題で、片方の桁だけを見る DFA を作るときに使います。

`count<T>(n)` は受理される長さ n の列の数を返し、`count_all<T>(n)` は長さ 0 から n までのそれぞれの数を返します。`dp(n, op, e, f, init)` は、長さ 0 の列に値 init を置き、列を 1 記号ずつ読むたびに値を `f(値, 記号, 位置)` で進め、同じ状態に来た値を op でまとめます。返すのは、受理される状態の値を op でまとめたものです。位置は 0 から数えます。op は可換なモノイドの演算で、e はその単位元です。まとめる順番は決まっていません。`dp_all` は長さ 0 から n までの答えを返します。長さを変えて何度も数えるときは、`count_all` か `dp_all` を 1 回呼ぶ方が速く済みます。記号の集まりは中で並べ替えるので、`f` に渡る記号の順番も並べ替えたあとのものです。

## 「N 以下」の DFA

`leq_digits(N)` は、上の桁から読んだ |N| 桁の数字の列のうち、10 進の文字列 N の表す数以下のものを受理します。記号は int の 0 から 9 です。`leq_digits(digits, base)` は、上の桁から並べた base 進の桁の列で N を受けます。`leq_digits(N, base, len)` は、整数 N を base 進の len 桁で受けます。どれも先頭の 0 を許し、長さが桁の数に等しい列についてだけ意味を持ちます。

`leq_digits_lsb(N, base)` は、下の桁から読んだ列のうち、値が N 以下のものを受理します。こちらは列の長さを問いません。`geq_digits` と `geq_digits_lsb` は、以下を以上にしたものです。

次の例は、N 以下の正の整数のうち、0 でない桁がちょうど K 個 (K ≥ 1) のものを数えます。

```cpp
std::vector alp= {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
DFA nonzero(alp, 0, [](int s, int c) { return s + (c != 0); }, [&](int s) { return s == K; }, K + 1);
long long ans= (leq_digits(N) & nonzero).count<long long>(N.size());
```

## 中身

構築では、初期状態から幅優先でたどって状態に番号を振り、状態の数と記号の数の積の大きさの表を作ります。番号を引く表は、状態の型が 64 bit に詰められるなら [ハッシュ表](../data_structure/HashMap.md) です。詰められるのは、整数や bool と、それらの pair、tuple、array で bit の和が 64 以下のものです。ほかの型は要素のハッシュを混ぜる std::unordered_map で、ハッシュを作れない型だけ std::map を使います。NFA では、NFA の状態の番号を先に振り、ε 遷移で閉じた移り先を前もって求めてから部分集合構成をします。NFA の状態が 64 個以下なら、部分集合を u64 の bit で持ちます。

DP は 1 歩ずつ進めます。成分が 1 つなら、その歩で来た状態の並びと値を持ち、表を引いて次の歩の並びを作ります。成分が 2 つ以上なら、成分ごとに、その歩で来られる状態の集合をその成分だけで求めて 0 から番号を振ります。来た組の置き場所は、成分ごとの番号を 2 の冪に切り上げた幅の bit に並べた値を添字にして、密な配列で引きます。成分ごとに次の歩の番号を bit の位置までずらした表を前もって作るので、移る先の添字は成分ごとの表の値の和です。ハッシュ表も待ち行列も使いません。成分の集合が前の歩と同じになったら、以後はその成分を求め直しません。

集合が止まった成分は bit の下位に、まだ変わる成分は上位に並べます。下位の成分どうしの表の和は、まとめた添字ごとの 1 つの表に前もってまとめます。止まった成分の表は変わらないので、まとめ直すのは止まる成分が増えたときだけです。上位の成分どうしも、組の数が少なければ毎歩 1 つの表にまとめます。そのため、移る先を 1 つ求めるのに引く表はたいてい 2 つです。

「N 以下」の DFA は位置を状態に持つので、状態の数は桁の数に比例します。それでも、ある歩で来られるのは上限と等しい状態と小さいと決まった状態の 2 つだけなので、積の配列はほかの成分の状態の数の 2 倍で済みます。手で書く桁 DP の dp[i][tight][...] と同じ大きさです。今の Library の Automaton は積を幅優先で表にしていたので、位置を状態に持つ DFA との積では、表の大きさが桁の数とほかの成分の状態の数の積になっていました。

添字の bit の和が 24 を超えたら、bit の多い 2 つの成分を 1 つにまとめ、最初からやり直します。まとめた成分は、来られる組だけを表にした積です。成分どうしが連動していて、組の大半には来られないときに起きます。

記録は algo-notes の notes/automaton-dp.md にあります。
