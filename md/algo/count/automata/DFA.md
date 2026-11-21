---
title: 有限オートマトンの上の DP
documentation_of: ../../neo/algo/count/automata/DFA.hpp
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
// 状態を見る口
template <class T> std::vector<std::pair<State, T>> count_states(int n) const;
template <class T, class A, class F> std::vector<std::pair<State, T>> dp_states(int n, A op, T e, F f, T init) const;
// 辞書順と k 番目、長さによらない数、長さが巨大なとき
template <class T> Ranker<T> ranker(int L) const;   // count(len)、count_less(w)、count_leq(w)、accepts(w)、kth(k, len)
bool is_finite() const;
template <class T> T count_any() const;
template <class T> T count_large(unsigned long long n) const;

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

## 状態を見る

dp の f は、4 つ目の引数に遷移の前後の状態 `const DFA<symbol_t>::Step& t` を受けてもかまいません。`t.to.get<S>(j)` は、`&` の左から j 番目の成分 (成分が 1 つなら 0) の、移った先の状態の値です。S は、その成分のコンストラクタに渡した状態の型です。`t.from` は移る前の状態で、`t.to[j]` は状態の値の代わりに成分の中での番号を返します。状態の値を引けるのは、遷移の関数で作った DFA の成分だけです (NFA と leq_digits などの成分は番号だけ)。Aho–Corasick のノードの得点や、KMP で一致した回数を足す問題に使います。

```cpp
// 長さ n の英小文字の列で、各位置で終わる語の得点 score[ノード] の合計の最大
DFA<char> ac(alp, 0, [&](int v, char c) { return go[v][c - 'a']; }, [](int) { return true; });
auto f= [&](long long x, char, int, const DFA<char>::Step& t) { return x + score[t.to.get<int>()]; };
long long best= ac.dp(n, [](long long a, long long b) { return std::max(a, b); }, LLONG_MIN, f, 0LL);
```

`count_states<T>(n)` と `dp_states(n, op, e, f, init)` は、長さ n の列で来られる状態ごとの数と dp の値を返します。受理されない状態も含み、`st.accepted()` で受理されるかが分かります。答えを最後の余りや LCS の長さごとに出す問題で、受理の判定を変えて何度も呼ばずに済みます。

```cpp
// N 以下の数を、桁和を 3 で割った余りごとに数える
DFA<int> sum3(alp, 0, [](int s, int c) { return (s + c) % 3; }, [](int) { return true; });
std::array<Mint, 3> by{};
for(auto& [st, c]: (leq_digits(N) & sum3).count_states<Mint>(N.size())) by[st.get<int>(1)]+= c;
```

状態を見る口は、積を初期状態から来られる組だけの 1 つの表にしてから進めます。位置を状態に持つ DFA と大きい成分の積では、表が長さと成分の状態の数の積の大きさになります。

## 辞書順と k 番目

`ranker<T>(L)` は、状態と残りの長さごとの受理される列の数を、後ろからの DP で前もって作った表を返します。大きさは (L + 1) × (積の来られる組の数) です。表からは、長さ L 以下の列について次のことを、長さ × 記号の数の手間で引けます。辞書順は、記号の集まりを並べ替えた順で比べます。

- `count(len)`: 受理される長さ len の列の数。
- `count_less(w)`: 受理される長さ |w| の列のうち、w より辞書順で小さいもの。`count_leq(w)` は w 以下のもの。
- `kth(k, len)`: 受理される長さ len の列のうち、辞書順で k 番目 (0 始まり) のもの。k が数以上なら空を返す。

T が整数型なら、数は型の最大で止めます。k 番目の数や、条件を満たす数のうち N 以上で最小のものを作る問題に使います。上限 N を変えて何度も数える問題では、leq_digits との積を毎回作る代わりに、ranker を 1 回作って N ごとに `count_leq` を呼びます。先頭の 0 を許した固定の長さの数字の列なら、辞書順は数の大小と同じです。

```cpp
// 4 と 7 だけからなる数のうち、7 の個数が偶数のものの、小さい方から k 番目 (k は 0 始まり、長さは len)
DFA<int> d({4, 7}, 0, [](int s, int c) { return s ^ (c == 7); }, [](int s) { return s == 0; });
auto R= d.ranker<unsigned long long>(18);
std::vector<int> w= R.kth(k, len);
```

## 長さによらない数と、長さが巨大なとき

`count_any<T>()` は、受理される列を長さによらず数えます。受理される列が有限個のとき (`is_finite()` が真) だけ使えます。部分列の自動機のように閉路の無い DFA で、長さごとに 1 歩ずつ進めずに、状態 × 記号の手間で数えます。

`count_large<T>(n)` は、受理される長さ n の列の数を、n が 10^18 のように大きくても求めます。T は体 (素数を法とする ZMod など) でなければなりません。最初の 2 × (積の来られる組の数) + 2 項を count_all で出し、Berlekamp–Massey で漸化式にしてから、n 項目を多項式の剰余の累乗で求めます。手間は、組の数を S として S^2 × 記号の数と S^2 log n なので、S が数百から千ほどまでに向きます。

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

構築では、初期状態から幅優先でたどって状態に番号を振り、状態の数と記号の数の積の大きさの表を作ります。状態の型が 64 bit に詰められるなら、番号は密な配列で引きます。詰めた値を整数の葉ごとの値に分け、葉ごとの値の範囲を 2 倍ずつ広げながら、葉の値を並べたものを添字にします。位置や余りを状態に持つ DFA は値が詰まっているので、ハッシュで散らすよりメモリを順に触れて速く済みます。配列が番号の数の 8 倍を超えるほどまばらなら [ハッシュ表](../data_structure/HashMap.md) で引き、番号の数が 2 倍になるたび、密な配列へ戻せるか見直します。詰められるのは、整数や bool と、それらの pair、tuple、array で bit の和が 64 以下のものです。ほかの型は要素のハッシュを混ぜる std::unordered_map で、ハッシュを作れない型だけ std::map を使います。NFA では、NFA の状態の番号を先に振り、ε 遷移で閉じた移り先を前もって求めてから部分集合構成をします。NFA の状態が 64 個以下なら、部分集合を u64 の bit で持ちます。

表、番号付けの配列、DP の置き場所の配列のうち 2 MB 以上のものは、Linux では 2 MB 境界に置いて huge page を頼みます。状態が数百万ある DFA の表は、倍々で伸ばすごとに新しい領域へ写され、そのページへ初めて触る手間が 4 KB のページでは大きいためです。reject との比べは、状態の型が 64 bit に詰められるなら詰めた値どうしで行います。

DP は 1 歩ずつ進めます。成分が 1 つなら、その歩で来た状態の並びと値を持ち、表を引いて次の歩の並びを作ります。成分が 2 つ以上なら、成分ごとに、その歩で来られる状態の集合をその成分だけで求めて 0 から番号を振ります。来た組の置き場所は、成分ごとの番号を 2 の冪に切り上げた幅の bit に並べた値を添字にして、密な配列で引きます。成分ごとに次の歩の番号を bit の位置までずらした表を前もって作るので、移る先の添字は成分ごとの表の値の和です。ハッシュ表も待ち行列も使いません。成分の集合が前の歩と同じになったら、以後はその成分を求め直しません。

集合が止まった成分は bit の下位に、まだ変わる成分は上位に並べます。下位の成分どうしの表の和は、まとめた添字ごとの 1 つの表に前もってまとめます。止まった成分の表は変わらないので、まとめ直すのは止まる成分が増えたときだけです。上位の成分どうしも、組の数が少なければ毎歩 1 つの表にまとめます。そのため、移る先を 1 つ求めるのに引く表はたいてい 2 つです。

「N 以下」の DFA は位置を状態に持つので、状態の数は桁の数に比例します。それでも、ある歩で来られるのは上限と等しい状態と小さいと決まった状態の 2 つだけなので、積の配列はほかの成分の状態の数の 2 倍で済みます。手で書く桁 DP の dp[i][tight][...] と同じ大きさです。今の Library の Automaton は積を幅優先で表にしていたので、位置を状態に持つ DFA との積では、表の大きさが桁の数とほかの成分の状態の数の積になっていました。

記号が 16 個以上の数え上げ (count と count_all) では、始めに各状態から同じ先へ移る記号をまとめてみます。移る先と記号の本数の組の数が遷移の数の 65% 以下なら、その組で進めます。同じ先へ多くの記号が集まると、その先の値への足し込みが 1 本の待ちの鎖になるためです。英小文字 26 個の列で、語が少ない Aho–Corasick の DFA の数え上げが 2〜7 倍速くなります。

添字の bit の和が 24 を超えたら、bit の多い 2 つの成分を 1 つにまとめ、最初からやり直します。まとめた成分は、来られる組だけを表にした積です。成分どうしが連動していて、組の大半には来られないときに起きます。

長さが 32 以上なら、成分を並べたまま進める前に、積を来られる組だけの表にできるか試します。積を幅優先で段ごとに作りながら、後ろ半分の段の組の増え方から、DP で表の 1 行を何歩使い回すかを見積もり、32 回より少なければやめます。「N 以下」の DFA と小さい成分の積では、小さいと決まった側の組が毎歩同じなので、使い回しが多くなります。表にすると、移る先を 1 つ求めるのに引く表が 1 つで済み、来た組の置き場所を引く配列も来られる組の数の大きさになります。位置を状態に持つ DFA のどの位置にも多くの組が並ぶ積は、組が毎歩新しくなるので表にしません。

procon-judge で、今の Library の Automaton を使う 24 問にこの DFA を使う提出を置いて比べました。次の時間は x64-gcc の最大のケースで、今の Library の値は同じ CPU モデルで測った別の回の記録です。

| 問題 | 中身 | CPU | 今の Library (ms) | DFA (ms) |
| --- | --- | --- | ---: | ---: |
| abc155-e | 位置と繰り下がりの組を状態に持つ 200 万状態の DFA で、最小を求める DP | EPYC 9V45 | 768 | 137 |
| abc235-f | 「N 以下」(10^4 桁) と使った数字の集合の積で、値の組の和 | Xeon 6975P-C | 620 | 337 |
| yuki-315 | 「N 以下」(20 万桁) と mod 3 と mod 8 の積 | Xeon 6975P-C | 578 | 331 |
| abc135-d | 位置を状態に持つ 10 万状態の DFA と mod 13 の積 | EPYC 9V74 | 277 | 36 |
| yuki-911 | 状態が長さ 20 万の vector を含む DFA と、「L 以上」と「R 以下」の積 | Xeon 6975P-C | 85 | 68 |
| yuki-362 | 長さを変えて何度も数える (count_all の 1 回にまとめた) | EPYC 9V45 | 97 | 9 |
| yuki-1740 | 「S より辞書順で小さい」と 'a' の数の積 | EPYC 9V74 | 142 | 95 |
| aoj-0570 | 「A 以上」と「B 以下」と mod M とジグザグの判定の積 | EPYC 9V45 | 65 | 43 |

積を表にする道と番号付けの密な配列を入れる前 (73d54cf) は、積の来られる組が 9000 ほどの aoj-0570 で、積を表にする今の Library の形と並んでいました (EPYC 9V45 で 62 ms)。clang の環境では、yuki-1740 と aoj-0570 が今の Library より 1 割から 3 割遅くなっていましたが、入れたあとはどちらも速くなりました。abc155-e は番号付けの分が大きく、入れる前は Xeon 6975P-C で 554 ms でした。ただし入れたあとも、x64-gcc では x64-clang の 1.4〜1.6 倍かかっていました。gcc は std::array の == (libstdc++ では memcmp) のために、遷移の先の状態を 4 byte ずつ書いた直後に 8 byte で読み直していて、x86 では読み込みが書き込みの完了を待つためです。この比べを詰めた値どうしにし、大きい配列に huge page を頼むようにしたあとは、EPYC 9V45 の x64-gcc で 366 ms から 137 ms になりました。記録は algo-notes の notes/automaton-dp.md にあります。
