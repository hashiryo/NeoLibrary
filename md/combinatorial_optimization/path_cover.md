---
title: 最小パス被覆と Dilworth の定理
documentation_of: ../../neo/combinatorial_optimization/path_cover.hpp
---

DAG の頂点を道で覆う問題を 2 つの意味で解きます。頂点を共有しない道で覆う最小パス被覆と、到達可能性の順序 (u から v へ道があれば u は v より前) を鎖に分ける最小の鎖の被覆です。後者の最小の本数は、最大反鎖の大きさと等しくなります (Dilworth の定理)。

```cpp
std::vector<std::vector<int>> min_path_cover(int n, const std::vector<std::array<int, 2>>& edges);
std::vector<std::vector<int>> min_chain_cover(int n, const std::vector<std::array<int, 2>>& edges);
std::vector<int> max_antichain(int n, const std::vector<std::array<int, 2>>& edges);
```

頂点は 0 から n - 1 の番号で、`edges` は有向辺 {u, v} (u から v) の列です。閉路と自己ループが無いことを前提にします。多重辺はあってもかまいません。

- `min_path_cover` は、頂点を共有しない道で全部の頂点を覆い、本数が最小の道の列を返します。道は頂点の列で、隣り合う頂点は辺で結ばれています。長さ 0 の道 (頂点 1 つ) も入ります。
- `min_chain_cover` は、到達可能性の順序で全部の頂点を鎖に分け、本数が最小の鎖の列を返します。鎖は頂点の列で、前の頂点から後の頂点へ道があります (辺で結ばれているとは限りません)。どの頂点もちょうど 1 つの鎖に入ります。
- `max_antichain` は、同じ順序の最大反鎖 (どの 2 つの間にも道が無い頂点の集合) を、番号の昇順で返します。大きさは `min_chain_cover` の本数と等しくなります。

道と鎖は、先頭の頂点の番号の順に並びます。

## どれを使うか

頂点を共有しない道が要る問題 (LOJ 6002、CF 1906I など) は `min_path_cover` です。道が頂点を共有してよい問題や、Dilworth の定理で反鎖に直す問題 (ABC 374 G、POJ 2594、洛谷 P4298 など) は `min_chain_cover` か `max_antichain` です。こちらは推移閉包を作らずに、元の DAG の辺のまま渡せます。

関係がもともと推移的 (u から v へ道があれば辺 {u, v} もある) なら、3 つの本数は同じで、`min_path_cover` の道はそのまま鎖になります。速さは渡す辺で決まります。比べられる組を全部辺で渡すなら、`min_path_cover` のほうが 2〜3 倍速くなります。Hasse 図の辺 (間に頂点を挟まない組) だけを渡せるなら、辺が少ないぶん `min_chain_cover` と `max_antichain` に渡すほうが速くなります。

## 中身

`min_path_cover` は、頂点 v を出る側と入る側に分けた二部グラフの最大マッチングを `BipartiteMatching` で求めます。辺 {u, v} は、u の出る側と v の入る側を結びます。u の出る側と w の入る側が組になっていれば、道の中で u の次を w にします。本数は n - (マッチングの大きさ) です。

`min_chain_cover` と `max_antichain` は、推移閉包の上の最大マッチングを、閉包を作らずに `MaxFlow` を 1 回流して求めます。網は次の形です。

- 頂点 v を、出る側 v' と入る側 v'' に分ける。
- s から v' へ容量 1、v'' から t へ容量 1 の辺を張る。
- 辺 {u, v} ごとに u' から v'' へ、頂点ごとに v'' から v' へ (v を通り抜ける) 容量 n の辺を張る。

u' から入って w'' から出る流れの 1 単位は、u から w への道をたどるので、閉包の上の組 (u, w) と見なせます。流れを 1 単位ずつ s から t へたどって組を取り出し、u の次を w にして鎖を作ります。最大反鎖は最小カットから取ります。v' が s の側、v'' が t の側にある頂点 v の集合です。この集合の 2 つの頂点の間に道があると、その道の上のどこかで容量 n の辺がカットをまたぎ、カットの容量が n 以上になります。鎖は 1 本以上あるので最大流は n - 1 以下で、最小カットではそうならず、集合は反鎖になります。

頂点を通り抜ける辺のため、鎖が多く深い疎な DAG では `MaxFlow` が遅くなります。手元の M2 で、頂点 20 万、辺 50 万の乱択の疎な DAG に 0.15〜2.5 s かかりました (同じ入力で `min_path_cover` は 30 ms ほど)。頂点が数千までの DAG なら数 ms で済み、推移閉包を作ってから `min_path_cover` へ渡すより速くなります (頂点 2000、辺 4 万で 3.5 ms。閉包を作る形は 30 ms)。記録は algo-notes の notes/network-flow.md にあります。
