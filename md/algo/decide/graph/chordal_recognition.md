---
title: 弦グラフの判定
documentation_of: ../../neo/algo/decide/graph/chordal_recognition.hpp
---

無向グラフが弦グラフ (長さ 4 以上の誘導閉路を持たない) かどうかを判定し、弦グラフなら完全除去順序 (PEO) を、そうでなければ長さ 4 以上の誘導閉路を返します。

```cpp
struct ChordalRecognition {
 bool chordal;
 std::vector<int> order;  // 弦グラフのときの PEO (長さ n)
 std::vector<int> cycle;  // 弦グラフでないときの誘導閉路 (巡回順、始点は繰り返さない)
};
ChordalRecognition chordal_recognition(int n, const std::vector<std::array<int, 2>>& edges);
```

頂点は 0 から n - 1 の番号で、`edges` は単純無向グラフの辺の列です (自己ループと多重辺は未定義)。最大基数探索 (MCS) を遅延バケットで行うので、頂点と辺の数の和に比例する時間で動きます。PEO の検査は、各頂点 x について「x より後ろの隣接のうち最も早いもの」1 本との比較で足ります。非弦のときの閉路は、MCS の証人 (x, z, y) から x と N(x)\{z} を避けた幅優先で作ります。

`order` は問題文の PEO の向き (そのまま出力できる順) です。`cycle` は巡回順で、始点を末尾に繰り返しません。
