# NeoLibrary のディレクトリ構成

ヘッダはすべて `neo/` の下に置きます（リポジトリ直下に `.hpp` を置かない）。その `neo/` を、次の 2 層に分けています。

- **`neo/base/`** … 問題の種類（最適化・数え上げなど）に依存しない部品
- **`neo/algo/`** … 問いの性質ごとの「終着点」になるアルゴリズム

説明用の Markdown は **`md/`** に、ヘッダと同じ相対パスで置きます（例: `neo/base/core/algebra/ZMod.hpp` と `md/base/core/algebra/ZMod.md`）。

## 依存の向き

**`algo/` から `base/` への `#include` だけ許可します。** `base/` のヘッダが `algo/` を読むことは禁止です。

## `#include` の書き方

- 形式: `#include "neo/base/..."` または `#include "neo/algo/..."`
- 探索パス: リポジトリルート（procon-judge では `-Ineolib`）

旧い `neo/data_structure/...` や `neo/combinatorial_optimization/...` といったトップカテゴリは使いません。

## `neo/base/`

| ディレクトリ | 入れるもの |
| --- | --- |
| **`core/`** | 数学的な対象・演算の部品。`algebra/`、`linalg/`（型本体）、`nt/`、`poly/`、`set_fps/`（集合冪級数の核）、`util/`（圧縮・ソートなど） |
| **`ds/`** | データ構造・索引。`range/`（セグ木・Wavelet など）、`set/`、`hash/`、`seq/`、`string/`、`tree/`（HLD・LCA など）、`2d/` |
| **`engine/`** | データ構造ではなく、処理の枠組み。例: 全方位 DP（`engine/tree/Rerooting.hpp`） |

`core/`、`ds/`、`engine/` のさらに下は、ドメイン名ではなく **扱う対象**（`math` に相当する `algebra` / `nt`、`graph` に相当する `tree` など）で切ります。ファイルが 1〜2 個しかないところでは、深くフォルダを増やさなくて構いません。

## `neo/algo/`

| ディレクトリ | 入れるもの |
| --- | --- |
| **`optimize/`** | 最大・最小、フロー、マッチング、カット、マトロイド、燃やす埋める系など |
| **`count/`** | 場合の数、数え上げ、Floor Sum、オートマトン DP、文字列の数え上げなど |
| **`decide/`** | 判定・構造解析（SCC、2-SAT、連結成分、二部性など） |
| **`compute/`** | 単体で答えを出す重量級の計算（因数分解、離散対数、行列演算の求解、多項式因数分解など） |

下には `flow/`、`burn_bury/`（燃やす埋める：`BinaryOptimization` など）、`cover/`（最小パス被覆・反鎖など）、`matching/` など、手法・対象ごとのサブディレクトリを置きます。

## 境界の例（迷ったとき）

- **セグ木・接尾辞配列・Rolling Hash** → `base/ds/`
- **Dijkstra・最大流・凸包チャート** → `algo/optimize/`（相当するサブディレクトリへ）
- **Z algorithm・回文前処理** → 索引として `base/ds/string/` に置くことが多い
- **集合冪級数**: ゼータ・畳み込みなどの核は `base/core/set_fps/`。グラフ上の閉路積やハフニアンなど、問題文に近い終着点は `algo/count/` や `algo/compute/` 側に置く

新規ヘッダを足すときは、上の表に当てはめ、依存が `algo → base` だけになるようにしてください。

## 旧構成からの対応

2026-11 に、旧 14 カテゴリ（`data_structure`、`graph`、`combinatorial_optimization` など）から上記ツリーへ移しました。機械的な旧→新の対応は `scripts/reorg_includes.py` の `MAPPING` にあります。