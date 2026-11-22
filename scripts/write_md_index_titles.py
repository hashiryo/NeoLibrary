#!/usr/bin/env python3
"""Create md/**/_index.md with Japanese sidebar titles."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent / "md"

# (relative path from md/, title, order)
ENTRIES: list[tuple[str, str, int]] = [
    ("base", "汎用基盤", 10),
    ("base/ds", "データ構造", 10),
    ("base/ds/2d", "二次元", 50),
    ("base/ds/hash", "ハッシュ", 40),
    ("base/ds/range", "区間・配列", 10),
    ("base/ds/range/segtree", "セグメント木の標準の構造体", 1),
    ("base/ds/seq", "列", 30),
    ("base/ds/set", "集合", 20),
    ("base/ds/string", "文字列", 60),
    ("base/ds/tree", "木", 70),
    ("base/core", "数学と代数", 20),
    ("base/core/algebra", "代数", 10),
    ("base/core/linalg", "線形代数", 20),
    ("base/core/nt", "数論", 30),
    ("base/core/poly", "多項式", 40),
    ("base/core/set_fps", "集合冪級数", 50),
    ("base/core/util", "ユーティリティ", 60),
    ("base/engine", "計算の枠組み", 30),
    ("base/engine/tree", "木", 10),
    ("algo", "問別アルゴリズム", 20),
    ("algo/optimize", "最適化", 10),
    ("algo/optimize/flow", "フロー", 10),
    ("algo/optimize/matching", "マッチング", 20),
    ("algo/optimize/cut", "カット", 30),
    ("algo/optimize/matroid", "マトロイド", 40),
    ("algo/optimize/burn_bury", "燃やす埋める", 50),
    ("algo/optimize/cover", "被覆", 55),
    ("algo/count", "数え上げ", 20),
    ("algo/count/automata", "オートマトン", 10),
    ("algo/count/graph", "グラフ", 20),
    ("algo/count/nt", "数論", 30),
    ("algo/count/set_fps", "集合冪級数", 40),
    ("algo/count/string", "文字列", 50),
    ("algo/decide", "判定", 30),
    ("algo/decide/graph", "グラフ", 10),
    ("algo/compute", "計算", 40),
    ("algo/compute/linalg", "線形代数", 10),
    ("algo/compute/nt", "数論", 20),
    ("algo/compute/poly", "多項式", 30),
    ("algo/compute/set_fps", "集合冪級数", 40),
]

TEMPLATE = """---
title: {title}
order: {order}
---
"""


def main() -> None:
    for rel, title, order in ENTRIES:
        path = ROOT / rel / "_index.md"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(TEMPLATE.format(title=title, order=order), encoding="utf-8")
        print("wrote", path.relative_to(ROOT.parent))


if __name__ == "__main__":
    main()