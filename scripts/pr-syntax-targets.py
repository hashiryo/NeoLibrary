#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# ///
"""PR で変わった neo/ のヘッダを 1 行 1 ファイルで出す。

PR の軽い検査 (.github/workflows/pr-check.yml) が、ここで出たヘッダを 1 環境で
-fsyntax-only でコンパイルする。テストデータもトークンも要らず、数分で終わる。
CE は拾えるが WA は拾えない。WA はマージ後に procon-judge が測って分かる。

ヘッダ単体のコンパイルはテンプレートを実体化しないので弱いが、構文と非テンプレート
の誤りは拾える。テンプレートまで実体化する検証は procon-judge の提出が担う
(test/ は持たない。Library が 2026-09-22 に test/ を消した形を引き継ぐ)。
"""

import argparse
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def changed_files(base: str, head: str) -> set[str]:
    """base から head までに変わったファイル (ROOT からの相対)。消したものは含めない。"""
    proc = subprocess.run(
        ["git", "diff", "--name-only", "--diff-filter=d", f"{base}...{head}"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    return {line.strip() for line in proc.stdout.splitlines() if line.strip()}


def select_targets(changed: set[str]) -> list[str]:
    """変わったファイルのうち neo/ のヘッダ。"""
    return sorted(f for f in changed if f.startswith("neo/") and f.endswith(".hpp"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("base")
    parser.add_argument("head", nargs="?", default="HEAD")
    args = parser.parse_args()
    for target in select_targets(changed_files(args.base, args.head)):
        print(target)


if __name__ == "__main__":
    main()
