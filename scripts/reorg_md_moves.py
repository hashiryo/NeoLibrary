#!/usr/bin/env python3
import importlib.util
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("reorg", ROOT / "scripts" / "reorg_includes.py")
mod = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(mod)


def tracked(p: Path) -> bool:
    return (
        subprocess.run(
            ["git", "ls-files", "--error-unmatch", str(p.relative_to(ROOT))],
            cwd=ROOT,
            capture_output=True,
        ).returncode
        == 0
    )


def main() -> None:
    ok = skip = 0
    for old, new in sorted(mod.MAPPING.items(), key=lambda x: len(x[0]), reverse=True):
        src = ROOT / "md" / old.replace(".hpp", ".md")
        dst = ROOT / "md" / new.replace(".hpp", ".md")
        if not src.exists():
            skip += 1
            continue
        if src.resolve() == dst.resolve():
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        rel_src = src.relative_to(ROOT)
        rel_dst = dst.relative_to(ROOT)
        if tracked(src):
            subprocess.run(["git", "mv", str(rel_src), str(rel_dst)], cwd=ROOT, check=True)
        else:
            shutil.move(str(src), str(dst))
        ok += 1
    print(f"moved {ok} skipped {skip}")


if __name__ == "__main__":
    main()