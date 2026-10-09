"""pr-syntax-targets.py の選び方。"""
import importlib.util
import unittest
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "pr_syntax_targets", Path(__file__).resolve().parent.parent / "pr-syntax-targets.py"
)
_mod = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_mod)
select_targets = _mod.select_targets


class TestSelectTargets(unittest.TestCase):
    def test_changed_headers_are_the_targets_in_path_order(self):
        self.assertEqual(
            select_targets({"neo/B.hpp", "neo/internal/X.hpp", "neo/A.hpp", "README.md"}),
            ["neo/A.hpp", "neo/B.hpp", "neo/internal/X.hpp"],
        )

    def test_unrelated_changes_select_nothing(self):
        self.assertEqual(
            select_targets({"scripts/x.ts", "README.md", "include/simde/simde-common.h", "md/algebra/A.md"}),
            [],
        )


if __name__ == "__main__":
    unittest.main()
