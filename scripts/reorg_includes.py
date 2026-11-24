#!/usr/bin/env python3
"""Replace neo/... include paths after base/algo reorg."""
from __future__ import annotations

import re
import sys
from pathlib import Path

# old path (without neo/) -> new path (without neo/)
MAPPING: dict[str, str] = {
    # base/core
    "algebra/ZMod.hpp": "base/core/algebra/ZMod.hpp",
    "algebra/GF2p64.hpp": "base/core/algebra/GF2p64.hpp",
    "algebra/semiring.hpp": "base/core/algebra/semiring.hpp",
    "sort/compress.hpp": "base/core/util/compress.hpp",
    "sort/argsort.hpp": "base/core/util/argsort.hpp",
    "sort/sort.hpp": "base/core/util/sort.hpp",
    "linalg/Matrix.hpp": "base/core/linalg/Matrix.hpp",
    "linalg/SquareMatrix.hpp": "base/core/linalg/SquareMatrix.hpp",
    "linalg/BitMatrix.hpp": "base/core/linalg/BitMatrix.hpp",
    "linalg/XorBasis.hpp": "base/core/linalg/XorBasis.hpp",
    "number_theory/gcd.hpp": "base/core/nt/gcd.hpp",
    "number_theory/inv_gcd.hpp": "base/core/nt/inv_gcd.hpp",
    "number_theory/is_prime.hpp": "base/core/nt/is_prime.hpp",
    "number_theory/pow_table.hpp": "base/core/nt/pow_table.hpp",
    "number_theory/enumerate_primes.hpp": "base/core/nt/enumerate_primes.hpp",
    "number_theory/sqrt_mod.hpp": "base/core/nt/sqrt_mod.hpp",
    "number_theory/stern_brocot.hpp": "base/core/nt/stern_brocot.hpp",
    "polynomial/convolve_naive.hpp": "base/core/convolve/naive.hpp",
    "polynomial/poly_divmod.hpp": "base/core/poly/poly_divmod.hpp",
    "polynomial/sparse_fps.hpp": "base/core/poly/sparse_fps.hpp",
    "fft/convolve_GF2p64.hpp": "base/core/convolve/GF2p64.hpp",
    "set_power_series/generic.hpp": "base/core/set_fps/generic.hpp",
    "set_power_series/ranked_xor.hpp": "base/core/set_fps/ranked_xor.hpp",
    "set_power_series/ranked_or64.hpp": "base/core/set_fps/ranked_or64.hpp",
    "set_power_series/semiring_simd.hpp": "base/core/set_fps/semiring_simd.hpp",
    "set_power_series/zeta_simd.hpp": "base/core/set_fps/zeta_simd.hpp",
    "set_power_series/digraph_simd.hpp": "base/core/set_fps/digraph_simd.hpp",
    "set_power_series/zeta.hpp": "base/core/set_fps/zeta.hpp",
    "set_power_series/exp.hpp": "base/core/set_fps/exp.hpp",
    "set_power_series/log.hpp": "base/core/set_fps/log.hpp",
    "set_power_series/inv.hpp": "base/core/set_fps/inv.hpp",
    "set_power_series/convolve.hpp": "base/core/set_fps/convolve.hpp",
    "set_power_series/composite.hpp": "base/core/set_fps/composite.hpp",
    "set_power_series/relaxed.hpp": "base/core/set_fps/relaxed.hpp",
    "set_power_series/graph.hpp": "algo/count/graph/graph.hpp",
    "set_power_series/hafnian.hpp": "algo/compute/set_fps/hafnian.hpp",
    "set_power_series/power_projection.hpp": "algo/count/set_fps/power_projection.hpp",
    # base/ds
    "data_structure/SegmentTree.hpp": "base/ds/range/SegmentTree.hpp",
    "data_structure/DynamicSegmentTree.hpp": "base/ds/range/DynamicSegmentTree.hpp",
    "data_structure/RangeMin.hpp": "base/ds/range/RangeMin.hpp",
    "data_structure/RangeLinearAddRangeMin.hpp": "base/ds/range/RangeLinearAddRangeMin.hpp",
    "data_structure/WaveletMatrix.hpp": "base/ds/range/WaveletMatrix.hpp",
    "data_structure/WaveletMatrixSum.hpp": "base/ds/range/WaveletMatrixSum.hpp",
    "data_structure/FastSet.hpp": "base/ds/set/FastSet.hpp",
    "data_structure/OrderedSet.hpp": "base/ds/set/OrderedSet.hpp",
    "data_structure/HashMap.hpp": "base/ds/hash/HashMap.hpp",
    "data_structure/FoldableQueue.hpp": "base/ds/seq/FoldableQueue.hpp",
    "data_structure/FoldableDeque.hpp": "base/ds/seq/FoldableDeque.hpp",
    "data_structure/WeightBalancedTree.hpp": "base/ds/seq/WeightBalancedTree.hpp",
    "data_structure/Sequence.hpp": "base/ds/seq/Sequence.hpp",
    "data_structure/SubsetSum.hpp": "base/ds/range/SubsetSum.hpp",
    "data_structure/PointAddRectangleSum.hpp": "base/ds/2d/PointAddRectangleSum.hpp",
    "data_structure/RectangleSum.hpp": "base/ds/2d/RectangleSum.hpp",
    "data_structure/RectangleCount.hpp": "base/ds/2d/RectangleCount.hpp",
    "graph/HeavyLightDecomposition.hpp": "base/ds/tree/HeavyLightDecomposition.hpp",
    "graph/LCA.hpp": "base/ds/tree/LCA.hpp",
    "graph/FunctionalGraph.hpp": "base/ds/tree/FunctionalGraph.hpp",
    "graph/PathProduct.hpp": "base/ds/tree/PathProduct.hpp",
    "string/RollingHash.hpp": "base/ds/string/RollingHash.hpp",
    "string/suffix_array.hpp": "base/ds/string/suffix_array.hpp",
    "string/prefix_function.hpp": "base/ds/string/prefix_function.hpp",
    "string/z_algorithm.hpp": "base/ds/string/z_algorithm.hpp",
    "string/palindrome_lengths.hpp": "base/ds/string/palindrome_lengths.hpp",
    "string/lyndon_factorization.hpp": "base/ds/string/lyndon_factorization.hpp",
    "string/aho_corasick.hpp": "base/ds/string/aho_corasick.hpp",
    "string/SuffixAutomaton.hpp": "base/ds/string/SuffixAutomaton.hpp",
    "string/Trie.hpp": "base/ds/string/Trie.hpp",
    "string/FMIndex.hpp": "base/ds/string/FMIndex.hpp",
    "string/SuffixLCP.hpp": "base/ds/string/SuffixLCP.hpp",
    "string/runs.hpp": "base/ds/string/runs.hpp",
    "string/eertree.hpp": "base/ds/string/eertree.hpp",
    "string/PalindromeDeque.hpp": "base/ds/string/PalindromeDeque.hpp",
    # segtree (longer prefix)
    "data_structure/segtree/AddMin.hpp": "base/ds/range/segtree/AddMin.hpp",
    "data_structure/segtree/AddSum.hpp": "base/ds/range/segtree/AddSum.hpp",
    "data_structure/segtree/Affine.hpp": "base/ds/range/segtree/Affine.hpp",
    "data_structure/segtree/Assign.hpp": "base/ds/range/segtree/Assign.hpp",
    "data_structure/segtree/AssignProd.hpp": "base/ds/range/segtree/AssignProd.hpp",
    "data_structure/segtree/ChminChmaxAddSum.hpp": "base/ds/range/segtree/ChminChmaxAddSum.hpp",
    "data_structure/segtree/Composite.hpp": "base/ds/range/segtree/Composite.hpp",
    "data_structure/segtree/Min.hpp": "base/ds/range/segtree/Min.hpp",
    "data_structure/segtree/Sum.hpp": "base/ds/range/segtree/Sum.hpp",
    # engine
    "graph/Rerooting.hpp": "base/engine/tree/Rerooting.hpp",
    # algo optimize
    "combinatorial_optimization/MaxFlow.hpp": "algo/optimize/flow/MaxFlow.hpp",
    "combinatorial_optimization/MinCostFlow.hpp": "algo/optimize/flow/MinCostFlow.hpp",
    "combinatorial_optimization/MaxFlowLowerBound.hpp": "algo/optimize/flow/MaxFlowLowerBound.hpp",
    "combinatorial_optimization/PrimalDual.hpp": "algo/optimize/flow/PrimalDual.hpp",
    "combinatorial_optimization/gomory_hu.hpp": "algo/optimize/flow/gomory_hu.hpp",
    "combinatorial_optimization/BipartiteMatching.hpp": "algo/optimize/matching/BipartiteMatching.hpp",
    "combinatorial_optimization/GeneralMatching.hpp": "algo/optimize/matching/GeneralMatching.hpp",
    "combinatorial_optimization/weighted_matching.hpp": "algo/optimize/matching/weighted_matching.hpp",
    "combinatorial_optimization/IncrementalBipartiteMatching.hpp": "algo/optimize/matching/IncrementalBipartiteMatching.hpp",
    "combinatorial_optimization/assignment.hpp": "algo/optimize/matching/assignment.hpp",
    "combinatorial_optimization/bipartite_edge_coloring.hpp": "algo/optimize/matching/bipartite_edge_coloring.hpp",
    "combinatorial_optimization/global_min_cut.hpp": "algo/optimize/cut/global_min_cut.hpp",
    "combinatorial_optimization/matroids.hpp": "algo/optimize/matroid/matroids.hpp",
    "combinatorial_optimization/matroid_intersection.hpp": "algo/optimize/matroid/matroid_intersection.hpp",
    "combinatorial_optimization/BinaryOptimization.hpp": "algo/optimize/burn_bury/BinaryOptimization.hpp",
    "combinatorial_optimization/KaryOptimization.hpp": "algo/optimize/burn_bury/KaryOptimization.hpp",
    "combinatorial_optimization/path_cover.hpp": "algo/optimize/cover/path_cover.hpp",
    # count
    "counting/DFA.hpp": "algo/count/automata/DFA.hpp",
    "number_theory/floor_sum.hpp": "algo/count/nt/floor_sum.hpp",
    "number_theory/mertens.hpp": "algo/count/nt/mertens.hpp",
    "number_theory/totient_sum.hpp": "algo/count/nt/totient_sum.hpp",
    "number_theory/prime_pi.hpp": "algo/count/nt/prime_pi.hpp",
    "string/distinct_subsequences.hpp": "algo/count/string/distinct_subsequences.hpp",
    "string/PrefixSubstringLCS.hpp": "algo/count/string/PrefixSubstringLCS.hpp",
    "string/find_all.hpp": "algo/count/string/find_all.hpp",
    # decide
    "graph/strongly_connected_components.hpp": "algo/decide/graph/strongly_connected_components.hpp",
    "graph/two_sat.hpp": "algo/decide/graph/two_sat.hpp",
    "graph/biconnected_components.hpp": "algo/decide/graph/biconnected_components.hpp",
    "graph/two_edge_connected_components.hpp": "algo/decide/graph/two_edge_connected_components.hpp",
    "graph/three_edge_connected_components.hpp": "algo/decide/graph/three_edge_connected_components.hpp",
    "graph/dominator_tree.hpp": "algo/decide/graph/dominator_tree.hpp",
    "graph/bipartite_coloring.hpp": "algo/decide/graph/bipartite_coloring.hpp",
    "graph/euler_split.hpp": "algo/decide/graph/euler_split.hpp",
    # compute
    "number_theory/factorize.hpp": "algo/compute/nt/factorize.hpp",
    "number_theory/discrete_log.hpp": "algo/compute/nt/discrete_log.hpp",
    "number_theory/DiscreteLogTable.hpp": "algo/compute/nt/DiscreteLogTable.hpp",
    "number_theory/kth_root_mod.hpp": "algo/compute/nt/kth_root_mod.hpp",
    "number_theory/primitive_root.hpp": "algo/compute/nt/primitive_root.hpp",
    "number_theory/kth_root.hpp": "algo/compute/nt/kth_root.hpp",
    "polynomial/lagrange.hpp": "algo/compute/poly/lagrange.hpp",
    "polynomial/linear_recurrence.hpp": "algo/compute/poly/linear_recurrence.hpp",
    "polynomial/poly_factorize.hpp": "algo/compute/poly/poly_factorize.hpp",
    "polynomial/poly_roots.hpp": "algo/compute/poly/poly_roots.hpp",
    "polynomial/berlekamp_massey.hpp": "algo/compute/poly/berlekamp_massey.hpp",
    "linalg/det.hpp": "algo/compute/linalg/det.hpp",
    "linalg/solve.hpp": "algo/compute/linalg/solve.hpp",
    "linalg/inverse.hpp": "algo/compute/linalg/inverse.hpp",
    "linalg/rank.hpp": "algo/compute/linalg/rank.hpp",
    "linalg/adjugate.hpp": "algo/compute/linalg/adjugate.hpp",
    "linalg/characteristic_polynomial.hpp": "algo/compute/linalg/characteristic_polynomial.hpp",
    "linalg/sparse_det.hpp": "algo/compute/linalg/sparse_det.hpp",
    "linalg/pfaffian.hpp": "algo/compute/linalg/pfaffian.hpp",
}

SKIP_SUFFIX = {
    ".png",
    ".jpg",
    ".gif",
    ".webp",
    ".woff",
    ".woff2",
    ".ico",
    ".pdf",
    ".zip",
    ".git",
}

SKIP_DIRS = {
    ".git",
    "node_modules",
    "site",
    "__pycache__",
}


def should_process(path: Path) -> bool:
    if path.suffix in SKIP_SUFFIX:
        return False
    for part in path.parts:
        if part in SKIP_DIRS:
            return False
    try:
        if path.is_file() and path.stat().st_size > 5_000_000:
            return False
    except OSError:
        return False
    return path.is_file()


def rewrite(text: str) -> tuple[str, int]:
    n = 0
    items = sorted(MAPPING.items(), key=lambda kv: len(kv[0]), reverse=True)
    for old, new in items:
        pat_old = f"neo/{old}"
        pat_new = f"neo/{new}"
        if pat_old not in text:
            continue
        count = text.count(pat_old)
        text = text.replace(pat_old, pat_new)
        n += count
    return text, n


def process_tree(root: Path) -> int:
    total = 0
    for path in root.rglob("*"):
        if not should_process(path):
            continue
        try:
            raw = path.read_bytes()
        except (OSError, PermissionError):
            continue
        if b"\0" in raw[:8192]:
            continue
        try:
            text = raw.decode("utf-8")
        except UnicodeDecodeError:
            continue
        new_text, n = rewrite(text)
        if n:
            path.write_text(new_text, encoding="utf-8")
            total += n
    return total


def main() -> None:
    roots = [Path(p) for p in sys.argv[1:]] or [
        Path(__file__).resolve().parent.parent,
        Path(__file__).resolve().parent.parent.parent / "procon-judge",
        Path(__file__).resolve().parent.parent.parent / "algo-notes",
    ]
    grand = 0
    for root in roots:
        if not root.exists():
            print(f"skip missing {root}")
            continue
        n = process_tree(root)
        print(f"{root}: {n} replacements")
        grand += n
    print(f"total: {grand}")


if __name__ == "__main__":
    main()