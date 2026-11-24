#!/usr/bin/env bash
# Step 2: git mv neo/* -> neo/base/* and neo/algo/*
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

mv_one() {
  local src="$1" dst="$2"
  if [[ ! -f "$src" ]]; then
    echo "missing: $src" >&2
    return 1
  fi
  mkdir -p "$(dirname "$dst")"
  git mv "$src" "$dst"
}

# --- base/core ---
mv_one neo/base/core/algebra/ZMod.hpp neo/base/core/algebra/ZMod.hpp
mv_one neo/base/core/algebra/GF2p64.hpp neo/base/core/algebra/GF2p64.hpp
mv_one neo/base/core/algebra/semiring.hpp neo/base/core/algebra/semiring.hpp
mv_one neo/base/core/util/compress.hpp neo/base/core/util/compress.hpp
mv_one neo/base/core/util/argsort.hpp neo/base/core/util/argsort.hpp
mv_one neo/base/core/util/sort.hpp neo/base/core/util/sort.hpp
mv_one neo/base/core/linalg/Matrix.hpp neo/base/core/linalg/Matrix.hpp
mv_one neo/base/core/linalg/SquareMatrix.hpp neo/base/core/linalg/SquareMatrix.hpp
mv_one neo/base/core/linalg/BitMatrix.hpp neo/base/core/linalg/BitMatrix.hpp
mv_one neo/base/core/linalg/XorBasis.hpp neo/base/core/linalg/XorBasis.hpp
mv_one neo/base/core/nt/gcd.hpp neo/base/core/nt/gcd.hpp
mv_one neo/base/core/nt/inv_gcd.hpp neo/base/core/nt/inv_gcd.hpp
mv_one neo/base/core/nt/is_prime.hpp neo/base/core/nt/is_prime.hpp
mv_one neo/base/core/nt/pow_table.hpp neo/base/core/nt/pow_table.hpp
mv_one neo/base/core/nt/enumerate_primes.hpp neo/base/core/nt/enumerate_primes.hpp
mv_one neo/base/core/nt/sqrt_mod.hpp neo/base/core/nt/sqrt_mod.hpp
mv_one neo/base/core/nt/stern_brocot.hpp neo/base/core/nt/stern_brocot.hpp
mv_one neo/base/core/convolve/naive.hpp neo/base/core/convolve/naive.hpp
mv_one neo/base/core/poly/poly_divmod.hpp neo/base/core/poly/poly_divmod.hpp
mv_one neo/base/core/poly/sparse_fps.hpp neo/base/core/poly/sparse_fps.hpp
mv_one neo/base/core/convolve/GF2p64.hpp neo/base/core/convolve/GF2p64.hpp
mv_one neo/base/core/set_fps/generic.hpp neo/base/core/set_fps/generic.hpp
mv_one neo/base/core/set_fps/ranked_xor.hpp neo/base/core/set_fps/ranked_xor.hpp
mv_one neo/base/core/set_fps/ranked_or64.hpp neo/base/core/set_fps/ranked_or64.hpp
mv_one neo/base/core/set_fps/semiring_simd.hpp neo/base/core/set_fps/semiring_simd.hpp
mv_one neo/base/core/set_fps/zeta_simd.hpp neo/base/core/set_fps/zeta_simd.hpp
mv_one neo/base/core/set_fps/digraph_simd.hpp neo/base/core/set_fps/digraph_simd.hpp
mv_one neo/base/core/set_fps/zeta.hpp neo/base/core/set_fps/zeta.hpp
mv_one neo/base/core/set_fps/exp.hpp neo/base/core/set_fps/exp.hpp
mv_one neo/base/core/set_fps/log.hpp neo/base/core/set_fps/log.hpp
mv_one neo/base/core/set_fps/inv.hpp neo/base/core/set_fps/inv.hpp
mv_one neo/base/core/set_fps/convolve.hpp neo/base/core/set_fps/convolve.hpp
mv_one neo/base/core/set_fps/composite.hpp neo/base/core/set_fps/composite.hpp
mv_one neo/base/core/set_fps/relaxed.hpp neo/base/core/set_fps/relaxed.hpp

# set_power_series 案B (algo)
mv_one neo/algo/count/graph/graph.hpp neo/algo/count/graph/graph.hpp
mv_one neo/algo/compute/set_fps/hafnian.hpp neo/algo/compute/set_fps/hafnian.hpp
mv_one neo/algo/count/set_fps/power_projection.hpp neo/algo/count/set_fps/power_projection.hpp

# --- base/ds ---
mv_one neo/base/ds/range/SegmentTree.hpp neo/base/ds/range/SegmentTree.hpp
mv_one neo/base/ds/range/DynamicSegmentTree.hpp neo/base/ds/range/DynamicSegmentTree.hpp
mv_one neo/base/ds/range/RangeMin.hpp neo/base/ds/range/RangeMin.hpp
mv_one neo/base/ds/range/RangeLinearAddRangeMin.hpp neo/base/ds/range/RangeLinearAddRangeMin.hpp
for f in AddMin AddSum Affine Assign AssignProd ChminChmaxAddSum Composite Min Sum; do
  mv_one "neo/data_structure/segtree/${f}.hpp" "neo/base/ds/range/segtree/${f}.hpp"
done
mv_one neo/base/ds/range/WaveletMatrix.hpp neo/base/ds/range/WaveletMatrix.hpp
mv_one neo/base/ds/range/WaveletMatrixSum.hpp neo/base/ds/range/WaveletMatrixSum.hpp
mv_one neo/base/ds/set/FastSet.hpp neo/base/ds/set/FastSet.hpp
mv_one neo/base/ds/set/OrderedSet.hpp neo/base/ds/set/OrderedSet.hpp
mv_one neo/base/ds/hash/HashMap.hpp neo/base/ds/hash/HashMap.hpp
mv_one neo/base/ds/seq/FoldableQueue.hpp neo/base/ds/seq/FoldableQueue.hpp
mv_one neo/base/ds/seq/FoldableDeque.hpp neo/base/ds/seq/FoldableDeque.hpp
mv_one neo/base/ds/seq/WeightBalancedTree.hpp neo/base/ds/seq/WeightBalancedTree.hpp
mv_one neo/base/ds/seq/Sequence.hpp neo/base/ds/seq/Sequence.hpp
mv_one neo/base/ds/range/SubsetSum.hpp neo/base/ds/range/SubsetSum.hpp
mv_one neo/base/ds/2d/PointAddRectangleSum.hpp neo/base/ds/2d/PointAddRectangleSum.hpp
mv_one neo/base/ds/2d/RectangleSum.hpp neo/base/ds/2d/RectangleSum.hpp
mv_one neo/base/ds/2d/RectangleCount.hpp neo/base/ds/2d/RectangleCount.hpp
mv_one neo/base/ds/tree/HeavyLightDecomposition.hpp neo/base/ds/tree/HeavyLightDecomposition.hpp
mv_one neo/base/ds/tree/LCA.hpp neo/base/ds/tree/LCA.hpp
mv_one neo/base/ds/tree/FunctionalGraph.hpp neo/base/ds/tree/FunctionalGraph.hpp
mv_one neo/base/ds/tree/PathProduct.hpp neo/base/ds/tree/PathProduct.hpp
for f in RollingHash suffix_array prefix_function z_algorithm palindrome_lengths lyndon_factorization \
  aho_corasick SuffixAutomaton Trie FMIndex SuffixLCP runs eertree PalindromeDeque; do
  mv_one "neo/string/${f}.hpp" "neo/base/ds/string/${f}.hpp"
done

# --- base/engine ---
mv_one neo/base/engine/tree/Rerooting.hpp neo/base/engine/tree/Rerooting.hpp

# --- algo/optimize ---
mv_one neo/algo/optimize/flow/MaxFlow.hpp neo/algo/optimize/flow/MaxFlow.hpp
mv_one neo/algo/optimize/flow/MinCostFlow.hpp neo/algo/optimize/flow/MinCostFlow.hpp
mv_one neo/algo/optimize/flow/MaxFlowLowerBound.hpp neo/algo/optimize/flow/MaxFlowLowerBound.hpp
mv_one neo/algo/optimize/flow/PrimalDual.hpp neo/algo/optimize/flow/PrimalDual.hpp
mv_one neo/algo/optimize/flow/gomory_hu.hpp neo/algo/optimize/flow/gomory_hu.hpp
mv_one neo/algo/optimize/matching/BipartiteMatching.hpp neo/algo/optimize/matching/BipartiteMatching.hpp
mv_one neo/algo/optimize/matching/GeneralMatching.hpp neo/algo/optimize/matching/GeneralMatching.hpp
mv_one neo/algo/optimize/matching/weighted_matching.hpp neo/algo/optimize/matching/weighted_matching.hpp
mv_one neo/algo/optimize/matching/IncrementalBipartiteMatching.hpp neo/algo/optimize/matching/IncrementalBipartiteMatching.hpp
mv_one neo/algo/optimize/matching/assignment.hpp neo/algo/optimize/matching/assignment.hpp
mv_one neo/algo/optimize/matching/bipartite_edge_coloring.hpp neo/algo/optimize/matching/bipartite_edge_coloring.hpp
mv_one neo/algo/optimize/cut/global_min_cut.hpp neo/algo/optimize/cut/global_min_cut.hpp
mv_one neo/algo/optimize/matroid/matroids.hpp neo/algo/optimize/matroid/matroids.hpp
mv_one neo/algo/optimize/matroid/matroid_intersection.hpp neo/algo/optimize/matroid/matroid_intersection.hpp
mv_one neo/algo/optimize/burn_bury/BinaryOptimization.hpp neo/algo/optimize/burn_bury/BinaryOptimization.hpp
mv_one neo/algo/optimize/burn_bury/KaryOptimization.hpp neo/algo/optimize/burn_bury/KaryOptimization.hpp
mv_one neo/algo/optimize/cover/path_cover.hpp neo/algo/optimize/cover/path_cover.hpp

# --- algo/count ---
mv_one neo/algo/count/automata/DFA.hpp neo/algo/count/automata/DFA.hpp
mv_one neo/algo/count/nt/floor_sum.hpp neo/algo/count/nt/floor_sum.hpp
mv_one neo/algo/count/nt/mertens.hpp neo/algo/count/nt/mertens.hpp
mv_one neo/algo/count/nt/totient_sum.hpp neo/algo/count/nt/totient_sum.hpp
mv_one neo/algo/count/nt/prime_pi.hpp neo/algo/count/nt/prime_pi.hpp
mv_one neo/algo/count/string/distinct_subsequences.hpp neo/algo/count/string/distinct_subsequences.hpp
mv_one neo/algo/count/string/PrefixSubstringLCS.hpp neo/algo/count/string/PrefixSubstringLCS.hpp
mv_one neo/algo/count/string/find_all.hpp neo/algo/count/string/find_all.hpp

# --- algo/decide ---
mv_one neo/algo/decide/graph/strongly_connected_components.hpp neo/algo/decide/graph/strongly_connected_components.hpp
mv_one neo/algo/decide/graph/two_sat.hpp neo/algo/decide/graph/two_sat.hpp
mv_one neo/algo/decide/graph/biconnected_components.hpp neo/algo/decide/graph/biconnected_components.hpp
mv_one neo/algo/decide/graph/two_edge_connected_components.hpp neo/algo/decide/graph/two_edge_connected_components.hpp
mv_one neo/algo/decide/graph/three_edge_connected_components.hpp neo/algo/decide/graph/three_edge_connected_components.hpp
mv_one neo/algo/decide/graph/dominator_tree.hpp neo/algo/decide/graph/dominator_tree.hpp
mv_one neo/algo/decide/graph/bipartite_coloring.hpp neo/algo/decide/graph/bipartite_coloring.hpp
mv_one neo/algo/decide/graph/euler_split.hpp neo/algo/decide/graph/euler_split.hpp

# --- algo/compute ---
mv_one neo/algo/compute/nt/factorize.hpp neo/algo/compute/nt/factorize.hpp
mv_one neo/algo/compute/nt/discrete_log.hpp neo/algo/compute/nt/discrete_log.hpp
mv_one neo/algo/compute/nt/DiscreteLogTable.hpp neo/algo/compute/nt/DiscreteLogTable.hpp
mv_one neo/algo/compute/nt/kth_root_mod.hpp neo/algo/compute/nt/kth_root_mod.hpp
mv_one neo/algo/compute/nt/primitive_root.hpp neo/algo/compute/nt/primitive_root.hpp
mv_one neo/algo/compute/nt/kth_root.hpp neo/algo/compute/nt/kth_root.hpp
mv_one neo/algo/compute/poly/lagrange.hpp neo/algo/compute/poly/lagrange.hpp
mv_one neo/algo/compute/poly/linear_recurrence.hpp neo/algo/compute/poly/linear_recurrence.hpp
mv_one neo/algo/compute/poly/poly_factorize.hpp neo/algo/compute/poly/poly_factorize.hpp
mv_one neo/algo/compute/poly/poly_roots.hpp neo/algo/compute/poly/poly_roots.hpp
mv_one neo/algo/compute/poly/berlekamp_massey.hpp neo/algo/compute/poly/berlekamp_massey.hpp
mv_one neo/algo/compute/linalg/det.hpp neo/algo/compute/linalg/det.hpp
mv_one neo/algo/compute/linalg/solve.hpp neo/algo/compute/linalg/solve.hpp
mv_one neo/algo/compute/linalg/inverse.hpp neo/algo/compute/linalg/inverse.hpp
mv_one neo/algo/compute/linalg/rank.hpp neo/algo/compute/linalg/rank.hpp
mv_one neo/algo/compute/linalg/adjugate.hpp neo/algo/compute/linalg/adjugate.hpp
mv_one neo/algo/compute/linalg/characteristic_polynomial.hpp neo/algo/compute/linalg/characteristic_polynomial.hpp
mv_one neo/algo/compute/linalg/sparse_det.hpp neo/algo/compute/linalg/sparse_det.hpp
mv_one neo/algo/compute/linalg/pfaffian.hpp neo/algo/compute/linalg/pfaffian.hpp

echo "neo hpp moves done"