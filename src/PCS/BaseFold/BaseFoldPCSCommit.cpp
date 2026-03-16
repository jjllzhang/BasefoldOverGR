#include "BaseFoldPCSInternal.hpp"

using NTL::LogicError;
using NTL::vec_ZZ_pE;

namespace basefold {

namespace {

ExtensionCommitRoundPrecomputation BuildExtensionCommitRoundPrecomputation(
    const FoldableCodeParams &params) {
  ExtensionCommitRoundPrecomputation precomputation;
  precomputation.levels.resize(static_cast<std::size_t>(params.d));
  if (params.d == 0) {
    return precomputation;
  }

  FieldElement one;
  NTL::set(one);
  const FieldElement zeta_minus_one = params.zeta - one;

  for (long level = 0; level < params.d; ++level) {
    const long n_i =
        basefold_pcs_internal::CodewordLengthAtLevelNoValidate(params, level);
    ExtensionCommitRoundLevelPrecomputation &level_cache =
        precomputation.levels[static_cast<std::size_t>(level)];
    if (n_i == 0) {
      continue;
    }

    const Oracle &diag = params.diag_T[static_cast<std::size_t>(level)];
    const FieldElement first_t = diag[0];
    bool all_equal = true;
    for (long j = 1; j < n_i; ++j) {
      if (diag[static_cast<std::size_t>(j)] != first_t) {
        all_equal = false;
        break;
      }
    }

    if (all_equal) {
      level_cache.inv_denoms.SetLength(1);
      const FieldElement denom = zeta_minus_one * first_t;
      if (!basefold_pcs_internal::TryInvertBaseUnit(level_cache.inv_denoms[0],
                                                    denom)) {
        LogicError(
            "BuildExtensionCommitRoundPrecomputation: denominator is not invertible");
      }
      continue;
    }

    std::vector<FieldElement> denoms(static_cast<std::size_t>(n_i));
    basefold_pcs_internal::ForEachIndexMaybeParallel(
        0, n_i, /*parallel_threshold=*/4096, [&](long j) {
          denoms[static_cast<std::size_t>(j)] =
              zeta_minus_one * diag[static_cast<std::size_t>(j)];
        });

    std::vector<FieldElement> inv_denoms;
    if (!basefold_pcs_internal::BatchInvertBaseUnits(inv_denoms, denoms)) {
      inv_denoms.resize(static_cast<std::size_t>(n_i));
      basefold_pcs_internal::ForEachIndexMaybeParallel(
          0, n_i, /*parallel_threshold=*/4096, [&](long j) {
            if (!basefold_pcs_internal::TryInvertBaseUnit(
                    inv_denoms[static_cast<std::size_t>(j)],
                    denoms[static_cast<std::size_t>(j)])) {
              LogicError(
                  "BuildExtensionCommitRoundPrecomputation: denominator is not invertible");
            }
          });
    }

    level_cache.inv_denoms.SetLength(n_i);
    for (long j = 0; j < n_i; ++j) {
      level_cache.inv_denoms[static_cast<std::size_t>(j)] =
          inv_denoms[static_cast<std::size_t>(j)];
    }
  }

  return precomputation;
}

}  // namespace

BaseFoldPCSCommitArtifacts BaseFoldPCSBuildCommitArtifacts(
    const vec_ZZ_pE &f_coeffs, const FoldableCodeParams &params) {
  basefold_pcs_internal::ValidateParamsOrThrow(params);
  if (f_coeffs.length() != MessageLength(params)) {
    LogicError("BaseFoldPCSBuildCommitArtifacts: f_coeffs has wrong length");
  }
  return BaseFoldPCSBuildCommitArtifactsUnchecked(f_coeffs, params);
}

BaseFoldPCSCommitArtifacts BaseFoldPCSBuildCommitArtifactsUnchecked(
    const vec_ZZ_pE &f_coeffs, const FoldableCodeParams &params) {
  BaseFoldPCSCommitArtifacts commit_artifacts;
  EncodeFoldableUnchecked(commit_artifacts.pi_d, f_coeffs, params);
  commit_artifacts.merkle_d = MerkleTree::Build(commit_artifacts.pi_d);
  commit_artifacts.root_d = commit_artifacts.merkle_d.Root();
  commit_artifacts.base_sumcheck_precomputation =
      BuildSumcheckMonomialPrecomputation(f_coeffs);
  commit_artifacts.extension_commit_precomputation =
      BuildExtensionCommitRoundPrecomputation(params);
  return commit_artifacts;
}

MerkleRoot BaseFoldPCSCommit(const vec_ZZ_pE &f_coeffs,
                             const FoldableCodeParams &params) {
  basefold_pcs_internal::ValidateParamsOrThrow(params);
  if (f_coeffs.length() != MessageLength(params)) {
    LogicError("BaseFoldPCSCommit: f_coeffs has wrong length");
  }

  Oracle pi_d;
  EncodeFoldableUnchecked(pi_d, f_coeffs, params);
  return MerkleTree::Build(pi_d).Root();
}

}  // namespace basefold
