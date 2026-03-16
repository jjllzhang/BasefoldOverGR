#include "BaseFoldPCSInternal.hpp"

using NTL::LogicError;
using NTL::vec_ZZ_pE;

namespace basefold {

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
