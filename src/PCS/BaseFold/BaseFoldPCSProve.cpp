#include "BaseFoldPCSInternal.hpp"

#include <string>
#include <vector>

#include "PCS/Common/Multilinear.hpp"

using NTL::LogicError;
using NTL::vec_ZZ_pE;

namespace basefold {
namespace {

SumcheckProver MakeBaseSumcheckProver(
    const vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
    const BaseFoldPCSCommitArtifacts &commit_artifacts) {
  if (commit_artifacts.base_sumcheck_precomputation.valid) {
    return SumcheckProver(commit_artifacts.base_sumcheck_precomputation, z);
  }
  return SumcheckProver(f_coeffs, z);
}

BaseFoldPCSEvalProof ProveEvalFromCommittedTopOracleUnchecked(
    const vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const FoldableCodeParams &params,
    const BaseFoldPCSCommitArtifacts &commit_artifacts) {
  basefold_pcs_internal::ValidateCommittedTopOracleArtifactsOrThrow(
      commit_artifacts,
      basefold_pcs_internal::CodewordLengthAtLevelNoValidate(params, params.d),
      "BaseFoldPCSProveEvalFromCommittedTopOracleUnchecked");

  BaseFoldPCSEvalProof proof;
  proof.commitments.roots_by_level.resize(
      static_cast<std::size_t>(params.d + 1));
  proof.h_by_level.resize(static_cast<std::size_t>(params.d));
  proof.commitments.roots_by_level[static_cast<std::size_t>(params.d)] =
      commit_artifacts.root_d;

  HashTranscript transcript = basefold_pcs_internal::MakeBaseFoldTranscript();
  basefold_pcs_internal::AbsorbPublicInput(transcript, commit_artifacts.root_d, z,
                                           claimed_y);

  if (params.d == 0) {
    proof.pi0_codeword = commit_artifacts.pi_d;
    return proof;
  }

  std::vector<Oracle> oracles(static_cast<std::size_t>(params.d));
  std::vector<MerkleTree> merkle(static_cast<std::size_t>(params.d));

  SumcheckProver sumcheck = MakeBaseSumcheckProver(f_coeffs, z, commit_artifacts);

  const QuadraticPoly h_d = sumcheck.CurrentPolynomial();
  proof.h_by_level[static_cast<std::size_t>(params.d - 1)] = h_d;
  AbsorbQuadraticPoly(transcript, h_d);

  std::vector<FieldElement> r_by_level(static_cast<std::size_t>(params.d));
  for (long i = params.d; i-- > 0;) {
    const FieldElement r_i =
        transcript.ChallengeFieldElement("r/" + std::to_string(i));
    r_by_level[static_cast<std::size_t>(i)] = r_i;

    const Oracle &upper_oracle =
        (i + 1 == params.d) ? commit_artifacts.pi_d
                            : oracles[static_cast<std::size_t>(i + 1)];
    basefold_pcs_internal::ProverCommitRoundNoValidate(
        oracles[static_cast<std::size_t>(i)], upper_oracle, r_i, i, params);

    merkle[static_cast<std::size_t>(i)] =
        MerkleTree::Build(oracles[static_cast<std::size_t>(i)]);
    const MerkleRoot root_i = merkle[static_cast<std::size_t>(i)].Root();
    proof.commitments.roots_by_level[static_cast<std::size_t>(i)] = root_i;
    transcript.AbsorbDigest(root_i);

    sumcheck.ReceiveChallenge(r_i);
    if (i > 0) {
      const QuadraticPoly h_i = sumcheck.CurrentPolynomial();
      proof.h_by_level[static_cast<std::size_t>(i - 1)] = h_i;
      AbsorbQuadraticPoly(transcript, h_i);
    }
  }

  proof.pi0_codeword = oracles[0];

  const long n_last =
      basefold_pcs_internal::CodewordLengthAtLevelNoValidate(params, params.d - 1);
  std::vector<IOPPQueryPlan> query_plans(static_cast<std::size_t>(num_queries));
  for (long q = 0; q < num_queries; ++q) {
    const long mu =
        transcript.ChallengeIndex("mu/" + std::to_string(q), n_last);
    query_plans[static_cast<std::size_t>(q)] =
        basefold_pcs_internal::MakeQueryPlanNoValidate(mu, params);
  }

  proof.query_multiproofs.resize(static_cast<std::size_t>(params.d + 1));
  const std::vector<std::vector<long>> requested_indices_by_tree =
      basefold_pcs_internal::CollectBaseQueryIndicesByTree(query_plans, params);
  for (long tree_level = 0; tree_level <= params.d; ++tree_level) {
    if (tree_level == params.d) {
      proof.query_multiproofs[static_cast<std::size_t>(tree_level)] =
          commit_artifacts.merkle_d.OpenMany(
              commit_artifacts.pi_d,
              requested_indices_by_tree[static_cast<std::size_t>(tree_level)]);
      continue;
    }
    proof.query_multiproofs[static_cast<std::size_t>(tree_level)] =
        merkle[static_cast<std::size_t>(tree_level)].OpenMany(
            oracles[static_cast<std::size_t>(tree_level)],
            requested_indices_by_tree[static_cast<std::size_t>(tree_level)]);
  }

  return proof;
}

}  // namespace

BaseFoldPCSEvalProof BaseFoldPCSProveEval(const vec_ZZ_pE &f_coeffs,
                                          const std::vector<FieldElement> &z,
                                          const FieldElement &claimed_y,
                                          long num_queries,
                                          const FoldableCodeParams &params) {
  basefold_pcs_internal::ValidateParamsOrThrow(params);
  const long kappa = basefold_pcs_internal::Log2ExactPowerOfTwoLong(params.k0);
  const long point_dim = params.d + kappa;
  if (static_cast<long>(z.size()) != point_dim) {
    LogicError("BaseFoldPCSProveEval: z has wrong dimension");
  }
  if (f_coeffs.length() != MessageLength(params)) {
    LogicError("BaseFoldPCSProveEval: f_coeffs has wrong length");
  }
  if (num_queries < 0) {
    LogicError("BaseFoldPCSProveEval: num_queries must be non-negative");
  }
  if (EvalMultilinearMonomialCoeffs(f_coeffs, z) != claimed_y) {
    LogicError("BaseFoldPCSProveEval: claimed_y != f(z)");
  }

  const BaseFoldPCSCommitArtifacts commit_artifacts =
      BaseFoldPCSBuildCommitArtifactsUnchecked(f_coeffs, params);
  return BaseFoldPCSProveEvalFromCommittedTopOracleUnchecked(
      f_coeffs, z, claimed_y, num_queries, params, commit_artifacts);
}

BaseFoldPCSEvalProof BaseFoldPCSProveEvalUnchecked(
    const vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const FoldableCodeParams &params) {
  const BaseFoldPCSCommitArtifacts commit_artifacts =
      BaseFoldPCSBuildCommitArtifactsUnchecked(f_coeffs, params);
  return BaseFoldPCSProveEvalFromCommittedTopOracleUnchecked(
      f_coeffs, z, claimed_y, num_queries, params, commit_artifacts);
}

BaseFoldPCSEvalProof BaseFoldPCSProveEvalFromCommittedTopOracle(
    const vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const FoldableCodeParams &params,
    const BaseFoldPCSCommitArtifacts &commit_artifacts) {
  basefold_pcs_internal::ValidateParamsOrThrow(params);
  const long kappa = basefold_pcs_internal::Log2ExactPowerOfTwoLong(params.k0);
  const long point_dim = params.d + kappa;
  if (static_cast<long>(z.size()) != point_dim) {
    LogicError("BaseFoldPCSProveEvalFromCommittedTopOracle: z has wrong dimension");
  }
  if (f_coeffs.length() != MessageLength(params)) {
    LogicError(
        "BaseFoldPCSProveEvalFromCommittedTopOracle: f_coeffs has wrong length");
  }
  if (num_queries < 0) {
    LogicError(
        "BaseFoldPCSProveEvalFromCommittedTopOracle: num_queries must be non-negative");
  }
  if (EvalMultilinearMonomialCoeffs(f_coeffs, z) != claimed_y) {
    LogicError(
        "BaseFoldPCSProveEvalFromCommittedTopOracle: claimed_y != f(z)");
  }

  return BaseFoldPCSProveEvalFromCommittedTopOracleUnchecked(
      f_coeffs, z, claimed_y, num_queries, params, commit_artifacts);
}

BaseFoldPCSEvalProof BaseFoldPCSProveEvalFromCommittedTopOracleUnchecked(
    const vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const FoldableCodeParams &params,
    const BaseFoldPCSCommitArtifacts &commit_artifacts) {
  return ProveEvalFromCommittedTopOracleUnchecked(
      f_coeffs, z, claimed_y, num_queries, params, commit_artifacts);
}

}  // namespace basefold
