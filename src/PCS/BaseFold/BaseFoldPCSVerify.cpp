#include "BaseFoldPCSInternal.hpp"

#include <string>
#include <vector>

#include "PCS/Common/MerkleMultiproofReplay.hpp"
#include "PCS/Common/Multilinear.hpp"

namespace basefold {

void ResetVerifierQueryParallelConfigFromEnv() {
  basefold_pcs_internal::MutableVerifierQueryParallelConfig() =
      basefold_pcs_internal::ParseVerifierQueryParallelConfigFromEnv();
}

void SetVerifierQueryParallelConfig(const VerifierQueryParallelConfig &cfg) {
  basefold_pcs_internal::MutableVerifierQueryParallelConfig() = cfg;
}

VerifierQueryParallelConfig GetVerifierQueryParallelConfig() {
  return basefold_pcs_internal::LoadVerifierQueryParallelConfig();
}

bool BaseFoldPCSVerifyEval(const MerkleRoot &commitment_C,
                           const std::vector<FieldElement> &z,
                           const FieldElement &claimed_y, long num_queries,
                           const BaseFoldPCSEvalProof &proof,
                           const FoldableCodeParams &params) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->pcs_verify_ns : nullptr,
                    prof ? &prof->pcs_verify_calls : nullptr);

  basefold_pcs_internal::ValidateParamsOrThrow(params);
  if (!basefold_pcs_internal::IsPowerOfTwoLong(params.k0)) {
    return false;
  }
  const long kappa = basefold_pcs_internal::Log2ExactPowerOfTwoLong(params.k0);
  const long point_dim = params.d + kappa;
  if (static_cast<long>(z.size()) != point_dim) {
    return false;
  }
  if (num_queries < 0) {
    return false;
  }
  if (static_cast<long>(proof.commitments.roots_by_level.size()) !=
      params.d + 1) {
    return false;
  }
  if (static_cast<long>(proof.h_by_level.size()) != params.d) {
    return false;
  }
  if (params.d > 0 &&
      static_cast<long>(proof.query_multiproofs.size()) != params.d + 1) {
    return false;
  }

  if (proof.commitments.roots_by_level[static_cast<std::size_t>(params.d)] !=
      commitment_C) {
    return false;
  }

  const long n0 = CodewordLengthAtLevel(params, 0);
  if (proof.pi0_codeword.length() != n0) {
    return false;
  }

  if (params.d == 0) {
    if (!proof.h_by_level.empty()) {
      return false;
    }
    if (!proof.query_multiproofs.empty()) {
      return false;
    }
    if (MerkleCommitOracle(proof.pi0_codeword) != commitment_C) {
      return false;
    }

    NTL::vec_ZZ_pE msg0;
    if (!DecodeC0(msg0, proof.pi0_codeword, params)) {
      return false;
    }
    if (msg0.length() != params.k0) {
      return false;
    }

    return EvalMultilinearMonomialCoeffs(msg0, z) == claimed_y;
  }

  if (MerkleCommitOracle(proof.pi0_codeword) !=
      proof.commitments.roots_by_level[0]) {
    return false;
  }

  HashTranscript transcript = basefold_pcs_internal::MakeBaseFoldTranscript();
  basefold_pcs_internal::AbsorbPublicInput(transcript, commitment_C, z, claimed_y);
  AbsorbQuadraticPoly(
      transcript, proof.h_by_level[static_cast<std::size_t>(params.d - 1)]);

  std::vector<FieldElement> r_by_level(static_cast<std::size_t>(params.d));
  for (long i = params.d; i-- > 0;) {
    const FieldElement r_i =
        transcript.ChallengeFieldElement("r/" + std::to_string(i));
    r_by_level[static_cast<std::size_t>(i)] = r_i;

    transcript.AbsorbDigest(
        proof.commitments.roots_by_level[static_cast<std::size_t>(i)]);
    if (i > 0) {
      AbsorbQuadraticPoly(
          transcript, proof.h_by_level[static_cast<std::size_t>(i - 1)]);
    }
  }

  if (!CheckSumcheckRelations(proof.h_by_level, r_by_level, claimed_y)) {
    return false;
  }

  const FieldElement r0 = r_by_level[0];
  const FieldElement h1_r0 = proof.h_by_level[0].Eval(r0);

  NTL::vec_ZZ_pE msg0;
  if (!DecodeC0(msg0, proof.pi0_codeword, params)) {
    return false;
  }
  if (msg0.length() != params.k0) {
    return false;
  }

  FieldElement suffix_eq;
  NTL::set(suffix_eq);
  for (long i = 0; i < params.d; ++i) {
    suffix_eq *= EqFactor(z[static_cast<std::size_t>(kappa + i)],
                          r_by_level[static_cast<std::size_t>(i)]);
  }

  NTL::vec_ZZ_pE f_eval = msg0;
  for (long bit = 0; bit < kappa; ++bit) {
    const long step = 1L << bit;
    for (long mask = 0; mask < f_eval.length(); ++mask) {
      if (mask & step) {
        f_eval[mask] += f_eval[mask ^ step];
      }
    }
  }

  NTL::vec_ZZ_pE prefix_eq;
  prefix_eq.SetLength(f_eval.length());
  if (kappa == 0) {
    prefix_eq[0] = FieldElement(1);
  } else {
    prefix_eq.SetLength(1L << kappa);
    prefix_eq[0] = FieldElement(1);
    for (long var = 0; var < kappa; ++var) {
      const long old = 1L << var;
      const FieldElement zi = z[static_cast<std::size_t>(var)];
      const FieldElement f0 = FieldElement(1) - zi;
      const FieldElement f1 = zi;
      for (long mask = 0; mask < old; ++mask) {
        const FieldElement base = prefix_eq[mask];
        prefix_eq[mask] = base * f0;
        prefix_eq[mask + old] = base * f1;
      }
    }
  }

  FieldElement sum = FieldElement(0);
  for (long mask = 0; mask < f_eval.length(); ++mask) {
    sum += f_eval[mask] * prefix_eq[mask];
  }

  if (suffix_eq * sum != h1_r0) {
    return false;
  }

  IOPPChallenges challenges;
  challenges.alphas = r_by_level;

  const long n_last = CodewordLengthAtLevel(params, params.d - 1);
  std::vector<IOPPQueryPlan> query_plans(static_cast<std::size_t>(num_queries));
  for (long q = 0; q < num_queries; ++q) {
    const long mu =
        transcript.ChallengeIndex("mu/" + std::to_string(q), n_last);
    query_plans[static_cast<std::size_t>(q)] = MakeQueryPlan(mu, params);
  }

  const std::vector<std::vector<long>> requested_indices_by_tree =
      basefold_pcs_internal::CollectBaseQueryIndicesByTree(query_plans, params);
  ScopedTimer query_timer(prof ? &prof->verify_query_merkle_ns : nullptr,
                          prof ? &prof->verify_query_merkle_calls : nullptr);

  for (long tree_level = 0; tree_level <= params.d; ++tree_level) {
    const MerkleMultiproof &multiproof =
        proof.query_multiproofs[static_cast<std::size_t>(tree_level)];
    const std::vector<long> &expected_indices =
        requested_indices_by_tree[static_cast<std::size_t>(tree_level)];
    if (static_cast<long>(multiproof.values.length()) !=
        static_cast<long>(expected_indices.size())) {
      return false;
    }
    const long leaf_count = CodewordLengthAtLevel(params, tree_level);
    if (!MerkleVerifyMultiproof(
            proof.commitments.roots_by_level[static_cast<std::size_t>(tree_level)],
            leaf_count, expected_indices, multiproof)) {
      return false;
    }
  }

  auto verify_one_query = [&](long q) -> bool {
    const IOPPQueryPlan &plan = query_plans[static_cast<std::size_t>(q)];

    for (long i = params.d; i-- > 0;) {
      const long mu = plan.mu_by_level[static_cast<std::size_t>(i)];
      const long n_i = CodewordLengthAtLevel(params, i);
      if (mu < 0 || mu >= n_i) {
        return false;
      }

      const MerkleMultiproof &upper_multiproof =
          proof.query_multiproofs[static_cast<std::size_t>(i + 1)];
      const MerkleMultiproof &folded_multiproof =
          proof.query_multiproofs[static_cast<std::size_t>(i)];
      const std::vector<long> &upper_indices =
          requested_indices_by_tree[static_cast<std::size_t>(i + 1)];
      const std::vector<long> &folded_indices =
          requested_indices_by_tree[static_cast<std::size_t>(i)];

      const FieldElement *left = multiproof_replay::FindMultiproofValue(
          upper_indices, upper_multiproof.values.length(), mu,
          [&](std::size_t pos) {
            return &upper_multiproof.values[static_cast<long>(pos)];
          });
      const FieldElement *right = multiproof_replay::FindMultiproofValue(
          upper_indices, upper_multiproof.values.length(), mu + n_i,
          [&](std::size_t pos) {
            return &upper_multiproof.values[static_cast<long>(pos)];
          });
      const FieldElement *folded = multiproof_replay::FindMultiproofValue(
          folded_indices, folded_multiproof.values.length(), mu,
          [&](std::size_t pos) {
            return &folded_multiproof.values[static_cast<long>(pos)];
          });
      if (left == nullptr || right == nullptr || folded == nullptr) {
        return false;
      }

      FieldElement x1;
      FieldElement x2;
      FoldingPoints(x1, x2, params, i, mu);
      const FieldElement expected =
          EvalLineAt(challenges.alphas[static_cast<std::size_t>(i)], x1, *left,
                     x2, *right);
      if (expected != *folded) {
        return false;
      }
    }

    return true;
  };

  return basefold_pcs_internal::VerifyQueriesMaybeParallel(num_queries, prof,
                                                           verify_one_query);
}

}  // namespace basefold
