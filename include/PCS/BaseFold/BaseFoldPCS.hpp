#ifndef BASEFOLD_BASEFOLDPCS_HPP_
#define BASEFOLD_BASEFOLDPCS_HPP_

#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/vec_ZZ_pE.h>

#include <vector>

#include "PCS/BaseFold/FoldableCode.hpp"
#include "PCS/BaseFold/IOPP.hpp"
#include "PCS/Common/Sumcheck.hpp"

namespace basefold {

// Degree-2 univariate polynomial over the outer challenge extension ring.
struct ExtensionQuadraticPoly {
  NTL::ZZ_pEX a0;
  NTL::ZZ_pEX a1;
  NTL::ZZ_pEX a2;
};

// A pruned Merkle multiproof for multiple leaves in the same extension oracle.
struct ExtensionMerkleMultiproof {
  std::vector<long> queried_indices;
  std::vector<NTL::ZZ_pEX> values;
  std::vector<Digest> sibling_hashes;
};

// Reusable top-level commitment artifacts for PCS eval proofs.
//
// These let callers precompute the top codeword/oracle pi_d and its Merkle tree,
// then build evaluation proofs without re-running the top encode/commit stage.
// `base_sumcheck_precomputation` is prover-local cache data and never enters the
// proof or verifier API.
struct ExtensionCommitRoundLevelPrecomputation {
  NTL::vec_ZZ_pE inv_denoms;
};

struct ExtensionCommitRoundPrecomputation {
  std::vector<ExtensionCommitRoundLevelPrecomputation> levels;
};

struct BaseFoldPCSCommitArtifacts {
  Oracle pi_d;
  MerkleTree merkle_d;
  MerkleRoot root_d;
  SumcheckMonomialPrecomputation base_sumcheck_precomputation;
  ExtensionCommitRoundPrecomputation extension_commit_precomputation;
  std::vector<NTL::ZZ_pEX> extension_lifted_pi_d;
};

// Additional proof payload used by the extension-challenge path.
//
// Commitment to π_d remains in the base ring (Merkle over ZZ_pE values), while
// sumcheck/folding arithmetic below runs in the extension ring E(U).
struct BaseFoldPCSExtensionProofData {
  bool has_extension_payload = false;

  // Merkle roots for extension oracles π_0..π_{d-1}.
  // Top-level π_d commitment is still in BaseFoldPCSEvalProof::commitments.
  std::vector<MerkleRoot> roots_by_level;

  // h_{i+1}(X) for i=0..d-1, represented in E(U).
  std::vector<ExtensionQuadraticPoly> h_by_level;

  // Optional explicit r_i for i=0..d-1, represented in E(U).
  // Verifier can recompute these from transcript; compact proofs may omit them.
  std::vector<NTL::ZZ_pEX> r_by_level;

  // Monomial coefficients of f(·, r_suffix) on the first κ variables.
  std::vector<NTL::ZZ_pEX> msg0_coeffs;

  // Full π_0 in E(U), derived from extension-ring folding.
  std::vector<NTL::ZZ_pEX> pi0_codeword;

  // Shared base openings into the top-level base commitment C = MerkleRoot(π_d).
  MerkleMultiproof base_top_query_multiproof;

  // Shared extension openings grouped by extension oracle level π_0..π_{d-1}.
  std::vector<ExtensionMerkleMultiproof> query_multiproofs;
};

// Non-interactive BaseFold PCS evaluation proof (Protocol 4 + Merkle+FS).
// Supports k0 = 2^κ (κ>=0) via the BaseFold paper's Remark 3 adaptation:
// the IOPP recursion has depth params.d, while the committed multilinear
// polynomial has dimension (params.d + κ).
struct BaseFoldPCSEvalProof {
  // Base commitment roots.
  // Base-challenge path stores roots for π_0..π_d.
  // Extension-challenge compact path may store only root(π_d) (or none).
  IOPPMerkleCommitments commitments;

  // Base sumcheck messages.
  // Base-challenge path: size == d.
  // Extension-challenge compact path: may be empty (use extension.h_by_level).
  std::vector<QuadraticPoly> h_by_level;

  // Base π_0 payload.
  // Base-challenge path: length n0 = c*k0.
  // Extension-challenge compact path: may be empty.
  Oracle pi0_codeword;

  // Shared base openings grouped by oracle level π_0..π_d.
  std::vector<MerkleMultiproof> query_multiproofs;

  // Optional extension-ring payload used when
  // BaseFoldPCSChallengeConfig::use_extension_challenges=true.
  BaseFoldPCSExtensionProofData extension;
};

// Challenge-domain configuration for opening/eval.
//
// Default behavior (`use_extension_challenges=false`) matches the current
// implementation: all Fiat-Shamir field challenges are sampled in the ambient
// ring/field represented by ZZ_pE.
//
// When `use_extension_challenges=true`, `challenge_extension_modulus` is
// expected to define the outer extension E(U) over the current ZZ_pE context.
struct BaseFoldPCSChallengeConfig {
  bool use_extension_challenges = false;
  NTL::ZZ_pEX challenge_extension_modulus;
};

// Runtime tuning knobs for verifier query-level parallelization.
// Defaults can also be provided via environment variables:
// - BASEFOLD_VERIFY_QUERY_QUERIES_PER_THREAD
// - BASEFOLD_VERIFY_QUERY_PARALLEL_THRESHOLD
// - BASEFOLD_VERIFY_QUERY_MAX_THREADS
struct VerifierQueryParallelConfig {
  long queries_per_thread = 1;
  long min_queries_for_parallelism = 2;
  int max_threads = 8;
};

// Loads verifier query parallel config from environment variables and applies it.
void ResetVerifierQueryParallelConfigFromEnv();

// Applies verifier query parallel config for the current process.
void SetVerifierQueryParallelConfig(const VerifierQueryParallelConfig &cfg);

// Returns the currently active verifier query parallel config.
VerifierQueryParallelConfig GetVerifierQueryParallelConfig();

// Runtime tuning knobs for prover commit-round element parallelism.
// Defaults can also be provided via environment variables:
// - BASEFOLD_PROVER_COMMIT_BASE_ELEMENTS_PER_THREAD
// - BASEFOLD_PROVER_COMMIT_EXT_ELEMENTS_PER_THREAD
//
// These values are passed as the `parallel_threshold` input to the existing
// round-local OpenMP helper, so smaller values allow threads to kick in at
// smaller round sizes. When `ext_elements_per_thread` is left at the default,
// the implementation may auto-tune it more aggressively on composite-modulus
// ring contexts; explicit env/CLI/process overrides disable that auto-tuning.
struct ProverCommitParallelConfig {
  long base_elements_per_thread = 4096;
  long ext_elements_per_thread = 128;
  bool base_elements_per_thread_overridden = false;
  bool ext_elements_per_thread_overridden = false;
};

// Loads prover commit-round parallel config from environment variables and
// applies it.
void ResetProverCommitParallelConfigFromEnv();

// Applies prover commit-round parallel config for the current process.
void SetProverCommitParallelConfig(const ProverCommitParallelConfig &cfg);

// Returns the currently active prover commit-round parallel config.
ProverCommitParallelConfig GetProverCommitParallelConfig();

// Computes the PCS commitment C := MerkleRoot(π_{f,d}) where π_{f,d} = Encd(f⃗).
// Preconditions:
// - f_coeffs.length() == MessageLength(params)
MerkleRoot BaseFoldPCSCommit(const NTL::vec_ZZ_pE &f_coeffs,
                             const FoldableCodeParams &params);

// Builds reusable top-level commitment artifacts for π_d.
//
// Preconditions:
// - f_coeffs.length() == MessageLength(params)
BaseFoldPCSCommitArtifacts BaseFoldPCSBuildCommitArtifacts(
    const NTL::vec_ZZ_pE &f_coeffs, const FoldableCodeParams &params);

// Unchecked variant of BaseFoldPCSBuildCommitArtifacts intended for hot paths.
BaseFoldPCSCommitArtifacts BaseFoldPCSBuildCommitArtifactsUnchecked(
    const NTL::vec_ZZ_pE &f_coeffs, const FoldableCodeParams &params);

// Proves an evaluation claim for a committed multilinear polynomial.
//
// Preconditions:
// - params.k0 is a power of two
// - z.size() == params.d + log2(params.k0)
// - f_coeffs.length() == MessageLength(params) == params.k0 * 2^params.d
// - claimed_y == f(z)
BaseFoldPCSEvalProof BaseFoldPCSProveEval(const NTL::vec_ZZ_pE &f_coeffs,
                                          const std::vector<FieldElement> &z,
                                          const FieldElement &claimed_y,
                                          long num_queries,
                                          const FoldableCodeParams &params);

// Unchecked variant of BaseFoldPCSProveEval intended for benchmarking/hot paths.
//
// Differences vs BaseFoldPCSProveEval:
// - Does NOT validate params (so it skips unit checks and dimension checks).
// - Does NOT check z/f_coeffs lengths.
// - Does NOT check claimed_y == f(z).
//
// Caller is responsible for satisfying the same preconditions as
// BaseFoldPCSProveEval.
BaseFoldPCSEvalProof BaseFoldPCSProveEvalUnchecked(
    const NTL::vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const FoldableCodeParams &params);

// Proves an evaluation claim assuming the top codeword/oracle π_d is already
// encoded and Merkle-committed via `commit_artifacts`.
BaseFoldPCSEvalProof BaseFoldPCSProveEvalFromCommittedTopOracle(
    const NTL::vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const FoldableCodeParams &params,
    const BaseFoldPCSCommitArtifacts &commit_artifacts);

// Unchecked variant of BaseFoldPCSProveEvalFromCommittedTopOracle intended for
// benchmarking/hot paths.
BaseFoldPCSEvalProof BaseFoldPCSProveEvalFromCommittedTopOracleUnchecked(
    const NTL::vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const FoldableCodeParams &params,
    const BaseFoldPCSCommitArtifacts &commit_artifacts);

// Verifies an evaluation proof for commitment `C` at point `z` with value `y`.
//
// Preconditions:
// - params.k0 is a power of two
// - z.size() == params.d + log2(params.k0)
bool BaseFoldPCSVerifyEval(const MerkleRoot &commitment_C,
                           const std::vector<FieldElement> &z,
                           const FieldElement &claimed_y,
                           long num_queries,
                           const BaseFoldPCSEvalProof &proof,
                           const FoldableCodeParams &params);

// Configurable entrypoints for challenge-domain experiments.
//
// - If `challenge_cfg.use_extension_challenges == false`, these forward to the
//   legacy BaseFoldPCSProveEval/BaseFoldPCSProveEvalUnchecked/VerifyEval.
// - If `challenge_cfg.use_extension_challenges == true`, they run the first
//   extension arithmetic path: challenges are sampled in E(U), and sumcheck /
//   folding checks run in E(U), while the committed codeword π_d stays in the
//   base ring.
BaseFoldPCSEvalProof BaseFoldPCSProveEvalWithChallengeConfig(
    const NTL::vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const FoldableCodeParams &params,
    const BaseFoldPCSChallengeConfig &challenge_cfg);

BaseFoldPCSEvalProof BaseFoldPCSProveEvalWithChallengeConfigUnchecked(
    const NTL::vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const FoldableCodeParams &params,
    const BaseFoldPCSChallengeConfig &challenge_cfg);

BaseFoldPCSEvalProof BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracle(
    const NTL::vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const FoldableCodeParams &params,
    const BaseFoldPCSCommitArtifacts &commit_artifacts,
    const BaseFoldPCSChallengeConfig &challenge_cfg);

BaseFoldPCSEvalProof
BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracleUnchecked(
    const NTL::vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const FoldableCodeParams &params,
    const BaseFoldPCSCommitArtifacts &commit_artifacts,
    const BaseFoldPCSChallengeConfig &challenge_cfg);

bool BaseFoldPCSVerifyEvalWithChallengeConfig(
    const MerkleRoot &commitment_C, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const BaseFoldPCSEvalProof &proof, const FoldableCodeParams &params,
    const BaseFoldPCSChallengeConfig &challenge_cfg);

// STIR-style protocol facade that separates setup/commit/prove/verify without
// changing the underlying proof system implementation.
struct BaseFoldPCSSetupInput {
  FoldableCodeParams params;
  BaseFoldPCSChallengeConfig challenge_config;
};

// Prover-local state produced by Commit and consumed by Prove.
struct BaseFoldPCSCommittedWitness {
  MerkleRoot commitment;
  NTL::vec_ZZ_pE f_coeffs;
  BaseFoldPCSCommitArtifacts commit_artifacts;
};

class BaseFoldPCSProver {
 public:
  explicit BaseFoldPCSProver(const FoldableCodeParams &params)
      : params_(params) {
    (void)MessageLength(params_);
  }

  explicit BaseFoldPCSProver(const BaseFoldPCSSetupInput &input)
      : params_(input.params), challenge_config_(input.challenge_config) {
    (void)MessageLength(params_);
  }

  BaseFoldPCSCommittedWitness Commit(const NTL::vec_ZZ_pE &f_coeffs) const {
    BaseFoldPCSCommittedWitness committed;
    committed.f_coeffs = f_coeffs;
    committed.commit_artifacts = BaseFoldPCSBuildCommitArtifacts(f_coeffs, params_);
    committed.commitment = committed.commit_artifacts.root_d;
    return committed;
  }

  BaseFoldPCSEvalProof Prove(const BaseFoldPCSCommittedWitness &committed,
                             const std::vector<FieldElement> &z,
                             const FieldElement &claimed_y,
                             long num_queries) const {
    if (committed.commitment != committed.commit_artifacts.root_d) {
      NTL::LogicError(
          "BaseFoldPCSProver::Prove: committed witness has inconsistent commitment");
    }
    return BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracle(
        committed.f_coeffs, z, claimed_y, num_queries, params_,
        committed.commit_artifacts, challenge_config_);
  }

  const FoldableCodeParams &params() const { return params_; }

  const BaseFoldPCSChallengeConfig &challenge_config() const {
    return challenge_config_;
  }

 private:
  FoldableCodeParams params_;
  BaseFoldPCSChallengeConfig challenge_config_;
};

class BaseFoldPCSVerifier {
 public:
  explicit BaseFoldPCSVerifier(const FoldableCodeParams &params)
      : params_(params) {
    (void)MessageLength(params_);
  }

  explicit BaseFoldPCSVerifier(const BaseFoldPCSSetupInput &input)
      : params_(input.params), challenge_config_(input.challenge_config) {
    (void)MessageLength(params_);
  }

  bool Verify(const MerkleRoot &commitment,
              const std::vector<FieldElement> &z,
              const FieldElement &claimed_y, long num_queries,
              const BaseFoldPCSEvalProof &proof) const {
    return BaseFoldPCSVerifyEvalWithChallengeConfig(
        commitment, z, claimed_y, num_queries, proof, params_,
        challenge_config_);
  }

  const FoldableCodeParams &params() const { return params_; }

  const BaseFoldPCSChallengeConfig &challenge_config() const {
    return challenge_config_;
  }

 private:
  FoldableCodeParams params_;
  BaseFoldPCSChallengeConfig challenge_config_;
};

struct BaseFoldPCSSetupOutput {
  BaseFoldPCSProver prover;
  BaseFoldPCSVerifier verifier;
};

inline BaseFoldPCSSetupOutput BaseFoldPCSSetup(
    const FoldableCodeParams &params) {
  return {BaseFoldPCSProver(params), BaseFoldPCSVerifier(params)};
}

inline BaseFoldPCSSetupOutput BaseFoldPCSSetup(
    const BaseFoldPCSSetupInput &input) {
  return {BaseFoldPCSProver(input), BaseFoldPCSVerifier(input)};
}

}  // namespace basefold

#endif  // BASEFOLD_BASEFOLDPCS_HPP_
