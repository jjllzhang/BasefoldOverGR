#ifndef BASEFOLD_BASEFOLDPCS_HPP_
#define BASEFOLD_BASEFOLDPCS_HPP_

#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/vec_ZZ_pE.h>

#include <vector>

#include "BaseFold/FoldableCode.hpp"
#include "BaseFold/IOPP.hpp"
#include "BaseFold/Sumcheck.hpp"

namespace basefold {

// A single IOPP.query repetition's Merkle openings (without duplicating π0).
struct BaseFoldPCSQueryProof {
  // Base-challenge path: left/right/folded each have size == params.d.
  //
  // Extension-challenge path (compact payload): this is used only for top-level
  // base openings into commitment C=MerkleRoot(π_d), so left/right have size 1
  // and folded is empty.
  std::vector<MerkleOpening> left;
  std::vector<MerkleOpening> right;
  std::vector<MerkleOpening> folded;
};

// Degree-2 univariate polynomial over the outer challenge extension ring.
struct ExtensionQuadraticPoly {
  NTL::ZZ_pEX a0;
  NTL::ZZ_pEX a1;
  NTL::ZZ_pEX a2;
};

// One opened leaf in an extension-ring Merkle-committed oracle.
struct ExtensionMerkleOpening {
  long index = 0;
  NTL::ZZ_pEX value;
  MerkleAuthPath auth_path;
};

// A pruned Merkle multiproof for multiple leaves in the same extension oracle.
struct ExtensionMerkleMultiproof {
  std::vector<long> queried_indices;
  std::vector<NTL::ZZ_pEX> values;
  std::vector<Digest> sibling_hashes;
};

// One query repetition's values along the extension-ring folding chain.
struct BaseFoldPCSQueryProofExtension {
  std::vector<ExtensionMerkleOpening> left;    // size == params.d
  std::vector<ExtensionMerkleOpening> right;   // size == params.d
  std::vector<ExtensionMerkleOpening> folded;  // size == params.d
};

// Additional proof payload used by the extension-challenge path.
//
// Commitment to π_d remains in the base ring (Merkle over ZZ_pE values), while
// sumcheck/folding arithmetic below runs in the extension ring E(U).
struct BaseFoldPCSExtensionProofData {
  bool enabled = false;

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
  std::vector<NTL::ZZ_pEX> pi0_full;

  // Shared base openings into the top-level base commitment C = MerkleRoot(π_d).
  // Compact multiproof layout may use this and leave BaseFoldPCSEvalProof::
  // query_proofs empty.
  MerkleMultiproof base_top_query_multiproof;

  // Shared extension openings grouped by extension oracle level π_0..π_{d-1}.
  // Compact multiproof layout may populate this and leave query_proofs empty.
  std::vector<ExtensionMerkleMultiproof> query_multiproofs;

  // Legacy per-query extension openings.
  std::vector<BaseFoldPCSQueryProofExtension> query_proofs;
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
  Oracle pi0_full;

  // Shared base openings grouped by oracle level π_0..π_d.
  // Base-challenge path may populate this and leave query_proofs empty.
  std::vector<MerkleMultiproof> query_multiproofs;

  // Base openings for ℓ independent IOPP.query repetitions.
  // Used by the extension-challenge compact path, and accepted as a legacy
  // base-challenge layout when query_multiproofs is empty.
  std::vector<BaseFoldPCSQueryProof> query_proofs;

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
  long parallel_query_threshold = 2;
  int max_threads = 8;
};

// Loads verifier query parallel config from environment variables and applies it.
void ResetVerifierQueryParallelConfigFromEnv();

// Applies verifier query parallel config for the current process.
void SetVerifierQueryParallelConfig(const VerifierQueryParallelConfig &cfg);

// Returns the currently active verifier query parallel config.
VerifierQueryParallelConfig GetVerifierQueryParallelConfig();

// Computes the PCS commitment C := MerkleRoot(π_{f,d}) where π_{f,d} = Encd(f⃗).
// Preconditions:
// - f_coeffs.length() == MessageLength(params)
MerkleRoot BaseFoldPCSCommit(const NTL::vec_ZZ_pE &f_coeffs,
                             const FoldableCodeParams &params);

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

bool BaseFoldPCSVerifyEvalWithChallengeConfig(
    const MerkleRoot &commitment_C, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const BaseFoldPCSEvalProof &proof, const FoldableCodeParams &params,
    const BaseFoldPCSChallengeConfig &challenge_cfg);

}  // namespace basefold

#endif  // BASEFOLD_BASEFOLDPCS_HPP_
