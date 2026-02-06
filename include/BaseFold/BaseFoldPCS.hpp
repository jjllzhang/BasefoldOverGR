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
  std::vector<MerkleOpening> left;    // size == params.d
  std::vector<MerkleOpening> right;   // size == params.d
  std::vector<MerkleOpening> folded;  // size == params.d
};

// Non-interactive BaseFold PCS evaluation proof (Protocol 4 + Merkle+FS).
// Supports k0 = 2^κ (κ>=0) via the BaseFold paper's Remark 3 adaptation:
// the IOPP recursion has depth params.d, while the committed multilinear
// polynomial has dimension (params.d + κ).
struct BaseFoldPCSEvalProof {
  // Merkle roots for π_0..π_d.
  IOPPMerkleCommitments commitments;

  // Sumcheck messages h_{i+1}(X) stored by level i (0<=i<d).
  // In particular, h_by_level.back() is h_d.
  std::vector<QuadraticPoly> h_by_level;

  // Full π_0 (length n0 = c*k0).
  Oracle pi0_full;

  // Openings for ℓ independent IOPP.query repetitions.
  std::vector<BaseFoldPCSQueryProof> query_proofs;
};

// Challenge-domain configuration for opening/eval.
//
// Default behavior (`use_extension_challenges=false`) matches the current
// implementation: all Fiat-Shamir field challenges are sampled in the ambient
// ring/field represented by ZZ_pE.
//
// When `use_extension_challenges=true`, `challenge_extension_modulus` is
// expected to define the outer extension E(U) over the current ZZ_pE context.
// The first API version only exposes this configuration and a Prove/Verify
// skeleton entrypoint.
struct BaseFoldPCSChallengeConfig {
  bool use_extension_challenges = false;
  NTL::ZZ_pEX challenge_extension_modulus;
};

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
//   extension-challenge skeleton path (to be completed in follow-up patches).
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
