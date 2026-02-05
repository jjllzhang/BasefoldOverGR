#ifndef BASEFOLD_BASEFOLDPCS_HPP_
#define BASEFOLD_BASEFOLDPCS_HPP_

#include <NTL/ZZ_pE.h>
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
//
// This is the k0==1 specialization from Basefold_over_GR.pdf (Protocol 4).
struct BaseFoldPCSEvalProof {
  // Merkle roots for π_0..π_d.
  IOPPMerkleCommitments commitments;

  // Sumcheck messages h_{i+1}(X) stored by level i (0<=i<d).
  // In particular, h_by_level.back() is h_d.
  std::vector<QuadraticPoly> h_by_level;

  // Full π_0 (n0 = c when k0==1).
  Oracle pi0_full;

  // Openings for ℓ independent IOPP.query repetitions.
  std::vector<BaseFoldPCSQueryProof> query_proofs;
};

// Computes the PCS commitment C := MerkleRoot(π_{f,d}) where π_{f,d} = Encd(f⃗).
// Preconditions:
// - f_coeffs.length() == MessageLength(params)
MerkleRoot BaseFoldPCSCommit(const NTL::vec_ZZ_pE &f_coeffs,
                             const FoldableCodeParams &params);

// Proves an evaluation claim for a committed multilinear polynomial.
//
// Preconditions:
// - params.k0 == 1
// - z.size() == params.d
// - f_coeffs.length() == MessageLength(params) == 2^params.d
// - claimed_y == f(z)
BaseFoldPCSEvalProof BaseFoldPCSProveEval(const NTL::vec_ZZ_pE &f_coeffs,
                                          const std::vector<FieldElement> &z,
                                          const FieldElement &claimed_y,
                                          long num_queries,
                                          const FoldableCodeParams &params);

// Verifies an evaluation proof for commitment `C` at point `z` with value `y`.
//
// Preconditions:
// - params.k0 == 1
// - z.size() == params.d
bool BaseFoldPCSVerifyEval(const MerkleRoot &commitment_C,
                           const std::vector<FieldElement> &z,
                           const FieldElement &claimed_y,
                           long num_queries,
                           const BaseFoldPCSEvalProof &proof,
                           const FoldableCodeParams &params);

}  // namespace basefold

#endif  // BASEFOLD_BASEFOLDPCS_HPP_

