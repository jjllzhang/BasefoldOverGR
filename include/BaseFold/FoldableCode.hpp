#ifndef BASEFOLD_FOLDABLECODE_HPP_
#define BASEFOLD_FOLDABLECODE_HPP_

#include <NTL/ZZ.h>
#include <NTL/ZZ_pE.h>
#include <NTL/mat_ZZ_pE.h>
#include <NTL/vec_ZZ_pE.h>

#include <vector>

namespace basefold {

// Parameters for a (c, k0, d)-foldable linear code over a commutative ring R
// represented by the current NTL contexts (ZZ_p / ZZ_pE).
//
// This matches the ζ-scaled foldable-code recursion in Basefold_over_GR.pdf
// (Definitions 7 and 8) and specializes to the binary-field construction in
// main.pdf (Algorithm 1).
//
// Preconditions (NTL contexts):
//   - For a field F_{p^r}: ZZ_p::init(p), ZZ_pE::init(F) (irreducible, deg=r).
//   - For a Galois ring GR(p^s, r): ZZ_p::init(p^s), ZZ_pE::init(F) where
//     F mod p is irreducible (deg=r).
//
// Shape invariants:
//   - k_i = k0 * 2^i
//   - n_i = c * k_i = c * k0 * 2^i
//   - G0 is k0 x n0
//   - diag_T[i] has length n_i and all entries are units in R
struct FoldableCodeParams {
  long c = 0;
  long k0 = 0;
  long d = 0;

  // Optional prime p used for unit checks in GR(p^s, r).
  // If p == 0, unit checks are skipped (non-zero checks still apply).
  NTL::ZZ p;

  // Fixed ζ ∈ R^× with (1-ζ) ∈ R^× (equivalently π(ζ) != 1 in the residue
  // field).
  NTL::ZZ_pE zeta;

  // Generator matrix of an [n0=ck0, k0]-MDS code.
  NTL::mat_ZZ_pE G0;

  // Diagonal entries of T0..T_{d-1}. Each vector must have length n_i.
  // (Index i corresponds to Ti in the paper.)
  std::vector<NTL::vec_ZZ_pE> diag_T;
};

// Returns k_d = k0 * 2^d.
long MessageLength(const FoldableCodeParams &params);

// Returns n_d = c * k_d.
long CodewordLength(const FoldableCodeParams &params);

// Implements the foldable-code encoding recursion Encd from Basefold_over_GR
// (and Algorithm 1 in main.pdf).
// Input: msg of length k_d
// Output: out of length n_d
void EncodeFoldable(NTL::vec_ZZ_pE &out, const NTL::vec_ZZ_pE &msg,
                    const FoldableCodeParams &params);

} // namespace basefold

#endif // BASEFOLD_FOLDABLECODE_HPP_
