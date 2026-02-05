#ifndef BASEFOLD_FOLDABLECODE_HPP_
#define BASEFOLD_FOLDABLECODE_HPP_

#include <NTL/ZZ_pE.h>
#include <NTL/mat_ZZ_pE.h>
#include <NTL/vec_ZZ_pE.h>

#include <vector>

namespace basefold {

// Parameters for a (c, k0, d)-foldable linear code over F_{2^s} as defined in
// main.pdf (Definition 6 / Algorithm 1).
//
// Preconditions (NTL contexts):
//   - ZZ_p::init(2) has been called
//   - ZZ_pE::init(F) has been called for an irreducible polynomial F of degree s
//
// Shape invariants:
//   - k_i = k0 * 2^i
//   - n_i = c * k_i = c * k0 * 2^i
//   - G0 is k0 x n0
//   - diag_T[i] has length n_i and all entries are non-zero
struct FoldableCodeParams {
  long c = 0;
  long k0 = 0;
  long d = 0;

  // Fixed ζ ∈ F^× with ζ != 1 (so 1-ζ is invertible).
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

// Implements Algorithm 1 (Encd) from main.pdf.
// Input: msg of length k_d
// Output: out of length n_d
void EncodeFoldable(NTL::vec_ZZ_pE &out, const NTL::vec_ZZ_pE &msg,
                    const FoldableCodeParams &params);

}  // namespace basefold

#endif  // BASEFOLD_FOLDABLECODE_HPP_

