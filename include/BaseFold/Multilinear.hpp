#ifndef BASEFOLD_MULTILINEAR_HPP_
#define BASEFOLD_MULTILINEAR_HPP_

#include <NTL/ZZ_pE.h>
#include <NTL/vec_ZZ_pE.h>

#include <vector>

namespace basefold {

using FieldElement = NTL::ZZ_pE;
using FieldVec = NTL::vec_ZZ_pE;

// Evaluates a multilinear polynomial given by its monomial-basis coefficient
// vector `coeffs` of length 2^d at the point `point` ∈ R^d.
//
// Coefficient ordering follows the paper's bit-decomposition convention:
// index v encodes the monomial ∏_{j=1}^d X_j^{bit(v)[j]}, where bit(v)[1] is the
// least-significant bit. This matches the split order in EncodeFoldable.
FieldElement EvalMultilinearMonomialCoeffs(
    const FieldVec &coeffs, const std::vector<FieldElement> &point);

// Returns the i-th factor of the multilinear equality polynomial:
//   eq_z(X) = ∏_i (z_i X_i + (1-z_i)(1-X_i))
FieldElement EqFactor(const FieldElement &z_i, const FieldElement &x_i);

// Evaluates eq_z(x) for z,x ∈ R^d.
FieldElement EqPolynomial(const std::vector<FieldElement> &z,
                          const std::vector<FieldElement> &x);

}  // namespace basefold

#endif  // BASEFOLD_MULTILINEAR_HPP_

