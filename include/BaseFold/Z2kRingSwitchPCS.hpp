#ifndef BASEFOLD_Z2K_RING_SWITCH_PCS_HPP_
#define BASEFOLD_Z2K_RING_SWITCH_PCS_HPP_

#include <NTL/ZZ.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vec_ZZ_pE.h>

#include <vector>

#include "BaseFold/Z2kPCSBackend.hpp"

namespace basefold {

enum class RingSwitchBasisKind {
  kPolynomial = 0,
};

struct RingSwitchBasisDescriptor {
  RingSwitchBasisKind kind = RingSwitchBasisKind::kPolynomial;
  long dimension = 0;
};

struct RingSwitchPCSSetupInput {
  long ell = 0;
  long kappa = 0;
  NTL::ZZ base_modulus;
  NTL::ZZ_pX extension_modulus;
  RingSwitchBasisDescriptor alpha_basis;
  RingSwitchBasisDescriptor beta_basis;
  Z2kPCSBackendHandle backend;
};

struct RingSwitchPCSParams {
  long ell = 0;
  long kappa = 0;
  long ell_prime = 0;
  NTL::ZZ base_modulus;
  NTL::ZZ_pX extension_modulus;
  RingSwitchBasisDescriptor alpha_basis;
  RingSwitchBasisDescriptor beta_basis;
  Z2kPCSBackendHandle backend;
};

struct RingSwitchComponentTensor {
  long basis_dimension = 0;
  long ell_prime = 0;
  // Row-major storage of A_{u||w}: row `u`, then column `w`.
  NTL::vec_ZZ_pE a_by_u_then_w;
  // Boolean-hypercube value table of the verifier-known polynomial r, with
  // linear index index(u || w) = u + 2^kappa * w.
  NTL::vec_ZZ_pE r_table;
  // Monomial-basis coefficients of r, cached for reuse with the current
  // backend core and existing multilinear helpers.
  NTL::vec_ZZ_pE r_monomial_coeffs;
};

RingSwitchBasisDescriptor ActivePolynomialBasisDescriptor();

void ValidateCurrentZ2kRingContextOrThrow(const NTL::ZZ &base_modulus,
                                          const NTL::ZZ_pX &extension_modulus,
                                          long kappa);

void ValidateRingSwitchPCSParamsOrThrow(const RingSwitchPCSParams &params);

RingSwitchPCSParams RingSwitchPCSSetup(const RingSwitchPCSSetupInput &input);

// Converts a length-2^d Boolean-hypercube value table (Lagrange-basis
// coefficients over {0,1}^d) into the corresponding monomial-basis
// coefficients.
NTL::vec_ZZ_pE BooleanHypercubeTableToMonomialCoeffs(
    const NTL::vec_ZZ_pE &table_values);

// Packs the paper-style Boolean-hypercube table of t over Z_{2^k} into the
// paper-style Boolean-hypercube table of t' over GR(2^k, 2^kappa).
NTL::vec_ZZ_pE PackZ2kCoeffsToGREvals(const RingSwitchPCSParams &params,
                                      const NTL::vec_ZZ_pE &t_table);

NTL::vec_ZZ_pE DecomposeGRElementToBaseCoeffsPolynomialBasis(
    const RingSwitchPCSParams &params, const NTL::ZZ_pE &element);

NTL::vec_ZZ_pE DecomposeGRElementToBaseCoeffs(
    const RingSwitchPCSParams &params, const NTL::ZZ_pE &element,
    const RingSwitchBasisDescriptor &basis);

RingSwitchComponentTensor BuildRingSwitchComponentTensor(
    const RingSwitchPCSParams &params,
    const std::vector<NTL::ZZ_pE> &r_suffix);

}  // namespace basefold

#endif  // BASEFOLD_Z2K_RING_SWITCH_PCS_HPP_
