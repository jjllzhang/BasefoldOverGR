#ifndef BASEFOLD_Z2K_RING_SWITCH_PCS_HPP_
#define BASEFOLD_Z2K_RING_SWITCH_PCS_HPP_

#include <NTL/ZZ.h>
#include <NTL/ZZ_pX.h>

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

RingSwitchBasisDescriptor ActivePolynomialBasisDescriptor();

void ValidateCurrentZ2kRingContextOrThrow(const NTL::ZZ &base_modulus,
                                          const NTL::ZZ_pX &extension_modulus,
                                          long kappa);

void ValidateRingSwitchPCSParamsOrThrow(const RingSwitchPCSParams &params);

RingSwitchPCSParams RingSwitchPCSSetup(const RingSwitchPCSSetupInput &input);

}  // namespace basefold

#endif  // BASEFOLD_Z2K_RING_SWITCH_PCS_HPP_
