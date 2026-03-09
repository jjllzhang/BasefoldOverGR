#include "BaseFold/Z2kRingSwitchPCS.hpp"

#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pXFactoring.h>

#include <string>

using NTL::LogicError;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pEBak;
using NTL::ZZ_pX;
using NTL::ZZ_pBak;

namespace basefold {
namespace {

long Pow2LongOrThrow(long exponent, const char *what) {
  if (exponent < 0) {
    LogicError(what);
  }
  if (exponent >= static_cast<long>(8 * sizeof(long) - 1)) {
    LogicError(what);
  }
  return 1L << exponent;
}

bool IsPowerOfTwoZZ(const ZZ &value) {
  if (value <= 0) {
    return false;
  }
  ZZ t = value;
  while ((t % 2) == 0) {
    t /= 2;
  }
  return t == 1;
}

ZZ NormalizeModNonNegative(const ZZ &value, const ZZ &modulus) {
  ZZ out = value % modulus;
  if (out < 0) {
    out += modulus;
  }
  return out;
}

ZZ_pX ReduceZZpXModPrime(const ZZ_pX &poly_over_pk, const ZZ &prime) {
  ZZ_pX out;
  NTL::clear(out);
  const long degree = NTL::deg(poly_over_pk);
  for (long i = 0; i <= degree; ++i) {
    ZZ_p coeff_mod_prime;
    NTL::conv(coeff_mod_prime,
              NormalizeModNonNegative(NTL::rep(NTL::coeff(poly_over_pk, i)),
                                      prime));
    if (coeff_mod_prime != 0) {
      NTL::SetCoeff(out, i, coeff_mod_prime);
    }
  }
  out.normalize();
  return out;
}

void ValidateBasisDescriptorOrThrow(const RingSwitchBasisDescriptor &basis,
                                    long expected_dimension,
                                    const char *label,
                                    const char *func_name) {
  if (basis.kind != RingSwitchBasisKind::kPolynomial) {
    LogicError((std::string(func_name) + ": " + label +
                " must be the active polynomial basis")
                   .c_str());
  }
  if (basis.dimension != expected_dimension) {
    LogicError((std::string(func_name) + ": " + label +
                " dimension does not match current ZZ_pE degree")
                   .c_str());
  }
}

void ValidateBasicIrreducibilityModTwoOrThrow(const ZZ_pX &extension_modulus,
                                              long expected_degree,
                                              const char *func_name) {
  ZZ_pBak modulus_bak;
  modulus_bak.save();
  ZZ_pEBak extension_bak;
  extension_bak.save();

  const ZZ two(2);
  ZZ_p::init(two);
  const ZZ_pX reduced = ReduceZZpXModPrime(extension_modulus, two);
  if (NTL::deg(reduced) != expected_degree) {
    LogicError((std::string(func_name) +
                ": extension_modulus must stay full-degree after mod-2 reduction")
                   .c_str());
  }
  if (NTL::IterIrredTest(reduced) != 1) {
    LogicError((std::string(func_name) +
                ": extension_modulus must be basic irreducible modulo 2")
                   .c_str());
  }
}

}  // namespace

RingSwitchBasisDescriptor ActivePolynomialBasisDescriptor() {
  if (!ZZ_pE::initialized()) {
    LogicError("ActivePolynomialBasisDescriptor: ZZ_pE context must be initialized");
  }
  const long degree = ZZ_pE::degree();
  if (degree <= 0) {
    LogicError(
        "ActivePolynomialBasisDescriptor: current ZZ_pE degree must be positive");
  }
  RingSwitchBasisDescriptor basis;
  basis.kind = RingSwitchBasisKind::kPolynomial;
  basis.dimension = degree;
  return basis;
}

void ValidateCurrentZ2kRingContextOrThrow(const ZZ &base_modulus,
                                          const ZZ_pX &extension_modulus,
                                          long kappa) {
  if (kappa < 1) {
    LogicError("ValidateCurrentZ2kRingContextOrThrow: kappa must be >= 1");
  }
  if (base_modulus <= 1) {
    LogicError(
        "ValidateCurrentZ2kRingContextOrThrow: base_modulus must be > 1");
  }
  if (!IsPowerOfTwoZZ(base_modulus)) {
    LogicError(
        "ValidateCurrentZ2kRingContextOrThrow: base_modulus must be a power of two");
  }
  if (ZZ_p::modulus() != base_modulus) {
    LogicError(
        "ValidateCurrentZ2kRingContextOrThrow: current ZZ_p modulus does not match base_modulus");
  }
  if (!ZZ_pE::initialized()) {
    LogicError(
        "ValidateCurrentZ2kRingContextOrThrow: ZZ_pE context must be initialized");
  }

  const long expected_degree = Pow2LongOrThrow(
      kappa,
      "ValidateCurrentZ2kRingContextOrThrow: kappa is too large for long");
  if (ZZ_pE::degree() != expected_degree) {
    LogicError(
        "ValidateCurrentZ2kRingContextOrThrow: current ZZ_pE degree must equal 2^kappa");
  }
  if (NTL::deg(extension_modulus) != expected_degree) {
    LogicError(
        "ValidateCurrentZ2kRingContextOrThrow: extension_modulus degree must equal 2^kappa");
  }
  if (ZZ_pE::modulus().val() != extension_modulus) {
    LogicError(
        "ValidateCurrentZ2kRingContextOrThrow: current ZZ_pE modulus does not match extension_modulus");
  }
  ValidateBasicIrreducibilityModTwoOrThrow(
      extension_modulus, expected_degree,
      "ValidateCurrentZ2kRingContextOrThrow");
}

void ValidateRingSwitchPCSParamsOrThrow(const RingSwitchPCSParams &params) {
  if (params.kappa < 1) {
    LogicError("ValidateRingSwitchPCSParamsOrThrow: kappa must be >= 1");
  }
  if (params.ell < params.kappa) {
    LogicError("ValidateRingSwitchPCSParamsOrThrow: ell must be >= kappa");
  }

  ValidateCurrentZ2kRingContextOrThrow(params.base_modulus,
                                       params.extension_modulus, params.kappa);
  Z2kPCSBackendValidateParamsOrThrow(params.backend);

  const long expected_degree = ZZ_pE::degree();
  ValidateBasisDescriptorOrThrow(params.alpha_basis, expected_degree, "alpha_basis",
                                 "ValidateRingSwitchPCSParamsOrThrow");
  ValidateBasisDescriptorOrThrow(params.beta_basis, expected_degree, "beta_basis",
                                 "ValidateRingSwitchPCSParamsOrThrow");

  const long ell_prime = params.ell - params.kappa;
  if (params.ell_prime != ell_prime) {
    LogicError(
        "ValidateRingSwitchPCSParamsOrThrow: ell_prime must equal ell-kappa");
  }

  const long expected_backend_message_length = Pow2LongOrThrow(
      ell_prime,
      "ValidateRingSwitchPCSParamsOrThrow: ell-kappa is too large for long");
  if (Z2kPCSBackendMessageLength(params.backend) !=
      expected_backend_message_length) {
    LogicError(
        "ValidateRingSwitchPCSParamsOrThrow: backend message length must equal 2^(ell-kappa)");
  }
  if (Z2kPCSBackendPointDimension(params.backend) != ell_prime) {
    LogicError(
        "ValidateRingSwitchPCSParamsOrThrow: backend point dimension must equal ell-kappa");
  }
}

RingSwitchPCSParams RingSwitchPCSSetup(const RingSwitchPCSSetupInput &input) {
  RingSwitchPCSParams params;
  params.ell = input.ell;
  params.kappa = input.kappa;
  params.ell_prime = input.ell - input.kappa;
  params.base_modulus = input.base_modulus;
  params.extension_modulus = input.extension_modulus;
  params.alpha_basis = input.alpha_basis;
  params.beta_basis = input.beta_basis;
  params.backend = input.backend;
  ValidateRingSwitchPCSParamsOrThrow(params);
  return params;
}

}  // namespace basefold
