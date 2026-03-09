#include "BaseFold/Z2kRingSwitchPCS.hpp"

#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pXFactoring.h>

#include <limits>
#include <string>

#include "BaseFold/Multilinear.hpp"

using NTL::LogicError;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pEBak;
using NTL::ZZ_pX;
using NTL::ZZ_pBak;
using NTL::vec_ZZ_pE;

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

bool IsPowerOfTwoLong(long value) {
  return value > 0 && (value & (value - 1)) == 0;
}

long Log2ExactPowerOfTwoLongOrThrow(long value, const char *what) {
  if (!IsPowerOfTwoLong(value)) {
    LogicError(what);
  }
  long out = 0;
  while (value > 1) {
    value >>= 1;
    ++out;
  }
  return out;
}

ZZ_pE BaseRingConstant(const ZZ_p &value) {
  ZZ_pX poly;
  if (value != 0) {
    NTL::SetCoeff(poly, 0, value);
  }
  ZZ_pE out;
  NTL::conv(out, poly);
  return out;
}

ZZ_pE BaseRingConstant(long value) {
  return BaseRingConstant(NTL::to_ZZ_p(value));
}

bool IsBaseRingConstant(const ZZ_pE &value) {
  const ZZ_pX poly = NTL::rep(value);
  const long degree = NTL::deg(poly);
  for (long i = 1; i <= degree; ++i) {
    if (NTL::coeff(poly, i) != 0) {
      return false;
    }
  }
  return true;
}

void ValidateBaseRingConstantOrThrow(const ZZ_pE &value, const char *label,
                                     long index, const char *func_name) {
  if (!IsBaseRingConstant(value)) {
    LogicError((std::string(func_name) + ": " + label + "[" +
                std::to_string(index) + "] must be a base-ring constant")
                   .c_str());
  }
}

void ValidateBaseRingVectorOrThrow(const vec_ZZ_pE &values, const char *label,
                                   const char *func_name) {
  for (long i = 0; i < values.length(); ++i) {
    ValidateBaseRingConstantOrThrow(values[i], label, i, func_name);
  }
}

std::vector<ZZ_pE> BooleanPointFromIndex(long index, long dimension) {
  std::vector<ZZ_pE> point(static_cast<std::size_t>(dimension));
  for (long i = 0; i < dimension; ++i) {
    point[static_cast<std::size_t>(i)] = BaseRingConstant((index >> i) & 1L);
  }
  return point;
}

long CheckedMultiplyLong(long a, long b, const char *what) {
  if (a < 0 || b < 0) {
    LogicError(what);
  }
  if (a == 0 || b == 0) {
    return 0;
  }
  if (a > std::numeric_limits<long>::max() / b) {
    LogicError(what);
  }
  return a * b;
}

vec_ZZ_pE DecomposeGRElementToBaseCoeffsPolynomialBasisUnchecked(
    long basis_dimension, const ZZ_pE &element) {
  vec_ZZ_pE coeffs;
  coeffs.SetLength(basis_dimension);
  const ZZ_pX poly = NTL::rep(element);
  for (long i = 0; i < basis_dimension; ++i) {
    coeffs[i] = BaseRingConstant(NTL::coeff(poly, i));
  }
  return coeffs;
}

struct PackedCommitInputs {
  vec_ZZ_pE t_packed_table;
  vec_ZZ_pE t_packed_monomial_coeffs;
};

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

PackedCommitInputs BuildPackedCommitInputs(const RingSwitchPCSParams &params,
                                           const vec_ZZ_pE &t_table) {
  PackedCommitInputs out;
  out.t_packed_table = PackZ2kCoeffsToGREvals(params, t_table);
  out.t_packed_monomial_coeffs =
      BooleanHypercubeTableToMonomialCoeffs(out.t_packed_table);
  return out;
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

vec_ZZ_pE BooleanHypercubeTableToMonomialCoeffs(const vec_ZZ_pE &table_values) {
  const long n = table_values.length();
  if (!IsPowerOfTwoLong(n)) {
    LogicError(
        "BooleanHypercubeTableToMonomialCoeffs: table length must be a power of two");
  }
  const long dimension = Log2ExactPowerOfTwoLongOrThrow(
      n,
      "BooleanHypercubeTableToMonomialCoeffs: table length must be a power of two");
  vec_ZZ_pE coeffs = table_values;
  for (long bit = 0; bit < dimension; ++bit) {
    const long step = 1L << bit;
    for (long mask = 0; mask < n; ++mask) {
      if (mask & step) {
        coeffs[mask] -= coeffs[mask ^ step];
      }
    }
  }
  return coeffs;
}

vec_ZZ_pE PackZ2kCoeffsToGREvals(const RingSwitchPCSParams &params,
                                 const vec_ZZ_pE &t_table) {
  ValidateRingSwitchPCSParamsOrThrow(params);
  const long expected_length = Pow2LongOrThrow(
      params.ell, "PackZ2kCoeffsToGREvals: ell is too large for long");
  if (t_table.length() != expected_length) {
    LogicError("PackZ2kCoeffsToGREvals: t_table length must equal 2^ell");
  }
  ValidateBaseRingVectorOrThrow(t_table, "t_table",
                                "PackZ2kCoeffsToGREvals");

  const long basis_dimension = params.beta_basis.dimension;
  const long packed_length = Pow2LongOrThrow(
      params.ell_prime,
      "PackZ2kCoeffsToGREvals: ell_prime is too large for long");
  vec_ZZ_pE packed;
  packed.SetLength(packed_length);

  if (params.beta_basis.kind != RingSwitchBasisKind::kPolynomial) {
    LogicError(
        "PackZ2kCoeffsToGREvals: beta_basis must be the active polynomial basis");
  }

  for (long w = 0; w < packed_length; ++w) {
    ZZ_pX poly;
    for (long v = 0; v < basis_dimension; ++v) {
      const long coeff_index = v + w * basis_dimension;
      const ZZ_p constant_coeff = NTL::coeff(NTL::rep(t_table[coeff_index]), 0);
      if (constant_coeff != 0) {
        NTL::SetCoeff(poly, v, constant_coeff);
      }
    }
    NTL::conv(packed[w], poly);
  }

  return packed;
}

MerkleRoot RingSwitchPCSCommit(const RingSwitchPCSParams &params,
                               const vec_ZZ_pE &t_table) {
  return RingSwitchPCSBuildCommitArtifacts(params, t_table).commitment;
}

RingSwitchPCSCommitArtifacts RingSwitchPCSBuildCommitArtifacts(
    const RingSwitchPCSParams &params, const vec_ZZ_pE &t_table) {
  ValidateRingSwitchPCSParamsOrThrow(params);
  const PackedCommitInputs packed = BuildPackedCommitInputs(params, t_table);

  RingSwitchPCSCommitArtifacts out;
  out.t_packed_table = packed.t_packed_table;
  out.t_packed_monomial_coeffs = packed.t_packed_monomial_coeffs;
  out.backend_commit_artifacts =
      Z2kPCSBackendBuildCommitArtifacts(params.backend,
                                        out.t_packed_monomial_coeffs);
  out.commitment = out.backend_commit_artifacts.commitment;
  return out;
}

vec_ZZ_pE DecomposeGRElementToBaseCoeffsPolynomialBasis(
    const RingSwitchPCSParams &params, const ZZ_pE &element) {
  ValidateRingSwitchPCSParamsOrThrow(params);
  return DecomposeGRElementToBaseCoeffsPolynomialBasisUnchecked(
      params.alpha_basis.dimension, element);
}

vec_ZZ_pE DecomposeGRElementToBaseCoeffs(
    const RingSwitchPCSParams &params, const ZZ_pE &element,
    const RingSwitchBasisDescriptor &basis) {
  ValidateRingSwitchPCSParamsOrThrow(params);
  ValidateBasisDescriptorOrThrow(basis, ZZ_pE::degree(), "basis",
                                 "DecomposeGRElementToBaseCoeffs");
  if (basis.kind != RingSwitchBasisKind::kPolynomial) {
    LogicError("DecomposeGRElementToBaseCoeffs: unsupported basis kind");
  }
  return DecomposeGRElementToBaseCoeffsPolynomialBasisUnchecked(
      basis.dimension, element);
}

RingSwitchComponentTensor BuildRingSwitchComponentTensor(
    const RingSwitchPCSParams &params, const std::vector<ZZ_pE> &r_suffix) {
  ValidateRingSwitchPCSParamsOrThrow(params);
  if (static_cast<long>(r_suffix.size()) != params.ell_prime) {
    LogicError(
        "BuildRingSwitchComponentTensor: r_suffix dimension must equal ell_prime");
  }

  const long basis_dimension = params.alpha_basis.dimension;
  const long num_w = Pow2LongOrThrow(
      params.ell_prime,
      "BuildRingSwitchComponentTensor: ell_prime is too large for long");
  const long total_coeffs = CheckedMultiplyLong(
      basis_dimension, num_w,
      "BuildRingSwitchComponentTensor: coefficient table is too large for long");

  RingSwitchComponentTensor tensor;
  tensor.basis_dimension = basis_dimension;
  tensor.ell_prime = params.ell_prime;
  tensor.a_by_u_then_w.SetLength(total_coeffs);
  tensor.r_table.SetLength(total_coeffs);

  for (long w = 0; w < num_w; ++w) {
    const std::vector<ZZ_pE> bool_point =
        BooleanPointFromIndex(w, params.ell_prime);
    const ZZ_pE eq_at_w = EqPolynomial(r_suffix, bool_point);
    const vec_ZZ_pE coeffs =
        DecomposeGRElementToBaseCoeffsPolynomialBasisUnchecked(basis_dimension,
                                                               eq_at_w);
    for (long u = 0; u < basis_dimension; ++u) {
      const long row_major_index = u * num_w + w;
      const long flat_index = u + w * basis_dimension;
      tensor.a_by_u_then_w[row_major_index] = coeffs[u];
      tensor.r_table[flat_index] = coeffs[u];
    }
  }
  tensor.r_monomial_coeffs =
      BooleanHypercubeTableToMonomialCoeffs(tensor.r_table);

  return tensor;
}

}  // namespace basefold
