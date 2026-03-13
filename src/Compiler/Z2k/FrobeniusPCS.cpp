#include "Compiler/Z2k/FrobeniusPCS.hpp"

#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

#include <limits>
#include <string>
#include <vector>

using NTL::LogicError;
using NTL::SetCoeff;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pX;
using std::size_t;
using std::string;
using std::vector;

namespace basefold {
namespace {

struct PackedCommitInputs {
  NTL::vec_ZZ_pE t_packed_table;
  NTL::vec_ZZ_pE t_packed_monomial_coeffs;
};

long Pow2LongOrThrow(long exponent, const char *what) {
  if (exponent < 0) {
    LogicError(what);
  }
  long out = 1;
  for (long i = 0; i < exponent; ++i) {
    if (out > std::numeric_limits<long>::max() / 2) {
      LogicError(what);
    }
    out *= 2;
  }
  return out;
}

long Log2ExactPowerOfTwoZZOrThrow(const ZZ &value, const char *what) {
  if (value <= 0) {
    LogicError(what);
  }
  ZZ tmp = value;
  long out = 0;
  while (tmp > 1) {
    if ((tmp % 2) != 0) {
      LogicError(what);
    }
    tmp /= 2;
    ++out;
  }
  return out;
}

ZZ_pE BaseRingConstant(const ZZ_p &value) {
  ZZ_pX poly;
  if (value != 0) {
    SetCoeff(poly, 0, value);
  }
  ZZ_pE out;
  NTL::conv(out, poly);
  return out;
}

void ValidateBaseRingConstantOrThrow(const ZZ_pE &value, const char *label,
                                     long index, const char *func_name) {
  if (!::IsBaseRingConstant(value)) {
    LogicError((string(func_name) + ": " + label + "[" +
                std::to_string(index) + "] must be a base-ring constant")
                   .c_str());
  }
}

void ValidateBaseRingVectorOrThrow(const NTL::vec_ZZ_pE &values,
                                   const char *label,
                                   const char *func_name) {
  for (long i = 0; i < values.length(); ++i) {
    ValidateBaseRingConstantOrThrow(values[i], label, i, func_name);
  }
}

NTL::vec_ZZ_pE BooleanHypercubeTableToMonomialCoeffsInternal(
    const NTL::vec_ZZ_pE &table_values) {
  const long n = table_values.length();
  if (n <= 0 || (n & (n - 1)) != 0) {
    LogicError("BooleanHypercubeTableToMonomialCoeffsInternal: table length must be a power of two");
  }
  long dimension = 0;
  for (long tmp = n; tmp > 1; tmp >>= 1) {
    ++dimension;
  }

  NTL::vec_ZZ_pE coeffs = table_values;
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

PackedCommitInputs BuildPackedCommitInputs(const FrobeniusPCSParams &params,
                                           const NTL::vec_ZZ_pE &t_table) {
  PackedCommitInputs out;
  out.t_packed_table = PackZ2kTableToFrobeniusGREvals(params, t_table);
  out.t_packed_monomial_coeffs =
      BooleanHypercubeTableToMonomialCoeffsInternal(out.t_packed_table);
  return out;
}

}  // namespace

void ValidateCurrentZ2kFrobeniusContextOrThrow(
    const ZZ &base_modulus, const ZZ_pX &extension_modulus, long kappa) {
  if (kappa < 1) {
    LogicError(
        "ValidateCurrentZ2kFrobeniusContextOrThrow: kappa must be >= 1");
  }
  if (base_modulus <= 1) {
    LogicError(
        "ValidateCurrentZ2kFrobeniusContextOrThrow: base_modulus must be > 1");
  }
  (void)Log2ExactPowerOfTwoZZOrThrow(
      base_modulus,
      "ValidateCurrentZ2kFrobeniusContextOrThrow: base_modulus must be a power of two");
  const long expected_degree = Pow2LongOrThrow(
      kappa,
      "ValidateCurrentZ2kFrobeniusContextOrThrow: kappa is too large for long");
  if (NTL::ZZ_p::modulus() != base_modulus) {
    LogicError(
        "ValidateCurrentZ2kFrobeniusContextOrThrow: current ZZ_p modulus does not match base_modulus");
  }
  if (!NTL::ZZ_pE::initialized()) {
    LogicError(
        "ValidateCurrentZ2kFrobeniusContextOrThrow: ZZ_pE context must be initialized");
  }
  if (NTL::ZZ_pE::degree() != expected_degree) {
    LogicError(
        "ValidateCurrentZ2kFrobeniusContextOrThrow: current ZZ_pE degree must equal 2^kappa");
  }
  if (NTL::deg(extension_modulus) != expected_degree) {
    LogicError(
        "ValidateCurrentZ2kFrobeniusContextOrThrow: extension_modulus degree must equal 2^kappa");
  }
  if (NTL::ZZ_pE::modulus() != extension_modulus) {
    LogicError(
        "ValidateCurrentZ2kFrobeniusContextOrThrow: current ZZ_pE modulus does not match extension_modulus");
  }
}

void ValidateFrobeniusPCSParamsOrThrow(const FrobeniusPCSParams &params) {
  if (params.kappa < 1) {
    LogicError("ValidateFrobeniusPCSParamsOrThrow: kappa must be >= 1");
  }
  if (params.ell < params.kappa) {
    LogicError("ValidateFrobeniusPCSParamsOrThrow: ell must be >= kappa");
  }

  ValidateCurrentZ2kFrobeniusContextOrThrow(params.base_modulus,
                                            params.extension_modulus,
                                            params.kappa);
  Z2kPCSBackendValidateParamsOrThrow(params.backend);

  const long base_log2 = Log2ExactPowerOfTwoZZOrThrow(
      params.base_modulus,
      "ValidateFrobeniusPCSParamsOrThrow: base_modulus must be a power of two");
  const long expected_degree = Pow2LongOrThrow(
      params.kappa,
      "ValidateFrobeniusPCSParamsOrThrow: kappa is too large for long");
  const long ell_prime = params.ell - params.kappa;
  if (params.ell_prime != ell_prime) {
    LogicError("ValidateFrobeniusPCSParamsOrThrow: ell_prime must equal ell-kappa");
  }

  if (params.basis_data.params.p != ZZ(2)) {
    LogicError("ValidateFrobeniusPCSParamsOrThrow: basis data must target p=2");
  }
  if (params.basis_data.params.k != base_log2) {
    LogicError(
        "ValidateFrobeniusPCSParamsOrThrow: basis data k does not match base_modulus");
  }
  if (params.basis_data.params.r != expected_degree) {
    LogicError(
        "ValidateFrobeniusPCSParamsOrThrow: basis data r does not match 2^kappa");
  }
  if (params.basis_data.modulus != params.base_modulus) {
    LogicError(
        "ValidateFrobeniusPCSParamsOrThrow: basis data modulus does not match base_modulus");
  }
  if (params.basis_data.extension_modulus != params.extension_modulus) {
    LogicError(
        "ValidateFrobeniusPCSParamsOrThrow: basis data extension modulus does not match params");
  }
  ::ValidateNormalBasisOrThrow(params.basis_data.normal_basis);

  const long basis_dimension =
      static_cast<long>(params.basis_data.normal_basis.beta.size());
  if (basis_dimension != expected_degree) {
    LogicError(
        "ValidateFrobeniusPCSParamsOrThrow: beta basis dimension must equal 2^kappa");
  }
  if (static_cast<long>(params.basis_data.normal_basis.alpha.size()) !=
      expected_degree) {
    LogicError(
        "ValidateFrobeniusPCSParamsOrThrow: alpha basis dimension must equal 2^kappa");
  }

  const long expected_backend_message_length = Pow2LongOrThrow(
      ell_prime,
      "ValidateFrobeniusPCSParamsOrThrow: ell-kappa is too large for long");
  if (Z2kPCSBackendMessageLength(params.backend) !=
      expected_backend_message_length) {
    LogicError(
        "ValidateFrobeniusPCSParamsOrThrow: backend message length must equal 2^(ell-kappa)");
  }
  if (Z2kPCSBackendPointDimension(params.backend) != ell_prime) {
    LogicError(
        "ValidateFrobeniusPCSParamsOrThrow: backend point dimension must equal ell-kappa");
  }
}

FrobeniusPCSParams FrobeniusPCSSetup(const FrobeniusPCSSetupInput &input) {
  if (input.kappa < 1) {
    LogicError("FrobeniusPCSSetup: kappa must be >= 1");
  }
  if (input.ell < input.kappa) {
    LogicError("FrobeniusPCSSetup: ell must be >= kappa");
  }

  const long base_log2 = Log2ExactPowerOfTwoZZOrThrow(
      input.base_modulus,
      "FrobeniusPCSSetup: base_modulus must be a power of two");
  const long expected_degree = Pow2LongOrThrow(
      input.kappa, "FrobeniusPCSSetup: kappa is too large for long");
  ValidateCurrentZ2kFrobeniusContextOrThrow(input.base_modulus,
                                            input.extension_modulus,
                                            input.kappa);
  Z2kPCSBackendValidateParamsOrThrow(input.backend);

  FrobeniusBasisParams basis_input;
  basis_input.p = ZZ(2);
  basis_input.k = base_log2;
  basis_input.r = expected_degree;
  basis_input.teichmuller_generator_max_trials =
      input.teichmuller_generator_max_trials;
  basis_input.affine_search_limit = input.affine_search_limit;

  FrobeniusPCSParams params;
  params.ell = input.ell;
  params.kappa = input.kappa;
  params.ell_prime = input.ell - input.kappa;
  params.base_modulus = input.base_modulus;
  params.extension_modulus = input.extension_modulus;
  params.basis_data =
      ::BuildFrobeniusBasisOrThrow(basis_input, input.extension_modulus);
  params.backend = input.backend;
  ValidateFrobeniusPCSParamsOrThrow(params);
  return params;
}

NTL::vec_ZZ_pE PackZ2kTableToFrobeniusGREvals(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table) {
  ValidateFrobeniusPCSParamsOrThrow(params);
  const long expected_length = Pow2LongOrThrow(
      params.ell, "PackZ2kTableToFrobeniusGREvals: ell is too large for long");
  if (t_table.length() != expected_length) {
    LogicError(
        "PackZ2kTableToFrobeniusGREvals: t_table length must equal 2^ell");
  }
  ValidateBaseRingVectorOrThrow(t_table, "t_table",
                                "PackZ2kTableToFrobeniusGREvals");

  const vector<ZZ_pE> &beta = params.basis_data.normal_basis.beta;
  const long basis_dimension = static_cast<long>(beta.size());
  const long packed_length = Pow2LongOrThrow(
      params.ell_prime,
      "PackZ2kTableToFrobeniusGREvals: ell_prime is too large for long");

  NTL::vec_ZZ_pE packed;
  packed.SetLength(packed_length);
  for (long w = 0; w < packed_length; ++w) {
    ZZ_pE acc;
    NTL::clear(acc);
    for (long v = 0; v < basis_dimension; ++v) {
      acc += beta[static_cast<size_t>(v)] *
             t_table[v + w * basis_dimension];
    }
    packed[w] = acc;
  }
  return packed;
}

NTL::vec_ZZ_pE DecomposeGRElementToBaseCoeffsFrobeniusBasis(
    const FrobeniusPCSParams &params, const ZZ_pE &element) {
  ValidateFrobeniusPCSParamsOrThrow(params);
  const vector<ZZ_p> coords =
      ::RecoverNormalBasisCoords(params.basis_data.normal_basis, element);
  NTL::vec_ZZ_pE out;
  out.SetLength(static_cast<long>(coords.size()));
  for (long i = 0; i < static_cast<long>(coords.size()); ++i) {
    out[i] = BaseRingConstant(coords[static_cast<size_t>(i)]);
  }
  return out;
}

MerkleRoot FrobeniusPCSCommit(const FrobeniusPCSParams &params,
                              const NTL::vec_ZZ_pE &t_table) {
  return FrobeniusPCSBuildCommitArtifacts(params, t_table).commitment;
}

FrobeniusPCSOuterCommitArtifacts FrobeniusPCSBuildOuterCommitArtifacts(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table) {
  ValidateFrobeniusPCSParamsOrThrow(params);
  const PackedCommitInputs packed = BuildPackedCommitInputs(params, t_table);

  FrobeniusPCSOuterCommitArtifacts out;
  out.t_packed_table = packed.t_packed_table;
  out.t_packed_monomial_coeffs = packed.t_packed_monomial_coeffs;
  return out;
}

FrobeniusPCSCommitArtifacts FrobeniusPCSBuildCommitArtifacts(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table) {
  const FrobeniusPCSOuterCommitArtifacts outer =
      FrobeniusPCSBuildOuterCommitArtifacts(params, t_table);

  FrobeniusPCSCommitArtifacts out;
  out.t_packed_table = outer.t_packed_table;
  out.t_packed_monomial_coeffs = outer.t_packed_monomial_coeffs;
  out.backend_commit_artifacts = Z2kPCSBackendBuildCommitArtifacts(
      params.backend, out.t_packed_monomial_coeffs);
  out.commitment = out.backend_commit_artifacts.commitment;
  return out;
}

}  // namespace basefold
