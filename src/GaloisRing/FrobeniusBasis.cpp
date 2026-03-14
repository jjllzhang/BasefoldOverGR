#include "GaloisRing/FrobeniusBasis.hpp"

#include <NTL/ZZ_pXFactoring.h>

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

#include "GaloisRing/PrimitiveElement.hpp"

using NTL::coeff;
using NTL::conv;
using NTL::deg;
using NTL::IsOne;
using NTL::LogicError;
using NTL::power;
using NTL::set;
using NTL::SetCoeff;
using NTL::to_ZZ;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pX;
using std::size_t;
using std::string;
using std::vector;

namespace {

struct BasisWithDualData {
  vector<ZZ_pE> basis;
  vector<ZZ_pE> dual;
};

struct PowerBasisFrobeniusContext {
  ZZ p;
  long p_long = 0;
  long r = 0;
  BasisWithDualData power_basis;
  vector<ZZ_pE> tau_images;
};

long PositiveZZToLongChecked(const ZZ &z, const char *label,
                             const char *func_name) {
  if (z <= 0) {
    LogicError((string(func_name) + ": " + label + " must be > 0").c_str());
  }

  const long max_long_bits = static_cast<long>(8 * sizeof(long) - 1);
  if (NTL::NumBits(z) > max_long_bits) {
    LogicError(
        (string(func_name) + ": " + label + " too large for long").c_str());
  }
  return NTL::to_long(z);
}

void ValidateExtensionPolynomialOrThrow(const ZZ_pX &extension_modulus, long r,
                                        const char *func_name) {
  if (deg(extension_modulus) != r) {
    LogicError(
        (string(func_name) + ": deg(extension_modulus) must equal r").c_str());
  }
  if (!IsOne(coeff(extension_modulus, r))) {
    LogicError(
        (string(func_name) + ": extension_modulus must be monic").c_str());
  }
}

void ValidatePositiveParamsOrThrow(const FrobeniusBasisParams &params,
                                   const char *func_name) {
  if (params.p <= 1) {
    LogicError((string(func_name) + ": p must be > 1").c_str());
  }
  if (params.k <= 0) {
    LogicError((string(func_name) + ": k must be > 0").c_str());
  }
  if (params.r <= 0) {
    LogicError((string(func_name) + ": r must be > 0").c_str());
  }
  if (params.teichmuller_generator_max_trials <= 0) {
    LogicError(
        (string(func_name) + ": teichmuller_generator_max_trials must be > 0")
            .c_str());
  }
  if (params.affine_search_limit < 0) {
    LogicError(
        (string(func_name) + ": affine_search_limit must be >= 0").c_str());
  }
}

ZZ_pE BaseRingConstant(const ZZ_p &value) {
  ZZ_pX poly;
  if (value != 0) {
    SetCoeff(poly, 0, value);
  }
  ZZ_pE out;
  conv(out, poly);
  return out;
}

ZZ_pE BaseRingConstant(long value) {
  return BaseRingConstant(NTL::to_ZZ_p(value));
}

unsigned long long PositiveZZToU64Checked(const ZZ &z, const char *label,
                                          const char *func_name) {
  if (z <= 0) {
    LogicError((string(func_name) + ": " + label + " must be > 0").c_str());
  }
  if (NTL::NumBits(z) > 64) {
    LogicError(
        (string(func_name) + ": " + label + " too large for u64").c_str());
  }
  std::ostringstream os;
  os << z;
  return std::stoull(os.str());
}

vector<ZZ> UniquePrimeFactorsU64(unsigned long long n) {
  vector<ZZ> factors;
  if (n <= 1) {
    return factors;
  }

  for (unsigned long long d = 2; d <= n / d; ++d) {
    if (n % d != 0) {
      continue;
    }
    factors.push_back(to_ZZ(static_cast<long>(d)));
    while (n % d == 0) {
      n /= d;
    }
  }
  if (n > 1) {
    std::ostringstream os;
    os << n;
    factors.push_back(to_ZZ(os.str().c_str()));
  }
  return factors;
}

bool HasExactOrder(const ZZ_pE &element, const ZZ &order,
                   const vector<ZZ> &prime_factors) {
  if (element == 0 || order <= 0) {
    return false;
  }
  ZZ_pE one;
  set(one);
  if (power(element, order) != one) {
    return false;
  }
  for (const ZZ &q : prime_factors) {
    if (power(element, order / q) == one) {
      return false;
    }
  }
  return true;
}

vector<ZZ_pE> BuildTeichmullerPowerBasis(const ZZ_pE &teichmuller_generator,
                                         long r) {
  vector<ZZ_pE> basis(static_cast<size_t>(r));
  basis[0] = BaseRingConstant(1);
  for (long i = 1; i < r; ++i) {
    basis[static_cast<size_t>(i)] =
        basis[static_cast<size_t>(i - 1)] * teichmuller_generator;
  }
  return basis;
}

PowerBasisFrobeniusContext
BuildPowerBasisFrobeniusContext(const FrobeniusBasisParams &params,
                                const ZZ_pE &teichmuller_generator,
                                const char *func_name) {
  PowerBasisFrobeniusContext context;
  context.p = params.p;
  context.p_long = PositiveZZToLongChecked(params.p, "p", func_name);
  context.r = params.r;
  context.power_basis.basis =
      BuildTeichmullerPowerBasis(teichmuller_generator, params.r);
  if (!IsBasisOverBaseRing(context.power_basis.basis)) {
    LogicError((string(func_name) +
                ": Teichmuller power basis does not span the extension")
                   .c_str());
  }
  context.power_basis.dual = BuildDualBasisOrThrow(context.power_basis.basis);

  context.tau_images.resize(static_cast<size_t>(params.r));
  for (long i = 0; i < params.r; ++i) {
    context.tau_images[static_cast<size_t>(i)] =
        power(teichmuller_generator, context.p_long * i);
  }
  return context;
}

ZZ_pE ApplyPowerBasisTau(const PowerBasisFrobeniusContext &context,
                         const ZZ_pE &element) {
  GaloisRingBasisData power_basis_data;
  power_basis_data.basis = context.power_basis.basis;
  power_basis_data.dual_basis = context.power_basis.dual;
  const vector<ZZ_p> coords = RecoverBasisCoordsOrThrow(
      power_basis_data, element, "ApplyPowerBasisTau");
  return ComposeFromBasisCoordsOrThrow(context.tau_images, coords,
                                       "ApplyPowerBasisTau");
}

bool HasFullOrbit(const vector<ZZ_pE> &orbit) {
  for (long i = 0; i < static_cast<long>(orbit.size()); ++i) {
    for (long j = i + 1; j < static_cast<long>(orbit.size()); ++j) {
      if (orbit[static_cast<size_t>(i)] == orbit[static_cast<size_t>(j)]) {
        return false;
      }
    }
  }
  return true;
}

bool TryPromoteNormalBasisCandidate(const PowerBasisFrobeniusContext &context,
                                    const ZZ_pE &candidate,
                                    NormalBasisData *normal_basis_out) {
  vector<ZZ_pE> orbit(static_cast<size_t>(context.r), ZZ_pE(0));
  orbit[0] = candidate;
  for (long i = 1; i < context.r; ++i) {
    orbit[static_cast<size_t>(i)] =
        ApplyPowerBasisTau(context, orbit[static_cast<size_t>(i - 1)]);
  }
  if (ApplyPowerBasisTau(context, orbit.back()) != orbit.front()) {
    return false;
  }
  if (!HasFullOrbit(orbit)) {
    return false;
  }
  if (!IsBasisOverBaseRing(orbit)) {
    return false;
  }

  normal_basis_out->beta = orbit;
  normal_basis_out->alpha = BuildDualBasisOrThrow(orbit);
  return true;
}

void ValidateTeichmullerGeneratorOrThrow(const FrobeniusBasisParams &params,
                                         const ZZ_pE &teichmuller_generator,
                                         const char *func_name) {
  const ZZ subgroup_order_zz = power(params.p, params.r) - ZZ(1);
  const unsigned long long subgroup_order_u64 =
      PositiveZZToU64Checked(subgroup_order_zz, "p^r-1", func_name);
  const vector<ZZ> subgroup_order_factors =
      UniquePrimeFactorsU64(subgroup_order_u64);
  if (!HasExactOrder(teichmuller_generator, subgroup_order_zz,
                     subgroup_order_factors)) {
    LogicError((string(func_name) +
                ": teichmuller_generator must have multiplicative order p^r-1")
                   .c_str());
  }
}

void ValidateNormalBasisAgainstTeichmullerOrThrow(
    const FrobeniusBasisParams &params, const NormalBasisData &normal_basis,
    const ZZ_pE &teichmuller_generator, const char *func_name) {
  ValidateTeichmullerGeneratorOrThrow(params, teichmuller_generator, func_name);
  const PowerBasisFrobeniusContext context =
      BuildPowerBasisFrobeniusContext(params, teichmuller_generator, func_name);
  const long n = static_cast<long>(normal_basis.beta.size());
  for (long i = 0; i < n; ++i) {
    const long next = (i + 1) % n;
    if (ApplyPowerBasisTau(context,
                           normal_basis.beta[static_cast<size_t>(i)]) !=
        normal_basis.beta[static_cast<size_t>(next)]) {
      LogicError((string(func_name) +
                  ": provided beta basis is not ordered as a Frobenius orbit")
                     .c_str());
    }
    if (ApplyPowerBasisTau(context,
                           normal_basis.alpha[static_cast<size_t>(i)]) !=
        normal_basis.alpha[static_cast<size_t>(next)]) {
      LogicError((string(func_name) +
                  ": provided alpha basis is not ordered as a Frobenius orbit")
                     .c_str());
    }
  }
}

} // namespace

void ValidateFrobeniusBasisContextsOrThrow(const FrobeniusBasisParams &params,
                                           const ZZ_pX &extension_modulus) {
  const char *const func_name = "ValidateFrobeniusBasisContextsOrThrow";
  ValidatePositiveParamsOrThrow(params, func_name);
  ValidateExtensionPolynomialOrThrow(extension_modulus, params.r, func_name);

  const ZZ expected_modulus = power(params.p, params.k);
  if (ZZ_p::modulus() != expected_modulus) {
    LogicError((string(func_name) + ": active ZZ_p modulus does not equal p^k")
                   .c_str());
  }
  if (!ZZ_pE::initialized()) {
    LogicError((string(func_name) + ": ZZ_pE must be initialized").c_str());
  }
  if (ZZ_pE::degree() != params.r) {
    LogicError(
        (string(func_name) + ": active ZZ_pE degree does not equal r").c_str());
  }
  if (ZZ_pE::modulus() != extension_modulus) {
    LogicError((string(func_name) +
                ": active ZZ_pE modulus does not match extension_modulus")
                   .c_str());
  }
}

FrobeniusBasisData
BuildFrobeniusBasisOrThrow(const FrobeniusBasisParams &params,
                           const ZZ_pX &extension_modulus) {
  ValidateFrobeniusBasisContextsOrThrow(params, extension_modulus);

  FrobeniusBasisData out;
  out.params = params;
  out.modulus = ZZ_p::modulus();
  out.extension_modulus = extension_modulus;
  out.teichmuller_generator =
      FindTeichmullerGenerator(params.p, params.k, params.r, extension_modulus,
                               params.teichmuller_generator_max_trials);
  out.normal_basis = FindNormalBasisOrThrow(params, out.teichmuller_generator);
  ValidateNormalBasisOrThrow(out.normal_basis);
  return out;
}

FrobeniusBasisData BuildFrobeniusBasisFromProvidedNormalBasisOrThrow(
    const FrobeniusBasisParams &params, const ZZ_pX &extension_modulus,
    const NormalBasisData &normal_basis, bool has_teichmuller_generator,
    const ZZ_pE &teichmuller_generator) {
  const char *const func_name =
      "BuildFrobeniusBasisFromProvidedNormalBasisOrThrow";
  ValidateFrobeniusBasisContextsOrThrow(params, extension_modulus);

  FrobeniusBasisData out;
  out.params = params;
  out.modulus = ZZ_p::modulus();
  out.extension_modulus = extension_modulus;
  out.normal_basis = normal_basis;
  ValidateNormalBasisOrThrow(out.normal_basis);

  if (has_teichmuller_generator) {
    out.teichmuller_generator = teichmuller_generator;
  } else {
    out.teichmuller_generator = FindTeichmullerGenerator(
        params.p, params.k, params.r, extension_modulus,
        params.teichmuller_generator_max_trials);
  }
  ValidateNormalBasisAgainstTeichmullerOrThrow(
      params, out.normal_basis, out.teichmuller_generator, func_name);
  return out;
}

NormalBasisData FindNormalBasisOrThrow(const FrobeniusBasisParams &params,
                                       const ZZ_pE &teichmuller_generator) {
  const char *const func_name = "FindNormalBasisOrThrow";
  ValidatePositiveParamsOrThrow(params, func_name);
  if (ZZ_pE::degree() != params.r) {
    LogicError(
        (string(func_name) + ": active ZZ_pE degree does not equal r").c_str());
  }

  const PowerBasisFrobeniusContext context =
      BuildPowerBasisFrobeniusContext(params, teichmuller_generator, func_name);
  const ZZ subgroup_order_zz = power(params.p, params.r) - ZZ(1);
  const long subgroup_order =
      PositiveZZToLongChecked(subgroup_order_zz, "p^r-1", func_name);

  NormalBasisData normal_basis;

  ZZ_pE teichmuller_power = BaseRingConstant(1);
  for (long exponent = 0; exponent < subgroup_order; ++exponent) {
    if (TryPromoteNormalBasisCandidate(context, teichmuller_power,
                                       &normal_basis)) {
      return normal_basis;
    }
    teichmuller_power *= teichmuller_generator;
  }

  for (long a0 = 0; a0 <= params.affine_search_limit; ++a0) {
    for (long a1 = 1; a1 <= params.affine_search_limit; ++a1) {
      if (a0 == 0 && a1 == 1) {
        continue;
      }
      teichmuller_power = BaseRingConstant(1);
      for (long exponent = 0; exponent < subgroup_order; ++exponent) {
        const ZZ_pE candidate =
            BaseRingConstant(a0) + BaseRingConstant(a1) * teichmuller_power;
        if (TryPromoteNormalBasisCandidate(context, candidate, &normal_basis)) {
          return normal_basis;
        }
        teichmuller_power *= teichmuller_generator;
      }
    }
  }

  LogicError((string(func_name) +
              ": failed to find a normal basis from deterministic search")
                 .c_str());
  return normal_basis;
}

void ValidateNormalBasisOrThrow(const NormalBasisData &normal_basis) {
  const char *const func_name = "ValidateNormalBasisOrThrow";
  GaloisRingBasisData basis_data;
  basis_data.basis = normal_basis.beta;
  basis_data.dual_basis = normal_basis.alpha;
  ValidateBasisDataOrThrow(basis_data, "normal_basis", func_name);
}

std::vector<ZZ_p> RecoverNormalBasisCoords(const NormalBasisData &normal_basis,
                                           const ZZ_pE &element) {
  GaloisRingBasisData basis_data;
  basis_data.basis = normal_basis.beta;
  basis_data.dual_basis = normal_basis.alpha;
  return RecoverBasisCoordsOrThrow(basis_data, element,
                                   "RecoverNormalBasisCoords");
}

ZZ_pE ComposeFromNormalBasisCoords(const NormalBasisData &normal_basis,
                                   const vector<ZZ_p> &coords) {
  return ComposeFromBasisCoordsOrThrow(normal_basis.beta, coords,
                                       "ComposeFromNormalBasisCoords");
}

ZZ_pE ApplyFrobeniusTau(const NormalBasisData &normal_basis,
                        const ZZ_pE &element) {
  const vector<ZZ_p> coords = RecoverNormalBasisCoords(normal_basis, element);
  vector<ZZ_p> rotated(coords.size(), ZZ_p(0));
  if (!coords.empty()) {
    rotated[0] = coords.back();
    for (long i = 1; i < static_cast<long>(coords.size()); ++i) {
      rotated[static_cast<size_t>(i)] = coords[static_cast<size_t>(i - 1)];
    }
  }
  return ComposeFromNormalBasisCoords(normal_basis, rotated);
}

ZZ_pE ApplyFrobeniusSigma(const NormalBasisData &normal_basis,
                          const ZZ_pE &element) {
  const vector<ZZ_p> coords = RecoverNormalBasisCoords(normal_basis, element);
  vector<ZZ_p> rotated(coords.size(), ZZ_p(0));
  if (!coords.empty()) {
    for (long i = 0; i + 1 < static_cast<long>(coords.size()); ++i) {
      rotated[static_cast<size_t>(i)] = coords[static_cast<size_t>(i + 1)];
    }
    rotated.back() = coords.front();
  }
  return ComposeFromNormalBasisCoords(normal_basis, rotated);
}
