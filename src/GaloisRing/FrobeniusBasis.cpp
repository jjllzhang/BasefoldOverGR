#include "GaloisRing/FrobeniusBasis.hpp"

#include <NTL/ZZ_pXFactoring.h>

#include <algorithm>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "GaloisRing/PrimitiveElement.hpp"

using NTL::GCD;
using NTL::IsOne;
using NTL::LogicError;
using NTL::SetCoeff;
using NTL::Vec;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pX;
using NTL::coeff;
using NTL::conv;
using NTL::deg;
using NTL::power;
using NTL::rep;
using NTL::set;
using NTL::to_ZZ;
using NTL::trace;
using std::size_t;
using std::string;
using std::swap;
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
    LogicError((string(func_name) + ": " + label + " too large for long")
                   .c_str());
  }
  return NTL::to_long(z);
}

void ValidateExtensionPolynomialOrThrow(const ZZ_pX &extension_modulus, long r,
                                        const char *func_name) {
  if (deg(extension_modulus) != r) {
    LogicError((string(func_name) + ": deg(extension_modulus) must equal r")
                   .c_str());
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
    LogicError((string(func_name) +
                ": teichmuller_generator_max_trials must be > 0")
                   .c_str());
  }
  if (params.affine_search_limit < 0) {
    LogicError((string(func_name) + ": affine_search_limit must be >= 0")
                   .c_str());
  }
}

void ValidateBasisShapeOrThrow(const vector<ZZ_pE> &basis, const char *label,
                               const char *func_name) {
  if (basis.empty()) {
    LogicError((string(func_name) + ": " + label + " must be non-empty")
                   .c_str());
  }
  if (static_cast<long>(basis.size()) != ZZ_pE::degree()) {
    LogicError((string(func_name) + ": " + label +
                ".size() must equal ZZ_pE::degree()")
                   .c_str());
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

bool IsUnit(const ZZ_p &value) {
  if (value == 0) {
    return false;
  }
  return GCD(rep(value), ZZ_p::modulus()) == 1;
}

vector<ZZ_pE> BuildPolynomialBasis(long r) {
  vector<ZZ_pE> basis(static_cast<size_t>(r));
  for (long i = 0; i < r; ++i) {
    ZZ_pX poly;
    SetCoeff(poly, i, 1);
    conv(basis[static_cast<size_t>(i)], poly);
  }
  return basis;
}

vector<vector<ZZ_p>> BuildCoordinateMatrix(const vector<ZZ_pE> &basis) {
  const long n = static_cast<long>(basis.size());
  vector<vector<ZZ_p>> matrix(static_cast<size_t>(n),
                              vector<ZZ_p>(static_cast<size_t>(n), ZZ_p(0)));
  for (long col = 0; col < n; ++col) {
    const ZZ_pX poly = rep(basis[static_cast<size_t>(col)]);
    for (long row = 0; row < n; ++row) {
      matrix[static_cast<size_t>(row)][static_cast<size_t>(col)] =
          coeff(poly, row);
    }
  }
  return matrix;
}

void SwapMatrixColumns(vector<vector<ZZ_p>> *matrix, long lhs, long rhs) {
  if (lhs == rhs) {
    return;
  }
  for (auto &row : *matrix) {
    swap(row[static_cast<size_t>(lhs)], row[static_cast<size_t>(rhs)]);
  }
}

bool TrySolveLinearSystemWithUnitPivots(const vector<vector<ZZ_p>> &matrix_in,
                                        const vector<ZZ_p> &rhs_in,
                                        vector<ZZ_p> *solution_out) {
  const long n = static_cast<long>(matrix_in.size());
  if (n == 0 || static_cast<long>(rhs_in.size()) != n) {
    return false;
  }
  for (const auto &row : matrix_in) {
    if (static_cast<long>(row.size()) != n) {
      return false;
    }
  }

  vector<vector<ZZ_p>> matrix = matrix_in;
  vector<ZZ_p> rhs = rhs_in;
  vector<long> column_perm(static_cast<size_t>(n));
  std::iota(column_perm.begin(), column_perm.end(), 0);

  for (long pivot_index = 0; pivot_index < n; ++pivot_index) {
    long pivot_row = -1;
    long pivot_col = -1;
    for (long row = pivot_index; row < n && pivot_row < 0; ++row) {
      for (long col = pivot_index; col < n; ++col) {
        if (IsUnit(matrix[static_cast<size_t>(row)]
                          [static_cast<size_t>(col)])) {
          pivot_row = row;
          pivot_col = col;
          break;
        }
      }
    }
    if (pivot_row < 0) {
      return false;
    }

    if (pivot_row != pivot_index) {
      swap(matrix[static_cast<size_t>(pivot_row)],
           matrix[static_cast<size_t>(pivot_index)]);
      swap(rhs[static_cast<size_t>(pivot_row)],
           rhs[static_cast<size_t>(pivot_index)]);
    }
    if (pivot_col != pivot_index) {
      SwapMatrixColumns(&matrix, pivot_col, pivot_index);
      swap(column_perm[static_cast<size_t>(pivot_col)],
           column_perm[static_cast<size_t>(pivot_index)]);
    }

    const ZZ_p pivot =
        matrix[static_cast<size_t>(pivot_index)]
              [static_cast<size_t>(pivot_index)];
    if (!IsUnit(pivot)) {
      return false;
    }
    const ZZ_p pivot_inverse = inv(pivot);
    for (long col = pivot_index; col < n; ++col) {
      matrix[static_cast<size_t>(pivot_index)][static_cast<size_t>(col)] *=
          pivot_inverse;
    }
    rhs[static_cast<size_t>(pivot_index)] *= pivot_inverse;

    for (long row = 0; row < n; ++row) {
      if (row == pivot_index) {
        continue;
      }
      const ZZ_p factor =
          matrix[static_cast<size_t>(row)][static_cast<size_t>(pivot_index)];
      if (factor == 0) {
        continue;
      }
      for (long col = pivot_index; col < n; ++col) {
        matrix[static_cast<size_t>(row)][static_cast<size_t>(col)] -=
            factor *
            matrix[static_cast<size_t>(pivot_index)][static_cast<size_t>(col)];
      }
      rhs[static_cast<size_t>(row)] -=
          factor * rhs[static_cast<size_t>(pivot_index)];
    }
  }

  vector<ZZ_p> solution(static_cast<size_t>(n), ZZ_p(0));
  for (long col = 0; col < n; ++col) {
    const long original_col = column_perm[static_cast<size_t>(col)];
    solution[static_cast<size_t>(original_col)] = rhs[static_cast<size_t>(col)];
  }
  *solution_out = solution;
  return true;
}

bool IsInvertibleOverBaseRing(const vector<vector<ZZ_p>> &matrix) {
  const long n = static_cast<long>(matrix.size());
  vector<ZZ_p> rhs(static_cast<size_t>(n), ZZ_p(0));
  vector<ZZ_p> solution;
  return TrySolveLinearSystemWithUnitPivots(matrix, rhs, &solution);
}

vector<ZZ_p> SolveLinearSystemWithUnitPivotsOrThrow(
    const vector<vector<ZZ_p>> &matrix, const vector<ZZ_p> &rhs,
    const char *func_name) {
  vector<ZZ_p> solution;
  if (!TrySolveLinearSystemWithUnitPivots(matrix, rhs, &solution)) {
    LogicError((string(func_name) +
                ": failed to solve linear system over the base ring")
                   .c_str());
  }
  return solution;
}

vector<ZZ_pE> BuildDualBasisByLinearSolveOrThrow(const vector<ZZ_pE> &basis,
                                                 const char *func_name) {
  ValidateBasisShapeOrThrow(basis, "basis", func_name);
  const long n = static_cast<long>(basis.size());
  const vector<ZZ_pE> polynomial_basis = BuildPolynomialBasis(n);

  vector<vector<ZZ_p>> trace_matrix(static_cast<size_t>(n),
                                    vector<ZZ_p>(static_cast<size_t>(n),
                                                 ZZ_p(0)));
  for (long row = 0; row < n; ++row) {
    for (long col = 0; col < n; ++col) {
      trace_matrix[static_cast<size_t>(row)][static_cast<size_t>(col)] =
          TraceToBaseRing(basis[static_cast<size_t>(row)] *
                          polynomial_basis[static_cast<size_t>(col)]);
    }
  }

  vector<ZZ_pE> dual_basis(static_cast<size_t>(n), ZZ_pE(0));
  for (long target = 0; target < n; ++target) {
    vector<ZZ_p> rhs(static_cast<size_t>(n), ZZ_p(0));
    rhs[static_cast<size_t>(target)] = NTL::to_ZZ_p(1);
    const vector<ZZ_p> solution =
        SolveLinearSystemWithUnitPivotsOrThrow(trace_matrix, rhs, func_name);

    ZZ_pE dual_element = ZZ_pE(0);
    for (long col = 0; col < n; ++col) {
      dual_element += BaseRingConstant(solution[static_cast<size_t>(col)]) *
                      polynomial_basis[static_cast<size_t>(col)];
    }
    dual_basis[static_cast<size_t>(target)] = dual_element;
  }
  return dual_basis;
}

vector<ZZ_p> RecoverCoordsWithDualBasis(const vector<ZZ_pE> &dual_basis,
                                        const ZZ_pE &element,
                                        const char *func_name) {
  ValidateBasisShapeOrThrow(dual_basis, "dual_basis", func_name);
  const long n = static_cast<long>(dual_basis.size());
  vector<ZZ_p> coords(static_cast<size_t>(n), ZZ_p(0));
  for (long i = 0; i < n; ++i) {
    coords[static_cast<size_t>(i)] =
        TraceToBaseRing(dual_basis[static_cast<size_t>(i)] * element);
  }
  return coords;
}

ZZ_pE ComposeFromBasisCoordsOrThrow(const vector<ZZ_pE> &basis,
                                    const vector<ZZ_p> &coords,
                                    const char *func_name) {
  ValidateBasisShapeOrThrow(basis, "basis", func_name);
  const long n = static_cast<long>(basis.size());
  if (static_cast<long>(coords.size()) != n) {
    LogicError((string(func_name) + ": coords.size() must equal basis.size()")
                   .c_str());
  }

  ZZ_pE out = ZZ_pE(0);
  for (long i = 0; i < n; ++i) {
    out += BaseRingConstant(coords[static_cast<size_t>(i)]) *
           basis[static_cast<size_t>(i)];
  }
  return out;
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

PowerBasisFrobeniusContext BuildPowerBasisFrobeniusContext(
    const FrobeniusBasisParams &params, const ZZ_pE &teichmuller_generator,
    const char *func_name) {
  PowerBasisFrobeniusContext context;
  context.p = params.p;
  context.p_long = PositiveZZToLongChecked(params.p, "p", func_name);
  context.r = params.r;
  context.power_basis.basis =
      BuildTeichmullerPowerBasis(teichmuller_generator, params.r);
  if (!IsInvertibleOverBaseRing(BuildCoordinateMatrix(context.power_basis.basis))) {
    LogicError((string(func_name) +
                ": Teichmuller power basis does not span the extension")
                   .c_str());
  }
  context.power_basis.dual = BuildDualBasisByLinearSolveOrThrow(
      context.power_basis.basis, func_name);

  context.tau_images.resize(static_cast<size_t>(params.r));
  for (long i = 0; i < params.r; ++i) {
    context.tau_images[static_cast<size_t>(i)] =
        power(teichmuller_generator, context.p_long * i);
  }
  return context;
}

ZZ_pE ApplyPowerBasisTau(const PowerBasisFrobeniusContext &context,
                         const ZZ_pE &element) {
  const vector<ZZ_p> coords = RecoverCoordsWithDualBasis(
      context.power_basis.dual, element, "ApplyPowerBasisTau");
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
  if (!IsInvertibleOverBaseRing(BuildCoordinateMatrix(orbit))) {
    return false;
  }

  normal_basis_out->beta = orbit;
  normal_basis_out->alpha =
      BuildDualBasisByLinearSolveOrThrow(orbit, "FindNormalBasisOrThrow");
  return true;
}

}  // namespace

void ValidateFrobeniusBasisContextsOrThrow(
    const FrobeniusBasisParams &params, const ZZ_pX &extension_modulus) {
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

FrobeniusBasisData BuildFrobeniusBasisOrThrow(
    const FrobeniusBasisParams &params, const ZZ_pX &extension_modulus) {
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

NormalBasisData FindNormalBasisOrThrow(
    const FrobeniusBasisParams &params,
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

std::vector<ZZ_pE> BuildDualBasisOrThrow(const vector<ZZ_pE> &basis) {
  return BuildDualBasisByLinearSolveOrThrow(basis, "BuildDualBasisOrThrow");
}

void ValidateNormalBasisOrThrow(const NormalBasisData &normal_basis) {
  const char *const func_name = "ValidateNormalBasisOrThrow";
  ValidateBasisShapeOrThrow(normal_basis.beta, "normal_basis.beta", func_name);
  ValidateBasisShapeOrThrow(normal_basis.alpha, "normal_basis.alpha",
                            func_name);
  if (!IsInvertibleOverBaseRing(BuildCoordinateMatrix(normal_basis.beta))) {
    LogicError((string(func_name) + ": normal_basis.beta is not a basis")
                   .c_str());
  }

  const long n = static_cast<long>(normal_basis.beta.size());
  for (long u = 0; u < n; ++u) {
    for (long v = 0; v < n; ++v) {
      const ZZ_p trace_value =
          TraceToBaseRing(normal_basis.alpha[static_cast<size_t>(u)] *
                          normal_basis.beta[static_cast<size_t>(v)]);
      const ZZ_p expected = (u == v) ? NTL::to_ZZ_p(1) : NTL::to_ZZ_p(0);
      if (trace_value != expected) {
        LogicError((string(func_name) +
                    ": dual-basis trace identity failed")
                       .c_str());
      }
    }
  }
}

bool IsBaseRingConstant(const ZZ_pE &value) {
  const ZZ_pX poly = rep(value);
  const long degree = deg(poly);
  for (long i = 1; i <= degree; ++i) {
    if (coeff(poly, i) != 0) {
      return false;
    }
  }
  return true;
}

ZZ_p TraceToBaseRing(const ZZ_pE &element) { return trace(element); }

std::vector<ZZ_p> RecoverNormalBasisCoords(
    const NormalBasisData &normal_basis, const ZZ_pE &element) {
  ValidateBasisShapeOrThrow(normal_basis.beta, "normal_basis.beta",
                            "RecoverNormalBasisCoords");
  ValidateBasisShapeOrThrow(normal_basis.alpha, "normal_basis.alpha",
                            "RecoverNormalBasisCoords");
  return RecoverCoordsWithDualBasis(normal_basis.alpha, element,
                                    "RecoverNormalBasisCoords");
}

ZZ_pE ComposeFromNormalBasisCoords(const NormalBasisData &normal_basis,
                                   const vector<ZZ_p> &coords) {
  return ComposeFromBasisCoordsOrThrow(normal_basis.beta, coords,
                                       "ComposeFromNormalBasisCoords");
}

ZZ_pE ApplyFrobeniusTau(const NormalBasisData &normal_basis,
                        const ZZ_pE &element) {
  const vector<ZZ_p> coords =
      RecoverNormalBasisCoords(normal_basis, element);
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
  const vector<ZZ_p> coords =
      RecoverNormalBasisCoords(normal_basis, element);
  vector<ZZ_p> rotated(coords.size(), ZZ_p(0));
  if (!coords.empty()) {
    for (long i = 0; i + 1 < static_cast<long>(coords.size()); ++i) {
      rotated[static_cast<size_t>(i)] = coords[static_cast<size_t>(i + 1)];
    }
    rotated.back() = coords.front();
  }
  return ComposeFromNormalBasisCoords(normal_basis, rotated);
}
