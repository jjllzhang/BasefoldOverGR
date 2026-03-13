#include "GaloisRing/Basis.hpp"

#include <NTL/ZZ_pX.h>

#include <numeric>
#include <string>
#include <utility>
#include <vector>

using NTL::GCD;
using NTL::LogicError;
using NTL::SetCoeff;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pX;
using NTL::coeff;
using NTL::conv;
using NTL::deg;
using NTL::inv;
using NTL::rep;
using NTL::trace;
using std::size_t;
using std::string;
using std::swap;
using std::vector;

namespace {

string QualifiedLabel(const char *label, const char *suffix) {
  return string(label) + suffix;
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

bool IsUnit(const ZZ_p &value) {
  if (value == 0) {
    return false;
  }
  return GCD(rep(value), ZZ_p::modulus()) == 1;
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

}  // namespace

std::vector<ZZ_pE> BuildPolynomialBasis(long dimension) {
  if (!ZZ_pE::initialized()) {
    LogicError("BuildPolynomialBasis: ZZ_pE context must be initialized");
  }
  if (dimension <= 0) {
    LogicError("BuildPolynomialBasis: dimension must be > 0");
  }

  vector<ZZ_pE> basis(static_cast<size_t>(dimension));
  for (long i = 0; i < dimension; ++i) {
    ZZ_pX poly;
    SetCoeff(poly, i, 1);
    conv(basis[static_cast<size_t>(i)], poly);
  }
  return basis;
}

void ValidateBasisShapeOrThrow(const vector<ZZ_pE> &basis, const char *label,
                               const char *func_name) {
  if (!ZZ_pE::initialized()) {
    LogicError((string(func_name) + ": ZZ_pE context must be initialized")
                   .c_str());
  }
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

void ValidateBasisDataOrThrow(const GaloisRingBasisData &basis,
                              const char *label, const char *func_name) {
  const string basis_label = QualifiedLabel(label, ".basis");
  const string dual_label = QualifiedLabel(label, ".dual_basis");
  ValidateBasisShapeOrThrow(basis.basis, basis_label.c_str(), func_name);
  ValidateBasisShapeOrThrow(basis.dual_basis, dual_label.c_str(), func_name);
  if (basis.basis.size() != basis.dual_basis.size()) {
    LogicError((string(func_name) + ": " + label +
                ".basis and .dual_basis must have the same size")
                   .c_str());
  }
  if (!IsBasisOverBaseRing(basis.basis)) {
    LogicError((string(func_name) + ": " + label + ".basis is not a basis")
                   .c_str());
  }
  if (!IsBasisOverBaseRing(basis.dual_basis)) {
    LogicError((string(func_name) + ": " + label +
                ".dual_basis is not a basis")
                   .c_str());
  }

  const long n = static_cast<long>(basis.basis.size());
  for (long u = 0; u < n; ++u) {
    for (long v = 0; v < n; ++v) {
      const ZZ_p trace_value =
          TraceToBaseRing(basis.dual_basis[static_cast<size_t>(u)] *
                          basis.basis[static_cast<size_t>(v)]);
      const ZZ_p expected = (u == v) ? NTL::to_ZZ_p(1) : NTL::to_ZZ_p(0);
      if (trace_value != expected) {
        LogicError((string(func_name) + ": " + label +
                    " failed the dual-basis trace identity")
                       .c_str());
      }
    }
  }
}

bool IsBasisOverBaseRing(const vector<ZZ_pE> &basis) {
  if (!ZZ_pE::initialized()) {
    return false;
  }
  if (basis.empty() || static_cast<long>(basis.size()) != ZZ_pE::degree()) {
    return false;
  }
  vector<ZZ_p> rhs(static_cast<size_t>(basis.size()), ZZ_p(0));
  vector<ZZ_p> solution;
  return TrySolveLinearSystemWithUnitPivots(BuildCoordinateMatrix(basis), rhs,
                                            &solution);
}

std::vector<ZZ_pE> BuildDualBasisOrThrow(const vector<ZZ_pE> &basis) {
  return BuildDualBasisByLinearSolveOrThrow(basis, "BuildDualBasisOrThrow");
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

std::vector<ZZ_p> RecoverBasisCoordsOrThrow(const GaloisRingBasisData &basis,
                                            const ZZ_pE &element,
                                            const char *func_name) {
  const string basis_label = QualifiedLabel("basis", ".basis");
  const string dual_label = QualifiedLabel("basis", ".dual_basis");
  ValidateBasisShapeOrThrow(basis.basis, basis_label.c_str(), func_name);
  ValidateBasisShapeOrThrow(basis.dual_basis, dual_label.c_str(), func_name);
  if (basis.basis.size() != basis.dual_basis.size()) {
    LogicError((string(func_name) +
                ": basis.basis and basis.dual_basis must have the same size")
                   .c_str());
  }
  return RecoverCoordsWithDualBasis(basis.dual_basis, element, func_name);
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
