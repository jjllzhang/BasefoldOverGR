#include "BaseFold/Multilinear.hpp"

#include <NTL/ZZ.h>

#include <cstddef>

using NTL::LogicError;
using NTL::vec_ZZ_pE;
using NTL::ZZ_pE;

namespace basefold {
namespace {

bool IsPowerOfTwoLong(long n) {
  return n > 0 && (n & (n - 1)) == 0;
}

long Log2ExactPowerOfTwoLong(long n) {
  if (!IsPowerOfTwoLong(n)) LogicError("Log2ExactPowerOfTwoLong: not a power of two");
  long d = 0;
  while (n > 1) {
    n >>= 1;
    ++d;
  }
  return d;
}

ZZ_pE One() {
  ZZ_pE one;
  NTL::set(one);
  return one;
}

}  // namespace

FieldElement EvalMultilinearMonomialCoeffs(
    const FieldVec &coeffs, const std::vector<FieldElement> &point) {
  const long n = coeffs.length();
  if (!IsPowerOfTwoLong(n)) LogicError("EvalMultilinearMonomialCoeffs: coeffs length must be 2^d");
  const long d = Log2ExactPowerOfTwoLong(n);
  if (point.size() != static_cast<std::size_t>(d)) {
    LogicError("EvalMultilinearMonomialCoeffs: point dimension mismatch");
  }

  vec_ZZ_pE cur = coeffs;
  long cur_len = n;
  for (long var = d; var-- > 0;) {
    const long half = cur_len / 2;
    for (long i = 0; i < half; ++i) {
      cur[i] = cur[i] + cur[i + half] * point[static_cast<std::size_t>(var)];
    }
    cur.SetLength(half);
    cur_len = half;
  }
  return cur[0];
}

FieldElement EqFactor(const FieldElement &z_i, const FieldElement &x_i) {
  const ZZ_pE one = One();
  return z_i * x_i + (one - z_i) * (one - x_i);
}

FieldElement EqPolynomial(const std::vector<FieldElement> &z,
                          const std::vector<FieldElement> &x) {
  if (z.size() != x.size()) LogicError("EqPolynomial: dimension mismatch");

  ZZ_pE acc = One();
  for (std::size_t i = 0; i < z.size(); ++i) {
    acc *= EqFactor(z[i], x[i]);
  }
  return acc;
}

}  // namespace basefold

