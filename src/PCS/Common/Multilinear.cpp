#include "PCS/Common/Multilinear.hpp"

#include <NTL/ZZ.h>

#include <cstddef>
#include <limits>

#include "PCS/Common/NtlParallel.hpp"

using NTL::LogicError;
using NTL::vec_ZZ_pE;
using NTL::ZZ_pE;

namespace basefold {
namespace {

constexpr long kEqualityTableParallelThreshold = 256;

bool IsPowerOfTwoLong(long n) { return n > 0 && (n & (n - 1)) == 0; }

long Log2ExactPowerOfTwoLong(long n) {
  if (!IsPowerOfTwoLong(n))
    LogicError("Log2ExactPowerOfTwoLong: not a power of two");
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

} // namespace

FieldElement
EvalMultilinearMonomialCoeffs(const FieldVec &coeffs,
                              const std::vector<FieldElement> &point) {
  const long n = coeffs.length();
  if (!IsPowerOfTwoLong(n))
    LogicError("EvalMultilinearMonomialCoeffs: coeffs length must be 2^d");
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
  if (z.size() != x.size())
    LogicError("EqPolynomial: dimension mismatch");

  ZZ_pE acc = One();
  for (std::size_t i = 0; i < z.size(); ++i) {
    acc *= EqFactor(z[i], x[i]);
  }
  return acc;
}

FieldVec EqualityTableFromPoint(const std::vector<FieldElement> &point) {
  long length = 1;
  for (std::size_t i = 0; i < point.size(); ++i) {
    if (length > std::numeric_limits<long>::max() / 2) {
      LogicError("EqualityTableFromPoint: dimension is too large for long");
    }
    length *= 2;
  }

  vec_ZZ_pE table;
  table.SetLength(length);
  table[0] = One();

  long filled = 1;
  const ZZ_pE one = One();
  for (std::size_t var = 0; var < point.size(); ++var) {
    const ZZ_pE zero_branch = one - point[var];
    const ZZ_pE one_branch = point[var];
    const long block_size = filled;
    const auto expand_one_entry = [&](long idx) {
      const ZZ_pE prev = table[idx];
      table[idx + block_size] = prev * one_branch;
      table[idx] = prev * zero_branch;
    };
    pcs_common_internal::ForEachIndexMaybeParallel(
        0, block_size, kEqualityTableParallelThreshold, expand_one_entry);
    filled *= 2;
  }

  return table;
}

} // namespace basefold
