#include "BaseFold/Sumcheck.hpp"

#include <NTL/ZZ.h>

#include <cstddef>
#include <vector>

#include "BaseFold/Multilinear.hpp"

using NTL::LogicError;
using NTL::vec_ZZ_pE;
using NTL::ZZ_pE;

namespace basefold {
namespace {

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

// In-place subset-sum (zeta) transform:
//   out[mask] = Σ_{S ⊆ mask} in[S]
vec_ZZ_pE BooleanEvalTableFromMonomialCoeffs(const vec_ZZ_pE &coeffs, long k) {
  const long n = coeffs.length();
  if (k < 0) LogicError("BooleanEvalTableFromMonomialCoeffs: negative k");
  if (n != (1L << k))
    LogicError("BooleanEvalTableFromMonomialCoeffs: length mismatch");

  vec_ZZ_pE eval = coeffs;
  for (long bit = 0; bit < k; ++bit) {
    const long step = 1L << bit;
    for (long mask = 0; mask < n; ++mask) {
      if (mask & step) {
        eval[mask] += eval[mask ^ step];
      }
    }
  }
  return eval;
}

std::vector<ZZ_pE> PrefixEqProducts(const std::vector<ZZ_pE> &z, long k_minus_1) {
  if (k_minus_1 < 0)
    LogicError("PrefixEqProducts: negative k_minus_1");
  if (k_minus_1 > static_cast<long>(z.size()))
    LogicError("PrefixEqProducts: k_minus_1 exceeds dimension");

  const ZZ_pE one = One();
  std::vector<ZZ_pE> prod;
  prod.resize(1);
  prod[0] = one;

  for (long var = 0; var < k_minus_1; ++var) {
    const ZZ_pE z_var = z[static_cast<std::size_t>(var)];
    const ZZ_pE f0 = one - z_var;  // factor(var, 0)
    const ZZ_pE f1 = z_var;        // factor(var, 1)

    const std::size_t old = prod.size();
    prod.resize(old * 2);
    for (std::size_t mask = 0; mask < old; ++mask) {
      const ZZ_pE base = prod[mask];
      prod[mask] = base * f0;
      prod[mask + old] = base * f1;
    }
  }

  return prod;
}

}  // namespace

FieldElement QuadraticPoly::Eval(const FieldElement &x) const {
  return a0 + a1 * x + a2 * x * x;
}

SumcheckProver::SumcheckProver(const FieldVec &f_coeffs,
                               const std::vector<FieldElement> &z)
    : z_(z) {
  const long n = f_coeffs.length();
  if (!IsPowerOfTwoLong(n))
    LogicError("SumcheckProver: f_coeffs length must be 2^d");
  d_ = Log2ExactPowerOfTwoLong(n);
  if (z_.size() != static_cast<std::size_t>(d_))
    LogicError("SumcheckProver: z dimension mismatch");

  cur_k_ = d_;
  coeffs_rem_ = f_coeffs;
  suffix_eq_prod_ = One();
}

QuadraticPoly SumcheckProver::CurrentPolynomial() const {
  if (cur_k_ <= 0)
    LogicError("SumcheckProver::CurrentPolynomial: no remaining variables");

  const long k = cur_k_;
  const long n = coeffs_rem_.length();
  if (n != (1L << k))
    LogicError("SumcheckProver::CurrentPolynomial: internal length mismatch");

  const vec_ZZ_pE eval = BooleanEvalTableFromMonomialCoeffs(coeffs_rem_, k);
  const long half = 1L << (k - 1);

  const std::vector<ZZ_pE> prefix = PrefixEqProducts(z_, k - 1);
  if (static_cast<long>(prefix.size()) != half)
    LogicError("SumcheckProver::CurrentPolynomial: prefix size mismatch");

  const ZZ_pE one = One();
  const ZZ_pE z_k = z_[static_cast<std::size_t>(k - 1)];
  const ZZ_pE factor0 = one - z_k;      // eq-factor at X_k=0
  const ZZ_pE factor1 = z_k;            // eq-factor at X_k=1
  const ZZ_pE delta_factor = factor1 - factor0;

  QuadraticPoly out;
  out.a0 = ZZ_pE(0);
  out.a1 = ZZ_pE(0);
  out.a2 = ZZ_pE(0);

  for (long mask = 0; mask < half; ++mask) {
    const ZZ_pE common =
        prefix[static_cast<std::size_t>(mask)] * suffix_eq_prod_;

    const ZZ_pE eq0 = common * factor0;
    const ZZ_pE delta_eq = common * delta_factor;

    const ZZ_pE f0 = eval[mask];
    const ZZ_pE f1 = eval[mask + half];
    const ZZ_pE delta_f = f1 - f0;

    out.a0 += f0 * eq0;
    out.a1 += f0 * delta_eq + delta_f * eq0;
    out.a2 += delta_f * delta_eq;
  }

  return out;
}

void SumcheckProver::ReceiveChallenge(const FieldElement &r_kminus1) {
  if (cur_k_ <= 0)
    LogicError("SumcheckProver::ReceiveChallenge: no remaining variables");

  const long k = cur_k_;
  const long n = coeffs_rem_.length();
  if (n != (1L << k))
    LogicError("SumcheckProver::ReceiveChallenge: internal length mismatch");

  suffix_eq_prod_ *= EqFactor(z_[static_cast<std::size_t>(k - 1)], r_kminus1);

  const long half = n / 2;
  for (long i = 0; i < half; ++i) {
    coeffs_rem_[i] = coeffs_rem_[i] + coeffs_rem_[i + half] * r_kminus1;
  }
  coeffs_rem_.SetLength(half);
  --cur_k_;
}

bool CheckSumcheckRelations(const std::vector<QuadraticPoly> &h_by_level,
                            const std::vector<FieldElement> &r,
                            const FieldElement &claimed_y) {
  const long d = static_cast<long>(h_by_level.size());
  if (d < 0)
    return false;
  if (static_cast<long>(r.size()) != d)
    return false;
  if (d == 0) {
    return claimed_y == claimed_y;
  }

  const ZZ_pE one = One();
  const ZZ_pE zero = ZZ_pE(0);

  const QuadraticPoly &h_d = h_by_level[static_cast<std::size_t>(d - 1)];
  if (h_d.Eval(zero) + h_d.Eval(one) != claimed_y)
    return false;

  for (long k = 1; k < d; ++k) {
    const QuadraticPoly &h_k = h_by_level[static_cast<std::size_t>(k - 1)];
    const QuadraticPoly &h_kp1 = h_by_level[static_cast<std::size_t>(k)];
    const ZZ_pE lhs = h_k.Eval(zero) + h_k.Eval(one);
    const ZZ_pE rhs = h_kp1.Eval(r[static_cast<std::size_t>(k)]);
    if (lhs != rhs)
      return false;
  }

  return true;
}

}  // namespace basefold

