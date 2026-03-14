#include "GaloisRing/PrimitiveElement.hpp"

#include <sstream>
#include <string>
#include <vector>

using NTL::coeff;
using NTL::conv;
using NTL::NumBits;
using NTL::power;
using NTL::random;
using NTL::rep;
using NTL::set;
using NTL::SetCoeff;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pBak;
using NTL::ZZ_pE;
using NTL::ZZ_pEBak;
using NTL::ZZ_pX;

namespace {

void ValidateExtensionPolynomial(const ZZ_pX &F, long s,
                                 const std::string &fn_name) {
  if (NTL::deg(F) != s) {
    const std::string msg = fn_name + ": deg(F) must equal s";
    NTL::LogicError(msg.c_str());
  }
  if (!NTL::IsOne(coeff(F, s))) {
    const std::string msg = fn_name + ": F must be monic";
    NTL::LogicError(msg.c_str());
  }
}

bool IsUnitByReductionModP(const ZZ_pE &a, const ZZ &p) {
  if (a == 0)
    return false;

  const ZZ_pX poly = rep(a);
  const long r = ZZ_pE::degree();
  for (long i = 0; i < r; ++i) {
    ZZ c = rep(coeff(poly, i));
    c %= p;
    if (c < 0)
      c += p;
    if (c != 0)
      return true;
  }
  return false;
}

unsigned long long PositiveZZToU64Checked(const ZZ &z, const std::string &label,
                                          const std::string &fn_name) {
  if (z <= 0) {
    const std::string msg = fn_name + ": " + label + " must be > 0";
    NTL::LogicError(msg.c_str());
  }
  if (NumBits(z) > 64) {
    const std::string msg = fn_name + ": " + label + " too large for u64";
    NTL::LogicError(msg.c_str());
  }
  std::ostringstream os;
  os << z;
  return std::stoull(os.str());
}

std::vector<ZZ> UniquePrimeFactorsU64(unsigned long long n) {
  std::vector<ZZ> factors;
  if (n <= 1)
    return factors;

  for (unsigned long long d = 2; d <= n / d; ++d) {
    if (n % d != 0)
      continue;
    factors.push_back(NTL::to_ZZ(static_cast<long>(d)));
    while (n % d == 0)
      n /= d;
  }
  if (n > 1) {
    std::ostringstream os;
    os << n;
    factors.push_back(NTL::to_ZZ(os.str().c_str()));
  }
  return factors;
}

bool HasExactOrder(const ZZ_pE &a, const ZZ &order,
                   const std::vector<ZZ> &prime_factors) {
  if (a == 0 || order <= 0)
    return false;
  ZZ_pE one;
  set(one);

  if (power(a, order) != one)
    return false;
  for (const ZZ &q : prime_factors) {
    if (power(a, order / q) == one)
      return false;
  }
  return true;
}

} // namespace

/*
    Returns a candidate "primitive element" over the Galois ring GR(p^k, s)
   (implementation-specific).

    IMPORTANT:
      - The extension/modulus polynomial F must be supplied by the caller and
        should match the intended (p, k, s) context. In particular, F must have
        degree s and be monic.
      - This function re-initializes NTL contexts (ZZ_p::init and ZZ_pE::init)
   internally, but restores the incoming contexts on return.

    Usage: ZZ_pE alpha = FindPrimitiveElement(p, k, s, F);
*/
ZZ_pE FindPrimitiveElement(ZZ p, long k, long s, const ZZ_pX &F) {
  if (p <= 1)
    NTL::LogicError("FindPrimitiveElement: p must be > 1");
  if (k <= 0)
    NTL::LogicError("FindPrimitiveElement: k must be > 0");
  if (s <= 0)
    NTL::LogicError("FindPrimitiveElement: s must be > 0");

  ZZ_pBak modulus_bak;
  modulus_bak.save();
  ZZ_pEBak extension_bak;
  extension_bak.save();

  const ZZ q3 = power(p, k);
  ZZ q6 = power(p, k - 1);

  long q5;
  conv(q5, q6);

  ZZ_p::init(q3);
  ValidateExtensionPolynomial(F, s, "FindPrimitiveElement");
  ZZ_pE::init(F);

  ZZ_pX H;
  SetCoeff(H, 1, 1);
  ZZ_pE b;
  conv(b, H);

  return power(b, q5);
}

ZZ_pE FindTeichmullerGenerator(ZZ p, long k, long s, const ZZ_pX &F,
                               long max_trials) {
  if (p <= 1)
    NTL::LogicError("FindTeichmullerGenerator: p must be > 1");
  if (k <= 0)
    NTL::LogicError("FindTeichmullerGenerator: k must be > 0");
  if (s <= 0)
    NTL::LogicError("FindTeichmullerGenerator: s must be > 0");
  if (max_trials <= 0)
    NTL::LogicError("FindTeichmullerGenerator: max_trials must be > 0");

  ZZ_pBak modulus_bak;
  modulus_bak.save();
  ZZ_pEBak extension_bak;
  extension_bak.save();

  const ZZ q3 = power(p, k);
  const ZZ projection_exp_zz = power(p, k - 1);
  const ZZ subgroup_order_zz = power(p, s) - ZZ(1);

  const unsigned long long subgroup_order_u64 = PositiveZZToU64Checked(
      subgroup_order_zz, "p^s-1", "FindTeichmullerGenerator");
  const std::vector<ZZ> subgroup_order_factors =
      UniquePrimeFactorsU64(subgroup_order_u64);

  ZZ_p::init(q3);
  ValidateExtensionPolynomial(F, s, "FindTeichmullerGenerator");
  ZZ_pE::init(F);

  ZZ_pX H;
  SetCoeff(H, 1, 1);
  ZZ_pE x;
  conv(x, H);

  const ZZ_pE deterministic = power(x, projection_exp_zz);
  if (HasExactOrder(deterministic, subgroup_order_zz, subgroup_order_factors)) {
    return deterministic;
  }

  for (long iter = 0; iter < max_trials; ++iter) {
    ZZ_pE unit_candidate;
    random(unit_candidate);
    if (!IsUnitByReductionModP(unit_candidate, p))
      continue;

    const ZZ_pE projected = power(unit_candidate, projection_exp_zz);
    if (HasExactOrder(projected, subgroup_order_zz, subgroup_order_factors)) {
      return projected;
    }
  }

  NTL::LogicError("FindTeichmullerGenerator: failed within max_trials");
  return ZZ_pE(0);
}
