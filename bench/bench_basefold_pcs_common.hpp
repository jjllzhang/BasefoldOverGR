#ifndef BASEFOLD_BENCH_PCS_COMMON_HPP_
#define BASEFOLD_BENCH_PCS_COMMON_HPP_

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pX.h>

#include <iostream>
#include <string>

#include "bench_common_helpers.hpp"
#include "GaloisRing/PrimitiveElement.hpp"
#include "PCS/BaseFold/BaseFoldPCS.hpp"
#include "PCS/BaseFold/ProofSize.hpp"
#include "PCS/Common/Hash.hpp"
#include "PCS/Common/Multilinear.hpp"
#include "PCS/Common/Profile.hpp"

namespace basefold_bench_pcs_common {

using NTL::LogicError;
using NTL::SetCoeff;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pEPush;
using NTL::ZZ_pEX;
using NTL::ZZ_pPush;
using NTL::ZZ_pX;
using NTL::conv;
using NTL::to_ZZ;
using NTL::vec_ZZ_pE;

using basefold_bench_common::BuildZZpE;
using basefold_bench_common::BuildZZpEX;
using basefold_bench_common::BuildZZpX;
using basefold_bench_common::ComputeStats;
using basefold_bench_common::DeduceBasePrimeAndExponent;
using basefold_bench_common::DeterministicResidue;
using basefold_bench_common::IsPowerOfTwoLong;
using basefold_bench_common::Log2ExactPowerOfTwoLong;
using basefold_bench_common::MakeDeterministicElement;
using basefold_bench_common::MakeDeterministicPoint;
using basefold_bench_common::MsSince;
using basefold_bench_common::NormalizeMod;
using basefold_bench_common::ParseCoeffList;
using basefold_bench_common::ParseInt;
using basefold_bench_common::ParseLong;
using basefold_bench_common::ParseNestedCoeffList;
using basefold_bench_common::ParseU64OrDie;
using basefold_bench_common::ParseZZ;
using basefold_bench_common::ParseZZString;
using basefold_bench_common::Pow2Checked;
using basefold_bench_common::SplitMix64;
using basefold_bench_common::Stats;
using basefold_bench_common::ValidateMonic;
using basefold_bench_common::ZZFromU64;

struct ContextSpec {
  std::string label;
  ZZ scalar_modulus = ZZ(0);
  ZZ base_prime = ZZ(0);
  std::vector<ZZ> F_coeffs;
  std::vector<ZZ> zeta_coeffs;
  std::vector<std::vector<ZZ>> challenge_ext_coeffs;
};

inline void DeduceBasePrimeAndExponent(const ContextSpec &spec, ZZ &p_out,
                                       long &k_out) {
  basefold_bench_common::DeduceBasePrimeAndExponent(spec, p_out, k_out);
}

inline basefold::FoldableCodeParams BuildParams_k0_1(long c, long d,
                                                     const ZZ &prime_p,
                                                     const ZZ_pE &zeta) {
  return basefold_bench_common::BuildFoldableParamsK0Eq1(
      c, d, prime_p, zeta, "BuildParams_k0_1");
}

inline NTL::mat_ZZ_pE BuildSystematicG0(long c, long k0) {
  if (c <= 0) {
    LogicError("BuildSystematicG0: c must be positive");
  }
  if (k0 <= 0) {
    LogicError("BuildSystematicG0: k0 must be positive");
  }
  if (c > std::numeric_limits<long>::max() / k0) {
    LogicError("BuildSystematicG0: overflow in n0");
  }
  const long n0 = c * k0;

  NTL::mat_ZZ_pE G0;
  G0.SetDims(k0, n0);

  const ZZ_pE one = ZZ_pE(1);
  for (long block = 0; block < c; ++block) {
    const long base = block * k0;
    for (long r = 0; r < k0; ++r) {
      G0[r][base + r] = one;
    }
  }
  return G0;
}

inline basefold::FoldableCodeParams BuildParams_k0_pow2(long c, long k0,
                                                        long d,
                                                        const ZZ &prime_p,
                                                        const ZZ_pE &zeta) {
  if (c <= 0) {
    LogicError("BuildParams_k0_pow2: c must be positive");
  }
  if (k0 <= 0) {
    LogicError("BuildParams_k0_pow2: k0 must be positive");
  }
  if (!IsPowerOfTwoLong(k0)) {
    LogicError("BuildParams_k0_pow2: k0 must be a power of two");
  }
  if (d < 0) {
    LogicError("BuildParams_k0_pow2: d must be non-negative");
  }
  if (c > std::numeric_limits<long>::max() / k0) {
    LogicError("BuildParams_k0_pow2: overflow in n0");
  }
  const long n0 = c * k0;

  basefold::FoldableCodeParams params;
  params.c = c;
  params.k0 = k0;
  params.d = d;
  params.p = prime_p;
  params.zeta = zeta;
  params.G0 = BuildSystematicG0(c, k0);

  params.diag_T.resize(static_cast<std::size_t>(d));
  for (long level = 0; level < d; ++level) {
    const long pow2 = Pow2Checked(level);
    if (n0 > std::numeric_limits<long>::max() / pow2) {
      LogicError("BuildParams_k0_pow2: overflow in n_i");
    }
    const long ni = n0 * pow2;
    params.diag_T[static_cast<std::size_t>(level)].SetLength(ni);
    for (long i = 0; i < ni; ++i) {
      params.diag_T[static_cast<std::size_t>(level)][i] = ZZ_pE(1);
    }
  }

  return params;
}

inline vec_ZZ_pE MakeDeterministicCoefficients(long coeff_count,
                                               std::uint64_t seed) {
  vec_ZZ_pE coeffs;
  coeffs.SetLength(coeff_count);
  for (long i = 0; i < coeff_count; ++i) {
    coeffs[i] = MakeDeterministicElement(seed ^ static_cast<std::uint64_t>(i));
  }
  return coeffs;
}

inline std::uint64_t ComputeProofSizeBytes(
    const basefold::BaseFoldPCSEvalProof &proof,
    bool use_extension_challenges, long challenge_ext_degree) {
  basefold::BaseFoldProofSizeOptions options;
  options.include_version_byte = true;
  if (use_extension_challenges || proof.extension.has_extension_payload) {
    options.challenge_ext_degree = challenge_ext_degree;
  }
  return basefold::BaseFoldPCSEvalProofSizeBytes(proof, options);
}

}  // namespace basefold_bench_pcs_common

#endif  // BASEFOLD_BENCH_PCS_COMMON_HPP_
