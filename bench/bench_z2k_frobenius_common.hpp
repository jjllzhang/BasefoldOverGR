#ifndef BASEFOLD_BENCH_Z2K_FROBENIUS_COMMON_HPP_
#define BASEFOLD_BENCH_Z2K_FROBENIUS_COMMON_HPP_

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

#include <iostream>
#include <string>

#include "bench_common_helpers.hpp"
#include "Compiler/Z2k/BaseFoldBackendAdapter.hpp"
#include "Compiler/Z2k/FrobeniusPCS.hpp"
#include "GaloisRing/PrimitiveElement.hpp"
#include "PCS/Common/Hash.hpp"
#include "PCS/Common/Multilinear.hpp"

namespace basefold_bench_z2k_frobenius_common {

using NTL::LogicError;
using NTL::SetCoeff;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pEPush;
using NTL::ZZ_pPush;
using NTL::ZZ_pX;
using NTL::conv;
using NTL::to_ZZ;
using NTL::vec_ZZ_pE;

using basefold_bench_common::BaseRingConstant;
using basefold_bench_common::BuildZZpE;
using basefold_bench_common::BuildZZpX;
using basefold_bench_common::ComputeStats;
using basefold_bench_common::DeduceBasePrimeAndExponent;
using basefold_bench_common::DeterministicResidue;
using basefold_bench_common::ExtensionElementBytesOrThrow;
using basefold_bench_common::FixedCoeffBytesOrThrow;
using basefold_bench_common::MakeDeterministicBaseRingTable;
using basefold_bench_common::MakeDeterministicElement;
using basefold_bench_common::MakeDeterministicPoint;
using basefold_bench_common::MsSince;
using basefold_bench_common::MulU64OrThrow;
using basefold_bench_common::NormalizeMod;
using basefold_bench_common::PackedVectorFixedBytesOrThrow;
using basefold_bench_common::ParseCoeffList;
using basefold_bench_common::ParseInt;
using basefold_bench_common::ParseLong;
using basefold_bench_common::ParseU64OrDie;
using basefold_bench_common::ParseZZ;
using basefold_bench_common::ParseZZString;
using basefold_bench_common::Pow2Checked;
using basefold_bench_common::SplitMix64;
using basefold_bench_common::Stats;
using basefold_bench_common::ValidateMonic;
using basefold_bench_common::ZZFromU64;

using ContextSpec = basefold_bench_common::BasicContextSpec;

inline std::vector<ZZ_pE> BooleanPointFromIndex(long index, long dimension) {
  std::vector<ZZ_pE> point(static_cast<std::size_t>(dimension));
  for (long i = 0; i < dimension; ++i) {
    point[static_cast<std::size_t>(i)] = BaseRingConstant(to_ZZ((index >> i) & 1L));
  }
  return point;
}

inline ZZ_pE EvalFromBooleanTable(const vec_ZZ_pE &table, long dimension,
                                  const std::vector<ZZ_pE> &point) {
  if (table.length() != Pow2Checked(dimension)) {
    LogicError("EvalFromBooleanTable: table length mismatch");
  }
  ZZ_pE acc = ZZ_pE(0);
  for (long idx = 0; idx < table.length(); ++idx) {
    acc += table[idx] * basefold::EqPolynomial(point,
                                               BooleanPointFromIndex(idx, dimension));
  }
  return acc;
}

inline basefold::FoldableCodeParams BuildBackendParams(long c, long d,
                                                       const ZZ &p,
                                                       const ZZ_pE &zeta) {
  return basefold_bench_common::BuildFoldableParamsK0Eq1(
      c, d, p, zeta, "BuildBackendParams");
}

inline basefold::FrobeniusPCSParams BuildFrobeniusParams(
    long c, long ell, long kappa, const ContextSpec &spec, const ZZ_pX &F,
    const ZZ_pE &zeta) {
  ZZ p_base;
  long k_base = 0;
  DeduceBasePrimeAndExponent(spec, p_base, k_base);

  basefold::FrobeniusPCSSetupInput input;
  input.ell = ell;
  input.kappa = kappa;
  input.base_modulus = spec.scalar_modulus;
  input.extension_modulus = F;
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(
      BuildBackendParams(c, ell - kappa, p_base, zeta));
  return basefold::FrobeniusPCSSetup(input);
}

}  // namespace basefold_bench_z2k_frobenius_common

#endif  // BASEFOLD_BENCH_Z2K_FROBENIUS_COMMON_HPP_
