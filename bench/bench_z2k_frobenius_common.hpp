#ifndef BASEFOLD_BENCH_Z2K_FROBENIUS_COMMON_HPP_
#define BASEFOLD_BENCH_Z2K_FROBENIUS_COMMON_HPP_

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <vector>

#include "bench_common_helpers.hpp"
#include "Compiler/Z2k/BaseFoldBackendAdapter.hpp"
#include "Compiler/Z2k/FrobeniusPCS.hpp"
#include "GaloisRing/FrobeniusBasis.hpp"
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

enum class FrobeniusBenchBasisOrigin {
  kRuntimeSearch,
  kHardcodedCommitLike,
  kHardcodedEvalLike,
};

struct FrobeniusBenchBasisSummary {
  FrobeniusBenchBasisOrigin origin =
      FrobeniusBenchBasisOrigin::kRuntimeSearch;
  bool has_preset_descriptor = false;
  long preset_a0 = 0;
  long preset_a1 = 0;
  long preset_exponent = 0;
  long candidate_count = 0;
  long calibrated_candidate_count = 0;
  std::uint64_t score = 0;
  std::uint64_t beta_weight = 0;
  std::uint64_t alpha_weight = 0;
  std::uint64_t tau_alpha_weight = 0;
  double calibration_total_ms = 0.0;
};

struct FrobeniusBenchSetupResult {
  basefold::FrobeniusPCSParams params;
  FrobeniusBenchBasisSummary basis_summary;
};

struct FrobeniusBenchBasisSelection {
  NormalBasisData normal_basis;
  ZZ_pE teichmuller_generator;
  FrobeniusBenchBasisSummary summary;
};

struct FrobeniusBenchCandidate {
  NormalBasisData normal_basis;
  FrobeniusBenchBasisSummary summary;
};

enum class FrobeniusBenchCalibrationMode {
  kCommit,
  kOuterCommit,
  kEval,
  kOuterProve,
  kOuterVerify,
};

enum class FrobeniusBenchPresetFamily {
  kCommitLike,
  kEvalLike,
};

struct FrobeniusBenchPresetDescriptor {
  long ring_power = 0;
  FrobeniusBenchPresetFamily family = FrobeniusBenchPresetFamily::kCommitLike;
  const char *teich_coeffs = nullptr;
  long a0 = 0;
  long a1 = 0;
  long exponent = 0;
};

inline const char *FrobeniusBenchPresetFamilyName(
    FrobeniusBenchPresetFamily family) {
  switch (family) {
    case FrobeniusBenchPresetFamily::kCommitLike:
      return "commit-like";
    case FrobeniusBenchPresetFamily::kEvalLike:
      return "eval-like";
  }
  return "unknown";
}

inline FrobeniusBenchPresetFamily NormalizeFrobeniusBenchPresetFamily(
    FrobeniusBenchCalibrationMode mode) {
  switch (mode) {
    case FrobeniusBenchCalibrationMode::kCommit:
    case FrobeniusBenchCalibrationMode::kOuterCommit:
      return FrobeniusBenchPresetFamily::kCommitLike;
    case FrobeniusBenchCalibrationMode::kEval:
    case FrobeniusBenchCalibrationMode::kOuterProve:
    case FrobeniusBenchCalibrationMode::kOuterVerify:
      return FrobeniusBenchPresetFamily::kEvalLike;
  }
  return FrobeniusBenchPresetFamily::kCommitLike;
}

inline std::vector<ZZ_pE> BooleanPointFromIndex(long index, long dimension) {
  std::vector<ZZ_pE> point(static_cast<std::size_t>(dimension));
  for (long i = 0; i < dimension; ++i) {
    point[static_cast<std::size_t>(i)] =
        BaseRingConstant(to_ZZ((index >> i) & 1L));
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
    acc += table[idx] *
           basefold::EqPolynomial(point, BooleanPointFromIndex(idx, dimension));
  }
  return acc;
}

inline basefold::FoldableCodeParams BuildBackendParams(long c, long d,
                                                       const ZZ &p,
                                                       const ZZ_pE &zeta) {
  return basefold_bench_common::BuildFoldableParamsK0Eq1(
      c, d, p, zeta, "BuildBackendParams");
}

inline ZZ_pE LiftBaseCoeff(const ZZ_p &value) {
  ZZ_pX poly;
  if (value != 0) {
    SetCoeff(poly, 0, value);
  }
  ZZ_pE out;
  conv(out, poly);
  return out;
}

inline long ElementCoeffWeight(const ZZ_pE &element) {
  const ZZ_pX poly = NTL::rep(element);
  long weight = 0;
  for (long i = 0; i <= NTL::deg(poly); ++i) {
    if (NTL::coeff(poly, i) != 0) {
      ++weight;
    }
  }
  return weight;
}

inline std::uint64_t BasisCoeffWeight(const std::vector<ZZ_pE> &basis) {
  std::uint64_t total = 0;
  for (const ZZ_pE &element : basis) {
    total += static_cast<std::uint64_t>(ElementCoeffWeight(element));
  }
  return total;
}

inline long RotateBasisIndex(long index, long offset, long dimension) {
  long out = (index + offset) % dimension;
  if (out < 0) {
    out += dimension;
  }
  return out;
}

inline ZZ_pE ComposeFromCoordsWithBasisRowOrThrow(
    const std::vector<ZZ_pE> &basis_row, const std::vector<ZZ_p> &coords,
    const char *func_name) {
  if (basis_row.size() != coords.size()) {
    LogicError((std::string(func_name) +
                ": basis row length must equal coordinate length")
                   .c_str());
  }
  ZZ_pE out = ZZ_pE(0);
  for (long i = 0; i < static_cast<long>(coords.size()); ++i) {
    out += LiftBaseCoeff(coords[static_cast<std::size_t>(i)]) *
           basis_row[static_cast<std::size_t>(i)];
  }
  return out;
}

inline std::uint64_t BuildTauAlphaCoeffWeightOrThrow(
    const NormalBasisData &normal_basis) {
  const char *const func_name = "BuildTauAlphaCoeffWeightOrThrow";
  const long basis_dimension =
      static_cast<long>(normal_basis.beta.size());
  if (basis_dimension <= 0 ||
      static_cast<long>(normal_basis.alpha.size()) != basis_dimension) {
    LogicError((std::string(func_name) + ": invalid normal basis shape")
                   .c_str());
  }

  std::uint64_t total = 0;
  for (long u = 0; u < basis_dimension; ++u) {
    const std::vector<ZZ_p> alpha_coords =
        RecoverNormalBasisCoords(normal_basis,
                                 normal_basis.alpha[static_cast<std::size_t>(u)]);
    for (long power = 0; power < basis_dimension; ++power) {
      std::vector<ZZ_pE> tau_row(static_cast<std::size_t>(basis_dimension));
      for (long j = 0; j < basis_dimension; ++j) {
        tau_row[static_cast<std::size_t>(j)] =
            normal_basis.beta[static_cast<std::size_t>(
                RotateBasisIndex(j, power, basis_dimension))];
      }
      total += static_cast<std::uint64_t>(ElementCoeffWeight(
          ComposeFromCoordsWithBasisRowOrThrow(
              tau_row, alpha_coords, "BuildTauAlphaCoeffWeightOrThrow")));
    }
  }
  return total;
}

inline std::uint64_t AddWeightedCostOrThrow(std::uint64_t total,
                                            std::uint64_t factor,
                                            std::uint64_t weight,
                                            const char *func_name) {
  if (weight == 0 || factor == 0) {
    return total;
  }
  const std::uint64_t term =
      MulU64OrThrow(factor, weight, func_name);
  if (term > std::numeric_limits<std::uint64_t>::max() - total) {
    LogicError((std::string(func_name) + ": score overflow").c_str());
  }
  return total + term;
}

inline FrobeniusBenchBasisSummary SummarizeFrobeniusBenchBasisOrThrow(
    long ell, long kappa, const NormalBasisData &normal_basis) {
  const char *const func_name = "SummarizeFrobeniusBenchBasisOrThrow";
  const long ell_prime = ell - kappa;
  if (ell_prime < 0) {
    LogicError((std::string(func_name) + ": ell must be >= kappa").c_str());
  }

  FrobeniusBenchBasisSummary out;
  out.beta_weight = BasisCoeffWeight(normal_basis.beta);
  out.alpha_weight = BasisCoeffWeight(normal_basis.alpha);
  out.tau_alpha_weight = BuildTauAlphaCoeffWeightOrThrow(normal_basis);

  const std::uint64_t packed_length =
      static_cast<std::uint64_t>(Pow2Checked(ell_prime));
  const std::uint64_t basis_dimension =
      static_cast<std::uint64_t>(normal_basis.beta.size());
  const std::uint64_t ell_prime_u64 =
      static_cast<std::uint64_t>(ell_prime);

  // Approximate the outer hot path by counting coefficient work in:
  // - packing t -> t' against beta,
  // - suffix orbit construction / indexed tau powers against beta,
  // - point/partial decomposition against alpha,
  // - precomputed tau(alpha_u) linear combinations.
  out.score = AddWeightedCostOrThrow(
      out.score, packed_length, out.beta_weight, func_name);
  out.score = AddWeightedCostOrThrow(
      out.score, basis_dimension * (ell_prime_u64 + 1ULL), out.beta_weight,
      func_name);
  out.score = AddWeightedCostOrThrow(
      out.score, ell_prime_u64 + basis_dimension, out.alpha_weight,
      func_name);
  out.score += out.tau_alpha_weight;
  return out;
}

inline std::string EncodeBasisElementForDedup(const ZZ_pE &element) {
  std::string out;
  const ZZ_pX poly = NTL::rep(element);
  const long degree = ZZ_pE::degree();
  for (long i = 0; i < degree; ++i) {
    if (!out.empty()) {
      out.push_back(',');
    }
    out += std::to_string(conv<long>(NTL::rep(NTL::coeff(poly, i))));
  }
  return out;
}

inline std::string EncodeBasisForDedup(const std::vector<ZZ_pE> &basis) {
  std::string out;
  for (const ZZ_pE &element : basis) {
    if (!out.empty()) {
      out.push_back(';');
    }
    out += EncodeBasisElementForDedup(element);
  }
  return out;
}

struct FrobeniusBenchPowerBasisContext {
  std::vector<ZZ_pE> basis;
  std::vector<ZZ_pE> dual_basis;
  std::vector<ZZ_pE> tau_images;
};

inline FrobeniusBenchPowerBasisContext BuildPowerBasisContextOrThrow(
    const ZZ &p_base, long basis_dimension, const ZZ_pE &teichmuller_generator) {
  const char *const func_name = "BuildPowerBasisContextOrThrow";
  const long p_long = conv<long>(p_base);
  if (basis_dimension <= 0) {
    LogicError((std::string(func_name) +
                ": basis dimension must be positive")
                   .c_str());
  }

  FrobeniusBenchPowerBasisContext out;
  out.basis.resize(static_cast<std::size_t>(basis_dimension));
  out.basis[0] = BaseRingConstant(ZZ(1));
  for (long i = 1; i < basis_dimension; ++i) {
    out.basis[static_cast<std::size_t>(i)] =
        out.basis[static_cast<std::size_t>(i - 1)] * teichmuller_generator;
  }
  out.dual_basis = BuildDualBasisOrThrow(out.basis);
  out.tau_images.resize(static_cast<std::size_t>(basis_dimension));
  for (long i = 0; i < basis_dimension; ++i) {
    out.tau_images[static_cast<std::size_t>(i)] =
        NTL::power(teichmuller_generator, p_long * i);
  }
  return out;
}

inline ZZ_pE ApplyPowerBasisTauOrThrow(
    const FrobeniusBenchPowerBasisContext &ctx, const ZZ_pE &element) {
  GaloisRingBasisData basis_data;
  basis_data.basis = ctx.basis;
  basis_data.dual_basis = ctx.dual_basis;
  const std::vector<ZZ_p> coords = RecoverBasisCoordsOrThrow(
      basis_data, element, "ApplyPowerBasisTauOrThrow");
  return ComposeFromBasisCoordsOrThrow(ctx.tau_images, coords,
                                       "ApplyPowerBasisTauOrThrow");
}

inline bool TryPromoteNormalBasisCandidateOrThrow(
    const FrobeniusBenchPowerBasisContext &ctx, long basis_dimension,
    const ZZ_pE &candidate, NormalBasisData *out) {
  std::vector<ZZ_pE> orbit(static_cast<std::size_t>(basis_dimension),
                           ZZ_pE(0));
  orbit[0] = candidate;
  for (long i = 1; i < basis_dimension; ++i) {
    orbit[static_cast<std::size_t>(i)] =
        ApplyPowerBasisTauOrThrow(ctx, orbit[static_cast<std::size_t>(i - 1)]);
  }
  if (ApplyPowerBasisTauOrThrow(ctx, orbit.back()) != orbit.front()) {
    return false;
  }
  for (long i = 0; i < basis_dimension; ++i) {
    for (long j = i + 1; j < basis_dimension; ++j) {
      if (orbit[static_cast<std::size_t>(i)] ==
          orbit[static_cast<std::size_t>(j)]) {
        return false;
      }
    }
  }
  if (!IsBasisOverBaseRing(orbit)) {
    return false;
  }

  out->beta = orbit;
  out->alpha = BuildDualBasisOrThrow(orbit);
  return true;
}

inline bool IsBetterFrobeniusBenchBasis(
    const FrobeniusBenchBasisSummary &candidate,
    const FrobeniusBenchBasisSummary &best_so_far) {
  if (candidate.score != best_so_far.score) {
    return candidate.score < best_so_far.score;
  }
  if (candidate.tau_alpha_weight != best_so_far.tau_alpha_weight) {
    return candidate.tau_alpha_weight < best_so_far.tau_alpha_weight;
  }
  if (candidate.alpha_weight != best_so_far.alpha_weight) {
    return candidate.alpha_weight < best_so_far.alpha_weight;
  }
  return candidate.beta_weight < best_so_far.beta_weight;
}

inline double CalibrateFrobeniusBenchCandidateMsOrThrow(
    long c, long ell, long kappa, const ContextSpec &spec, const ZZ_pX &F,
    const ZZ_pE &zeta, const ZZ_pE &teichmuller_generator,
    const NormalBasisData &normal_basis,
    FrobeniusBenchCalibrationMode mode) {
  const int calibration_timed_reps = (NTL::deg(F) >= 64) ? 1 : 2;
  ZZ p_base;
  long k_base = 0;
  DeduceBasePrimeAndExponent(spec, p_base, k_base);

  basefold::FrobeniusPCSSetupInput input;
  input.ell = ell;
  input.kappa = kappa;
  input.base_modulus = spec.scalar_modulus;
  input.extension_modulus = F;
  input.use_provided_basis = true;
  input.provided_basis.normal_basis = normal_basis;
  input.provided_basis.has_teichmuller_generator = true;
  input.provided_basis.teichmuller_generator = teichmuller_generator;
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(
      BuildBackendParams(c, ell - kappa, p_base, zeta));
  const basefold::FrobeniusPCSParams params = basefold::FrobeniusPCSSetup(input);

  const vec_ZZ_pE t_table = MakeDeterministicBaseRingTable(
      Pow2Checked(ell), 0x51eed1234ULL);
  const std::vector<ZZ_pE> z = MakeDeterministicPoint(
      ell, 0x51eed5678ULL);
  const ZZ_pE claimed_s = EvalFromBooleanTable(t_table, ell, z);

  switch (mode) {
    case FrobeniusBenchCalibrationMode::kCommit: {
      (void)basefold::FrobeniusPCSBuildCommitArtifactsUnchecked(params,
                                                                t_table);
      double total_ms = 0.0;
      for (int rep = 0; rep < calibration_timed_reps; ++rep) {
        const auto t0 = std::chrono::steady_clock::now();
        (void)basefold::FrobeniusPCSBuildCommitArtifactsUnchecked(params,
                                                                  t_table);
        const auto t1 = std::chrono::steady_clock::now();
        total_ms += MsSince(t0, t1);
      }
      return total_ms / static_cast<double>(calibration_timed_reps);
    }
    case FrobeniusBenchCalibrationMode::kOuterCommit: {
      (void)basefold::FrobeniusPCSBuildOuterCommitArtifactsUnchecked(params,
                                                                     t_table);
      double total_ms = 0.0;
      for (int rep = 0; rep < calibration_timed_reps; ++rep) {
        const auto t0 = std::chrono::steady_clock::now();
        (void)basefold::FrobeniusPCSBuildOuterCommitArtifactsUnchecked(params,
                                                                       t_table);
        const auto t1 = std::chrono::steady_clock::now();
        total_ms += MsSince(t0, t1);
      }
      return total_ms / static_cast<double>(calibration_timed_reps);
    }
    case FrobeniusBenchCalibrationMode::kOuterProve: {
      const basefold::FrobeniusPCSOuterCommitArtifacts warmup_artifacts =
          basefold::FrobeniusPCSBuildOuterCommitArtifactsUnchecked(params,
                                                                   t_table);
      const basefold::MerkleRoot commitment = basefold::Z2kPCSBackendCommit(
          params.backend, warmup_artifacts.t_packed_monomial_coeffs);
      (void)basefold::FrobeniusPCSProveOuterEvalFromCommitArtifactsUnchecked(
              params, t_table, commitment, z, claimed_s, /*num_queries=*/2,
              warmup_artifacts);

      const basefold::FrobeniusPCSOuterCommitArtifacts outer_artifacts =
          basefold::FrobeniusPCSBuildOuterCommitArtifactsUnchecked(params,
                                                                   t_table);
      const basefold::MerkleRoot timed_commitment = basefold::Z2kPCSBackendCommit(
          params.backend, outer_artifacts.t_packed_monomial_coeffs);
      double total_ms = 0.0;
      for (int rep = 0; rep < calibration_timed_reps; ++rep) {
        const auto t0 = std::chrono::steady_clock::now();
        (void)basefold::FrobeniusPCSProveOuterEvalFromCommitArtifactsUnchecked(
                params, t_table, timed_commitment, z, claimed_s,
                /*num_queries=*/2, outer_artifacts);
        const auto t1 = std::chrono::steady_clock::now();
        total_ms += MsSince(t0, t1);
      }
      return total_ms / static_cast<double>(calibration_timed_reps);
    }
    case FrobeniusBenchCalibrationMode::kEval: {
      const basefold::FrobeniusPCSCommitArtifacts warmup_artifacts =
          basefold::FrobeniusPCSBuildCommitArtifactsUnchecked(params, t_table);
      const basefold::FrobeniusPCSEvalProof warmup_proof =
          basefold::FrobeniusPCSProveEvalFromCommitArtifactsUnchecked(
              params, t_table, z, claimed_s, /*num_queries=*/2,
              warmup_artifacts);
      (void)basefold::FrobeniusPCSVerifyEvalUnchecked(
          params, warmup_artifacts.commitment, z, claimed_s,
          /*num_queries=*/2, warmup_proof);

      double total_ms = 0.0;
      for (int rep = 0; rep < calibration_timed_reps; ++rep) {
        const basefold::FrobeniusPCSCommitArtifacts commit_artifacts =
            basefold::FrobeniusPCSBuildCommitArtifactsUnchecked(params, t_table);
        const auto t0 = std::chrono::steady_clock::now();
        const basefold::FrobeniusPCSEvalProof proof =
            basefold::FrobeniusPCSProveEvalFromCommitArtifactsUnchecked(
                params, t_table, z, claimed_s, /*num_queries=*/2,
                commit_artifacts);
        const auto t1 = std::chrono::steady_clock::now();
        (void)basefold::FrobeniusPCSVerifyEvalUnchecked(
            params, commit_artifacts.commitment, z, claimed_s,
            /*num_queries=*/2, proof);
        const auto t2 = std::chrono::steady_clock::now();
        total_ms += MsSince(t0, t1) + MsSince(t1, t2);
      }
      return total_ms / static_cast<double>(calibration_timed_reps);
    }
    case FrobeniusBenchCalibrationMode::kOuterVerify: {
      const basefold::FrobeniusPCSOuterCommitArtifacts outer_artifacts =
          basefold::FrobeniusPCSBuildOuterCommitArtifactsUnchecked(params,
                                                                   t_table);
      const basefold::MerkleRoot commitment = basefold::Z2kPCSBackendCommit(
          params.backend, outer_artifacts.t_packed_monomial_coeffs);
      const basefold::FrobeniusPCSOuterEvalProof proof =
          basefold::FrobeniusPCSProveOuterEvalFromCommitArtifactsUnchecked(
              params, t_table, commitment, z, claimed_s, /*num_queries=*/2,
              outer_artifacts);
      (void)basefold::FrobeniusPCSVerifyOuterEvalUnchecked(
          params, commitment, z, claimed_s, /*num_queries=*/2, proof);

      double total_ms = 0.0;
      for (int rep = 0; rep < calibration_timed_reps; ++rep) {
        const auto t0 = std::chrono::steady_clock::now();
        (void)basefold::FrobeniusPCSVerifyOuterEvalUnchecked(
            params, commitment, z, claimed_s, /*num_queries=*/2, proof);
        const auto t1 = std::chrono::steady_clock::now();
        total_ms += MsSince(t0, t1);
      }
      return total_ms / static_cast<double>(calibration_timed_reps);
    }
  }

  LogicError("CalibrateFrobeniusBenchCandidateMsOrThrow: unknown mode");
  return 0.0;
}

inline bool CoeffListsEqual(const std::vector<ZZ> &lhs,
                            const std::vector<ZZ> &rhs) {
  return lhs == rhs;
}

inline const std::vector<ZZ> &DefaultFrobeniusBenchRingF64Coefficients() {
  static const std::vector<ZZ> coeffs = ParseCoeffList(
      "1,1,1,0,0,1,1,1,0,1,1,1,1,0,1,1,0,0,1,0,1,0,0,0,1,0,1,0,0,0,0,0,"
      "1,1,0,1,0,0,0,0,0,0,0,1,0,0,1,1,1,0,1,0,1,0,0,1,1,0,0,0,0,0,1,0,1");
  return coeffs;
}

inline const std::vector<ZZ> &DefaultFrobeniusBenchRingZetaCoefficients() {
  static const std::vector<ZZ> coeffs = {to_ZZ(0), to_ZZ(1)};
  return coeffs;
}

inline const char *HardcodedFrobeniusBenchTeichGR2P16Commit() {
  return "45633,60894,6171,28286,31313,54020,30093,30737,42881,36340,64837,"
         "10874,27619,58618,8028,17647,21434,53557,13020,57960,49637,44951,"
         "9873,28419,33538,15230,24487,2314,5551,8593,43431,49621,30857,"
         "25470,24931,22200,4914,47650,15694,42543,47182,64000,34107,61286,"
         "7196,58692,17583,37262,22841,40475,46028,38084,50565,50864,7868,"
         "2156,23139,19317,57310,37588,52734,61140,59947,46582";
}

inline const char *HardcodedFrobeniusBenchTeichGR2P16Eval() {
  return "824,23993,58432,19654,58935,61882,35611,47037,59638,46673,16824,"
         "53297,52055,35257,12010,13220,26646,21751,51060,4642,58681,23759,"
         "14885,22615,54937,44100,51915,38636,313,30789,63199,63468,59859,"
         "9730,62376,53109,57346,40025,33393,27873,47240,53074,35565,48343,"
         "56427,39607,55074,57385,10807,29495,21381,59797,33909,45650,47393,"
         "53293,26918,64826,23047,30952,41635,30549,47575,37377";
}

inline const char *HardcodedFrobeniusBenchTeichGR2P32Commit() {
  return "3974052975,2754983006,4233417734,4286062292,2647968052,1152099070,"
         "3191733157,2197650923,2678686011,35917027,2621736083,3137249822,"
         "1490288372,1657281415,3802663348,2006581548,842067576,2887479654,"
         "367457937,3675219403,2433415119,2721732469,807700258,510556651,"
         "361712481,2970006975,485594793,471070752,2038679100,3666484067,"
         "4226130582,344323174,1937451277,2267160503,2867304793,1897316502,"
         "1693747775,2327826984,811285306,824014658,2819084680,3238314736,"
         "3476670004,1556945612,3841388351,1575770472,482152851,4149052931,"
         "834474169,1210037668,1178654359,2055291126,1584885930,2922721310,"
         "1212809260,4243758229,759418010,3522998699,3122852354,1723278475,"
         "213087313,4189060713,509320022,2555454484";
}

inline const char *HardcodedFrobeniusBenchTeichGR2P32Eval() {
  return "3575023388,324388685,1809141884,2440694764,1409980591,1310425013,"
         "2999281438,3414419292,1333184447,1459765452,774575863,420527729,"
         "3389464615,3600163478,2260671799,2904513476,1169385683,1964449836,"
         "4157797162,1022178604,2826109917,2346694200,1878387841,745805553,"
         "734286914,4123059607,3868246688,2054144111,1390739461,1934700618,"
         "2508764913,3244671549,3712221732,697619988,2368185306,3846405786,"
         "66765973,3787904848,3104258881,1239097115,2385014714,2831196982,"
         "466370487,614734239,899164129,3435731300,1074857077,4149637875,"
         "2426538900,1431661860,3389178296,1737305073,1818253879,2885546600,"
         "2520120042,3051710406,1886865358,361209921,604932529,1434100238,"
         "2998135094,2436902134,1621589318,4001084751";
}

inline const char *HardcodedFrobeniusBenchTeichGR2P64Commit() {
  return "3186157747636644455,6664373155789465109,-7467126605502528218,"
         "4738866194998687222,-6664898714321795214,4667173809883564059,"
         "7590691955035782485,-3213916289553817147,8868041461640204220,"
         "-4995930518861564922,-5874748727437011985,4548636185917099587,"
         "-6366416747481845307,6651394545241464163,6132756885820962452,"
         "-175396245523218901,1777624117269406575,932930460914344793,"
         "1222687882271927860,8898247032323614832,616518313758600897,"
         "-2375940611803283274,-1918104367781710431,-8037251819158249629,"
         "8433409356695854537,-4674411131257593124,60421177534402699,"
         "5318494741392306051,5382096564327236148,1073233429360821614,"
         "8533206948932758995,8760318816168552553,4246651848110607857,"
         "-2433016345851728173,2402485392567355292,3715347568670939551,"
         "901484155740385596,6386101317100216599,6186800466014839144,"
         "-7915374939124227648,5988552236457976066,4408079761145728006,"
         "-6140925610082486349,-1581714811113774628,-2995382121648495298,"
         "8283707250997976083,283944844832099598,183110046522961328,"
         "7935894059075804765,3320605857055605124,-3627218936484973742,"
         "5355618174986832003,-3856045703331435675,4275303121674936542,"
         "1465102657704972778,-3928901260856851930,-6336161739815395477,"
         "-7857253276191863916,5214267507850507125,7311097195948604319,"
         "5075575028024542693,1126906802058995161,856590484498440379,"
         "4138971116645268018";
}

inline const char *HardcodedFrobeniusBenchTeichGR2P64Eval() {
  return "82526573900118045,-2282837098308641719,-3306656548453409127,"
         "9022936343877642449,7068320003879135530,-2508813252541016339,"
         "5848795030410004195,-28596047853599206,2438342501074330054,"
         "-2516774141905409618,2947392116476571262,1645195605597459089,"
         "-3020815675360013600,-5310403043169701789,-3347433393807343710,"
         "-6257552943425278090,-2610385849123687788,6284470502933593051,"
         "6396925008995730551,-8206270101101774244,3085221556086522228,"
         "-985625458521267734,-1101593912545722449,-7617933044929969264,"
         "2782951515653038602,-2716079648649794869,-2869734189683392300,"
         "6204155532645004211,879184742927851371,813151942551041681,"
         "1524569054680098370,-8461672504676547194,-4102250126696170492,"
         "6879922801319668910,3904026284100670168,-3758321073789612513,"
         "-8432030617279528313,-7361219622093522768,-7696368562015306382,"
         "-5828020373698145880,4062697247319967260,-8742097471527109335,"
         "8857351818695337275,2937075236925714112,-2701092953530215611,"
         "-4551653721096530753,-3316084553416752398,2997777575058573407,"
         "7696214705588831498,8607039367301565458,2647778563341458604,"
         "1922474644614051097,-3532057046123138559,-633956283470510545,"
         "6226459519943423439,2492195671290678114,5414716148253397351,"
         "1077238954522780210,506963380075380767,-2664313219180469033,"
         "-164469715334574457,-5111078476474036964,6713763085052410424,"
         "-992480368530444889";
}

inline bool LookupHardcodedFrobeniusBenchPresetDescriptor(
    long ring_power, FrobeniusBenchPresetFamily family,
    FrobeniusBenchPresetDescriptor *out) {
  if (out == nullptr) {
    return false;
  }

  switch (family) {
    case FrobeniusBenchPresetFamily::kCommitLike:
      switch (ring_power) {
        case 16:
          *out = FrobeniusBenchPresetDescriptor{
              16, family, HardcodedFrobeniusBenchTeichGR2P16Commit(), 1, 1, 14};
          return true;
        case 32:
          *out = FrobeniusBenchPresetDescriptor{
              32, family, HardcodedFrobeniusBenchTeichGR2P32Commit(), 1, 1, 10};
          return true;
        case 64:
          *out = FrobeniusBenchPresetDescriptor{
              64, family, HardcodedFrobeniusBenchTeichGR2P64Commit(), 1, 1, 15};
          return true;
        default:
          return false;
      }
    case FrobeniusBenchPresetFamily::kEvalLike:
      switch (ring_power) {
        case 16:
          *out = FrobeniusBenchPresetDescriptor{
              16, family, HardcodedFrobeniusBenchTeichGR2P16Eval(), 1, 1, 15};
          return true;
        case 32:
          *out = FrobeniusBenchPresetDescriptor{
              32, family, HardcodedFrobeniusBenchTeichGR2P32Eval(), 1, 1, 8};
          return true;
        case 64:
          *out = FrobeniusBenchPresetDescriptor{
              64, family, HardcodedFrobeniusBenchTeichGR2P64Eval(), 0, 1, 2};
          return true;
        default:
          return false;
      }
  }

  return false;
}

inline bool TrySelectHardcodedFrobeniusBenchBasisOrThrow(
    long ell, long kappa, const ContextSpec &spec, const ZZ_pX &F,
    const ZZ_pE &zeta, FrobeniusBenchCalibrationMode mode,
    FrobeniusBenchBasisSelection *selection) {
  const char *const func_name = "TrySelectHardcodedFrobeniusBenchBasisOrThrow";
  if (selection == nullptr) {
    LogicError((std::string(func_name) + ": selection must not be null")
                   .c_str());
  }
  if (kappa != 6 || NTL::deg(F) != 64 || spec.base_prime != to_ZZ(2)) {
    return false;
  }
  if (!CoeffListsEqual(spec.F_coeffs, DefaultFrobeniusBenchRingF64Coefficients())) {
    return false;
  }
  if (zeta != BuildZZpE(DefaultFrobeniusBenchRingZetaCoefficients())) {
    return false;
  }

  long ring_power = 0;
  if (spec.scalar_modulus == NTL::power(to_ZZ(2), 16)) {
    ring_power = 16;
  } else if (spec.scalar_modulus == NTL::power(to_ZZ(2), 32)) {
    ring_power = 32;
  } else if (spec.scalar_modulus == NTL::power(to_ZZ(2), 64)) {
    ring_power = 64;
  } else {
    return false;
  }

  FrobeniusBenchPresetDescriptor preset;
  if (!LookupHardcodedFrobeniusBenchPresetDescriptor(
          ring_power, NormalizeFrobeniusBenchPresetFamily(mode), &preset)) {
    return false;
  }

  const ZZ_pE teichmuller_generator = BuildZZpE(ParseCoeffList(preset.teich_coeffs));
  const FrobeniusBenchPowerBasisContext power_ctx =
      BuildPowerBasisContextOrThrow(spec.base_prime, NTL::deg(F),
                                    teichmuller_generator);

  const ZZ_pE candidate =
      BaseRingConstant(to_ZZ(preset.a0)) +
      BaseRingConstant(to_ZZ(preset.a1)) *
          NTL::power(teichmuller_generator, preset.exponent);

  NormalBasisData normal_basis;
  if (!TryPromoteNormalBasisCandidateOrThrow(power_ctx, NTL::deg(F), candidate,
                                             &normal_basis)) {
    LogicError((std::string(func_name) +
                ": hardcoded preset candidate did not promote to a normal basis")
                   .c_str());
  }

  selection->normal_basis = normal_basis;
  selection->teichmuller_generator = teichmuller_generator;
  selection->summary =
      SummarizeFrobeniusBenchBasisOrThrow(ell, kappa, normal_basis);
  selection->summary.origin =
      (preset.family == FrobeniusBenchPresetFamily::kCommitLike)
          ? FrobeniusBenchBasisOrigin::kHardcodedCommitLike
          : FrobeniusBenchBasisOrigin::kHardcodedEvalLike;
  selection->summary.has_preset_descriptor = true;
  selection->summary.preset_a0 = preset.a0;
  selection->summary.preset_a1 = preset.a1;
  selection->summary.preset_exponent = preset.exponent;
  return true;
}

inline FrobeniusBenchBasisSelection SelectPreferredFrobeniusBenchBasisOrThrow(
    long c, long ell, long kappa, const ContextSpec &spec, const ZZ_pX &F,
    const ZZ_pE &zeta, FrobeniusBenchCalibrationMode mode) {
  const char *const func_name = "SelectPreferredFrobeniusBenchBasisOrThrow";
  ZZ p_base;
  long k_base = 0;
  DeduceBasePrimeAndExponent(spec, p_base, k_base);

  const long basis_dimension = NTL::deg(F);
  if (basis_dimension <= 0) {
    LogicError((std::string(func_name) +
                ": extension degree must be positive")
                   .c_str());
  }

  const long kDistinctCandidateLimit = (basis_dimension >= 64) ? 16 : 64;
  const long kPerFamilyAttemptLimit = (basis_dimension >= 64) ? 16 : 64;
  const long kAffineSearchLimit = (basis_dimension >= 64) ? 2 : 4;
  const long kCalibrationFinalistLimit = kDistinctCandidateLimit;

  FrobeniusBenchBasisSelection selection;
  std::set<std::string> seen;
  std::vector<FrobeniusBenchCandidate> candidates;

  const ZZ_pE teichmuller_generator =
      FindTeichmullerGenerator(p_base, k_base, basis_dimension, F);
  const FrobeniusBenchPowerBasisContext power_ctx =
      BuildPowerBasisContextOrThrow(p_base, basis_dimension,
                                    teichmuller_generator);

  auto consider_candidate = [&](const ZZ_pE &candidate) {
    if (static_cast<long>(candidates.size()) >= kDistinctCandidateLimit) {
      return;
    }
    NormalBasisData normal_basis;
    if (!TryPromoteNormalBasisCandidateOrThrow(power_ctx, basis_dimension,
                                               candidate, &normal_basis)) {
      return;
    }
    const std::string key = EncodeBasisForDedup(normal_basis.beta);
    if (!seen.insert(key).second) {
      return;
    }
    FrobeniusBenchBasisSummary summary =
        SummarizeFrobeniusBenchBasisOrThrow(ell, kappa, normal_basis);
    summary.candidate_count = static_cast<long>(seen.size());
    candidates.push_back(
        FrobeniusBenchCandidate{std::move(normal_basis), summary});
  };

  ZZ_pE teichmuller_power = BaseRingConstant(ZZ(1));
  for (long exponent = 0; exponent < kPerFamilyAttemptLimit &&
                          static_cast<long>(candidates.size()) <
                              kDistinctCandidateLimit;
       ++exponent) {
    consider_candidate(teichmuller_power);
    teichmuller_power *= teichmuller_generator;
  }

  for (long a0 = 0; a0 <= kAffineSearchLimit &&
                   static_cast<long>(candidates.size()) < kDistinctCandidateLimit;
       ++a0) {
    for (long a1 = 1; a1 <= kAffineSearchLimit &&
                     static_cast<long>(candidates.size()) <
                         kDistinctCandidateLimit;
         ++a1) {
      if (a0 == 0 && a1 == 1) {
        continue;
      }
      teichmuller_power = BaseRingConstant(ZZ(1));
      for (long exponent = 0; exponent < kPerFamilyAttemptLimit &&
                              static_cast<long>(candidates.size()) <
                                  kDistinctCandidateLimit;
           ++exponent) {
        consider_candidate(BaseRingConstant(to_ZZ(a0)) +
                           BaseRingConstant(to_ZZ(a1)) * teichmuller_power);
        teichmuller_power *= teichmuller_generator;
      }
    }
  }

  if (!candidates.empty()) {
    std::sort(candidates.begin(), candidates.end(),
              [](const FrobeniusBenchCandidate &lhs,
                 const FrobeniusBenchCandidate &rhs) {
                return IsBetterFrobeniusBenchBasis(lhs.summary, rhs.summary);
              });
    const long finalist_count = std::min<long>(
        static_cast<long>(candidates.size()), kCalibrationFinalistLimit);
    double best_calibration_ms = std::numeric_limits<double>::infinity();
    long best_index = 0;
    for (long i = 0; i < finalist_count; ++i) {
      const double calibration_ms = CalibrateFrobeniusBenchCandidateMsOrThrow(
          c, ell, kappa, spec, F, zeta, teichmuller_generator,
          candidates[static_cast<std::size_t>(i)].normal_basis, mode);
      candidates[static_cast<std::size_t>(i)].summary.calibration_total_ms =
          calibration_ms;
      candidates[static_cast<std::size_t>(i)].summary.calibrated_candidate_count =
          finalist_count;
      if (calibration_ms < best_calibration_ms) {
        best_calibration_ms = calibration_ms;
        best_index = i;
      } else if (calibration_ms == best_calibration_ms &&
                 IsBetterFrobeniusBenchBasis(
                     candidates[static_cast<std::size_t>(i)].summary,
                     candidates[static_cast<std::size_t>(best_index)].summary)) {
        best_index = i;
      }
    }

    selection.normal_basis =
        candidates[static_cast<std::size_t>(best_index)].normal_basis;
    selection.summary = candidates[static_cast<std::size_t>(best_index)].summary;
    selection.summary.candidate_count = static_cast<long>(candidates.size());
    selection.teichmuller_generator = teichmuller_generator;
    return selection;
  }

  FrobeniusBasisParams basis_input;
  basis_input.p = p_base;
  basis_input.k = k_base;
  basis_input.r = basis_dimension;
  basis_input.teichmuller_generator_max_trials = 1024;
  basis_input.affine_search_limit = kAffineSearchLimit;
  const FrobeniusBasisData fallback =
      BuildFrobeniusBasisOrThrow(basis_input, F);

  selection.normal_basis = fallback.normal_basis;
  selection.teichmuller_generator = fallback.teichmuller_generator;
  selection.summary =
      SummarizeFrobeniusBenchBasisOrThrow(ell, kappa, selection.normal_basis);
  selection.summary.candidate_count = 1;
  selection.summary.calibrated_candidate_count = 1;
  selection.summary.calibration_total_ms =
      CalibrateFrobeniusBenchCandidateMsOrThrow(
          c, ell, kappa, spec, F, zeta, selection.teichmuller_generator,
          selection.normal_basis, mode);
  return selection;
}

inline FrobeniusBenchSetupResult BuildFrobeniusBenchSetupResult(
    FrobeniusBenchCalibrationMode mode, long c, long ell, long kappa,
    const ContextSpec &spec, const ZZ_pX &F, const ZZ_pE &zeta) {
  ZZ p_base;
  long k_base = 0;
  DeduceBasePrimeAndExponent(spec, p_base, k_base);

  FrobeniusBenchBasisSelection selection;
  if (!TrySelectHardcodedFrobeniusBenchBasisOrThrow(ell, kappa, spec, F, zeta,
                                                    mode, &selection)) {
    selection = SelectPreferredFrobeniusBenchBasisOrThrow(c, ell, kappa, spec,
                                                          F, zeta, mode);
  }

  basefold::FrobeniusPCSSetupInput input;
  input.ell = ell;
  input.kappa = kappa;
  input.base_modulus = spec.scalar_modulus;
  input.extension_modulus = F;
  input.use_provided_basis = true;
  input.provided_basis.normal_basis = selection.normal_basis;
  input.provided_basis.has_teichmuller_generator = true;
  input.provided_basis.teichmuller_generator =
      selection.teichmuller_generator;
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(
      BuildBackendParams(c, ell - kappa, p_base, zeta));

  FrobeniusBenchSetupResult out;
  out.params = basefold::FrobeniusPCSSetup(input);
  out.basis_summary = selection.summary;
  return out;
}

inline basefold::FrobeniusPCSParams BuildFrobeniusParams(
    long c, long ell, long kappa, const ContextSpec &spec, const ZZ_pX &F,
    const ZZ_pE &zeta) {
  return BuildFrobeniusBenchSetupResult(FrobeniusBenchCalibrationMode::kEval,
                                        c, ell, kappa, spec, F, zeta)
      .params;
}

inline void PrintFrobeniusBasisSummary(
    std::ostream &os, const FrobeniusBenchBasisSummary &summary) {
  if (summary.origin != FrobeniusBenchBasisOrigin::kRuntimeSearch) {
    const FrobeniusBenchPresetFamily family =
        (summary.origin == FrobeniusBenchBasisOrigin::kHardcodedCommitLike)
            ? FrobeniusBenchPresetFamily::kCommitLike
            : FrobeniusBenchPresetFamily::kEvalLike;
    os << "  basis mode bench-hardcoded preferred-normal\n"
       << "  basis preset " << FrobeniusBenchPresetFamilyName(family);
    if (summary.has_preset_descriptor) {
      os << " (a0=" << summary.preset_a0 << ", a1=" << summary.preset_a1
         << ", exponent=" << summary.preset_exponent << ")";
    }
    os << "\n"
       << "  basis score " << summary.score << " (beta weight "
       << summary.beta_weight << ", alpha weight " << summary.alpha_weight
       << ", tau(alpha) weight " << summary.tau_alpha_weight << ")\n";
    return;
  }

  os << "  basis mode bench-fixed preferred-normal\n"
     << "  basis search distinct candidates " << summary.candidate_count
     << " (calibrated " << summary.calibrated_candidate_count << ")\n"
     << "  basis score " << summary.score << " (beta weight "
     << summary.beta_weight << ", alpha weight " << summary.alpha_weight
     << ", tau(alpha) weight " << summary.tau_alpha_weight << ")\n"
     << "  basis calibration total " << summary.calibration_total_ms
     << " ms\n";
}

}  // namespace basefold_bench_z2k_frobenius_common

#endif  // BASEFOLD_BENCH_Z2K_FROBENIUS_COMMON_HPP_
