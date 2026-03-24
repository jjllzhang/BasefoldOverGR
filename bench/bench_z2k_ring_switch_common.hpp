#ifndef BASEFOLD_BENCH_Z2K_RING_SWITCH_COMMON_HPP_
#define BASEFOLD_BENCH_Z2K_RING_SWITCH_COMMON_HPP_

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "bench_common_helpers.hpp"
#include "Compiler/Z2k/BaseFoldBackendAdapter.hpp"
#include "Compiler/Z2k/RingSwitchPCS.hpp"
#include "Compiler/Z2k/RingSwitchProofSerialize.hpp"
#include "GaloisRing/PrimitiveElement.hpp"
#include "PCS/Common/Hash.hpp"
#include "PCS/Common/Multilinear.hpp"

namespace basefold_bench_z2k_ring_switch_common {

using NTL::LogicError;
using NTL::SetCoeff;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pEPush;
using NTL::ZZ_pPush;
using NTL::ZZ_pX;
using NTL::conv;
using NTL::power2_ZZ;
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

enum class BasisMode {
  kDefault,
  kProvided,
};

enum class BenchDefaultBasisPreset {
  kCommitLike,
  kEvalLike,
};

enum class BenchShiftDirection {
  kLower,
  kUpper,
};

struct BenchShiftBasisSpec {
  const char *name = "poly";
  BenchShiftDirection direction = BenchShiftDirection::kLower;
  unsigned shift_mask = 0;
  bool use_dual_as_primal = false;
};

struct BasisCliData {
  bool has_basis = false;
  std::vector<ZZ_pE> basis;
  bool has_dual_basis = false;
  std::vector<ZZ_pE> dual_basis;
};

struct EncodedBasisCliData {
  bool has_basis = false;
  std::string basis_spec;
  bool has_dual_basis = false;
  std::string dual_basis_spec;
};

struct RingSwitchBenchCliArgs {
  ContextSpec context;
  BasisMode basis_mode = BasisMode::kDefault;
  EncodedBasisCliData alpha;
  EncodedBasisCliData beta;
};

struct RingSwitchBenchCliConfig {
  ContextSpec context;
  BasisMode basis_mode = BasisMode::kDefault;
  bool has_default_basis_preset = false;
  BenchDefaultBasisPreset default_basis_preset =
      BenchDefaultBasisPreset::kEvalLike;
  std::string default_alpha_preset_name;
  std::string default_beta_preset_name;
  BasisCliData alpha;
  BasisCliData beta;
};


inline std::vector<ZZ_pE> ParseBasisElementList(const std::string &s) {
  std::vector<ZZ_pE> out;
  std::size_t pos = 0;
  while (pos < s.size()) {
    const std::size_t sep = s.find(';', pos);
    const std::size_t end = (sep == std::string::npos) ? s.size() : sep;
    std::string token = s.substr(pos, end - pos);
    const std::size_t first = token.find_first_not_of(" \t");
    const std::size_t last = token.find_last_not_of(" \t");
    if (first == std::string::npos) {
      LogicError("ParseBasisElementList: empty basis element");
    }
    token = token.substr(first, last - first + 1);
    out.push_back(BuildZZpE(ParseCoeffList(token)));
    pos = (sep == std::string::npos) ? s.size() : (sep + 1);
  }
  if (out.empty()) {
    LogicError("ParseBasisElementList: empty basis list");
  }
  return out;
}

inline bool ParseBasisMode(const char *s, BasisMode &out) {
  if (s == nullptr) {
    return false;
  }
  const std::string mode(s);
  if (mode == "default") {
    out = BasisMode::kDefault;
    return true;
  }
  if (mode == "provided") {
    out = BasisMode::kProvided;
    return true;
  }
  return false;
}

inline const char *BasisModeName(BasisMode mode) {
  return mode == BasisMode::kProvided ? "provided" : "default";
}

inline const char *BenchDefaultBasisPresetName(
    BenchDefaultBasisPreset preset) {
  return preset == BenchDefaultBasisPreset::kCommitLike ? "commit-like"
                                                        : "eval-like";
}

inline void ApplyBenchDefaultBasisPresetOrThrow(
    RingSwitchBenchCliConfig &config, const ZZ_pX &F,
    BenchDefaultBasisPreset preset);

inline const char *NeedValueOrExit(int &i, int argc, char **argv,
                                   const char *flag) {
  if (i + 1 >= argc) {
    std::cerr << "Missing value for " << flag << "\n";
    std::exit(2);
  }
  return argv[++i];
}

inline void SetEncodedBasisCliData(EncodedBasisCliData &out,
                                   const std::string &encoded,
                                   bool is_dual_basis) {
  if (is_dual_basis) {
    out.has_dual_basis = true;
    out.dual_basis_spec = encoded;
  } else {
    out.has_basis = true;
    out.basis_spec = encoded;
  }
}

inline bool TryParseBasisCliArg(const std::string &arg, int &i, int argc,
                                char **argv, RingSwitchBenchCliArgs &args) {
  if (arg == "--basis-mode") {
    if (!ParseBasisMode(NeedValueOrExit(i, argc, argv, "--basis-mode"),
                        args.basis_mode)) {
      std::cerr << "Invalid --basis-mode\n";
      std::exit(2);
    }
    return true;
  }
  if (arg == "--alpha-basis") {
    SetEncodedBasisCliData(
        args.alpha, NeedValueOrExit(i, argc, argv, "--alpha-basis"), false);
    return true;
  }
  if (arg == "--beta-basis") {
    SetEncodedBasisCliData(
        args.beta, NeedValueOrExit(i, argc, argv, "--beta-basis"), false);
    return true;
  }
  if (arg == "--alpha-dual-basis") {
    SetEncodedBasisCliData(args.alpha,
                           NeedValueOrExit(i, argc, argv, "--alpha-dual-basis"),
                           true);
    return true;
  }
  if (arg == "--beta-dual-basis") {
    SetEncodedBasisCliData(args.beta,
                           NeedValueOrExit(i, argc, argv, "--beta-dual-basis"),
                           true);
    return true;
  }
  return false;
}

inline void ValidateBasisCliDataOrThrow(const BasisCliData &data,
                                        long expected_dimension,
                                        const char *label,
                                        const char *func_name) {
  if (data.has_basis &&
      static_cast<long>(data.basis.size()) != expected_dimension) {
    LogicError((std::string(func_name) + ": " + label +
                ".basis size must equal extension degree")
                   .c_str());
  }
  if (data.has_dual_basis &&
      static_cast<long>(data.dual_basis.size()) != expected_dimension) {
    LogicError((std::string(func_name) + ": " + label +
                ".dual_basis size must equal extension degree")
                   .c_str());
  }
  if (data.has_dual_basis && !data.has_basis) {
    LogicError((std::string(func_name) + ": " + label +
                ".dual_basis requires the corresponding basis")
                   .c_str());
  }
}

inline BasisCliData DecodeBasisCliDataOrThrow(const EncodedBasisCliData &encoded,
                                              long expected_dimension,
                                              const char *label,
                                              const char *func_name) {
  BasisCliData out;
  out.has_basis = encoded.has_basis;
  out.has_dual_basis = encoded.has_dual_basis;
  if (encoded.has_basis) {
    out.basis = ParseBasisElementList(encoded.basis_spec);
  }
  if (encoded.has_dual_basis) {
    out.dual_basis = ParseBasisElementList(encoded.dual_basis_spec);
  }
  ValidateBasisCliDataOrThrow(out, expected_dimension, label, func_name);
  return out;
}

inline RingSwitchBenchCliConfig DecodeRingSwitchBenchCliConfigOrThrow(
    const RingSwitchBenchCliArgs &args, const ZZ_pX &F) {
  RingSwitchBenchCliConfig out;
  out.context = args.context;
  out.basis_mode = args.basis_mode;

  const long basis_dimension = NTL::deg(F);
  if (basis_dimension <= 0) {
    LogicError(
        "DecodeRingSwitchBenchCliConfigOrThrow: extension degree must be positive");
  }

  if (args.basis_mode == BasisMode::kDefault) {
    if (args.alpha.has_basis || args.alpha.has_dual_basis || args.beta.has_basis ||
        args.beta.has_dual_basis) {
      LogicError(
          "DecodeRingSwitchBenchCliConfigOrThrow: basis flags require --basis-mode provided");
    }
    return out;
  }

  out.alpha = DecodeBasisCliDataOrThrow(args.alpha, basis_dimension, "alpha",
                                        "DecodeRingSwitchBenchCliConfigOrThrow");
  out.beta = DecodeBasisCliDataOrThrow(args.beta, basis_dimension, "beta",
                                       "DecodeRingSwitchBenchCliConfigOrThrow");
  if (!out.alpha.has_basis || !out.beta.has_basis) {
    LogicError(
        "DecodeRingSwitchBenchCliConfigOrThrow: provided basis mode requires both alpha and beta");
  }
  return out;
}

inline RingSwitchBenchCliConfig DecodeRingSwitchBenchCliConfigOrThrow(
    const RingSwitchBenchCliArgs &args, const ZZ_pX &F,
    BenchDefaultBasisPreset default_basis_preset) {
  RingSwitchBenchCliConfig out =
      DecodeRingSwitchBenchCliConfigOrThrow(args, F);
  if (out.basis_mode == BasisMode::kDefault) {
    ApplyBenchDefaultBasisPresetOrThrow(out, F, default_basis_preset);
  }
  return out;
}

inline basefold::FoldableCodeParams BuildBackendParams(long c, long d,
                                                       const ZZ &p,
                                                       const ZZ_pE &zeta) {
  return basefold_bench_common::BuildFoldableParamsK0Eq1(
      c, d, p, zeta, "BuildBackendParams");
}

inline GaloisRingBasisData BuildBenchFixedPolynomialBasisDataOrThrow(
    long basis_dimension, const char *func_name) {
  GaloisRingBasisData basis_data;
  if (basis_dimension <= 0) {
    LogicError((std::string(func_name) +
                ": extension degree must be positive")
                   .c_str());
  }
  basis_data.basis = BuildPolynomialBasis(basis_dimension);
  basis_data.dual_basis = BuildDualBasisOrThrow(basis_data.basis);
  return basis_data;
}

inline std::vector<long> BuildBenchShiftListFromMask(unsigned shift_mask) {
  static const long kShiftOptions[] = {1, 2, 4, 8, 16, 32, 64};
  std::vector<long> shifts;
  for (std::size_t i = 0; i < std::size(kShiftOptions); ++i) {
    if ((shift_mask & (1u << i)) != 0u) {
      shifts.push_back(kShiftOptions[i]);
    }
  }
  return shifts;
}

inline std::vector<ZZ_pE> ApplyBenchShiftTransformToPolynomialBasisOrThrow(
    const std::vector<ZZ_pE> &polynomial_basis, BenchShiftDirection direction,
    unsigned shift_mask, const char *func_name) {
  ValidateBasisShapeOrThrow(polynomial_basis, "polynomial_basis", func_name);
  if (shift_mask == 0u) {
    return polynomial_basis;
  }
  const std::vector<long> shifts = BuildBenchShiftListFromMask(shift_mask);
  const long basis_dimension = static_cast<long>(polynomial_basis.size());
  std::vector<ZZ_pE> out(static_cast<std::size_t>(basis_dimension), ZZ_pE(0));
  for (long i = 0; i < basis_dimension; ++i) {
    ZZ_pE element = polynomial_basis[static_cast<std::size_t>(i)];
    for (const long shift : shifts) {
      const long source = direction == BenchShiftDirection::kLower ? i - shift
                                                                   : i + shift;
      if (source >= 0 && source < basis_dimension) {
        element += polynomial_basis[static_cast<std::size_t>(source)];
      }
    }
    out[static_cast<std::size_t>(i)] = element;
  }
  return out;
}

inline GaloisRingBasisData BuildBenchShiftBasisDataOrThrow(
    long basis_dimension, const BenchShiftBasisSpec &spec,
    const char *func_name) {
  GaloisRingBasisData basis_data =
      BuildBenchFixedPolynomialBasisDataOrThrow(basis_dimension, func_name);
  basis_data.basis = ApplyBenchShiftTransformToPolynomialBasisOrThrow(
      basis_data.basis, spec.direction, spec.shift_mask, func_name);
  basis_data.dual_basis = BuildDualBasisOrThrow(basis_data.basis);
  if (spec.use_dual_as_primal) {
    std::swap(basis_data.basis, basis_data.dual_basis);
  }
  return basis_data;
}

inline bool MatchesBenchRingSwitchF64DefaultContext(const ContextSpec &context,
                                                    const ZZ_pX &F) {
  if (context.base_prime != ZZ(2) || context.scalar_modulus <= 0 ||
      NTL::deg(F) != 64) {
    return false;
  }
  const ZZ_pX expected_F = BuildZZpX(ParseCoeffList(
      "1,1,1,0,0,1,1,1,0,1,1,1,1,0,1,1,0,0,1,0,1,0,0,0,1,0,1,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,1,0,0,1,1,1,0,1,0,1,0,0,1,1,0,0,0,0,0,1,0,1"));
  if (F != expected_F) {
    return false;
  }
  return context.scalar_modulus == power2_ZZ(16) ||
         context.scalar_modulus == power2_ZZ(32) ||
         context.scalar_modulus == power2_ZZ(64);
}

inline bool MatchesBenchRingSwitchF128DefaultContext(const ContextSpec &context,
                                                     const ZZ_pX &F) {
  if (context.base_prime != ZZ(2) || context.scalar_modulus <= 0 ||
      NTL::deg(F) != 128) {
    return false;
  }
  const ZZ_pX expected_F = BuildZZpX(ParseCoeffList(
      "1,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"
      "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"
      "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"
      "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"
      "1"));
  return F == expected_F && context.scalar_modulus == power2_ZZ(16);
}

inline void SelectBenchDefaultBasisPresetSpecsOrThrow(
    const ContextSpec &context, const ZZ_pX &F, BenchDefaultBasisPreset preset,
    BenchShiftBasisSpec &alpha_spec, BenchShiftBasisSpec &beta_spec) {
  alpha_spec = BenchShiftBasisSpec();
  beta_spec = BenchShiftBasisSpec();
  if (MatchesBenchRingSwitchF128DefaultContext(context, F)) {
    beta_spec.name = "lower_64";
    beta_spec.shift_mask = 1u << 6;
    return;
  }
  if (!MatchesBenchRingSwitchF64DefaultContext(context, F)) {
    return;
  }
  if (context.scalar_modulus == power2_ZZ(16)) {
    if (preset == BenchDefaultBasisPreset::kCommitLike) {
      beta_spec.name = "lower_16";
      beta_spec.shift_mask = 1u << 4;
    } else {
      beta_spec.name = "lower_32";
      beta_spec.shift_mask = 1u << 5;
    }
    return;
  }
  if (context.scalar_modulus == power2_ZZ(32)) {
    if (preset == BenchDefaultBasisPreset::kCommitLike) {
      beta_spec.name = "lower_32";
      beta_spec.shift_mask = 1u << 5;
    } else {
      alpha_spec.name = "poly_dual";
      alpha_spec.use_dual_as_primal = true;
      beta_spec.name = "lower_16";
      beta_spec.shift_mask = 1u << 4;
    }
    return;
  }
}

inline void ApplyBenchDefaultBasisPresetOrThrow(
    RingSwitchBenchCliConfig &config, const ZZ_pX &F,
    BenchDefaultBasisPreset preset) {
  const long basis_dimension = NTL::deg(F);
  if (basis_dimension <= 0) {
    LogicError(
        "ApplyBenchDefaultBasisPresetOrThrow: extension degree must be positive");
  }
  BenchShiftBasisSpec alpha_spec;
  BenchShiftBasisSpec beta_spec;
  SelectBenchDefaultBasisPresetSpecsOrThrow(config.context, F, preset,
                                            alpha_spec, beta_spec);
  config.has_default_basis_preset = true;
  config.default_basis_preset = preset;
  config.default_alpha_preset_name = alpha_spec.name;
  config.default_beta_preset_name = beta_spec.name;
  const GaloisRingBasisData alpha_basis =
      BuildBenchShiftBasisDataOrThrow(basis_dimension, alpha_spec,
                                      "ApplyBenchDefaultBasisPresetOrThrow");
  const GaloisRingBasisData beta_basis =
      BuildBenchShiftBasisDataOrThrow(basis_dimension, beta_spec,
                                      "ApplyBenchDefaultBasisPresetOrThrow");
  config.alpha.has_basis = true;
  config.alpha.basis = alpha_basis.basis;
  config.alpha.has_dual_basis = true;
  config.alpha.dual_basis = alpha_basis.dual_basis;
  config.beta.has_basis = true;
  config.beta.basis = beta_basis.basis;
  config.beta.has_dual_basis = true;
  config.beta.dual_basis = beta_basis.dual_basis;
}

inline basefold::RingSwitchPCSSetupInput BuildRingSwitchSetupInput(
    long c, long ell, long kappa, const RingSwitchBenchCliConfig &config,
    const ZZ_pX &F, const ZZ_pE &zeta) {
  ZZ p_base;
  long k_base = 0;
  DeduceBasePrimeAndExponent(config.context, p_base, k_base);
  (void)k_base;

  const long basis_dimension = NTL::deg(F);
  if (basis_dimension <= 0) {
    LogicError("BuildRingSwitchSetupInput: extension degree must be positive");
  }
  ValidateBasisCliDataOrThrow(config.alpha, basis_dimension, "alpha",
                              "BuildRingSwitchSetupInput");
  ValidateBasisCliDataOrThrow(config.beta, basis_dimension, "beta",
                              "BuildRingSwitchSetupInput");

  basefold::RingSwitchPCSSetupInput input;
  input.ell = ell;
  input.kappa = kappa;
  input.base_modulus = config.context.scalar_modulus;
  input.extension_modulus = F;
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(
      BuildBackendParams(c, ell - kappa, p_base, zeta));
  if (config.basis_mode == BasisMode::kProvided ||
      (config.alpha.has_basis && config.beta.has_basis)) {
    if (!config.alpha.has_basis || !config.beta.has_basis) {
      LogicError(
          "BuildRingSwitchSetupInput: active bench basis preset requires both alpha and beta");
    }
    input.use_provided_basis = true;
    input.provided_basis.has_alpha_basis = true;
    input.provided_basis.alpha_basis.basis = config.alpha.basis;
    if (config.alpha.has_dual_basis) {
      input.provided_basis.alpha_basis.dual_basis = config.alpha.dual_basis;
    }
    input.provided_basis.has_beta_basis = true;
    input.provided_basis.beta_basis.basis = config.beta.basis;
    if (config.beta.has_dual_basis) {
      input.provided_basis.beta_basis.dual_basis = config.beta.dual_basis;
    }
  } else {
    input.use_provided_basis = true;
    input.provided_basis.has_alpha_basis = true;
    input.provided_basis.alpha_basis =
        BuildBenchFixedPolynomialBasisDataOrThrow(basis_dimension,
                                                  "BuildRingSwitchSetupInput");
    input.provided_basis.has_beta_basis = true;
    input.provided_basis.beta_basis =
        BuildBenchFixedPolynomialBasisDataOrThrow(basis_dimension,
                                                  "BuildRingSwitchSetupInput");
  }
  return input;
}

inline basefold::RingSwitchPCSParams BuildRingSwitchParams(
    long c, long ell, long kappa, const ContextSpec &spec, const ZZ_pX &F,
    const ZZ_pE &zeta) {
  RingSwitchBenchCliConfig config;
  config.context = spec;
  return basefold::RingSwitchPCSSetup(
      BuildRingSwitchSetupInput(c, ell, kappa, config, F, zeta));
}

inline basefold::RingSwitchPCSParams BuildRingSwitchParams(
    long c, long ell, long kappa, const RingSwitchBenchCliConfig &config,
    const ZZ_pX &F, const ZZ_pE &zeta) {
  return basefold::RingSwitchPCSSetup(
      BuildRingSwitchSetupInput(c, ell, kappa, config, F, zeta));
}

inline void PrintProvidedBasisFlagHelp(std::ostream &os, const char *indent) {
  os << indent << "[--basis-mode <default|provided>]\n"
     << indent << "[--alpha-basis <elem0;elem1;...>] [--beta-basis <elem0;elem1;...>]\n"
     << indent
     << "[--alpha-dual-basis <elem0;elem1;...>] [--beta-dual-basis <elem0;elem1;...>]\n";
}

inline void PrintBasisCliNotes(std::ostream &os, const char *indent) {
  os << indent
     << "Default mode uses hardcoded bench presets for selected GR(2^k,r) contexts and polynomial fallback elsewhere.\n"
     << indent
     << "In provided mode, basis elements are entered in polynomial-basis coordinates.\n"
     << indent
     << "If a provided dual basis is omitted, setup derives it automatically.\n";
}

inline const char *DualBasisSourceName(const BasisCliData &basis) {
  return basis.has_dual_basis ? "supplied" : "derived";
}

inline void PrintBasisModeSummary(std::ostream &os,
                                  const RingSwitchBenchCliConfig &config,
                                  const basefold::RingSwitchPCSParams &params) {
  if (config.basis_mode == BasisMode::kProvided) {
    os << "  basis mode provided\n";
  } else {
    os << "  basis mode default (bench-fixed "
       << BenchDefaultBasisPresetName(config.default_basis_preset)
       << " preset)\n";
    if (!config.default_alpha_preset_name.empty()) {
      os << "  default alpha preset " << config.default_alpha_preset_name
         << "\n";
    }
    if (!config.default_beta_preset_name.empty()) {
      os << "  default beta preset " << config.default_beta_preset_name
         << "\n";
    }
  }
  if (config.alpha.has_basis) {
    os << "  alpha basis dim " << params.alpha_basis.basis.size() << " (dual "
       << DualBasisSourceName(config.alpha) << ")\n";
  }
  if (config.beta.has_basis) {
    os << "  beta basis dim " << params.beta_basis.basis.size() << " (dual "
       << DualBasisSourceName(config.beta) << ")\n";
  }
}

}  // namespace basefold_bench_z2k_ring_switch_common

#endif  // BASEFOLD_BENCH_Z2K_RING_SWITCH_COMMON_HPP_
