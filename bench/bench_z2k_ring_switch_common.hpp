#ifndef BASEFOLD_BENCH_Z2K_RING_SWITCH_COMMON_HPP_
#define BASEFOLD_BENCH_Z2K_RING_SWITCH_COMMON_HPP_

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

#include <iostream>
#include <string>

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

inline basefold::FoldableCodeParams BuildBackendParams(long c, long d,
                                                       const ZZ &p,
                                                       const ZZ_pE &zeta) {
  return basefold_bench_common::BuildFoldableParamsK0Eq1(
      c, d, p, zeta, "BuildBackendParams");
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
  if (config.basis_mode == BasisMode::kProvided) {
    if (!config.alpha.has_basis || !config.beta.has_basis) {
      LogicError(
          "BuildRingSwitchSetupInput: provided basis mode requires both alpha and beta");
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
     << "Default mode uses the active polynomial alpha/beta basis.\n"
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
  os << "  basis mode " << BasisModeName(config.basis_mode) << "\n";
  if (config.basis_mode != BasisMode::kProvided) {
    return;
  }
  os << "  alpha basis dim " << params.alpha_basis.basis.size() << " (dual "
     << DualBasisSourceName(config.alpha) << ")\n";
  os << "  beta basis dim " << params.beta_basis.basis.size() << " (dual "
     << DualBasisSourceName(config.beta) << ")\n";
}

}  // namespace basefold_bench_z2k_ring_switch_common

#endif  // BASEFOLD_BENCH_Z2K_RING_SWITCH_COMMON_HPP_
