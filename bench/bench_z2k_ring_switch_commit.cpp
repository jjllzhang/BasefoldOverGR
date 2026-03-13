#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "bench_z2k_ring_switch_common.hpp"

namespace {

using namespace basefold_bench_z2k_ring_switch_common;

struct BenchResult {
  Stats packing;
  Stats commit;
  std::uint64_t anti_opt_checksum = 0;
};

BenchResult RunCommitBenchmark(const basefold::RingSwitchPCSParams &params,
                               const vec_ZZ_pE &t_table, int warmup,
                               int reps) {
  if (warmup < 0) {
    LogicError("RunCommitBenchmark: warmup must be >= 0");
  }
  if (reps <= 0) {
    LogicError("RunCommitBenchmark: reps must be > 0");
  }

  std::vector<double> packing_ms;
  std::vector<double> commit_ms;
  packing_ms.reserve(static_cast<std::size_t>(reps));
  commit_ms.reserve(static_cast<std::size_t>(reps));

  std::uint64_t anti_opt_checksum = 0;
  for (int iter = -warmup; iter < reps; ++iter) {
    const auto t0 = std::chrono::steady_clock::now();
    const vec_ZZ_pE packed_table =
        basefold::PackZ2kCoeffsToGREvals(params, t_table);
    const vec_ZZ_pE packed_monomial =
        basefold::BooleanHypercubeTableToMonomialCoeffs(packed_table);
    const auto t1 = std::chrono::steady_clock::now();

    const auto t2 = std::chrono::steady_clock::now();
    const basefold::MerkleRoot commitment =
        basefold::RingSwitchPCSCommit(params, t_table);
    const auto t3 = std::chrono::steady_clock::now();

    const basefold::MerkleRoot direct_commitment =
        basefold::Z2kPCSBackendCommit(params.backend, packed_monomial);
    if (commitment != direct_commitment) {
      LogicError("RunCommitBenchmark: compiler commitment mismatch");
    }

    if (packed_table.length() > 0) {
      anti_opt_checksum ^=
          static_cast<std::uint64_t>(packed_table[0] != ZZ_pE(0));
    }
    if (!commitment.empty()) {
      anti_opt_checksum ^= static_cast<std::uint64_t>(commitment[0]);
    }

    if (iter >= 0) {
      packing_ms.push_back(MsSince(t0, t1));
      commit_ms.push_back(MsSince(t2, t3));
    }
  }

  BenchResult out;
  out.packing = ComputeStats(packing_ms);
  out.commit = ComputeStats(commit_ms);
  out.anti_opt_checksum = anti_opt_checksum;
  return out;
}

void PrintResult(long c, long ell, long kappa, int warmup, int reps,
                 const BenchResult &result) {
  std::cout << "\n[ring-switch commit]"
            << " c=" << c << " ell=" << ell << " kappa=" << kappa
            << " ell'=" << (ell - kappa) << " warmup=" << warmup
            << " reps=" << reps << "\n";
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "  hash backend " << basefold::SelectedHashBackendName() << "\n";
  std::cout << "  packing mean " << result.packing.mean_ms << " ms  (min "
            << result.packing.min_ms << ", max " << result.packing.max_ms
            << ")\n";
  std::cout << "  commit  mean " << result.commit.mean_ms << " ms  (min "
            << result.commit.min_ms << ", max " << result.commit.max_ms
            << ")\n";
  std::cout << "  anti-opt checksum " << result.anti_opt_checksum << "\n";
}

void PrintHelp() {
  std::cout
      << "bench_z2k_ring_switch_commit\n\n"
      << "Usage:\n"
      << "  bench_z2k_ring_switch_commit [--c <int>] [--ell <int>] [--kappa <int>]\n"
      << "                               [--warmup <int>] [--reps <int>] [--seed <u64>]\n"
      << "                               [--auto-zeta teich]\n"
      << "                               [--ring-mod <decimal-int>] [--ring-p <decimal-int>]\n"
      << "                               [--ring-F <a0,a1,...>] [--ring-zeta <b0,b1,...>]\n\n"
      << "Notes:\n"
      << "  Headline packing time measures t -> t' packing plus Boolean-table to monomial conversion.\n"
      << "  Headline commit time measures the full RingSwitchPCSCommit path.\n";
  PrintCurrentBasisModeNotes(std::cout, "  ");
}

}  // namespace

int main(int argc, char **argv) {
  long c = 2;
  long ell = 3;
  long kappa = 1;
  int warmup = 1;
  int reps = 3;
  std::uint64_t seed = 0x5eed1234ULL;
  bool auto_zeta_teich = false;

  ContextSpec spec;
  spec.scalar_modulus = to_ZZ(4);
  spec.base_prime = to_ZZ(2);
  spec.F_coeffs = {to_ZZ(1), to_ZZ(1), to_ZZ(1)};
  spec.zeta_coeffs = {to_ZZ(0), to_ZZ(1)};

  for (int i = 1; i < argc; ++i) {
    const std::string arg(argv[i]);
    const auto need_value = [&](const char *flag) -> const char * {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << flag << "\n";
        std::exit(2);
      }
      return argv[++i];
    };

    if (arg == "--help" || arg == "-h") {
      PrintHelp();
      return 0;
    } else if (arg == "--c") {
      if (!ParseLong(need_value("--c"), c)) {
        std::cerr << "Invalid --c\n";
        return 2;
      }
    } else if (arg == "--ell") {
      if (!ParseLong(need_value("--ell"), ell)) {
        std::cerr << "Invalid --ell\n";
        return 2;
      }
    } else if (arg == "--kappa") {
      if (!ParseLong(need_value("--kappa"), kappa)) {
        std::cerr << "Invalid --kappa\n";
        return 2;
      }
    } else if (arg == "--warmup") {
      if (!ParseInt(need_value("--warmup"), warmup)) {
        std::cerr << "Invalid --warmup\n";
        return 2;
      }
    } else if (arg == "--reps") {
      if (!ParseInt(need_value("--reps"), reps)) {
        std::cerr << "Invalid --reps\n";
        return 2;
      }
    } else if (arg == "--seed") {
      seed = ParseU64OrDie(need_value("--seed"), "--seed");
    } else if (arg == "--auto-zeta") {
      const std::string mode = need_value("--auto-zeta");
      if (mode != "teich") {
        std::cerr << "Unsupported --auto-zeta mode\n";
        return 2;
      }
      auto_zeta_teich = true;
    } else if (arg == "--ring-mod") {
      if (!ParseZZ(need_value("--ring-mod"), spec.scalar_modulus)) {
        std::cerr << "Invalid --ring-mod\n";
        return 2;
      }
    } else if (arg == "--ring-p") {
      if (!ParseZZ(need_value("--ring-p"), spec.base_prime)) {
        std::cerr << "Invalid --ring-p\n";
        return 2;
      }
    } else if (arg == "--ring-F") {
      spec.F_coeffs = ParseCoeffList(need_value("--ring-F"));
    } else if (arg == "--ring-zeta") {
      spec.zeta_coeffs = ParseCoeffList(need_value("--ring-zeta"));
    } else {
      std::cerr << "Unknown argument: " << arg << "\n";
      return 2;
    }
  }

  try {
    if (ell < kappa) {
      LogicError("main: ell must be >= kappa");
    }
    if (c <= 0) {
      LogicError("main: c must be positive");
    }

    ZZ_pPush mod_push(spec.scalar_modulus);
    ValidateMonic(spec.F_coeffs, spec.scalar_modulus, "ring F");
    const ZZ_pX F = BuildZZpX(spec.F_coeffs);
    ZZ_pEPush ext_push(F);

    ZZ p_base;
    long k_base = 0;
    DeduceBasePrimeAndExponent(spec, p_base, k_base);
    ZZ_pE zeta;
    if (auto_zeta_teich) {
      zeta = FindTeichmullerGenerator(p_base, k_base, NTL::deg(F), F);
    } else {
      zeta = BuildZZpE(spec.zeta_coeffs);
    }

    const basefold::RingSwitchPCSParams params =
        BuildRingSwitchParams(c, ell, kappa, spec, F, zeta);
    const vec_ZZ_pE t_table =
        MakeDeterministicBaseRingTable(Pow2Checked(ell), seed);

    const BenchResult result =
        RunCommitBenchmark(params, t_table, warmup, reps);
    PrintResult(c, ell, kappa, warmup, reps, result);
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 2;
  }
}
