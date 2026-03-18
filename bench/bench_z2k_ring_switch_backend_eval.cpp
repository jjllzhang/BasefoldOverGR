#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "PCS/Common/Multilinear.hpp"
#include "bench_z2k_ring_switch_common.hpp"

namespace {

using namespace basefold_bench_z2k_ring_switch_common;

struct BenchResult {
  Stats commit;
  Stats prove;
  Stats verify;
  std::uint64_t proof_size_bytes = 0;
  double proof_size_kb = 0.0;
  std::uint64_t anti_opt_checksum = 0;
};

BenchResult RunBackendEvalBenchmark(const basefold::Z2kPCSBackendHandle &backend,
                                    const vec_ZZ_pE &packed_monomial_coeffs,
                                    const std::vector<ZZ_pE> &z_backend,
                                    const ZZ_pE &claimed_y, long num_queries,
                                    int warmup, int reps) {
  if (warmup < 0) {
    LogicError("RunBackendEvalBenchmark: warmup must be >= 0");
  }
  if (reps <= 0) {
    LogicError("RunBackendEvalBenchmark: reps must be > 0");
  }
  if (num_queries < 0) {
    LogicError("RunBackendEvalBenchmark: num_queries must be >= 0");
  }

  std::vector<double> commit_ms;
  std::vector<double> prove_ms;
  std::vector<double> verify_ms;
  commit_ms.reserve(static_cast<std::size_t>(reps));
  prove_ms.reserve(static_cast<std::size_t>(reps));
  verify_ms.reserve(static_cast<std::size_t>(reps));

  std::uint64_t anti_opt_checksum = 0;
  std::uint64_t proof_size_bytes = 0;
  for (int iter = -warmup; iter < reps; ++iter) {
    const auto t0 = std::chrono::steady_clock::now();
    const basefold::Z2kPCSBackendCommitArtifacts commit_artifacts =
        basefold::Z2kPCSBackendBuildCommitArtifacts(backend, packed_monomial_coeffs);
    const auto t1 = std::chrono::steady_clock::now();

    const auto t2 = std::chrono::steady_clock::now();
    const basefold::Z2kPCSBackendEvalProof proof = basefold::Z2kPCSBackendProveEval(
        backend, packed_monomial_coeffs, z_backend, claimed_y, num_queries,
        &commit_artifacts);
    const auto t3 = std::chrono::steady_clock::now();

    proof_size_bytes = basefold::Z2kPCSBackendEvalProofSizeBytes(backend, proof);

    const auto t4 = std::chrono::steady_clock::now();
    const bool ok = basefold::Z2kPCSBackendVerifyEval(
        backend, commit_artifacts.commitment, z_backend, claimed_y, num_queries,
        proof);
    const auto t5 = std::chrono::steady_clock::now();
    if (!ok) {
      LogicError("RunBackendEvalBenchmark: verification failed");
    }

    anti_opt_checksum ^= proof_size_bytes;
    anti_opt_checksum ^= static_cast<std::uint64_t>(ok);
    if (!commit_artifacts.commitment.empty()) {
      anti_opt_checksum ^= static_cast<std::uint64_t>(commit_artifacts.commitment[0]);
    }

    if (iter >= 0) {
      commit_ms.push_back(MsSince(t0, t1));
      prove_ms.push_back(MsSince(t2, t3));
      verify_ms.push_back(MsSince(t4, t5));
    }
  }

  BenchResult out;
  out.commit = ComputeStats(commit_ms);
  out.prove = ComputeStats(prove_ms);
  out.verify = ComputeStats(verify_ms);
  out.proof_size_bytes = proof_size_bytes;
  out.proof_size_kb = static_cast<double>(proof_size_bytes) / 1024.0;
  out.anti_opt_checksum = anti_opt_checksum;
  return out;
}

void PrintResult(long c, long ell, long kappa, long queries, int warmup,
                 int reps, const RingSwitchBenchCliConfig &config,
                 const basefold::RingSwitchPCSParams &params,
                 const BenchResult &result) {
  std::cout << "\n[ring-switch backend eval]"
            << " c=" << c << " ell=" << ell << " kappa=" << kappa
            << " ell'=" << (ell - kappa) << " queries=" << queries
            << " warmup=" << warmup << " reps=" << reps << "\n";
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "  hash backend " << basefold::SelectedHashBackendName() << "\n";
  PrintBasisModeSummary(std::cout, config, params);
  std::cout << "  backend commit mean " << result.commit.mean_ms << " ms  (min "
            << result.commit.min_ms << ", max " << result.commit.max_ms
            << ")\n";
  std::cout << "  prove-phase mean " << result.prove.mean_ms << " ms  (min "
            << result.prove.min_ms << ", max " << result.prove.max_ms
            << ")\n";
  std::cout << "  verifier mean " << result.verify.mean_ms << " ms  (min "
            << result.verify.min_ms << ", max " << result.verify.max_ms
            << ")\n";
  std::cout << "  proof size  " << result.proof_size_kb << " KB  ("
            << result.proof_size_bytes << " B)\n";
  std::cout << "  anti-opt checksum " << result.anti_opt_checksum << "\n";
}

void PrintHelp() {
  std::cout
      << "bench_z2k_ring_switch_backend_eval\n\n"
      << "Usage:\n"
      << "  bench_z2k_ring_switch_backend_eval [--c <int>] [--ell <int>] [--kappa <int>]\n"
      << "                                     [--queries <int>] [--warmup <int>] [--reps <int>]\n"
      << "                                     [--seed <u64>] [--auto-zeta teich]\n"
      << "                                     [--ring-mod <decimal-int>] [--ring-p <decimal-int>]\n"
      << "                                     [--ring-F <a0,a1,...>] [--ring-zeta <b0,b1,...>]\n";
  PrintProvidedBasisFlagHelp(std::cout, "                                     ");
  std::cout << "\nNotes:\n"
      << "  This bench measures only the packed backend PCS over GR(2^k,2^kappa) with ell' = ell-kappa.\n"
      << "  Packing and compiler outer proof work are excluded from all timed regions.\n";
  PrintBasisCliNotes(std::cout, "  ");
}

}  // namespace

int main(int argc, char **argv) {
  long c = 4;
  long ell = 3;
  long kappa = 1;
  long queries = 2;
  int warmup = 1;
  int reps = 3;
  std::uint64_t seed = 0x5eed5678ULL;
  bool auto_zeta_teich = false;

  RingSwitchBenchCliArgs cli;
  cli.context.scalar_modulus = to_ZZ(4);
  cli.context.base_prime = to_ZZ(2);
  cli.context.F_coeffs = {to_ZZ(1), to_ZZ(1), to_ZZ(1)};
  cli.context.zeta_coeffs = {to_ZZ(0), to_ZZ(1)};

  for (int i = 1; i < argc; ++i) {
    const std::string arg(argv[i]);

    if (arg == "--help" || arg == "-h") {
      PrintHelp();
      return 0;
    } else if (arg == "--c") {
      if (!ParseLong(NeedValueOrExit(i, argc, argv, "--c"), c)) {
        std::cerr << "Invalid --c\n";
        return 2;
      }
    } else if (arg == "--ell") {
      if (!ParseLong(NeedValueOrExit(i, argc, argv, "--ell"), ell)) {
        std::cerr << "Invalid --ell\n";
        return 2;
      }
    } else if (arg == "--kappa") {
      if (!ParseLong(NeedValueOrExit(i, argc, argv, "--kappa"), kappa)) {
        std::cerr << "Invalid --kappa\n";
        return 2;
      }
    } else if (arg == "--queries") {
      if (!ParseLong(NeedValueOrExit(i, argc, argv, "--queries"), queries)) {
        std::cerr << "Invalid --queries\n";
        return 2;
      }
    } else if (arg == "--warmup") {
      if (!ParseInt(NeedValueOrExit(i, argc, argv, "--warmup"), warmup)) {
        std::cerr << "Invalid --warmup\n";
        return 2;
      }
    } else if (arg == "--reps") {
      if (!ParseInt(NeedValueOrExit(i, argc, argv, "--reps"), reps)) {
        std::cerr << "Invalid --reps\n";
        return 2;
      }
    } else if (arg == "--seed") {
      seed = ParseU64OrDie(NeedValueOrExit(i, argc, argv, "--seed"), "--seed");
    } else if (arg == "--auto-zeta") {
      const std::string mode = NeedValueOrExit(i, argc, argv, "--auto-zeta");
      if (mode != "teich") {
        std::cerr << "Unsupported --auto-zeta mode\n";
        return 2;
      }
      auto_zeta_teich = true;
    } else if (arg == "--ring-mod") {
      if (!ParseZZ(NeedValueOrExit(i, argc, argv, "--ring-mod"),
                   cli.context.scalar_modulus)) {
        std::cerr << "Invalid --ring-mod\n";
        return 2;
      }
    } else if (arg == "--ring-p") {
      if (!ParseZZ(NeedValueOrExit(i, argc, argv, "--ring-p"),
                   cli.context.base_prime)) {
        std::cerr << "Invalid --ring-p\n";
        return 2;
      }
    } else if (arg == "--ring-F") {
      cli.context.F_coeffs = ParseCoeffList(
          NeedValueOrExit(i, argc, argv, "--ring-F"));
    } else if (arg == "--ring-zeta") {
      cli.context.zeta_coeffs = ParseCoeffList(
          NeedValueOrExit(i, argc, argv, "--ring-zeta"));
    } else if (TryParseBasisCliArg(arg, i, argc, argv, cli)) {
      continue;
    } else {
      std::cerr << "Unknown argument: " << arg << "\n";
      return 2;
    }
  }

  try {
    if (ell < kappa) {
      LogicError("main: ell must be >= kappa");
    }
    if (queries < 0) {
      LogicError("main: queries must be non-negative");
    }
    if (c <= 0) {
      LogicError("main: c must be positive");
    }

    ZZ_pPush mod_push(cli.context.scalar_modulus);
    ValidateMonic(cli.context.F_coeffs, cli.context.scalar_modulus, "ring F");
    const ZZ_pX F = BuildZZpX(cli.context.F_coeffs);
    ZZ_pEPush ext_push(F);

    const RingSwitchBenchCliConfig config =
        DecodeRingSwitchBenchCliConfigOrThrow(
            cli, F, BenchDefaultBasisPreset::kEvalLike);
    ZZ p_base;
    long k_base = 0;
    DeduceBasePrimeAndExponent(cli.context, p_base, k_base);
    ZZ_pE zeta;
    if (auto_zeta_teich) {
      zeta = FindTeichmullerGenerator(p_base, k_base, NTL::deg(F), F);
    } else {
      zeta = BuildZZpE(cli.context.zeta_coeffs);
    }

    const basefold::RingSwitchPCSParams params =
        BuildRingSwitchParams(c, ell, kappa, config, F, zeta);
    const vec_ZZ_pE t_table =
        MakeDeterministicBaseRingTable(Pow2Checked(ell), seed);
    const basefold::RingSwitchPCSOuterCommitArtifacts outer_commit_artifacts =
        basefold::RingSwitchPCSBuildOuterCommitArtifacts(params, t_table);
    const std::vector<ZZ_pE> z_backend =
        MakeDeterministicPoint(params.ell_prime, seed ^ 0xabcddcbaULL);
    const ZZ_pE claimed_y = basefold::EvalMultilinearMonomialCoeffs(
        outer_commit_artifacts.t_packed_monomial_coeffs, z_backend);

    const BenchResult result = RunBackendEvalBenchmark(
        params.backend, outer_commit_artifacts.t_packed_monomial_coeffs,
        z_backend, claimed_y, queries, warmup, reps);
    PrintResult(c, ell, kappa, queries, warmup, reps, config, params, result);
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 2;
  }
}
