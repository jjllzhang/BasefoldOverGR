#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "Compiler/Z2k/FrobeniusProofSerialize.hpp"
#include "bench_z2k_frobenius_common.hpp"

using NTL::ZZ;
using NTL::ZZ_pE;
using NTL::ZZ_pEPush;
using NTL::ZZ_pPush;
using NTL::ZZ_pX;
using NTL::to_ZZ;
using NTL::vec_ZZ_pE;

namespace {

using basefold_bench_z2k_frobenius_common::BuildFrobeniusParams;
using basefold_bench_z2k_frobenius_common::BuildZZpE;
using basefold_bench_z2k_frobenius_common::BuildZZpX;
using basefold_bench_z2k_frobenius_common::ComputeStats;
using basefold_bench_z2k_frobenius_common::ContextSpec;
using basefold_bench_z2k_frobenius_common::EvalFromBooleanTable;
using basefold_bench_z2k_frobenius_common::MakeDeterministicBaseRingTable;
using basefold_bench_z2k_frobenius_common::MakeDeterministicPoint;
using basefold_bench_z2k_frobenius_common::MsSince;
using basefold_bench_z2k_frobenius_common::ParseCoeffList;
using basefold_bench_z2k_frobenius_common::ParseInt;
using basefold_bench_z2k_frobenius_common::ParseLong;
using basefold_bench_z2k_frobenius_common::ParseU64OrDie;
using basefold_bench_z2k_frobenius_common::ParseZZ;
using basefold_bench_z2k_frobenius_common::Pow2Checked;
using basefold_bench_z2k_frobenius_common::Stats;
using basefold_bench_z2k_frobenius_common::ValidateMonic;

struct BenchResult {
  Stats verify;
  std::uint64_t proof_size_bytes = 0;
  double proof_size_kb = 0.0;
  std::uint64_t anti_opt_checksum = 0;
};

BenchResult RunVerifyBenchmark(const basefold::FrobeniusPCSParams &params,
                               const vec_ZZ_pE &t_table,
                               const std::vector<ZZ_pE> &z,
                               const ZZ_pE &claimed_s, long num_queries,
                               int warmup, int reps) {
  if (warmup < 0) {
    NTL::LogicError("RunVerifyBenchmark: warmup must be >= 0");
  }
  if (reps <= 0) {
    NTL::LogicError("RunVerifyBenchmark: reps must be > 0");
  }
  if (num_queries < 0) {
    NTL::LogicError("RunVerifyBenchmark: num_queries must be >= 0");
  }

  const basefold::FrobeniusPCSOuterCommitArtifacts outer_commit_artifacts =
      basefold::FrobeniusPCSBuildOuterCommitArtifacts(params, t_table);
  const basefold::MerkleRoot commitment = basefold::Z2kPCSBackendCommit(
      params.backend, outer_commit_artifacts.t_packed_monomial_coeffs);
  const basefold::FrobeniusPCSOuterEvalProof proof =
      basefold::FrobeniusPCSProveOuterEvalFromCommitArtifacts(
          params, t_table, commitment, z, claimed_s, num_queries,
          outer_commit_artifacts);
  const std::uint64_t proof_size_bytes =
      basefold::FrobeniusPCSOuterProofSizeBytes(params, proof);

  std::vector<double> verify_ms;
  verify_ms.reserve(static_cast<std::size_t>(reps));
  std::uint64_t anti_opt_checksum = proof_size_bytes;

  for (int iter = -warmup; iter < reps; ++iter) {
    const auto t0 = std::chrono::steady_clock::now();
    const bool ok = basefold::FrobeniusPCSVerifyOuterEval(
        params, commitment, z, claimed_s, num_queries, proof);
    const auto t1 = std::chrono::steady_clock::now();

    if (!ok) {
      NTL::LogicError("RunVerifyBenchmark: verification failed");
    }

    anti_opt_checksum ^= static_cast<std::uint64_t>(ok);
    if (!commitment.empty()) {
      anti_opt_checksum ^= static_cast<std::uint64_t>(commitment[0]);
    }

    if (iter >= 0) {
      verify_ms.push_back(MsSince(t0, t1));
    }
  }

  BenchResult out;
  out.verify = ComputeStats(verify_ms);
  out.proof_size_bytes = proof_size_bytes;
  out.proof_size_kb = static_cast<double>(proof_size_bytes) / 1024.0;
  out.anti_opt_checksum = anti_opt_checksum;
  return out;
}

void PrintResult(long c, long ell, long kappa, long queries, int warmup,
                 int reps, const BenchResult &result) {
  std::cout << "\n[frobenius outer verify]"
            << " c=" << c << " ell=" << ell << " kappa=" << kappa
            << " ell'=" << (ell - kappa) << " queries=" << queries
            << " warmup=" << warmup << " reps=" << reps << "\n";
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "  hash backend " << basefold::SelectedHashBackendName() << "\n";
  std::cout << "  verifier mean " << result.verify.mean_ms << " ms  (min "
            << result.verify.min_ms << ", max " << result.verify.max_ms
            << ")\n";
  std::cout << "  input proof size " << result.proof_size_kb << " KB  ("
            << result.proof_size_bytes << " B)\n";
  std::cout << "  anti-opt checksum " << result.anti_opt_checksum << "\n";
}

void PrintHelp() {
  std::cout
      << "bench_z2k_frobenius_outer_verify\n\n"
      << "Usage:\n"
      << "  bench_z2k_frobenius_outer_verify [--c <int>] [--ell <int>] [--kappa <int>]\n"
      << "                                   [--queries <int>] [--warmup <int>] [--reps <int>]\n"
      << "                                   [--seed <u64>] [--auto-zeta teich]\n"
      << "                                   [--ring-mod <decimal-int>] [--ring-p <decimal-int>]\n"
      << "                                   [--ring-F <a0,a1,...>] [--ring-zeta <b0,b1,...>]\n\n"
      << "Notes:\n"
      << "  Commitment and outer proof are prebuilt outside timed verify.\n"
      << "  Timed verify only measures the compiler-side outer verification.\n";
}

}  // namespace

int main(int argc, char **argv) {
  long c = 2;
  long ell = 3;
  long kappa = 1;
  long queries = 2;
  int warmup = 1;
  int reps = 3;
  std::uint64_t seed = 0x5eed5678ULL;
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
    } else if (arg == "--queries") {
      if (!ParseLong(need_value("--queries"), queries)) {
        std::cerr << "Invalid --queries\n";
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

  if (c <= 0 || ell <= 0 || kappa < 0 || kappa > ell || queries < 0 ||
      warmup < 0 || reps <= 0) {
    std::cerr << "Invalid benchmark parameters\n";
    return 2;
  }
  ValidateMonic(spec.F_coeffs, spec.scalar_modulus, "--ring-F");

  try {
    ZZ p_base = spec.base_prime;
    if (p_base <= 1) {
      p_base = spec.scalar_modulus;
    }

    ZZ_pPush mod_push(spec.scalar_modulus);
    const ZZ_pX F = BuildZZpX(spec.F_coeffs);
    ZZ_pEPush ext_push(F);

    ZZ_pE zeta;
    if (auto_zeta_teich) {
      long k_base = 0;
      basefold_bench_z2k_frobenius_common::DeduceBasePrimeAndExponent(
          spec, p_base, k_base);
      zeta = FindTeichmullerGenerator(p_base, k_base, NTL::deg(F), F);
    } else {
      zeta = BuildZZpE(spec.zeta_coeffs);
    }

    const basefold::FrobeniusPCSParams params =
        BuildFrobeniusParams(c, ell, kappa, spec, F, zeta);
    const vec_ZZ_pE t_table =
        MakeDeterministicBaseRingTable(Pow2Checked(ell), seed);
    const std::vector<ZZ_pE> z =
        MakeDeterministicPoint(ell, seed ^ 0xabcddcbaULL);
    const ZZ_pE claimed_s = EvalFromBooleanTable(t_table, ell, z);

    const BenchResult result = RunVerifyBenchmark(params, t_table, z, claimed_s,
                                                  queries, warmup, reps);
    PrintResult(c, ell, kappa, queries, warmup, reps, result);
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  } catch (...) {
    std::cerr << "Error: unknown exception\n";
    return 1;
  }

  return 0;
}
