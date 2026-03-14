#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "PCS/Common/Profile.hpp"
#include "Compiler/Z2k/FrobeniusProofSerialize.hpp"
#include "bench_z2k_frobenius_common.hpp"

using NTL::SetCoeff;
using NTL::ZZ;
using NTL::ZZ_pE;
using NTL::ZZ_pEPush;
using NTL::ZZ_pPush;
using NTL::ZZ_pX;
using NTL::conv;
using NTL::to_ZZ;
using NTL::vec_ZZ_pE;

namespace {

using basefold_bench_z2k_frobenius_common::BuildFrobeniusBenchSetupResult;
using basefold_bench_z2k_frobenius_common::BuildZZpE;
using basefold_bench_z2k_frobenius_common::BuildZZpX;
using basefold_bench_z2k_frobenius_common::ComputeStats;
using basefold_bench_z2k_frobenius_common::ContextSpec;
using basefold_bench_z2k_frobenius_common::EvalFromBooleanTable;
using basefold_bench_z2k_frobenius_common::FrobeniusBenchCalibrationMode;
using basefold_bench_z2k_frobenius_common::MakeDeterministicBaseRingTable;
using basefold_bench_z2k_frobenius_common::MakeDeterministicPoint;
using basefold_bench_z2k_frobenius_common::MsSince;
using basefold_bench_z2k_frobenius_common::ParseCoeffList;
using basefold_bench_z2k_frobenius_common::ParseInt;
using basefold_bench_z2k_frobenius_common::ParseLong;
using basefold_bench_z2k_frobenius_common::ParseU64OrDie;
using basefold_bench_z2k_frobenius_common::ParseZZ;
using basefold_bench_z2k_frobenius_common::PrintFrobeniusBasisSummary;
using basefold_bench_z2k_frobenius_common::Stats;
using basefold_bench_z2k_frobenius_common::ValidateMonic;

struct BenchResult {
  Stats prove_total;
  Stats prove_outer;
  Stats prove_backend;
  Stats verify_total;
  Stats verify_outer;
  Stats verify_backend;
  std::uint64_t outer_proof_size_bytes = 0;
  std::uint64_t proof_size_bytes = 0;
  double outer_proof_size_kb = 0.0;
  double proof_size_kb = 0.0;
  std::uint64_t anti_opt_checksum = 0;
};

BenchResult RunEvalBenchmark(const basefold::FrobeniusPCSParams &params,
                             const vec_ZZ_pE &t_table,
                             const std::vector<ZZ_pE> &z,
                             const ZZ_pE &claimed_s, long num_queries,
                             int warmup, int reps) {
  if (warmup < 0) {
    NTL::LogicError("RunEvalBenchmark: warmup must be >= 0");
  }
  if (reps <= 0) {
    NTL::LogicError("RunEvalBenchmark: reps must be > 0");
  }
  if (num_queries < 0) {
    NTL::LogicError("RunEvalBenchmark: num_queries must be >= 0");
  }

  std::vector<double> prove_total_ms;
  std::vector<double> prove_outer_ms;
  std::vector<double> prove_backend_ms;
  std::vector<double> verify_total_ms;
  std::vector<double> verify_outer_ms;
  std::vector<double> verify_backend_ms;
  prove_total_ms.reserve(static_cast<std::size_t>(reps));
  prove_outer_ms.reserve(static_cast<std::size_t>(reps));
  prove_backend_ms.reserve(static_cast<std::size_t>(reps));
  verify_total_ms.reserve(static_cast<std::size_t>(reps));
  verify_outer_ms.reserve(static_cast<std::size_t>(reps));
  verify_backend_ms.reserve(static_cast<std::size_t>(reps));

  std::uint64_t anti_opt_checksum = 0;
  std::uint64_t outer_proof_size_bytes = 0;
  std::uint64_t proof_size_bytes = 0;

  for (int iter = -warmup; iter < reps; ++iter) {
    const basefold::FrobeniusPCSCommitArtifacts commit_artifacts =
        basefold::FrobeniusPCSBuildCommitArtifacts(params, t_table);

    basefold::Profile prover_prof;
    basefold::Profile verifier_prof;
    basefold::ResetProfile(prover_prof);
    basefold::ResetProfile(verifier_prof);

    const auto t0 = std::chrono::steady_clock::now();
    basefold::FrobeniusPCSEvalProof proof;
    {
      basefold::ProfileGuard guard(&prover_prof);
      proof = basefold::FrobeniusPCSProveEvalFromCommitArtifacts(
          params, t_table, z, claimed_s, num_queries, commit_artifacts);
    }
    const auto t1 = std::chrono::steady_clock::now();

    outer_proof_size_bytes =
        basefold::FrobeniusPCSOuterProofSizeBytes(params, proof);
    proof_size_bytes = basefold::FrobeniusPCSEvalProofSizeBytes(params, proof);

    const auto t2 = std::chrono::steady_clock::now();
    bool ok = false;
    {
      basefold::ProfileGuard guard(&verifier_prof);
      ok = basefold::FrobeniusPCSVerifyEval(params, commit_artifacts.commitment,
                                            z, claimed_s, num_queries, proof);
    }
    const auto t3 = std::chrono::steady_clock::now();
    if (!ok) {
      NTL::LogicError("RunEvalBenchmark: verification failed");
    }

    const double prove_total = MsSince(t0, t1);
    const double prove_backend =
        basefold::NsToMs(prover_prof.z2k_backend_prove_ns);
    const double verify_total = MsSince(t2, t3);
    const double verify_backend =
        basefold::NsToMs(verifier_prof.z2k_backend_verify_ns);
    const double prove_outer = std::max(0.0, prove_total - prove_backend);
    const double verify_outer = std::max(0.0, verify_total - verify_backend);

    anti_opt_checksum ^= static_cast<std::uint64_t>(ok);
    anti_opt_checksum ^= proof_size_bytes;
    if (!commit_artifacts.commitment.empty()) {
      anti_opt_checksum ^=
          static_cast<std::uint64_t>(commit_artifacts.commitment[0]);
    }

    if (iter >= 0) {
      prove_total_ms.push_back(prove_total);
      prove_outer_ms.push_back(prove_outer);
      prove_backend_ms.push_back(prove_backend);
      verify_total_ms.push_back(verify_total);
      verify_outer_ms.push_back(verify_outer);
      verify_backend_ms.push_back(verify_backend);
    }
  }

  BenchResult out;
  out.prove_total = ComputeStats(prove_total_ms);
  out.prove_outer = ComputeStats(prove_outer_ms);
  out.prove_backend = ComputeStats(prove_backend_ms);
  out.verify_total = ComputeStats(verify_total_ms);
  out.verify_outer = ComputeStats(verify_outer_ms);
  out.verify_backend = ComputeStats(verify_backend_ms);
  out.outer_proof_size_bytes = outer_proof_size_bytes;
  out.proof_size_bytes = proof_size_bytes;
  out.outer_proof_size_kb =
      static_cast<double>(outer_proof_size_bytes) / 1024.0;
  out.proof_size_kb = static_cast<double>(proof_size_bytes) / 1024.0;
  out.anti_opt_checksum = anti_opt_checksum;
  return out;
}

void PrintResult(long c, long ell, long kappa, long queries, int warmup,
                 int reps,
                 const basefold_bench_z2k_frobenius_common::FrobeniusBenchBasisSummary
                     &basis_summary,
                 const BenchResult &result) {
  std::cout << "\n[frobenius eval]"
            << " c=" << c << " ell=" << ell << " kappa=" << kappa
            << " ell'=" << (ell - kappa) << " queries=" << queries
            << " warmup=" << warmup << " reps=" << reps << "\n";
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "  hash backend " << basefold::SelectedHashBackendName() << "\n";
  PrintFrobeniusBasisSummary(std::cout, basis_summary);
  std::cout << "  prove-phase mean " << result.prove_total.mean_ms << " ms  (min "
            << result.prove_total.min_ms << ", max "
            << result.prove_total.max_ms << ")\n";
  std::cout << "  outer prover mean " << result.prove_outer.mean_ms << " ms  (min "
            << result.prove_outer.min_ms << ", max "
            << result.prove_outer.max_ms << ")\n";
  std::cout << "  backend prover mean " << result.prove_backend.mean_ms
            << " ms  (min " << result.prove_backend.min_ms << ", max "
            << result.prove_backend.max_ms << ")\n";
  std::cout << "  verifier mean " << result.verify_total.mean_ms << " ms  (min "
            << result.verify_total.min_ms << ", max "
            << result.verify_total.max_ms << ")\n";
  std::cout << "  outer verifier mean " << result.verify_outer.mean_ms
            << " ms  (min " << result.verify_outer.min_ms << ", max "
            << result.verify_outer.max_ms << ")\n";
  std::cout << "  backend verifier mean " << result.verify_backend.mean_ms
            << " ms  (min " << result.verify_backend.min_ms << ", max "
            << result.verify_backend.max_ms << ")\n";
  std::cout << "  outer proof size " << result.outer_proof_size_kb << " KB  ("
            << result.outer_proof_size_bytes << " B)\n";
  std::cout << "  proof size  " << result.proof_size_kb << " KB  ("
            << result.proof_size_bytes << " B)\n";
  std::cout << "  anti-opt checksum " << result.anti_opt_checksum << "\n";
}

void PrintHelp() {
  std::cout
      << "bench_z2k_frobenius_eval\n\n"
      << "Usage:\n"
      << "  bench_z2k_frobenius_eval [--c <int>] [--ell <int>] [--kappa <int>]\n"
      << "                           [--queries <int>] [--warmup <int>] [--reps <int>]\n"
      << "                           [--seed <u64>] [--auto-zeta teich]\n"
      << "                           [--ring-mod <decimal-int>] [--ring-p <decimal-int>]\n"
      << "                           [--ring-F <a0,a1,...>] [--ring-zeta <b0,b1,...>]\n\n"
      << "Notes:\n"
      << "  Commit/artifact construction happens before timed prove, matching current bench_basefold_pcs_eval semantics.\n"
      << "  Outer prover/verifier times are total times with the backend prove/verify subcall removed.\n"
      << "  Proof size reports exact serializer-backed bytes for the public FrobeniusPCSEvalProof.\n";
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

    const auto setup = BuildFrobeniusBenchSetupResult(
        FrobeniusBenchCalibrationMode::kEval, c, ell, kappa, spec, F, zeta);
    const basefold::FrobeniusPCSParams &params = setup.params;
    const vec_ZZ_pE t_table =
        MakeDeterministicBaseRingTable(1L << ell, seed);
    const std::vector<ZZ_pE> z =
        MakeDeterministicPoint(ell, seed ^ 0xabcddcbaULL);
    const ZZ_pE claimed_s = EvalFromBooleanTable(t_table, ell, z);

    const BenchResult result =
        RunEvalBenchmark(params, t_table, z, claimed_s, queries, warmup, reps);
    PrintResult(c, ell, kappa, queries, warmup, reps, setup.basis_summary,
                result);
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  } catch (...) {
    std::cerr << "Error: unknown exception\n";
    return 1;
  }

  return 0;
}
