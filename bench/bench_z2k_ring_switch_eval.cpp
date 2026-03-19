#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "PCS/Common/Multilinear.hpp"
#include "PCS/Common/Profile.hpp"
#include "bench_z2k_ring_switch_common.hpp"

namespace {

using namespace basefold_bench_z2k_ring_switch_common;

struct BenchResult {
  Stats commit_outer;
  Stats commit_backend;
  Stats commit_total;
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

BenchResult RunEvalBenchmark(const basefold::RingSwitchPCSParams &params,
                             const vec_ZZ_pE &t_table,
                             const std::vector<ZZ_pE> &z,
                             const ZZ_pE &claimed_s, long num_queries,
                             bool use_checked_prover_path, int warmup,
                             int reps) {
  if (warmup < 0) {
    LogicError("RunEvalBenchmark: warmup must be >= 0");
  }
  if (reps <= 0) {
    LogicError("RunEvalBenchmark: reps must be > 0");
  }
  if (num_queries < 0) {
    LogicError("RunEvalBenchmark: num_queries must be >= 0");
  }

  std::vector<double> prove_total_ms;
  std::vector<double> prove_outer_ms;
  std::vector<double> prove_backend_ms;
  std::vector<double> commit_outer_ms;
  std::vector<double> commit_backend_ms;
  std::vector<double> commit_total_ms;
  std::vector<double> verify_total_ms;
  std::vector<double> verify_outer_ms;
  std::vector<double> verify_backend_ms;
  commit_outer_ms.reserve(static_cast<std::size_t>(reps));
  commit_backend_ms.reserve(static_cast<std::size_t>(reps));
  commit_total_ms.reserve(static_cast<std::size_t>(reps));
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
    const auto c0 = std::chrono::steady_clock::now();
    const basefold::RingSwitchPCSOuterCommitArtifacts outer_commit_artifacts =
        basefold::RingSwitchPCSBuildOuterCommitArtifacts(params, t_table);
    const auto c1 = std::chrono::steady_clock::now();

    const auto c2 = std::chrono::steady_clock::now();
    const basefold::Z2kPCSBackendCommitArtifacts backend_commit_artifacts =
        basefold::Z2kPCSBackendBuildCommitArtifacts(
            params.backend, outer_commit_artifacts.t_packed_monomial_coeffs);
    const auto c3 = std::chrono::steady_clock::now();

    basefold::RingSwitchPCSCommitArtifacts commit_artifacts;
    commit_artifacts.t_packed_table = outer_commit_artifacts.t_packed_table;
    commit_artifacts.t_packed_monomial_coeffs =
        outer_commit_artifacts.t_packed_monomial_coeffs;
    commit_artifacts.backend_commit_artifacts = backend_commit_artifacts;
    commit_artifacts.commitment = backend_commit_artifacts.commitment;

    basefold::Profile prover_prof;
    basefold::Profile verifier_prof;
    basefold::ResetProfile(prover_prof);
    basefold::ResetProfile(verifier_prof);

    const auto t0 = std::chrono::steady_clock::now();
    basefold::RingSwitchPCSEvalProof proof;
    {
      basefold::ProfileGuard guard(&prover_prof);
      if (use_checked_prover_path) {
        proof = basefold::RingSwitchPCSProveEvalFromCommitArtifacts(
            params, t_table, z, claimed_s, num_queries, commit_artifacts);
      } else {
        proof = basefold::RingSwitchPCSProveEvalFromCommitArtifactsUnchecked(
            params, t_table, z, claimed_s, num_queries, commit_artifacts);
      }
    }
    const auto t1 = std::chrono::steady_clock::now();

    outer_proof_size_bytes =
        basefold::RingSwitchPCSOuterProofSizeBytes(params, proof);
    proof_size_bytes = basefold::RingSwitchPCSEvalProofSizeBytes(params, proof);

    const auto t2 = std::chrono::steady_clock::now();
    bool ok = false;
    {
      basefold::ProfileGuard guard(&verifier_prof);
      ok = basefold::RingSwitchPCSVerifyEvalUnchecked(
          params, commit_artifacts.commitment, z, claimed_s, num_queries,
          proof);
    }
    const auto t3 = std::chrono::steady_clock::now();
    if (!ok) {
      LogicError("RunEvalBenchmark: verification failed");
    }

    const double prove_total = MsSince(t0, t1);
    const double prove_backend =
        basefold::NsToMs(prover_prof.z2k_backend_prove_ns);
    const double commit_outer = MsSince(c0, c1);
    const double commit_backend = MsSince(c2, c3);
    const double commit_total = commit_outer + commit_backend;
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
      commit_outer_ms.push_back(commit_outer);
      commit_backend_ms.push_back(commit_backend);
      commit_total_ms.push_back(commit_total);
      prove_total_ms.push_back(prove_total);
      prove_outer_ms.push_back(prove_outer);
      prove_backend_ms.push_back(prove_backend);
      verify_total_ms.push_back(verify_total);
      verify_outer_ms.push_back(verify_outer);
      verify_backend_ms.push_back(verify_backend);
    }
  }

  BenchResult out;
  out.commit_outer = ComputeStats(commit_outer_ms);
  out.commit_backend = ComputeStats(commit_backend_ms);
  out.commit_total = ComputeStats(commit_total_ms);
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
                 int reps, const RingSwitchBenchCliConfig &config,
                 const basefold::RingSwitchPCSParams &params,
                 const BenchResult &result) {
  std::cout << "\n[ring-switch eval]"
            << " c=" << c << " ell=" << ell << " kappa=" << kappa
            << " ell'=" << (ell - kappa) << " queries=" << queries
            << " warmup=" << warmup << " reps=" << reps << "\n";
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "  hash backend " << basefold::SelectedHashBackendName() << "\n";
  PrintBasisModeSummary(std::cout, config, params);
  std::cout << "  outer commit mean " << result.commit_outer.mean_ms
            << " ms  (min " << result.commit_outer.min_ms << ", max "
            << result.commit_outer.max_ms << ")\n";
  std::cout << "  backend commit mean " << result.commit_backend.mean_ms
            << " ms  (min " << result.commit_backend.min_ms << ", max "
            << result.commit_backend.max_ms << ")\n";
  std::cout << "  commit total mean " << result.commit_total.mean_ms
            << " ms  (min " << result.commit_total.min_ms << ", max "
            << result.commit_total.max_ms << ")\n";
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
      << "bench_z2k_ring_switch_eval\n\n"
      << "Usage:\n"
      << "  bench_z2k_ring_switch_eval [--c <int>] [--ell <int>] [--kappa <int>]\n"
      << "                             [--queries <int>] [--warmup <int>] [--reps <int>] [--checked]\n"
      << "                             [--seed <u64>] [--auto-zeta teich]\n"
      << "                             [--ring-mod <decimal-int>] [--ring-p <decimal-int>]\n"
      << "                             [--ring-F <a0,a1,...>] [--ring-zeta <b0,b1,...>]\n";
  PrintProvidedBasisFlagHelp(std::cout, "                             ");
  std::cout << "\nNotes:\n"
      << "  This bench measures one full compiler eval run and reports commit split as outer commit + backend commit.\n"
      << "  Timed prove defaults to the unchecked prover hot path; pass --checked to include prover-side input and honest-witness self-checks.\n"
      << "  Timed verify defaults to the unchecked verifier hot path with trusted params.\n"
      << "  Outer prover/verifier times are total times with the backend prove/verify subcall removed.\n"
      << "  Proof size reports exact serializer-backed bytes for the public RingSwitchPCSEvalProof.\n";
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
  bool use_checked_prover_path = false;

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
    } else if (arg == "--checked") {
      use_checked_prover_path = true;
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
    const std::vector<ZZ_pE> z =
        MakeDeterministicPoint(ell, seed ^ 0xabcddcbaULL);
    const ZZ_pE claimed_s = basefold::EvalMultilinearMonomialCoeffs(
        basefold::BooleanHypercubeTableToMonomialCoeffs(t_table), z);

    const BenchResult result =
        RunEvalBenchmark(params, t_table, z, claimed_s, queries,
                         use_checked_prover_path, warmup, reps);
    PrintResult(c, ell, kappa, queries, warmup, reps, config, params, result);
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 2;
  }
}
