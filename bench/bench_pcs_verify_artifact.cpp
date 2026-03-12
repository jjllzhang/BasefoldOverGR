#include <iomanip>
#include <iostream>
#include <string>

#include "bench_pcs_artifact_common.hpp"

using namespace basefold_bench_pcs_artifact;
namespace fs = std::filesystem;

namespace {

void PrintHelp() {
  std::cout
      << "bench_pcs_verify_artifact (PCS verify benchmark from one dumped artifact case)\n\n"
      << "Usage:\n"
      << "  bench_pcs_verify_artifact --artifact-root <dir> --artifact-id <id>\n"
      << "                            [--warmup <int>] [--reps <int>] [--profile]\n"
      << "                            [--verifier-query-per-thread <int>]\n"
      << "                            [--verifier-query-threshold <int>]\n"
      << "                            [--verifier-query-max-threads <int>]\n\n"
      << "Notes:\n"
      << "  File IO and proof/public-input deserialization are excluded from verifier mean.\n"
      << "  Merkle build thread flags are intentionally unsupported here because artifact verify\n"
      << "  does not build Merkle trees inside the measured path.\n";
}

void PrintResult(const LoadedArtifactCase &artifact,
                 const ArtifactVerifyBenchResult &result, int warmup, int reps) {
  std::cout << "\n[artifact verify]"
            << " artifact_id=" << artifact.metadata.artifact_id
            << " mode=" << artifact.metadata.mode << " c=" << artifact.metadata.c
            << " k0=" << artifact.metadata.k0 << " d=" << artifact.metadata.d
            << "  mod=" << artifact.metadata.scalar_modulus
            << "  poly_dim=" << artifact.metadata.poly_dim
            << "  queries=" << artifact.metadata.queries
            << "  warmup=" << warmup << " reps=" << reps;
  if (artifact.metadata.use_extension_challenges) {
    std::cout << "  ext_challenges=on"
              << "  ext_deg=" << artifact.challenge_ext_degree;
  } else {
    std::cout << "  ext_challenges=off";
  }
  std::cout << "\n";
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "  display_key " << artifact.metadata.display_key << "\n";
  std::cout << "  hash backend " << artifact.metadata.hash_backend << "\n";
  std::cout << "  artifact load wall time " << artifact.load_wall_ms
            << " ms  (excluded from verifier mean)\n";
  std::cout << "  artifact deserialize wall time "
            << artifact.deserialize_wall_ms
            << " ms  (excluded from verifier mean)\n";
  std::cout << "  verifier mean " << result.verifier.mean_ms << " ms  (min "
            << result.verifier.min_ms << ", max " << result.verifier.max_ms
            << ")\n";
  std::cout << "  input proof size " << result.proof_size_kb << " KB  ("
            << result.proof_size_bytes << " B)\n";
  std::cout << "  anti-opt checksum " << result.anti_opt_checksum << "\n";
  if (result.has_profile) {
    basefold::PrintProfile(std::cout, result.verifier_profile);
  }
}

}  // namespace

int main(int argc, char **argv) {
  fs::path artifact_root;
  std::string artifact_id;
  bool enable_profile = false;
  int warmup = 1;
  int reps = 3;

  basefold::ResetVerifierQueryParallelConfigFromEnv();
  basefold::VerifierQueryParallelConfig verifier_query_cfg =
      basefold::GetVerifierQueryParallelConfig();

  for (int i = 1; i < argc; ++i) {
    const std::string arg(argv[i]);
    auto RequireNextArgValue = [&](const char *flag) -> const char * {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << flag << "\n";
        std::exit(2);
      }
      return argv[++i];
    };

    if (arg == "--artifact-root") {
      artifact_root = RequireNextArgValue("--artifact-root");
    } else if (arg == "--artifact-id") {
      artifact_id = RequireNextArgValue("--artifact-id");
    } else if (arg == "--profile") {
      enable_profile = true;
    } else if (arg == "--warmup") {
      if (!basefold_bench_pcs_common::ParseInt(
              RequireNextArgValue("--warmup"), warmup) ||
          warmup < 0) {
        std::cerr << "Invalid --warmup\n";
        return 2;
      }
    } else if (arg == "--reps") {
      if (!basefold_bench_pcs_common::ParseInt(
              RequireNextArgValue("--reps"), reps) ||
          reps <= 0) {
        std::cerr << "Invalid --reps\n";
        return 2;
      }
    } else if (arg == "--verifier-query-per-thread") {
      if (!basefold_bench_pcs_common::ParseLong(
              RequireNextArgValue("--verifier-query-per-thread"),
              verifier_query_cfg.queries_per_thread) ||
          verifier_query_cfg.queries_per_thread <= 0) {
        std::cerr << "Invalid --verifier-query-per-thread\n";
        return 2;
      }
    } else if (arg == "--verifier-query-threshold") {
      if (!basefold_bench_pcs_common::ParseLong(
              RequireNextArgValue("--verifier-query-threshold"),
              verifier_query_cfg.min_queries_for_parallelism) ||
          verifier_query_cfg.min_queries_for_parallelism <= 0) {
        std::cerr << "Invalid --verifier-query-threshold\n";
        return 2;
      }
    } else if (arg == "--verifier-query-max-threads") {
      if (!basefold_bench_pcs_common::ParseInt(
              RequireNextArgValue("--verifier-query-max-threads"),
              verifier_query_cfg.max_threads) ||
          verifier_query_cfg.max_threads <= 0) {
        std::cerr << "Invalid --verifier-query-max-threads\n";
        return 2;
      }
    } else if (arg == "--help" || arg == "-h") {
      PrintHelp();
      return 0;
    } else if (arg == "--merkle-leaves-per-thread" ||
               arg == "--merkle-level-threshold" ||
               arg == "--merkle-max-threads") {
      std::cerr << "Unsupported " << arg
                << " (artifact verify does not build Merkle trees)\n";
      return 2;
    } else {
      std::cerr << "Unknown argument: " << arg << "\n";
      PrintHelp();
      return 2;
    }
  }

  if (artifact_root.empty()) {
    std::cerr << "Missing required --artifact-root\n";
    return 2;
  }
  if (artifact_id.empty()) {
    std::cerr << "Missing required --artifact-id\n";
    return 2;
  }

  try {
    basefold::SetVerifierQueryParallelConfig(verifier_query_cfg);
    const LoadedArtifactCase artifact =
        LoadArtifactCaseForVerify(artifact_root, artifact_id);
    const ArtifactVerifyBenchResult result =
        RunArtifactVerifyBenchmark(artifact, enable_profile, warmup, reps);
    PrintResult(artifact, result, warmup, reps);
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 2;
  }

  return 0;
}
