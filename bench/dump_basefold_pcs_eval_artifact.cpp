#include <iostream>
#include <string>

#include "bench_basefold_pcs_artifact_common.hpp"

using namespace basefold_bench_pcs_common;
using namespace basefold_bench_pcs_artifact;
namespace fs = std::filesystem;

namespace {

void PrintHelp() {
  std::cout
      << "dump_basefold_pcs_eval_artifact (persist one BaseFold eval-proof artifact case)\n\n"
      << "Usage:\n"
      << "  dump_basefold_pcs_eval_artifact --artifact-root <dir> [--artifact-id <id>]\n"
      << "                         [--context-id <id>] [--context-label <label>]\n"
      << "                         [--lambda <text>] [--gamma <text>]\n"
      << "                         --mode field|ring\n"
      << "                         [--c <int>] [--k0 <int>] [--d <int>] [--queries <int>]\n"
      << "                         [--checked] [--seed <u64>]\n"
      << "                         [--use-extension-challenges]\n"
      << "                         [--field-mod <decimal-int>] [--field-F <a0,a1,...>] [--field-zeta <b0,b1,...>]\n"
      << "                         [--ring-mod <decimal-int>]  [--ring-p <decimal-int>] [--ring-F <a0,a1,...>] [--ring-zeta <b0,b1,...>]\n"
      << "                         [--field-challenge-ext <a0;a1;...>] [--ring-challenge-ext <a0;a1;...>]\n"
      << "                         [--auto-zeta teich]\n\n"
      << "Notes:\n"
      << "  Exactly one case is dumped per invocation.\n"
      << "  Existing artifacts are never overwritten.\n"
      << "  Artifact files are written under manifest.jsonl + objects/<artifact_id>/.\n";
}

}  // namespace

int main(int argc, char **argv) {
  DumpArtifactRequest request;
  fs::path artifact_root;

  ContextSpec field;
  field.label = "Field";
  field.scalar_modulus = to_ZZ(2);
  field.base_prime = ZZ(0);
  field.F_coeffs = {to_ZZ(1), to_ZZ(1), to_ZZ(1)};
  field.zeta_coeffs = {to_ZZ(0), to_ZZ(1)};

  ContextSpec ring;
  ring.label = "Ring";
  ring.scalar_modulus = to_ZZ(4);
  ring.base_prime = to_ZZ(2);
  ring.F_coeffs = {to_ZZ(1), to_ZZ(1), to_ZZ(1)};
  ring.zeta_coeffs = {to_ZZ(0), to_ZZ(1)};

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
      request.artifact_id = RequireNextArgValue("--artifact-id");
    } else if (arg == "--context-id") {
      request.context_id = RequireNextArgValue("--context-id");
    } else if (arg == "--context-label") {
      request.context_label = RequireNextArgValue("--context-label");
    } else if (arg == "--lambda") {
      request.lambda = RequireNextArgValue("--lambda");
    } else if (arg == "--gamma") {
      request.gamma = RequireNextArgValue("--gamma");
    } else if (arg == "--mode") {
      request.mode = RequireNextArgValue("--mode");
      if (request.mode == "both") {
        std::cerr << "Invalid --mode (dump tool supports field|ring only)\n";
        return 2;
      }
      if (request.mode != "field" && request.mode != "ring") {
        std::cerr << "Invalid --mode (expected field|ring)\n";
        return 2;
      }
    } else if (arg == "--c") {
      if (!ParseLong(RequireNextArgValue("--c"), request.c) || request.c <= 0) {
        std::cerr << "Invalid --c\n";
        return 2;
      }
    } else if (arg == "--k0") {
      if (!ParseLong(RequireNextArgValue("--k0"), request.k0) ||
          request.k0 <= 0) {
        std::cerr << "Invalid --k0\n";
        return 2;
      }
    } else if (arg == "--d") {
      if (!ParseLong(RequireNextArgValue("--d"), request.d) || request.d < 0) {
        std::cerr << "Invalid --d\n";
        return 2;
      }
    } else if (arg == "--queries") {
      if (!ParseLong(RequireNextArgValue("--queries"), request.queries) ||
          request.queries < 0) {
        std::cerr << "Invalid --queries\n";
        return 2;
      }
    } else if (arg == "--checked") {
      request.use_checked_prover_path = true;
    } else if (arg == "--seed") {
      request.seed = ParseU64OrDie(RequireNextArgValue("--seed"), "--seed");
    } else if (arg == "--use-extension-challenges") {
      request.use_extension_challenges = true;
    } else if (arg == "--auto-zeta") {
      const std::string mode = RequireNextArgValue("--auto-zeta");
      if (mode != "teich") {
        std::cerr << "Invalid --auto-zeta (expected teich)\n";
        return 2;
      }
      request.auto_zeta_teich = true;
    } else if (arg == "--field-mod") {
      if (!ParseZZ(RequireNextArgValue("--field-mod"), field.scalar_modulus) ||
          field.scalar_modulus <= 1) {
        std::cerr << "Invalid --field-mod\n";
        return 2;
      }
    } else if (arg == "--field-F") {
      field.F_coeffs = ParseCoeffList(RequireNextArgValue("--field-F"));
    } else if (arg == "--field-zeta") {
      field.zeta_coeffs = ParseCoeffList(RequireNextArgValue("--field-zeta"));
    } else if (arg == "--field-challenge-ext") {
      field.challenge_ext_coeffs =
          ParseNestedCoeffList(RequireNextArgValue("--field-challenge-ext"));
    } else if (arg == "--ring-mod") {
      if (!ParseZZ(RequireNextArgValue("--ring-mod"), ring.scalar_modulus) ||
          ring.scalar_modulus <= 1) {
        std::cerr << "Invalid --ring-mod\n";
        return 2;
      }
    } else if (arg == "--ring-p") {
      if (!ParseZZ(RequireNextArgValue("--ring-p"), ring.base_prime) ||
          ring.base_prime <= 1) {
        std::cerr << "Invalid --ring-p\n";
        return 2;
      }
    } else if (arg == "--ring-F") {
      ring.F_coeffs = ParseCoeffList(RequireNextArgValue("--ring-F"));
    } else if (arg == "--ring-zeta") {
      ring.zeta_coeffs = ParseCoeffList(RequireNextArgValue("--ring-zeta"));
    } else if (arg == "--ring-challenge-ext") {
      ring.challenge_ext_coeffs =
          ParseNestedCoeffList(RequireNextArgValue("--ring-challenge-ext"));
    } else if (arg == "--help" || arg == "-h") {
      PrintHelp();
      return 0;
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
  if (request.mode.empty()) {
    std::cerr << "Missing required --mode\n";
    return 2;
  }

  request.artifact_root = artifact_root;
  const ContextSpec &spec = (request.mode == "field") ? field : ring;

  try {
    const DumpArtifactResult result = DumpEvalArtifact(spec, request);
    std::cout << "artifact_id " << result.metadata.artifact_id << "\n";
    std::cout << "case_path " << result.object_dir.string() << "\n";
    std::cout << "proof size " << result.metadata.proof_size_bytes << " B\n";
    std::cout << "display_key " << result.metadata.display_key << "\n";
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 2;
  }

  return 0;
}
