#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pX.h>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "bench_common_helpers.hpp"
#include "GaloisRing/PrimitiveElement.hpp"
#include "PCS/BaseFold/BaseFoldPCS.hpp"
#include "PCS/BaseFold/ProofSize.hpp"
#include "PCS/Common/Hash.hpp"
#include "PCS/Common/Multilinear.hpp"
#include "PCS/Common/Profile.hpp"

using NTL::conv;
using NTL::LogicError;
using NTL::SetCoeff;
using NTL::to_ZZ;
using NTL::vec_ZZ_pE;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pEX;
using NTL::ZZ_pEPush;
using NTL::ZZ_pPush;
using NTL::ZZ_pX;

namespace {

using basefold_bench_common::BuildZZpE;
using basefold_bench_common::BuildZZpEX;
using basefold_bench_common::BuildZZpX;
using basefold_bench_common::ComputeStats;
using basefold_bench_common::DeduceBasePrimeAndExponent;
using basefold_bench_common::IsPowerOfTwoLong;
using basefold_bench_common::Log2ExactPowerOfTwoLong;
using basefold_bench_common::MakeDeterministicCoefficients;
using basefold_bench_common::MakeDeterministicPoint;
using basefold_bench_common::MsSince;
using basefold_bench_common::NormalizeMod;
using basefold_bench_common::ParseCoeffList;
using basefold_bench_common::ParseInt;
using basefold_bench_common::ParseLong;
using basefold_bench_common::ParseNestedCoeffList;
using basefold_bench_common::ParseU64OrDie;
using basefold_bench_common::ParseZZ;
using basefold_bench_common::Pow2Checked;
using basefold_bench_common::Stats;
using basefold_bench_common::ValidateMonic;

struct ContextSpec {
  std::string label;
  ZZ scalar_modulus = ZZ(0);  // ZZ_p modulus (p for fields, p^s for rings)
  ZZ base_prime = ZZ(0);      // optional: the prime p (only used by unit checks)
  std::vector<ZZ> F_coeffs;     // extension modulus polynomial coefficients
  std::vector<ZZ> zeta_coeffs;  // ζ element coefficients
  // Coefficients for challenge extension modulus E(U), represented as
  // "a0;a1;...;ad", where each ai is a ZZ_pE element written "c0,c1,...".
  // Empty means "use default E(U) = U^2 + U + zeta".
  std::vector<std::vector<ZZ>> challenge_ext_coeffs;
};

basefold::FoldableCodeParams BuildParams_k0_1(long c, long d, const ZZ &prime_p,
                                             const ZZ_pE &zeta) {
  return basefold_bench_common::BuildFoldableParamsK0Eq1(
      c, d, prime_p, zeta, "BuildParams_k0_1");
}

basefold::FoldableCodeParams BuildParams_k0_pow2(long c, long k0, long d,
                                                 const ZZ &prime_p,
                                                 const ZZ_pE &zeta) {
  return basefold_bench_common::BuildFoldableParamsK0Pow2(
      c, k0, d, prime_p, zeta, "BuildParams_k0_pow2");
}

struct BenchResult {
  Stats prove_phase;
  Stats verifier;
  std::uint64_t anti_opt_checksum = 0;
  std::uint64_t proof_size_bytes = 0;
  double proof_size_kb = 0.0;
  basefold::Profile prover_profile;
  basefold::Profile verifier_profile;
  bool has_profile = false;
};

std::uint64_t ComputeProofSizeBytes(
    const basefold::BaseFoldPCSEvalProof &proof,
    bool use_extension_challenges, long challenge_ext_degree) {
  basefold::BaseFoldProofSizeOptions options;
  options.include_version_byte = true;
  if (use_extension_challenges || proof.extension.has_extension_payload) {
    options.challenge_ext_degree = challenge_ext_degree;
  }
  return basefold::BaseFoldPCSEvalProofSizeBytes(proof, options);
}

BenchResult RunEvalBenchmark(const vec_ZZ_pE &f_coeffs,
                             const std::vector<ZZ_pE> &z,
                             const ZZ_pE &y,
                             long num_queries,
                             const basefold::FoldableCodeParams &params,
                             const basefold::BaseFoldPCSChallengeConfig *challenge_cfg,
                             long challenge_ext_degree,
                             bool use_checked_prover_path,
                             bool enable_profile,
                             int warmup, int reps) {
  if (warmup < 0) LogicError("RunEvalBenchmark: warmup must be >= 0");
  if (reps <= 0) LogicError("RunEvalBenchmark: reps must be > 0");
  if (num_queries < 0) LogicError("RunEvalBenchmark: num_queries must be >= 0");

  std::vector<double> prove_phase_ms;
  std::vector<double> verifier_ms;
  prove_phase_ms.reserve(static_cast<std::size_t>(reps));
  verifier_ms.reserve(static_cast<std::size_t>(reps));

  basefold::Profile prover_prof;
  basefold::Profile verifier_prof;
  basefold::ResetProfile(prover_prof);
  basefold::ResetProfile(verifier_prof);

  std::uint64_t anti_opt_checksum = 0;
  std::uint64_t proof_size_bytes_last = 0;
  double proof_size_kb_last = 0.0;

  for (int iter = -warmup; iter < reps; ++iter) {
    const basefold::BaseFoldPCSCommitArtifacts commit_artifacts =
        basefold::BaseFoldPCSBuildCommitArtifactsUnchecked(f_coeffs, params);
    const basefold::MerkleRoot &commitment_root = commit_artifacts.root_d;

    const auto t0 = std::chrono::steady_clock::now();
    const basefold::BaseFoldPCSEvalProof proof = [&] {
      if (enable_profile && iter >= 0) {
        basefold::ProfileGuard guard(&prover_prof);
        basefold::ScopedTimer timer(&prover_prof.pcs_prove_ns,
                                    &prover_prof.pcs_prove_calls);
        if (challenge_cfg != nullptr) {
          return use_checked_prover_path
                     ? basefold::BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracle(
                           f_coeffs, z, y, num_queries, params,
                           commit_artifacts, *challenge_cfg)
                     : basefold::BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracleUnchecked(
                           f_coeffs, z, y, num_queries, params,
                           commit_artifacts, *challenge_cfg);
        }
        return use_checked_prover_path
                   ? basefold::BaseFoldPCSProveEvalFromCommittedTopOracle(
                         f_coeffs, z, y, num_queries, params, commit_artifacts)
                   : basefold::BaseFoldPCSProveEvalFromCommittedTopOracleUnchecked(
                         f_coeffs, z, y, num_queries, params, commit_artifacts);
      }
      if (challenge_cfg != nullptr) {
        return use_checked_prover_path
                   ? basefold::BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracle(
                         f_coeffs, z, y, num_queries, params,
                         commit_artifacts, *challenge_cfg)
                   : basefold::BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracleUnchecked(
                         f_coeffs, z, y, num_queries, params,
                         commit_artifacts, *challenge_cfg);
      }
      return use_checked_prover_path
                 ? basefold::BaseFoldPCSProveEvalFromCommittedTopOracle(
                       f_coeffs, z, y, num_queries, params, commit_artifacts)
                 : basefold::BaseFoldPCSProveEvalFromCommittedTopOracleUnchecked(
                       f_coeffs, z, y, num_queries, params, commit_artifacts);
    }();
    const auto t1 = std::chrono::steady_clock::now();

    const std::uint64_t proof_size_bytes = ComputeProofSizeBytes(
        proof, challenge_cfg != nullptr, challenge_ext_degree);
    const double proof_size_kb =
        static_cast<double>(proof_size_bytes) / 1024.0;

    const auto t2 = std::chrono::steady_clock::now();
    bool ok = false;
    if (enable_profile && iter >= 0) {
      basefold::ProfileGuard guard(&verifier_prof);
      ok = (challenge_cfg != nullptr)
               ? basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
                     commitment_root, z, y, num_queries, proof, params,
                     *challenge_cfg)
               : basefold::BaseFoldPCSVerifyEval(commitment_root, z, y,
                                                 num_queries, proof, params);
    } else {
      ok = (challenge_cfg != nullptr)
               ? basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
                     commitment_root, z, y, num_queries, proof, params,
                     *challenge_cfg)
               : basefold::BaseFoldPCSVerifyEval(commitment_root, z, y,
                                                 num_queries, proof, params);
    }
    const auto t3 = std::chrono::steady_clock::now();

    if (!ok) {
      LogicError("RunEvalBenchmark: verification failed");
    }

    // Prevent over-optimization: fold in a few bytes.
    if (!commitment_root.empty())
      anti_opt_checksum ^= static_cast<std::uint64_t>(commitment_root[0]);
    anti_opt_checksum ^= static_cast<std::uint64_t>(ok);

    if (iter >= 0) {
      prove_phase_ms.push_back(MsSince(t0, t1));
      verifier_ms.push_back(MsSince(t2, t3));
      proof_size_bytes_last = proof_size_bytes;
      proof_size_kb_last = proof_size_kb;
    }
  }

  BenchResult out;
  out.prove_phase = ComputeStats(prove_phase_ms);
  out.verifier = ComputeStats(verifier_ms);
  out.anti_opt_checksum = anti_opt_checksum;
  out.proof_size_bytes = proof_size_bytes_last;
  out.proof_size_kb = proof_size_kb_last;
  out.prover_profile = prover_prof;
  out.verifier_profile = verifier_prof;
  out.has_profile = enable_profile;
  return out;
}

void PrintResult(const std::string &label, const ZZ &scalar_modulus, long c,
                 long d, long k0, long num_queries, int warmup, int reps,
                 const BenchResult &r, bool use_extension_challenges,
                 long challenge_ext_degree) {
  const long pow2_d = Pow2Checked(d);
  if (k0 > std::numeric_limits<long>::max() / pow2_d) {
    LogicError("PrintResult: overflow in k_d");
  }
  const long k_d = k0 * pow2_d;
  if (c > std::numeric_limits<long>::max() / k_d) {
    LogicError("PrintResult: overflow in n_d");
  }
  const long n_d = c * k_d;

  std::cout << "\n[" << label << "] c=" << c << " k0=" << k0 << " d=" << d
            << "  mod=" << scalar_modulus << "  k_d=" << k_d
            << "  n_d=" << n_d
            << "  queries=" << num_queries << "  warmup=" << warmup
            << " reps=" << reps;
  if (use_extension_challenges) {
    std::cout << "  ext_challenges=on"
              << "  ext_deg=" << challenge_ext_degree;
  } else {
    std::cout << "  ext_challenges=off";
  }
  std::cout << "\n";
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "  hash backend " << basefold::SelectedHashBackendName() << "\n";
  std::cout << "  prove-phase mean " << r.prove_phase.mean_ms << " ms  (min "
            << r.prove_phase.min_ms << ", max " << r.prove_phase.max_ms
            << ")\n";
  std::cout << "  verifier mean " << r.verifier.mean_ms << " ms  (min "
            << r.verifier.min_ms << ", max " << r.verifier.max_ms << ")\n";
  std::cout << "  proof size  " << r.proof_size_kb << " KB  ("
            << r.proof_size_bytes << " B)\n";

  std::cout << "  anti-opt checksum " << r.anti_opt_checksum << "\n";
  if (r.has_profile) {
    basefold::PrintProfile(std::cout, r.prover_profile);
    basefold::PrintProfile(std::cout, r.verifier_profile);
  }
}

void PrintHelp() {
  std::cout
      << "bench_basefold_pcs_eval (PCS prove+verify; prover excludes top encode+commit)\n\n"
      << "Usage:\n"
      << "  bench_basefold_pcs_eval [--mode field|ring|both] [--c <int>] [--k0 <int>] [--d <int>]\n"
      << "               [--queries <int>] [--checked] [--profile] [--warmup <int>] [--reps <int>] [--seed <u64>]\n"
      << "               [--merkle-leaves-per-thread <int>] [--merkle-level-threshold <int>] [--merkle-max-threads <int>]\n"
      << "               [--verifier-query-per-thread <int>] [--verifier-query-threshold <int>] [--verifier-query-max-threads <int>]\n"
      << "               [--prover-commit-base-per-thread <int>] [--prover-commit-ext-per-thread <int>]\n"
      << "               [--use-extension-challenges]\n"
      << "               [--field-challenge-ext <a0;a1;...>] [--ring-challenge-ext <a0;a1;...>]\n"
      << "               [--auto-zeta teich]\n"
      << "               [--field-mod <decimal-int>] [--field-F <a0,a1,...>] [--field-zeta <b0,b1,...>]\n"
      << "               [--ring-mod <decimal-int>]  [--ring-p <decimal-int>] [--ring-F <a0,a1,...>] [--ring-zeta <b0,b1,...>]\n\n"
      << "Notes:\n"
      << "  By default, prover uses BaseFoldPCSProveEvalFromCommittedTopOracleUnchecked\n"
      << "  (skips validation and claimed_y check), so headline prover time excludes\n"
      << "  the top-level EncodeFoldable + Merkle commit stage.\n\n"
      << "  With --auto-zeta teich, zeta is derived as a Teichmuller generator from (p,k,F);\n"
      << "  then --field-zeta/--ring-zeta are ignored.\n\n"
      << "  With --use-extension-challenges, prove/verify uses the extension-challenge\n"
      << "  path in BaseFoldPCSChallengeConfig.\n"
      << "  --field-challenge-ext / --ring-challenge-ext use ';' to separate ZZ_pE\n"
      << "  coefficients and ',' for each ZZ_pE coefficient polynomial.\n"
      << "  Example: '0,1;1;1' means E(U)=x + U + U^2.\n\n"
      << "  Exact fixed-width proof payload size is computed from the generated proof\n"
      << "  with transcript-recoverable indices/challenges omitted, and printed\n"
      << "  together with prover/verifier timing.\n\n"
      << "  If --*-challenge-ext is omitted, default is E(U)=zeta + U + U^2.\n\n"
      << "  Merkle build parallel tuning can be configured by env vars:\n"
      << "    BASEFOLD_MERKLE_LEAVES_PER_THREAD\n"
      << "    BASEFOLD_MERKLE_PARALLEL_LEVEL_THRESHOLD\n"
      << "    BASEFOLD_MERKLE_MAX_THREADS\n"
      << "  Verifier query parallel tuning can be configured by env vars:\n"
      << "    BASEFOLD_VERIFY_QUERY_QUERIES_PER_THREAD\n"
      << "    BASEFOLD_VERIFY_QUERY_PARALLEL_THRESHOLD\n"
      << "    BASEFOLD_VERIFY_QUERY_MAX_THREADS\n"
      << "  Prover commit-round parallel tuning can be configured by env vars:\n"
      << "    BASEFOLD_PROVER_COMMIT_BASE_ELEMENTS_PER_THREAD\n"
      << "    BASEFOLD_PROVER_COMMIT_EXT_ELEMENTS_PER_THREAD\n"
      << "  and overridden by the CLI flags above.\n\n"
      << "  PCS Eval supports k0 = 2^κ. The multilinear point dimension is (d + κ).\n\n"
      << "Examples:\n"
      << "  # GF(2^2) with F(x)=x^2+x+1 and zeta=x\n"
      << "  bench_basefold_pcs_eval --mode field --field-mod 2 --field-F 1,1,1 --field-zeta 0,1 --d 16 --queries 4\n"
      << "  # GR(4,2) with the same extension polynomial and zeta=x\n"
      << "  bench_basefold_pcs_eval --mode ring  --ring-mod 4 --ring-p 2 --ring-F 1,1,1 --ring-zeta 0,1 --d 16 --queries 4\n"
      << "  # Auto zeta from Teichmuller subgroup generator\n"
      << "  bench_basefold_pcs_eval --mode ring  --ring-mod 4 --ring-p 2 --ring-F 1,1,1 --auto-zeta teich --d 16 --queries 4\n"
      << "  # Extension-challenge path (quote ';' argument)\n"
      << "  bench_basefold_pcs_eval --mode field --field-mod 2 --field-F 1,1,1 --field-zeta 0,1 --use-extension-challenges --field-challenge-ext '0,1;1;1' --d 16 --queries 4\n"
      << "  bench_basefold_pcs_eval --mode ring  --ring-mod 4 --ring-p 2 --ring-F 1,1,1 --ring-zeta 0,1 --use-extension-challenges --ring-challenge-ext '0,1;1;1' --d 16 --queries 4\n";
}

void RunOneContext(const ContextSpec &spec, long c, long k0, long d,
                   long num_queries,
                   bool use_extension_challenges,
                   bool use_checked_prover_path, bool enable_profile, int warmup,
                   int reps, bool auto_zeta_teich,
                   std::uint64_t seed) {
  if (spec.scalar_modulus <= 1) LogicError("RunOneContext: modulus must be > 1");
  if (k0 <= 0) LogicError("RunOneContext: k0 must be positive");
  if (!IsPowerOfTwoLong(k0)) LogicError("RunOneContext: k0 must be a power of two");

  const ZZ modulus = spec.scalar_modulus;
  ZZ_pPush mod_push(modulus);

  ValidateMonic(spec.F_coeffs, spec.scalar_modulus, "F");
  const ZZ_pX F = BuildZZpX(spec.F_coeffs);
  ZZ_pEPush e_push(F);

  ZZ_pE zeta;
  if (auto_zeta_teich) {
    ZZ p_base;
    long k_base = 0;
    DeduceBasePrimeAndExponent(spec, p_base, k_base);
    const long s = NTL::deg(F);
    zeta = FindTeichmullerGenerator(p_base, k_base, s, F);
  } else {
    const ZZ_pX zpoly = BuildZZpX(spec.zeta_coeffs);
    conv(zeta, zpoly);
  }

  const basefold::FoldableCodeParams params = (k0 == 1)
                                                  ? BuildParams_k0_1(
                                                        c, d,
                                                        (spec.base_prime > 1)
                                                            ? spec.base_prime
                                                            : spec.scalar_modulus,
                                                        zeta)
                                                  : BuildParams_k0_pow2(
                                                        c, k0, d,
                                                        (spec.base_prime > 1)
                                                            ? spec.base_prime
                                                            : spec.scalar_modulus,
                                                        zeta);

  const long pow2_d = Pow2Checked(d);
  if (k0 > std::numeric_limits<long>::max() / pow2_d) {
    LogicError("RunOneContext: overflow in k_d");
  }
  const long k_d = k0 * pow2_d;
  const vec_ZZ_pE f_coeffs = MakeDeterministicCoefficients(k_d, seed);
  const long point_dim = d + Log2ExactPowerOfTwoLong(k0);
  const std::vector<ZZ_pE> z =
      MakeDeterministicPoint(point_dim, seed ^ 0xdeadbeefULL);
  const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);

  basefold::BaseFoldPCSChallengeConfig challenge_cfg;
  long challenge_ext_degree = 0;
  const basefold::BaseFoldPCSChallengeConfig *challenge_cfg_ptr = nullptr;
  if (use_extension_challenges) {
    ZZ_pEX challenge_modulus;
    if (spec.challenge_ext_coeffs.empty()) {
      NTL::clear(challenge_modulus);
      NTL::SetCoeff(challenge_modulus, 0, zeta);
      NTL::SetCoeff(challenge_modulus, 1, ZZ_pE(1));
      NTL::SetCoeff(challenge_modulus, 2, ZZ_pE(1));
    } else {
      challenge_modulus = BuildZZpEX(spec.challenge_ext_coeffs);
    }
    challenge_cfg.use_extension_challenges = true;
    challenge_cfg.challenge_extension_modulus = challenge_modulus;
    challenge_ext_degree = NTL::deg(challenge_modulus);
    challenge_cfg_ptr = &challenge_cfg;
  }

  const BenchResult r = RunEvalBenchmark(f_coeffs, z, y, num_queries, params,
                                         challenge_cfg_ptr, challenge_ext_degree,
                                         use_checked_prover_path,
                                         enable_profile, warmup, reps);
  PrintResult(spec.label, spec.scalar_modulus, c, d, k0, num_queries, warmup,
              reps, r, use_extension_challenges, challenge_ext_degree);
}

}  // namespace

int main(int argc, char **argv) {
  long d = 16;
  long c = 2;
  long k0 = 1;
  long num_queries = 4;
  bool use_extension_challenges = false;
  bool use_checked_prover_path = false;
  bool enable_profile = false;
  int warmup = 1;
  int reps = 3;
  std::uint64_t seed = 0;
  bool auto_zeta_teich = false;
  basefold::ResetMerkleBuildParallelConfigFromEnv();
  basefold::MerkleBuildParallelConfig merkle_cfg =
      basefold::GetMerkleBuildParallelConfig();
  basefold::ResetVerifierQueryParallelConfigFromEnv();
  basefold::VerifierQueryParallelConfig verifier_query_cfg =
      basefold::GetVerifierQueryParallelConfig();
  basefold::ResetProverCommitParallelConfigFromEnv();
  basefold::ProverCommitParallelConfig prover_commit_cfg =
      basefold::GetProverCommitParallelConfig();

  bool do_field = true;
  bool do_ring = true;

  ContextSpec field;
  field.label = "Field";
  field.scalar_modulus = to_ZZ(2);
  field.base_prime = ZZ(0);
  field.F_coeffs = {to_ZZ(1), to_ZZ(1), to_ZZ(1)};   // x^2 + x + 1
  field.zeta_coeffs = {to_ZZ(0), to_ZZ(1)};   // x

  ContextSpec ring;
  ring.label = "Ring";
  ring.scalar_modulus = to_ZZ(4);
  ring.base_prime = to_ZZ(2);
  ring.F_coeffs = {to_ZZ(1), to_ZZ(1), to_ZZ(1)};    // x^2 + x + 1
  ring.zeta_coeffs = {to_ZZ(0), to_ZZ(1)};    // x

  for (int i = 1; i < argc; ++i) {
    const std::string arg(argv[i]);
    auto RequireNextArgValue = [&](const char *flag) -> const char * {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << flag << "\n";
        std::exit(2);
      }
      return argv[++i];
    };

    if (arg == "--mode") {
      const std::string m = RequireNextArgValue("--mode");
      if (m == "field") {
        do_field = true;
        do_ring = false;
      } else if (m == "ring") {
        do_field = false;
        do_ring = true;
      } else if (m == "both") {
        do_field = true;
        do_ring = true;
      } else {
        std::cerr << "Invalid --mode (expected field|ring|both)\n";
        return 2;
      }
    } else if (arg == "--d") {
      if (!ParseLong(RequireNextArgValue("--d"), d) || d < 0) {
        std::cerr << "Invalid --d\n";
        return 2;
      }
    } else if (arg == "--c") {
      if (!ParseLong(RequireNextArgValue("--c"), c) || c <= 0) {
        std::cerr << "Invalid --c\n";
        return 2;
      }
    } else if (arg == "--k0") {
      if (!ParseLong(RequireNextArgValue("--k0"), k0) || k0 <= 0) {
        std::cerr << "Invalid --k0\n";
        return 2;
      }
    } else if (arg == "--queries") {
      if (!ParseLong(RequireNextArgValue("--queries"), num_queries) ||
          num_queries < 0) {
        std::cerr << "Invalid --queries\n";
        return 2;
      }
    } else if (arg == "--checked") {
      use_checked_prover_path = true;
    } else if (arg == "--use-extension-challenges") {
      use_extension_challenges = true;
    } else if (arg == "--profile") {
      enable_profile = true;
    } else if (arg == "--warmup") {
      if (!ParseInt(RequireNextArgValue("--warmup"), warmup) || warmup < 0) {
        std::cerr << "Invalid --warmup\n";
        return 2;
      }
    } else if (arg == "--merkle-leaves-per-thread") {
      if (!ParseLong(RequireNextArgValue("--merkle-leaves-per-thread"),
                     merkle_cfg.leaves_per_thread) ||
          merkle_cfg.leaves_per_thread <= 0) {
        std::cerr << "Invalid --merkle-leaves-per-thread\n";
        return 2;
      }
    } else if (arg == "--merkle-level-threshold") {
      if (!ParseLong(RequireNextArgValue("--merkle-level-threshold"),
                     merkle_cfg.parallel_level_threshold) ||
          merkle_cfg.parallel_level_threshold <= 0) {
        std::cerr << "Invalid --merkle-level-threshold\n";
        return 2;
      }
    } else if (arg == "--merkle-max-threads") {
      if (!ParseInt(RequireNextArgValue("--merkle-max-threads"),
                    merkle_cfg.max_threads) ||
          merkle_cfg.max_threads <= 0) {
        std::cerr << "Invalid --merkle-max-threads\n";
        return 2;
      }
    } else if (arg == "--verifier-query-per-thread") {
      if (!ParseLong(RequireNextArgValue("--verifier-query-per-thread"),
                     verifier_query_cfg.queries_per_thread) ||
          verifier_query_cfg.queries_per_thread <= 0) {
        std::cerr << "Invalid --verifier-query-per-thread\n";
        return 2;
      }
    } else if (arg == "--verifier-query-threshold") {
      if (!ParseLong(RequireNextArgValue("--verifier-query-threshold"),
                     verifier_query_cfg.min_queries_for_parallelism) ||
          verifier_query_cfg.min_queries_for_parallelism <= 0) {
        std::cerr << "Invalid --verifier-query-threshold\n";
        return 2;
      }
    } else if (arg == "--verifier-query-max-threads") {
      if (!ParseInt(RequireNextArgValue("--verifier-query-max-threads"),
                    verifier_query_cfg.max_threads) ||
          verifier_query_cfg.max_threads <= 0) {
        std::cerr << "Invalid --verifier-query-max-threads\n";
        return 2;
      }
    } else if (arg == "--prover-commit-base-per-thread") {
      if (!ParseLong(RequireNextArgValue("--prover-commit-base-per-thread"),
                     prover_commit_cfg.base_elements_per_thread) ||
          prover_commit_cfg.base_elements_per_thread <= 0) {
        std::cerr << "Invalid --prover-commit-base-per-thread\n";
        return 2;
      }
      prover_commit_cfg.base_elements_per_thread_overridden = true;
    } else if (arg == "--prover-commit-ext-per-thread") {
      if (!ParseLong(RequireNextArgValue("--prover-commit-ext-per-thread"),
                     prover_commit_cfg.ext_elements_per_thread) ||
          prover_commit_cfg.ext_elements_per_thread <= 0) {
        std::cerr << "Invalid --prover-commit-ext-per-thread\n";
        return 2;
      }
      prover_commit_cfg.ext_elements_per_thread_overridden = true;
    } else if (arg == "--reps") {
      if (!ParseInt(RequireNextArgValue("--reps"), reps) || reps <= 0) {
        std::cerr << "Invalid --reps\n";
        return 2;
      }
    } else if (arg == "--seed") {
      seed = ParseU64OrDie(RequireNextArgValue("--seed"), "--seed");
    } else if (arg == "--auto-zeta") {
      const std::string mode = RequireNextArgValue("--auto-zeta");
      if (mode == "teich") {
        auto_zeta_teich = true;
      } else {
        std::cerr << "Invalid --auto-zeta (expected teich)\n";
        return 2;
      }
    } else if (arg == "--field-mod") {
      if (!ParseZZ(RequireNextArgValue("--field-mod"),
                   field.scalar_modulus) ||
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
      if (!ParseZZ(RequireNextArgValue("--ring-mod"),
                   ring.scalar_modulus) ||
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
      std::cerr << "Unknown arg: " << arg << "\n";
      return 2;
    }
  }

  try {
    basefold::SetMerkleBuildParallelConfig(merkle_cfg);
    basefold::SetVerifierQueryParallelConfig(verifier_query_cfg);
    basefold::SetProverCommitParallelConfig(prover_commit_cfg);
    if (!do_field && !do_ring) {
      std::cerr << "Nothing to do: --mode disabled both field and ring\n";
      return 2;
    }
    if (do_field)
      RunOneContext(field, c, k0, d, num_queries, use_extension_challenges,
                    use_checked_prover_path, enable_profile,
                    warmup, reps, auto_zeta_teich, seed);
    if (do_ring)
      RunOneContext(ring, c, k0, d, num_queries, use_extension_challenges,
                    use_checked_prover_path, enable_profile,
                    warmup, reps, auto_zeta_teich, seed);
  } catch (const std::exception &e) {
    std::cerr << "Unhandled std::exception: " << e.what() << "\n";
    return 2;
  } catch (...) {
    std::cerr << "Unhandled non-std exception\n";
    return 2;
  }

  return 0;
}
