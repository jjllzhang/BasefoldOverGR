#include "BaseFoldPCSInternal.hpp"

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

#include "GaloisRing/Inverse.hpp"

using NTL::LogicError;
using NTL::ZZ;
using NTL::ZZ_pE;
using NTL::ZZ_pX;

namespace basefold {
namespace basefold_pcs_internal {
namespace {

void SortAndUniqueIndices(std::vector<long> &indices) {
  std::sort(indices.begin(), indices.end());
  indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
}

}  // namespace

long ParsePositiveEnvLong(const char *name, long fallback) {
  const char *raw = std::getenv(name);
  if (raw == nullptr || *raw == '\0') {
    return fallback;
  }

  char *end = nullptr;
  const long value = std::strtol(raw, &end, 10);
  if (end == raw || *end != '\0' || value <= 0) {
    return fallback;
  }
  return value;
}

int ParsePositiveEnvInt(const char *name, int fallback) {
  const char *raw = std::getenv(name);
  if (raw == nullptr || *raw == '\0') {
    return fallback;
  }

  char *end = nullptr;
  const long value = std::strtol(raw, &end, 10);
  if (end == raw || *end != '\0' || value <= 0 ||
      value > static_cast<long>(std::numeric_limits<int>::max())) {
    return fallback;
  }
  return static_cast<int>(value);
}

VerifierQueryParallelConfig ParseVerifierQueryParallelConfigFromEnv() {
  VerifierQueryParallelConfig cfg;
  cfg.queries_per_thread = ParsePositiveEnvLong(
      "BASEFOLD_VERIFY_QUERY_QUERIES_PER_THREAD", cfg.queries_per_thread);
  cfg.min_queries_for_parallelism =
      ParsePositiveEnvLong("BASEFOLD_VERIFY_QUERY_PARALLEL_THRESHOLD",
                           cfg.min_queries_for_parallelism);
  cfg.max_threads =
      ParsePositiveEnvInt("BASEFOLD_VERIFY_QUERY_MAX_THREADS", cfg.max_threads);
  return cfg;
}

ProverCommitParallelConfig ParseProverCommitParallelConfigFromEnv() {
  ProverCommitParallelConfig cfg;
  if (const char *base_env =
          std::getenv("BASEFOLD_PROVER_COMMIT_BASE_ELEMENTS_PER_THREAD")) {
    cfg.base_elements_per_thread =
        ParsePositiveEnvLong("BASEFOLD_PROVER_COMMIT_BASE_ELEMENTS_PER_THREAD",
                             cfg.base_elements_per_thread);
    cfg.base_elements_per_thread_overridden = (base_env[0] != '\0');
  }
  if (const char *ext_env =
          std::getenv("BASEFOLD_PROVER_COMMIT_EXT_ELEMENTS_PER_THREAD")) {
    cfg.ext_elements_per_thread =
        ParsePositiveEnvLong("BASEFOLD_PROVER_COMMIT_EXT_ELEMENTS_PER_THREAD",
                             cfg.ext_elements_per_thread);
    cfg.ext_elements_per_thread_overridden = (ext_env[0] != '\0');
  }
  return cfg;
}

VerifierQueryParallelConfig &MutableVerifierQueryParallelConfig() {
  static VerifierQueryParallelConfig cfg =
      ParseVerifierQueryParallelConfigFromEnv();
  return cfg;
}

ProverCommitParallelConfig &MutableProverCommitParallelConfig() {
  static ProverCommitParallelConfig cfg =
      ParseProverCommitParallelConfigFromEnv();
  return cfg;
}

VerifierQueryParallelConfig LoadVerifierQueryParallelConfig() {
  VerifierQueryParallelConfig cfg = MutableVerifierQueryParallelConfig();
  if (cfg.queries_per_thread <= 0) {
    cfg.queries_per_thread = 1;
  }
  if (cfg.min_queries_for_parallelism <= 0) {
    cfg.min_queries_for_parallelism = 2;
  }
  if (cfg.max_threads <= 0) {
    cfg.max_threads = 8;
  }
  return cfg;
}

ProverCommitParallelConfig LoadProverCommitParallelConfig() {
  ProverCommitParallelConfig cfg = MutableProverCommitParallelConfig();
  if (cfg.base_elements_per_thread <= 0) {
    cfg.base_elements_per_thread = 4096;
  }
  if (cfg.ext_elements_per_thread <= 0) {
    cfg.ext_elements_per_thread = 128;
  }
  return cfg;
}

long LoadEffectiveExtElementsPerThread() {
  const ProverCommitParallelConfig cfg = LoadProverCommitParallelConfig();
  if (!cfg.ext_elements_per_thread_overridden &&
      cfg.ext_elements_per_thread == 128) {
    struct CachedPrimeStatus {
      bool initialized = false;
      ZZ modulus;
      bool is_prime = false;
    };
    thread_local CachedPrimeStatus cache;
    const ZZ current_modulus = NTL::ZZ_p::modulus();
    if (!cache.initialized || cache.modulus != current_modulus) {
      cache.modulus = current_modulus;
      cache.is_prime = NTL::ProbPrime(current_modulus);
      cache.initialized = true;
    }
    if (!cache.is_prime) {
      return 64;
    }
  }
  return cfg.ext_elements_per_thread;
}

int ChooseQueryVerifyThreads(long num_queries) {
#if defined(BASEFOLD_USE_OPENMP)
  const VerifierQueryParallelConfig cfg = LoadVerifierQueryParallelConfig();
  if (num_queries < cfg.min_queries_for_parallelism) {
    return 1;
  }

  const long blocks =
      (num_queries + cfg.queries_per_thread - 1) / cfg.queries_per_thread;
  int threads_to_use = static_cast<int>(blocks);
  if (threads_to_use > cfg.max_threads) {
    threads_to_use = cfg.max_threads;
  }
  const int max_threads = omp_get_max_threads();
  if (threads_to_use > max_threads) {
    threads_to_use = max_threads;
  }
  if (threads_to_use < 1) {
    threads_to_use = 1;
  }
  return threads_to_use;
#else
  (void)num_queries;
  return 1;
#endif
}

void ValidateParamsOrThrow(const FoldableCodeParams &params) {
  (void)MessageLength(params);
}

long Pow2Checked(long e) {
  if (e < 0) {
    LogicError("Pow2Checked: negative exponent");
  }
  if (e >= static_cast<long>(8 * sizeof(long) - 1)) {
    LogicError("Pow2Checked: exponent too large for long");
  }
  return 1L << e;
}

bool IsPowerOfTwoLong(long n) {
  return n > 0 && (n & (n - 1)) == 0;
}

long Log2ExactPowerOfTwoLong(long n) {
  if (!IsPowerOfTwoLong(n)) {
    LogicError("Log2ExactPowerOfTwoLong: not a power of two");
  }
  long out = 0;
  while (n > 1) {
    n >>= 1;
    ++out;
  }
  return out;
}

long CodewordLengthAtLevelNoValidate(const FoldableCodeParams &params,
                                     long level) {
  if (level < 0 || level > params.d) {
    LogicError("CodewordLengthAtLevelNoValidate: level out of range");
  }
  const long pow2 = Pow2Checked(level);
  if (params.k0 <= 0 || params.c <= 0) {
    LogicError("CodewordLengthAtLevelNoValidate: invalid c/k0");
  }
  if (params.k0 > std::numeric_limits<long>::max() / pow2) {
    LogicError("CodewordLengthAtLevelNoValidate: overflow");
  }
  const long k = params.k0 * pow2;
  if (params.c > std::numeric_limits<long>::max() / k) {
    LogicError("CodewordLengthAtLevelNoValidate: overflow");
  }
  return params.c * k;
}

IOPPQueryPlan MakeQueryPlanNoValidate(long initial_mu,
                                      const FoldableCodeParams &params) {
  IOPPQueryPlan plan;
  plan.initial_mu = initial_mu;
  plan.mu_by_level.resize(static_cast<std::size_t>(params.d));

  if (params.d == 0) {
    return plan;
  }

  const long n_last = CodewordLengthAtLevelNoValidate(params, params.d - 1);
  if (initial_mu < 0 || initial_mu >= n_last) {
    LogicError("MakeQueryPlanNoValidate: initial_mu out of range");
  }

  long mu = initial_mu;
  for (long i = params.d; i-- > 0;) {
    plan.mu_by_level[static_cast<std::size_t>(i)] = mu;
    if (i > 0) {
      const long n_prev = CodewordLengthAtLevelNoValidate(params, i - 1);
      if (mu >= n_prev) {
        mu -= n_prev;
      }
    }
  }

  return plan;
}

void ValidateCommittedTopOracleArtifactsOrThrow(
    const BaseFoldPCSCommitArtifacts &commit_artifacts, long expected_length,
    const char *func_name) {
  if (commit_artifacts.pi_d.length() != expected_length) {
    LogicError((std::string(func_name) +
                ": commit_artifacts.pi_d has wrong length")
                   .c_str());
  }
  if (commit_artifacts.merkle_d.Root() != commit_artifacts.root_d) {
    LogicError((std::string(func_name) + ": commit_artifacts root mismatch")
                   .c_str());
  }
}

std::vector<std::vector<long>> CollectBaseQueryIndicesByTree(
    const std::vector<IOPPQueryPlan> &query_plans,
    const FoldableCodeParams &params) {
  std::vector<std::vector<long>> requested(
      static_cast<std::size_t>(params.d + 1));
  for (const IOPPQueryPlan &plan : query_plans) {
    for (long i = 0; i < params.d; ++i) {
      const long mu_i = plan.mu_by_level[static_cast<std::size_t>(i)];
      const long n_i = CodewordLengthAtLevelNoValidate(params, i);
      requested[static_cast<std::size_t>(i)].push_back(mu_i);
      requested[static_cast<std::size_t>(i + 1)].push_back(mu_i);
      requested[static_cast<std::size_t>(i + 1)].push_back(mu_i + n_i);
    }
  }
  for (std::vector<long> &indices : requested) {
    SortAndUniqueIndices(indices);
  }
  return requested;
}

bool HasMerkleMultiproofPayload(const MerkleMultiproof &proof) {
  return !proof.queried_indices.empty() || proof.values.length() != 0 ||
         !proof.sibling_hashes.empty();
}

HashTranscript MakeBaseFoldTranscript() {
  HashTranscriptConfig config;
  config.domain_separator = "BaseFoldPCS/v1";
  config.byte_order = TranscriptByteOrder::kBigEndian;
  config.error_prefix = "BaseFoldHashTranscript";
  return HashTranscript(config);
}

void AbsorbPublicInput(HashTranscript &transcript,
                       const MerkleRoot &commitment,
                       const std::vector<FieldElement> &z,
                       const FieldElement &y) {
  transcript.AbsorbDigest(commitment);
  for (const FieldElement &zi : z) {
    transcript.AbsorbFieldElement(zi);
  }
  transcript.AbsorbFieldElement(y);
}

bool TryInvertBaseUnit(FieldElement &inv_out, const FieldElement &a) {
  if (a == 0) {
    return false;
  }

  const ZZ modulus = NTL::ZZ_p::modulus();
  if (modulus <= 1) {
    LogicError("TryInvertBaseUnit: invalid base modulus");
  }

  if (NTL::ProbPrime(modulus)) {
    try {
      inv_out = NTL::inv(a);
      return true;
    } catch (...) {
      return false;
    }
  }

  const long r = ZZ_pE::degree();
  if (r <= 0) {
    LogicError("TryInvertBaseUnit: invalid extension degree");
  }

  if (r == 1) {
    const ZZ a_rep = NTL::rep(NTL::coeff(NTL::rep(a), 0));
    ZZ inv_rep;
    if (NTL::InvModStatus(inv_rep, a_rep, modulus) != 0) {
      return false;
    }
    NTL::ZZ_p inv_base;
    NTL::conv(inv_base, inv_rep);
    ZZ_pX poly;
    NTL::clear(poly);
    NTL::SetCoeff(poly, 0, inv_base);
    NTL::conv(inv_out, poly);
    return true;
  }

  inv_out = Inv(a, r);
  if (inv_out == 0) {
    return false;
  }

  FieldElement one;
  NTL::set(one);
  return a * inv_out == one;
}

bool BatchInvertBaseUnits(std::vector<FieldElement> &inverses,
                          const std::vector<FieldElement> &values) {
  const long n = static_cast<long>(values.size());
  inverses.resize(static_cast<std::size_t>(n));
  if (n == 0) {
    return true;
  }

  std::vector<FieldElement> prefix(static_cast<std::size_t>(n + 1));
  FieldElement one;
  NTL::set(one);
  prefix[0] = one;

  for (long i = 0; i < n; ++i) {
    const FieldElement &value = values[static_cast<std::size_t>(i)];
    if (value == 0) {
      return false;
    }
    prefix[static_cast<std::size_t>(i + 1)] =
        prefix[static_cast<std::size_t>(i)] * value;
  }

  FieldElement inv_total;
  if (!TryInvertBaseUnit(inv_total, prefix[static_cast<std::size_t>(n)])) {
    return false;
  }

  FieldElement suffix = inv_total;
  for (long i = n; i-- > 0;) {
    inverses[static_cast<std::size_t>(i)] =
        suffix * prefix[static_cast<std::size_t>(i)];
    suffix *= values[static_cast<std::size_t>(i)];
  }
  return true;
}

FieldElement EvalLineAtWithInvDenom(const FieldElement &x,
                                    const FieldElement &x1,
                                    const FieldElement &y1,
                                    const FieldElement &y2,
                                    const FieldElement &inv_denom) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->eval_line_at_ns : nullptr,
                    prof ? &prof->eval_line_at_calls : nullptr);
  return y1 + (x - x1) * (y2 - y1) * inv_denom;
}

void ProverCommitRoundNoValidate(Oracle &pi_i, const Oracle &pi_ip1,
                                 const FieldElement &alpha_i, long level_i,
                                 const FoldableCodeParams &params) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->prover_commit_round_ns : nullptr,
                    prof ? &prof->prover_commit_round_calls : nullptr);

  const long n_i = CodewordLengthAtLevelNoValidate(params, level_i);
  pi_i.SetLength(n_i);

  const Oracle &diag = params.diag_T[static_cast<std::size_t>(level_i)];
  const ProverCommitParallelConfig cfg = LoadProverCommitParallelConfig();
  const long parallel_threshold = cfg.base_elements_per_thread;

  if (n_i > 0) {
    const FieldElement first_x1 = diag[0];
    bool all_equal = true;
    for (long j = 1; j < n_i; ++j) {
      if (diag[static_cast<std::size_t>(j)] != first_x1) {
        all_equal = false;
        break;
      }
    }
    if (all_equal) {
      const FieldElement denom = (params.zeta * first_x1) - first_x1;
      if (denom == 0) {
        LogicError("ProverCommitRoundNoValidate: x1 must not equal x2");
      }
      FieldElement inv_denom;
      if (!TryInvertBaseUnit(inv_denom, denom)) {
        LogicError("ProverCommitRoundNoValidate: x2-x1 must be a unit");
      }
      ForEachIndexMaybeParallel(0, n_i, parallel_threshold, [&](long j) {
        pi_i[static_cast<std::size_t>(j)] = EvalLineAtWithInvDenom(
            alpha_i, first_x1, pi_ip1[static_cast<std::size_t>(j)],
            pi_ip1[static_cast<std::size_t>(j + n_i)], inv_denom);
      });
      return;
    }
  }

  std::vector<FieldElement> denoms(static_cast<std::size_t>(n_i));
  ForEachIndexMaybeParallel(0, n_i, parallel_threshold, [&](long j) {
    const FieldElement &x1 = diag[static_cast<std::size_t>(j)];
    denoms[static_cast<std::size_t>(j)] = (params.zeta * x1) - x1;
  });

  std::vector<FieldElement> inv_denoms;
  if (!BatchInvertBaseUnits(inv_denoms, denoms)) {
    inv_denoms.resize(static_cast<std::size_t>(n_i));
    ForEachIndexMaybeParallel(0, n_i, parallel_threshold, [&](long j) {
      const FieldElement &denom = denoms[static_cast<std::size_t>(j)];
      if (denom == 0) {
        LogicError("ProverCommitRoundNoValidate: x1 must not equal x2");
      }
      if (!TryInvertBaseUnit(inv_denoms[static_cast<std::size_t>(j)], denom)) {
        LogicError("ProverCommitRoundNoValidate: x2-x1 must be a unit");
      }
    });
  }

  ForEachIndexMaybeParallel(0, n_i, parallel_threshold, [&](long j) {
    const FieldElement &x1 = diag[static_cast<std::size_t>(j)];
    pi_i[static_cast<std::size_t>(j)] = EvalLineAtWithInvDenom(
        alpha_i, x1, pi_ip1[static_cast<std::size_t>(j)],
        pi_ip1[static_cast<std::size_t>(j + n_i)],
        inv_denoms[static_cast<std::size_t>(j)]);
  });
}

}  // namespace basefold_pcs_internal
}  // namespace basefold
