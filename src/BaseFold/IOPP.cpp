#include "BaseFold/IOPP.hpp"
#include "BaseFold/Hash.hpp"
#include "BaseFold/Profile.hpp"

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/mat_ZZ_pE.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "GaloisRing/Inverse.hpp"

#if defined(BASEFOLD_USE_OPENMP)
#include <omp.h>
#endif

using NTL::BytesFromZZ;
using NTL::coeff;
using NTL::LogicError;
using NTL::mat_ZZ_pE;
using NTL::mul;
using NTL::NumBytes;
using NTL::rep;
using NTL::vec_ZZ_pE;
using NTL::ZZ;
using NTL::ZZ_pE;
using NTL::ZZ_pX;

namespace basefold {

thread_local Profile *g_active_profile = nullptr;

namespace {

constexpr std::size_t kDigestBytes = Digest{}.size();

long Pow2Checked(long e) {
  if (e < 0)
    LogicError("Pow2Checked: negative exponent");
  if (e >= static_cast<long>(8 * sizeof(long) - 1)) {
    LogicError("Pow2Checked: exponent too large for long");
  }
  return 1L << e;
}

long ParsePositiveEnvLong(const char *name, long fallback) {
  const char *raw = std::getenv(name);
  if (raw == nullptr || *raw == '\0') return fallback;
  char *end = nullptr;
  const long v = std::strtol(raw, &end, 10);
  if (end == raw || *end != '\0' || v <= 0) return fallback;
  return v;
}

int ParsePositiveEnvInt(const char *name, int fallback) {
  const char *raw = std::getenv(name);
  if (raw == nullptr || *raw == '\0') return fallback;
  char *end = nullptr;
  const long v = std::strtol(raw, &end, 10);
  if (end == raw || *end != '\0' || v <= 0 ||
      v > static_cast<long>(std::numeric_limits<int>::max())) {
    return fallback;
  }
  return static_cast<int>(v);
}

MerkleBuildParallelConfig ParseMerkleBuildParallelConfigFromEnv() {
  MerkleBuildParallelConfig cfg;
  cfg.leafs_per_thread = ParsePositiveEnvLong(
      "BASEFOLD_MERKLE_LEAFS_PER_THREAD", cfg.leafs_per_thread);
  cfg.parallel_level_threshold = ParsePositiveEnvLong(
      "BASEFOLD_MERKLE_PARALLEL_LEVEL_THRESHOLD",
      cfg.parallel_level_threshold);
  cfg.max_threads =
      ParsePositiveEnvInt("BASEFOLD_MERKLE_MAX_THREADS", cfg.max_threads);
  return cfg;
}

MerkleBuildParallelConfig &MutableMerkleBuildParallelConfig() {
  static MerkleBuildParallelConfig cfg = ParseMerkleBuildParallelConfigFromEnv();
  return cfg;
}

MerkleBuildParallelConfig LoadMerkleBuildParallelConfig() {
  MerkleBuildParallelConfig cfg = MutableMerkleBuildParallelConfig();
  if (cfg.leafs_per_thread <= 0) cfg.leafs_per_thread = 32768;
  if (cfg.parallel_level_threshold <= 0) cfg.parallel_level_threshold = 8192;
  if (cfg.max_threads <= 0) cfg.max_threads = 8;
  return cfg;
}

int ChooseMerkleBuildThreads(long leaf_count,
                             const MerkleBuildParallelConfig &cfg) {
#if defined(BASEFOLD_USE_OPENMP)
  if (leaf_count < cfg.leafs_per_thread) return 1;
  const int max_threads = omp_get_max_threads();
  int threads_to_use = static_cast<int>(leaf_count / cfg.leafs_per_thread);
  if (threads_to_use > cfg.max_threads) threads_to_use = cfg.max_threads;
  if (threads_to_use > max_threads) threads_to_use = max_threads;
  if (threads_to_use < 1) threads_to_use = 1;
  return threads_to_use;
#else
  (void)leaf_count;
  (void)cfg;
  return 1;
#endif
}

void ValidateParamsOrThrow(const FoldableCodeParams &params) {
  (void)MessageLength(params);
}

void WriteU64BE(Byte *out, std::uint64_t v) {
  for (int i = 7; i >= 0; --i) {
    out[7 - i] = static_cast<Byte>((v >> (8 * i)) & 0xff);
  }
}

void WriteUIntBE(Byte *out, std::uint64_t v, std::size_t bytes) {
  for (std::size_t i = 0; i < bytes; ++i) {
    out[bytes - 1 - i] = static_cast<Byte>((v >> (8 * i)) & 0xff);
  }
}

struct FieldEncodingInfo {
  long degree = 0;        // extension degree r
  long coeff_bytes = 0;   // fixed big-endian bytes per ZZ_p coefficient
};

FieldEncodingInfo CurrentFieldEncodingInfo() {
  const long r = ZZ_pE::degree();
  if (r <= 0) LogicError("CurrentFieldEncodingInfo: invalid extension degree");

  static thread_local ZZ cached_modulus;
  static thread_local long cached_r = 0;
  static thread_local long cached_coeff_bytes = 0;
  static thread_local bool has_cache = false;

  const ZZ &modulus = NTL::ZZ_p::modulus();
  if (!has_cache || r != cached_r || modulus != cached_modulus) {
    cached_modulus = modulus;
    cached_r = r;
    const ZZ max_coeff = modulus - 1;
    const long bytes = NumBytes(max_coeff);
    cached_coeff_bytes = (bytes > 0) ? bytes : 1;
    has_cache = true;
  }

  FieldEncodingInfo info;
  info.degree = cached_r;
  info.coeff_bytes = cached_coeff_bytes;
  return info;
}

bool IsUnitInBaseRing(const NTL::ZZ_p &a) {
  if (a == 0) return false;
  NTL::ZZ g;
  NTL::GCD(g, NTL::rep(a), NTL::ZZ_p::modulus());
  return g == 1;
}

bool BaseModulusIsPrime() {
  static thread_local ZZ cached_modulus;
  static thread_local bool cached_is_prime = false;
  static thread_local bool has_cache = false;

  const ZZ &modulus = NTL::ZZ_p::modulus();
  if (!has_cache || modulus != cached_modulus) {
    cached_modulus = modulus;
    cached_is_prime = NTL::ProbPrime(modulus);
    has_cache = true;
  }
  return cached_is_prime;
}

bool IsUnit(const ZZ_pE &a) {
  if (a == 0) return false;
  const long r = ZZ_pE::degree();
  if (r <= 0) LogicError("IsUnit: invalid extension degree");
  if (r == 1) {
    return IsUnitInBaseRing(coeff(rep(a), 0));
  }
  const NTL::ZZ_pX &poly = rep(a);
  for (long i = 0; i < r; ++i) {
    if (IsUnitInBaseRing(coeff(poly, i))) return true;
  }
  return false;
}

bool TryInvertUnit(ZZ_pE &inv_out, const ZZ_pE &a) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->try_invert_unit_ns : nullptr,
                    prof ? &prof->try_invert_unit_calls : nullptr);

  if (a == 0) return false;
  const long r = ZZ_pE::degree();
  if (r <= 0) LogicError("TryInvertUnit: invalid extension degree");

  // Context-aware single-entry cache: inversions are often repeated (e.g. when
  // diag_T is constant, all folding denominators match).
  static thread_local bool cache_valid = false;
  static thread_local ZZ cache_modulus;
  static thread_local long cache_r = 0;
  static thread_local ZZ_pE cache_a;
  static thread_local ZZ_pE cache_inv;

  const ZZ modulus = NTL::ZZ_p::modulus();
  if (cache_valid && cache_r == r && modulus == cache_modulus && a == cache_a) {
    if (prof != nullptr) ++prof->try_invert_unit_cache_hits;
    inv_out = cache_inv;
    return true;
  }

  ZZ_pE inv_candidate;
  bool ok = false;

  if (BaseModulusIsPrime()) {
    try {
      inv_candidate = NTL::inv(a);
      ok = true;
    } catch (...) {
      ok = false;
    }
    if (!ok) return false;
  } else {
    bool is_unit = false;
    if (prof != nullptr) {
      ScopedTimer is_unit_timer(&prof->is_unit_ns, &prof->is_unit_calls);
      is_unit = IsUnit(a);
    } else {
      is_unit = IsUnit(a);
    }
    if (!is_unit) return false;

    if (r == 1) {
      if (prof != nullptr) {
        ScopedTimer inv_timer(&prof->inv_fallback_ns,
                              &prof->inv_fallback_calls);
        const NTL::ZZ a_rep = NTL::rep(coeff(rep(a), 0));
        NTL::ZZ inv_rep;
        if (NTL::InvModStatus(inv_rep, a_rep, NTL::ZZ_p::modulus()) != 0)
          return false;
        NTL::ZZ_p inv_base;
        NTL::conv(inv_base, inv_rep);
        NTL::ZZ_pX poly;
        NTL::clear(poly);
        NTL::SetCoeff(poly, 0, inv_base);
        NTL::conv(inv_candidate, poly);
        ok = true;
      } else {
        const NTL::ZZ a_rep = NTL::rep(coeff(rep(a), 0));
        NTL::ZZ inv_rep;
        if (NTL::InvModStatus(inv_rep, a_rep, NTL::ZZ_p::modulus()) != 0)
          return false;
        NTL::ZZ_p inv_base;
        NTL::conv(inv_base, inv_rep);
        NTL::ZZ_pX poly;
        NTL::clear(poly);
        NTL::SetCoeff(poly, 0, inv_base);
        NTL::conv(inv_candidate, poly);
        ok = true;
      }
    } else {
      if (prof != nullptr) {
        ScopedTimer inv_timer(&prof->inv_fallback_ns,
                              &prof->inv_fallback_calls);
        inv_candidate = Inv(a, r);
      } else {
        inv_candidate = Inv(a, r);
      }
      if (inv_candidate == 0) return false;
      ZZ_pE one;
      NTL::set(one);
      ok = (a * inv_candidate == one);
    }

    if (!ok) return false;
  }

  inv_out = inv_candidate;

  cache_valid = true;
  cache_modulus = modulus;
  cache_r = r;
  cache_a = a;
  cache_inv = inv_candidate;

  return true;
}

Digest HashWithPrefix(Byte prefix) {
  return HashDigest(&prefix, 1, "HashWithPrefix");
}

Digest HashNode(const Digest &left, const Digest &right) {
  std::array<Byte, 1 + 2 * kDigestBytes> in;
  in[0] = static_cast<Byte>(0x01);
  std::memcpy(in.data() + 1, left.data(), left.size());
  std::memcpy(in.data() + 1 + left.size(), right.data(), right.size());
  return HashDigest(in.data(), in.size(), "HashNode");
}

Digest HashLeaf(long index, const FieldElement &value) {
  const FieldEncodingInfo enc = CurrentFieldEncodingInfo();
  const std::size_t coeff_bytes = static_cast<std::size_t>(enc.coeff_bytes);
  const std::size_t payload_bytes =
      1 + 8 + static_cast<std::size_t>(enc.degree) * coeff_bytes;

  static thread_local std::vector<Byte> in;
  if (in.size() != payload_bytes) {
    in.resize(payload_bytes);
  }
  in[0] = static_cast<Byte>(0x00);
  WriteU64BE(in.data() + 1, static_cast<std::uint64_t>(index));

  const ZZ_pX &poly = rep(value);
  std::size_t off = 1 + 8;
  for (long i = 0; i < enc.degree; ++i) {
    if (enc.coeff_bytes <= 8) {
      unsigned long c_ul = 0;
      NTL::conv(c_ul, rep(coeff(poly, i)));
      WriteUIntBE(in.data() + off, static_cast<std::uint64_t>(c_ul),
                  coeff_bytes);
    } else {
      const ZZ &c = rep(coeff(poly, i));
      BytesFromZZ(reinterpret_cast<unsigned char *>(in.data() + off), c,
                  enc.coeff_bytes);
    }
    off += coeff_bytes;
  }

  return HashDigest(in.data(), in.size(), "HashLeaf");
}

Digest HashRootWithCount(long leaf_count, const Digest &raw_root) {
  std::array<Byte, 1 + 8 + kDigestBytes> in;
  in[0] = static_cast<Byte>(0x03);
  WriteU64BE(in.data() + 1, static_cast<std::uint64_t>(leaf_count));
  std::memcpy(in.data() + 1 + 8, raw_root.data(), raw_root.size());
  return HashDigest(in.data(), in.size(), "HashRootWithCount");
}

std::size_t ExpectedMerkleHeight(long leaf_count) {
  if (leaf_count <= 0)
    return 0;
  std::size_t height = 0;
  long n = leaf_count;
  while (n > 1) {
    if (n & 1L)
      n += 1;
    n /= 2;
    ++height;
  }
  return height;
}

struct MerkleMultiproofPlanLevel {
  std::vector<long> node_indices;
  std::vector<long> sibling_indices;
  std::vector<long> parent_indices;
};

struct MerkleMultiproofPlan {
  long tree_leaf_count = 0;
  long padded_leaf_count = 0;
  std::vector<long> queried_indices;
  std::vector<MerkleMultiproofPlanLevel> levels;
  MerkleMultiproofStats stats;
};

long NextPowerOfTwoLong(long x) {
  if (x <= 1)
    return 1;
  long value = 1;
  while (value < x) {
    if (value > std::numeric_limits<long>::max() / 2) {
      LogicError("NextPowerOfTwoLong: overflow");
    }
    value *= 2;
  }
  return value;
}

bool IsSortedUniqueMerkleIndicesInRange(long leaf_count,
                                        const std::vector<long> &indices) {
  if (leaf_count < 0) {
    return false;
  }
  long prev = -1;
  for (long index : indices) {
    if (index < 0 || index >= leaf_count) {
      return false;
    }
    if (prev >= 0 && index <= prev) {
      return false;
    }
    prev = index;
  }
  return true;
}

std::vector<long> SortAndValidateMerkleIndicesOrThrow(
    long leaf_count, const std::vector<long> &queried_indices,
    const char *func_name) {
  if (leaf_count < 0) {
    const std::string msg = std::string(func_name) + ": invalid leaf count";
    LogicError(msg.c_str());
  }
  std::vector<long> unique = queried_indices;
  std::sort(unique.begin(), unique.end());
  unique.erase(std::unique(unique.begin(), unique.end()), unique.end());
  for (long index : unique) {
    if (index < 0 || index >= leaf_count) {
      const std::string msg = std::string(func_name) + ": index out of range";
      LogicError(msg.c_str());
    }
  }
  return unique;
}

MerkleMultiproofPlan BuildMerkleMultiproofPlanFromSortedUnique(
    long leaf_count, const std::vector<long> &queried_indices) {
  MerkleMultiproofPlan plan;
  plan.tree_leaf_count = leaf_count;
  if (leaf_count <= 0) {
    return plan;
  }

  plan.padded_leaf_count = NextPowerOfTwoLong(leaf_count);
  plan.queried_indices = queried_indices;
  plan.stats.opened_leaf_count =
      static_cast<std::uint64_t>(queried_indices.size());

  if (queried_indices.empty()) {
    return plan;
  }

  std::vector<long> current = queried_indices;
  long level_width = plan.padded_leaf_count;
  while (level_width > 1) {
    MerkleMultiproofPlanLevel level;
    level.node_indices = current;
    level.sibling_indices.reserve((current.size() + 1U) / 2U);
    level.parent_indices.reserve((current.size() + 1U) / 2U);

    for (std::size_t i = 0; i < current.size();) {
      const long node = current[i];
      const bool has_paired_sibling =
          (i + 1U < current.size()) && ((node & 1L) == 0L) &&
          (current[i + 1U] == node + 1L);
      if (!has_paired_sibling) {
        level.sibling_indices.push_back(node ^ 1L);
      }
      level.parent_indices.push_back(node / 2L);
      i += has_paired_sibling ? 2U : 1U;
    }

    plan.stats.unique_sibling_count +=
        static_cast<std::uint64_t>(level.sibling_indices.size());
    current = level.parent_indices;
    plan.levels.push_back(std::move(level));
    level_width /= 2L;
  }

  plan.stats.verifier_hashes = plan.stats.unique_sibling_count;
  return plan;
}

Digest MerkleRootRaw(std::vector<Digest> level) {
  if (level.empty())
    return HashWithPrefix(static_cast<Byte>(0x04));

  while (level.size() > 1) {
    if (level.size() % 2 == 1)
      level.push_back(level.back());
    std::vector<Digest> next;
    next.reserve(level.size() / 2);
    for (std::size_t i = 0; i < level.size(); i += 2) {
      next.push_back(HashNode(level[i], level[i + 1]));
    }
    level = std::move(next);
  }
  return level[0];
}

bool SolveLinearSystemRref(vec_ZZ_pE &x_out, mat_ZZ_pE &aug) {
  const long m = aug.NumRows();
  const long n_plus_1 = aug.NumCols();
  if (n_plus_1 <= 0)
    LogicError("SolveLinearSystemRref: empty system");
  const long n = n_plus_1 - 1;

  x_out.SetLength(n);
  for (long i = 0; i < n; ++i)
    x_out[i] = ZZ_pE(0);

  std::vector<long> pivot_col_orig;
  pivot_col_orig.resize(static_cast<std::size_t>(m), -1);

  std::vector<long> col_perm;
  col_perm.resize(static_cast<std::size_t>(n));
  for (long c = 0; c < n; ++c)
    col_perm[static_cast<std::size_t>(c)] = c;

  long row = 0;
  for (long col = 0; col < n && row < m; ++col) {
    long pivot_row = -1;
    long pivot_col = -1;
    ZZ_pE inv_pivot;

    for (long c = col; c < n; ++c) {
      for (long r = row; r < m; ++r) {
        if (TryInvertUnit(inv_pivot, aug[r][c])) {
          pivot_row = r;
          pivot_col = c;
          break;
        }
      }
      if (pivot_row >= 0)
        break;
    }

    if (pivot_row < 0) {
      break;
    }

    if (pivot_col != col) {
      for (long r = 0; r < m; ++r) {
        std::swap(aug[r][col], aug[r][pivot_col]);
      }
      std::swap(col_perm[static_cast<std::size_t>(col)],
                col_perm[static_cast<std::size_t>(pivot_col)]);
    }

    if (pivot_row != row) {
      for (long c = 0; c < n_plus_1; ++c) {
        std::swap(aug[row][c], aug[pivot_row][c]);
      }
    }

    for (long c = col; c < n_plus_1; ++c) {
      aug[row][c] *= inv_pivot;
    }

    for (long r = 0; r < m; ++r) {
      if (r == row)
        continue;
      if (aug[r][col] == 0)
        continue;
      const ZZ_pE factor = aug[r][col];
      for (long c = col; c < n_plus_1; ++c) {
        aug[r][c] -= factor * aug[row][c];
      }
    }

    pivot_col_orig[static_cast<std::size_t>(row)] =
        col_perm[static_cast<std::size_t>(col)];
    ++row;
  }

  for (long r = 0; r < m; ++r) {
    bool all_zero = true;
    for (long c = 0; c < n; ++c) {
      if (aug[r][c] != 0) {
        all_zero = false;
        break;
      }
    }
    if (all_zero && aug[r][n] != 0)
      return false;
  }

  for (long r = 0; r < m; ++r) {
    const long pc_orig = pivot_col_orig[static_cast<std::size_t>(r)];
    if (pc_orig >= 0)
      x_out[pc_orig] = aug[r][n];
  }

  return true;
}

} // namespace

void ResetMerkleBuildParallelConfigFromEnv() {
  MutableMerkleBuildParallelConfig() = ParseMerkleBuildParallelConfigFromEnv();
}

void SetMerkleBuildParallelConfig(const MerkleBuildParallelConfig &cfg) {
  MutableMerkleBuildParallelConfig() = cfg;
}

MerkleBuildParallelConfig GetMerkleBuildParallelConfig() {
  return LoadMerkleBuildParallelConfig();
}

long MessageLengthAtLevel(const FoldableCodeParams &params, long level) {
  ValidateParamsOrThrow(params);
  if (level < 0 || level > params.d) {
    LogicError("MessageLengthAtLevel: level out of range");
  }
  const long pow2 = Pow2Checked(level);
  if (params.k0 > std::numeric_limits<long>::max() / pow2) {
    LogicError("MessageLengthAtLevel: overflow");
  }
  return params.k0 * pow2;
}

long CodewordLengthAtLevel(const FoldableCodeParams &params, long level) {
  ValidateParamsOrThrow(params);
  const long k = MessageLengthAtLevel(params, level);
  if (params.c > std::numeric_limits<long>::max() / k) {
    LogicError("CodewordLengthAtLevel: overflow");
  }
  return params.c * k;
}

void FoldingPoints(FieldElement &x_left, FieldElement &x_right,
                   const FoldableCodeParams &params, long level_i, long j) {
  ValidateParamsOrThrow(params);
  if (level_i < 0 || level_i >= params.d) {
    LogicError("FoldingPoints: level out of range");
  }
  const long n_i = CodewordLengthAtLevel(params, level_i);
  if (j < 0 || j >= n_i)
    LogicError("FoldingPoints: index out of range");

  const FieldElement &t = params.diag_T[static_cast<std::size_t>(level_i)][j];
  x_left = t;
  x_right = params.zeta * t;
}

FieldElement EvalLineAt(const FieldElement &x, const FieldElement &x1,
                        const FieldElement &y1, const FieldElement &x2,
                        const FieldElement &y2) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->eval_line_at_ns : nullptr,
                    prof ? &prof->eval_line_at_calls : nullptr);

  const FieldElement denom = x2 - x1;
  if (denom == 0)
    LogicError("EvalLineAt: x1 must not equal x2");
  ZZ_pE inv_denom;
  if (!TryInvertUnit(inv_denom, denom)) {
    LogicError("EvalLineAt: x2-x1 must be a unit");
  }
  return y1 + (x - x1) * (y2 - y1) * inv_denom;
}

void ProverCommitRound(Oracle &pi_i, const Oracle &pi_ip1,
                       const FieldElement &alpha_i, long level_i,
                       const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  if (level_i < 0 || level_i >= params.d) {
    LogicError("ProverCommitRound: level out of range");
  }

  const long n_i = CodewordLengthAtLevel(params, level_i);
  if (pi_ip1.length() != 2 * n_i) {
    LogicError("ProverCommitRound: pi_{i+1} has wrong length");
  }

  pi_i.SetLength(n_i);
  for (long j = 0; j < n_i; ++j) {
    FieldElement x1, x2;
    FoldingPoints(x1, x2, params, level_i, j);
    const FieldElement &y1 = pi_ip1[j];
    const FieldElement &y2 = pi_ip1[j + n_i];
    pi_i[j] = EvalLineAt(alpha_i, x1, y1, x2, y2);
  }
}

void ProverCommitAll(IOPPOracles &oracles, const Oracle &pi_d,
                     const IOPPChallenges &challenges,
                     const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  if (static_cast<long>(challenges.alphas.size()) != params.d) {
    LogicError("ProverCommitAll: challenges.alphas has wrong size");
  }

  const long n_d = CodewordLengthAtLevel(params, params.d);
  if (pi_d.length() != n_d)
    LogicError("ProverCommitAll: pi_d has wrong length");

  oracles.pi.resize(static_cast<std::size_t>(params.d + 1));
  oracles.pi[static_cast<std::size_t>(params.d)] = pi_d;

  for (long i = params.d; i-- > 0;) {
    ProverCommitRound(oracles.pi[static_cast<std::size_t>(i)],
                      oracles.pi[static_cast<std::size_t>(i + 1)],
                      challenges.alphas[static_cast<std::size_t>(i)], i,
                      params);
  }
}

IOPPQueryPlan MakeQueryPlan(long initial_mu, const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  IOPPQueryPlan plan;
  plan.initial_mu = initial_mu;
  plan.mu_by_level.resize(static_cast<std::size_t>(params.d));

  if (params.d == 0)
    return plan;

  const long n_last = CodewordLengthAtLevel(params, params.d - 1);
  if (initial_mu < 0 || initial_mu >= n_last) {
    LogicError("MakeQueryPlan: initial_mu out of range");
  }

  long mu = initial_mu;
  for (long i = params.d; i-- > 0;) {
    plan.mu_by_level[static_cast<std::size_t>(i)] = mu;
    if (i > 0) {
      const long n_prev = CodewordLengthAtLevel(params, i - 1);
      if (mu >= n_prev)
        mu -= n_prev;
    }
  }

  return plan;
}

bool VerifyQueryFromOpenings(const IOPPQueryPlan &plan,
                             const IOPPChallenges &challenges,
                             const IOPPQueryOpenings &openings,
                             const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  if (static_cast<long>(challenges.alphas.size()) != params.d)
    return false;
  if (static_cast<long>(plan.mu_by_level.size()) != params.d)
    return false;
  if (static_cast<long>(openings.left.size()) != params.d)
    return false;
  if (static_cast<long>(openings.right.size()) != params.d)
    return false;
  if (static_cast<long>(openings.folded.size()) != params.d)
    return false;

  for (long i = params.d; i-- > 0;) {
    const long mu = plan.mu_by_level[static_cast<std::size_t>(i)];
    const long n_i = CodewordLengthAtLevel(params, i);
    if (mu < 0 || mu >= n_i)
      return false;

    FieldElement x1, x2;
    FoldingPoints(x1, x2, params, i, mu);
    const FieldElement &y1 = openings.left[static_cast<std::size_t>(i)];
    const FieldElement &y2 = openings.right[static_cast<std::size_t>(i)];
    const FieldElement expected = EvalLineAt(
        challenges.alphas[static_cast<std::size_t>(i)], x1, y1, x2, y2);
    if (expected != openings.folded[static_cast<std::size_t>(i)])
      return false;
  }

  return IsCodewordC0(openings.pi0_full, params);
}

bool VerifyQueryFromOracles(const IOPPQueryPlan &plan,
                            const IOPPChallenges &challenges,
                            const IOPPOracles &oracles,
                            const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  if (static_cast<long>(challenges.alphas.size()) != params.d)
    return false;
  if (static_cast<long>(plan.mu_by_level.size()) != params.d)
    return false;
  if (static_cast<long>(oracles.pi.size()) != params.d + 1)
    return false;

  for (long i = 0; i <= params.d; ++i) {
    if (oracles.pi[static_cast<std::size_t>(i)].length() !=
        CodewordLengthAtLevel(params, i)) {
      return false;
    }
  }

  for (long i = params.d; i-- > 0;) {
    const long mu = plan.mu_by_level[static_cast<std::size_t>(i)];
    const long n_i = CodewordLengthAtLevel(params, i);
    if (mu < 0 || mu >= n_i)
      return false;

    FieldElement x1, x2;
    FoldingPoints(x1, x2, params, i, mu);
    const FieldElement &y1 = oracles.pi[static_cast<std::size_t>(i + 1)][mu];
    const FieldElement &y2 =
        oracles.pi[static_cast<std::size_t>(i + 1)][mu + n_i];
    const FieldElement expected = EvalLineAt(
        challenges.alphas[static_cast<std::size_t>(i)], x1, y1, x2, y2);
    if (expected != oracles.pi[static_cast<std::size_t>(i)][mu])
      return false;
  }

  return IsCodewordC0(oracles.pi[0], params);
}

MerkleRoot MerkleCommitOracle(const Oracle &oracle) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->merkle_commit_oracle_ns : nullptr,
                    prof ? &prof->merkle_commit_oracle_calls : nullptr);

  const long leaf_count = oracle.length();
  if (leaf_count < 0)
    LogicError("MerkleCommitOracle: invalid leaf count");
  if (leaf_count == 0)
    return HashRootWithCount(0, HashWithPrefix(static_cast<Byte>(0x04)));

  std::vector<Digest> leaf_hashes;
  leaf_hashes.reserve(static_cast<std::size_t>(leaf_count));
  for (long i = 0; i < leaf_count; ++i) {
    leaf_hashes.push_back(HashLeaf(i, oracle[i]));
  }
  const Digest raw_root = MerkleRootRaw(std::move(leaf_hashes));
  return HashRootWithCount(leaf_count, raw_root);
}

MerkleOpening MerkleOpenOracle(const Oracle &oracle, long index) {
  const long leaf_count = oracle.length();
  if (index < 0 || index >= leaf_count) {
    LogicError("MerkleOpenOracle: index out of range");
  }

  std::vector<Digest> level;
  level.reserve(static_cast<std::size_t>(leaf_count));
  for (long i = 0; i < leaf_count; ++i) {
    level.push_back(HashLeaf(i, oracle[i]));
  }

  MerkleAuthPath path;
  const std::size_t height = ExpectedMerkleHeight(leaf_count);
  path.sibling_hashes.reserve(height);

  long idx = index;
  while (level.size() > 1) {
    if (level.size() % 2 == 1)
      level.push_back(level.back());

    const long sibling = (idx % 2 == 0) ? (idx + 1) : (idx - 1);
    path.sibling_hashes.push_back(level[static_cast<std::size_t>(sibling)]);

    std::vector<Digest> next;
    next.reserve(level.size() / 2);
    for (std::size_t i = 0; i < level.size(); i += 2) {
      next.push_back(HashNode(level[i], level[i + 1]));
    }
    idx /= 2;
    level = std::move(next);
  }

  MerkleOpening opening;
  opening.index = index;
  opening.value = oracle[index];
  opening.auth_path = std::move(path);
  return opening;
}

MerkleMultiproof MerkleOpenOracleMany(const Oracle &oracle,
                                      const std::vector<long> &queried_indices) {
  const MerkleTree tree = MerkleTree::Build(oracle);
  return tree.OpenMany(oracle, queried_indices);
}

bool MerkleVerifyOpening(const MerkleRoot &root, long leaf_count,
                         const MerkleOpening &opening) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->merkle_verify_opening_ns : nullptr,
                    prof ? &prof->merkle_verify_opening_calls : nullptr);

  if (leaf_count < 0)
    return false;
  if (opening.index < 0 || opening.index >= leaf_count)
    return false;

  const std::size_t expected_height = ExpectedMerkleHeight(leaf_count);
  if (opening.auth_path.sibling_hashes.size() != expected_height)
    return false;

  Digest cur = HashLeaf(opening.index, opening.value);
  long idx = opening.index;
  for (std::size_t i = 0; i < expected_height; ++i) {
    const Digest &sib = opening.auth_path.sibling_hashes[i];
    if (idx % 2 == 1) {
      cur = HashNode(sib, cur);
    } else {
      cur = HashNode(cur, sib);
    }
    idx /= 2;
  }

  const Digest expected_root = HashRootWithCount(leaf_count, cur);
  return expected_root == root;
}

MerkleMultiproofStats PlanMerkleMultiproof(
    long leaf_count, const std::vector<long> &queried_indices) {
  const std::vector<long> unique = SortAndValidateMerkleIndicesOrThrow(
      leaf_count, queried_indices, "PlanMerkleMultiproof");
  return BuildMerkleMultiproofPlanFromSortedUnique(leaf_count, unique).stats;
}

bool MerkleVerifyMultiproof(const MerkleRoot &root, long leaf_count,
                            const MerkleMultiproof &proof) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->merkle_verify_opening_ns : nullptr,
                    prof ? &prof->merkle_verify_opening_calls : nullptr);

  if (leaf_count < 0) {
    return false;
  }
  if (static_cast<long>(proof.values.length()) !=
      static_cast<long>(proof.queried_indices.size())) {
    return false;
  }

  if (leaf_count == 0) {
    if (!proof.queried_indices.empty() || proof.values.length() != 0 ||
        !proof.sibling_hashes.empty()) {
      return false;
    }
    return root == HashRootWithCount(0, HashWithPrefix(static_cast<Byte>(0x04)));
  }

  if (!IsSortedUniqueMerkleIndicesInRange(leaf_count, proof.queried_indices)) {
    return false;
  }

  const MerkleMultiproofPlan plan =
      BuildMerkleMultiproofPlanFromSortedUnique(leaf_count, proof.queried_indices);
  if (proof.sibling_hashes.size() !=
      static_cast<std::size_t>(plan.stats.unique_sibling_count)) {
    return false;
  }

  if (proof.queried_indices.empty()) {
    return proof.sibling_hashes.empty();
  }

  std::vector<std::pair<long, Digest>> current;
  current.resize(proof.queried_indices.size());
  for (std::size_t i = 0; i < proof.queried_indices.size(); ++i) {
    current[i] = {proof.queried_indices[i],
                  HashLeaf(proof.queried_indices[i], proof.values[static_cast<long>(i)])};
  }

  std::size_t sibling_cursor = 0;
  for (const MerkleMultiproofPlanLevel &level : plan.levels) {
    if (current.size() != level.node_indices.size()) {
      return false;
    }
    for (std::size_t i = 0; i < current.size(); ++i) {
      if (current[i].first != level.node_indices[i]) {
        return false;
      }
    }

    std::vector<std::pair<long, Digest>> parents;
    parents.resize(level.parent_indices.size());
    std::size_t parent_count = 0;
    for (std::size_t i = 0; i < current.size();) {
      const long index = current[i].first;
      const Digest &hash = current[i].second;
      const bool has_adjacent =
          (i + 1U < current.size()) && ((index & 1L) == 0L) &&
          (current[i + 1U].first == index + 1L);
      if (parent_count >= parents.size()) {
        return false;
      }

      const long parent_index = index / 2L;
      if (parent_index != level.parent_indices[parent_count]) {
        return false;
      }

      Digest parent_hash;
      if (has_adjacent) {
        parent_hash = HashNode(hash, current[i + 1U].second);
        i += 2U;
      } else {
        if (sibling_cursor >= proof.sibling_hashes.size()) {
          return false;
        }
        const Digest &sibling = proof.sibling_hashes[sibling_cursor++];
        parent_hash = ((index & 1L) == 0L) ? HashNode(hash, sibling)
                                           : HashNode(sibling, hash);
        ++i;
      }

      parents[parent_count] = {parent_index, parent_hash};
      ++parent_count;
    }
    if (parent_count != parents.size()) {
      return false;
    }
    current = std::move(parents);
  }

  if (sibling_cursor != proof.sibling_hashes.size() || current.size() != 1U) {
    return false;
  }
  return HashRootWithCount(leaf_count, current.front().second) == root;
}

MerkleTree MerkleTree::Build(const Oracle &oracle) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->merkle_tree_build_ns : nullptr,
                    prof ? &prof->merkle_tree_build_calls : nullptr);

  const long leaf_count = oracle.length();
  if (leaf_count < 0)
    LogicError("MerkleTree::Build: invalid leaf count");

  MerkleTree t;
  t.leaf_count_ = leaf_count;
  t.levels_.clear();
  t.levels_.reserve(ExpectedMerkleHeight(leaf_count));
  t.raw_root_ = Digest{};

  if (leaf_count == 0) {
    t.raw_root_ = HashWithPrefix(static_cast<Byte>(0x04));
    return t;
  }

  std::vector<Digest> level;
  level.resize(static_cast<std::size_t>(leaf_count));
  const MerkleBuildParallelConfig merkle_cfg = LoadMerkleBuildParallelConfig();
#if defined(BASEFOLD_USE_OPENMP)
  const int threads_to_use = ChooseMerkleBuildThreads(leaf_count, merkle_cfg);
  if (threads_to_use >= 2) {
    const ZZ base_modulus = NTL::ZZ_p::modulus();
    const ZZ_pX extension_modulus = NTL::ZZ_pE::modulus().val();
    std::vector<Digest> next;
    long next_size = 0;
    bool done = false;
    bool parallel_level = false;
    const long parallel_level_threshold = merkle_cfg.parallel_level_threshold;

#pragma omp parallel num_threads(threads_to_use) shared(level, next, next_size, done, parallel_level, t)
    {
      NTL::ZZ_p::init(base_modulus);
      NTL::ZZ_pE::init(extension_modulus);

#pragma omp for schedule(static)
      for (long i = 0; i < leaf_count; ++i) {
        level[static_cast<std::size_t>(i)] = HashLeaf(i, oracle[i]);
      }

      while (true) {
#pragma omp single
        {
          done = (level.size() <= 1);
          if (!done) {
            if (level.size() % 2 == 1) {
              level.push_back(level.back());
            }
            next_size = static_cast<long>(level.size() / 2);
            next.resize(static_cast<std::size_t>(next_size));
            parallel_level = (next_size >= parallel_level_threshold);
          }
        }
        if (done) break;

        if (parallel_level) {
#pragma omp for schedule(static)
          for (long i = 0; i < next_size; ++i) {
            const std::size_t j = static_cast<std::size_t>(2 * i);
            next[static_cast<std::size_t>(i)] = HashNode(level[j], level[j + 1]);
          }
        } else {
#pragma omp single
          {
            for (long i = 0; i < next_size; ++i) {
              const std::size_t j = static_cast<std::size_t>(2 * i);
              next[static_cast<std::size_t>(i)] = HashNode(level[j], level[j + 1]);
            }
          }
        }

#pragma omp single
        {
          t.levels_.push_back(std::move(level));
          level = std::move(next);
        }
      }
    }
  } else
#endif
  {
    for (long i = 0; i < leaf_count; ++i) {
      level[static_cast<std::size_t>(i)] = HashLeaf(i, oracle[i]);
    }
    while (level.size() > 1) {
      if (level.size() % 2 == 1)
        level.push_back(level.back());

      std::vector<Digest> next;
      const long next_size = static_cast<long>(level.size() / 2);
      next.resize(static_cast<std::size_t>(next_size));
      for (long i = 0; i < next_size; ++i) {
        const std::size_t j = static_cast<std::size_t>(2 * i);
        next[static_cast<std::size_t>(i)] = HashNode(level[j], level[j + 1]);
      }
      t.levels_.push_back(std::move(level));
      level = std::move(next);
    }
  }

  t.raw_root_ = level[0];
  return t;
}

MerkleRoot MerkleTree::Root() const {
  return HashRootWithCount(leaf_count_, raw_root_);
}

MerkleOpening MerkleTree::Open(const Oracle &oracle, long index) const {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->merkle_tree_open_ns : nullptr,
                    prof ? &prof->merkle_tree_open_calls : nullptr);

  if (oracle.length() != leaf_count_) {
    LogicError("MerkleTree::Open: oracle length mismatch");
  }
  if (index < 0 || index >= leaf_count_) {
    LogicError("MerkleTree::Open: index out of range");
  }

  MerkleAuthPath path;
  const std::size_t height = ExpectedMerkleHeight(leaf_count_);
  path.sibling_hashes.reserve(height);

  long idx = index;
  for (std::size_t h = 0; h < height; ++h) {
    const std::vector<Digest> &level = levels_[h];
    const long sibling = (idx % 2 == 0) ? (idx + 1) : (idx - 1);
    path.sibling_hashes.push_back(level[static_cast<std::size_t>(sibling)]);
    idx /= 2;
  }

  MerkleOpening opening;
  opening.index = index;
  opening.value = oracle[index];
  opening.auth_path = std::move(path);
  return opening;
}

MerkleMultiproof MerkleTree::OpenMany(
    const Oracle &oracle, const std::vector<long> &queried_indices) const {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->merkle_tree_open_ns : nullptr,
                    prof ? &prof->merkle_tree_open_calls : nullptr);

  if (oracle.length() != leaf_count_) {
    LogicError("MerkleTree::OpenMany: oracle length mismatch");
  }
  const std::vector<long> unique = SortAndValidateMerkleIndicesOrThrow(
      leaf_count_, queried_indices, "MerkleTree::OpenMany");

  MerkleMultiproof proof;
  proof.queried_indices = unique;
  proof.values.SetLength(static_cast<long>(unique.size()));
  if (unique.empty()) {
    return proof;
  }

  const MerkleMultiproofPlan plan =
      BuildMerkleMultiproofPlanFromSortedUnique(leaf_count_, unique);
  for (std::size_t i = 0; i < unique.size(); ++i) {
    proof.values[static_cast<long>(i)] =
        oracle[static_cast<long>(unique[i])];
  }
  proof.sibling_hashes.reserve(
      static_cast<std::size_t>(plan.stats.unique_sibling_count));
  for (std::size_t level_index = 0; level_index < plan.levels.size();
       ++level_index) {
    const std::vector<Digest> &level = levels_[level_index];
    for (long sibling : plan.levels[level_index].sibling_indices) {
      proof.sibling_hashes.push_back(level[static_cast<std::size_t>(sibling)]);
    }
  }
  return proof;
}

IOPPChallenges
FiatShamirDeriveChallenges(FiatShamirTranscript &transcript,
                           const IOPPMerkleCommitments &commitments,
                           const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  if (static_cast<long>(commitments.roots_by_level.size()) != params.d + 1) {
    LogicError("FiatShamirDeriveChallenges: roots_by_level has wrong size");
  }

  IOPPChallenges out;
  out.alphas.resize(static_cast<std::size_t>(params.d));
  if (params.d == 0)
    return out;

  {
    const MerkleRoot &root_d =
        commitments.roots_by_level[static_cast<std::size_t>(params.d)];
    transcript.AbsorbBytes(root_d.data(), root_d.size());
  }
  for (long i = params.d; i-- > 0;) {
    out.alphas[static_cast<std::size_t>(i)] =
        transcript.ChallengeFieldElement("alpha/" + std::to_string(i));
    const MerkleRoot &root_i =
        commitments.roots_by_level[static_cast<std::size_t>(i)];
    transcript.AbsorbBytes(root_i.data(), root_i.size());
  }

  return out;
}

std::vector<IOPPQueryPlan>
FiatShamirDeriveQueryPlans(FiatShamirTranscript &transcript, long num_queries,
                           const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  if (num_queries < 0)
    LogicError("FiatShamirDeriveQueryPlans: negative count");

  std::vector<IOPPQueryPlan> plans;
  plans.reserve(static_cast<std::size_t>(num_queries));

  if (params.d == 0) {
    for (long q = 0; q < num_queries; ++q)
      plans.push_back(MakeQueryPlan(0, params));
    return plans;
  }

  const long n_last = CodewordLengthAtLevel(params, params.d - 1);
  for (long q = 0; q < num_queries; ++q) {
    const long mu =
        transcript.ChallengeIndex("mu/" + std::to_string(q), n_last);
    plans.push_back(MakeQueryPlan(mu, params));
  }

  return plans;
}

bool VerifyQueryFromMerkleOpenings(const IOPPQueryPlan &plan,
                                   const IOPPChallenges &challenges,
                                   const IOPPQueryMerkleOpenings &openings,
                                   const IOPPMerkleCommitments &commitments,
                                   const FoldableCodeParams &params) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->verify_query_merkle_ns : nullptr,
                    prof ? &prof->verify_query_merkle_calls : nullptr);

  ValidateParamsOrThrow(params);
  if (static_cast<long>(commitments.roots_by_level.size()) != params.d + 1) {
    return false;
  }
  if (static_cast<long>(challenges.alphas.size()) != params.d)
    return false;
  if (static_cast<long>(plan.mu_by_level.size()) != params.d)
    return false;
  if (static_cast<long>(openings.left.size()) != params.d)
    return false;
  if (static_cast<long>(openings.right.size()) != params.d)
    return false;
  if (static_cast<long>(openings.folded.size()) != params.d)
    return false;

  if (MerkleCommitOracle(openings.pi0_full) != commitments.roots_by_level[0]) {
    return false;
  }

  for (long i = params.d; i-- > 0;) {
    const long mu = plan.mu_by_level[static_cast<std::size_t>(i)];
    const long n_i = CodewordLengthAtLevel(params, i);
    if (mu < 0 || mu >= n_i)
      return false;

    const MerkleOpening &left = openings.left[static_cast<std::size_t>(i)];
    const MerkleOpening &right = openings.right[static_cast<std::size_t>(i)];
    const MerkleOpening &folded = openings.folded[static_cast<std::size_t>(i)];

    if (left.index != mu)
      return false;
    if (right.index != mu + n_i)
      return false;
    if (folded.index != mu)
      return false;

    const long n_ip1 = 2 * n_i;
    if (!MerkleVerifyOpening(
            commitments.roots_by_level[static_cast<std::size_t>(i + 1)], n_ip1,
            left)) {
      return false;
    }
    if (!MerkleVerifyOpening(
            commitments.roots_by_level[static_cast<std::size_t>(i + 1)], n_ip1,
            right)) {
      return false;
    }
    if (!MerkleVerifyOpening(
            commitments.roots_by_level[static_cast<std::size_t>(i)], n_i,
            folded)) {
      return false;
    }

    FieldElement x1, x2;
    FoldingPoints(x1, x2, params, i, mu);
    const FieldElement expected =
        EvalLineAt(challenges.alphas[static_cast<std::size_t>(i)], x1,
                   left.value, x2, right.value);
    if (expected != folded.value)
      return false;
  }

  return IsCodewordC0(openings.pi0_full, params);
}

bool IsCodewordC0(const Oracle &pi0, const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  const long n0 = CodewordLengthAtLevel(params, 0);
  if (pi0.length() != n0) return false;

  vec_ZZ_pE msg0;
  return DecodeC0(msg0, pi0, params);
}

bool DecodeC0(vec_ZZ_pE &msg0_out, const Oracle &pi0,
              const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  const long n0 = CodewordLengthAtLevel(params, 0);
  if (pi0.length() != n0) return false;

  if (params.k0 == 1) {
    ZZ_pE inv_g;
    long pivot = -1;
    for (long j = 0; j < n0; ++j) {
      if (TryInvertUnit(inv_g, params.G0[0][j])) {
        pivot = j;
        break;
      }
    }
    if (pivot < 0) return false;

    msg0_out.SetLength(1);
    msg0_out[0] = pi0[pivot] * inv_g;

    vec_ZZ_pE rec;
    mul(rec, msg0_out, params.G0);
    return rec == pi0;
  }

  mat_ZZ_pE a;
  a.SetDims(n0, params.k0 + 1);

  for (long row = 0; row < n0; ++row) {
    for (long col = 0; col < params.k0; ++col) {
      a[row][col] = params.G0[col][row];
    }
    a[row][params.k0] = pi0[row];
  }

  if (!SolveLinearSystemRref(msg0_out, a)) return false;

  vec_ZZ_pE rec;
  mul(rec, msg0_out, params.G0);
  return rec == pi0;
}

} // namespace basefold
