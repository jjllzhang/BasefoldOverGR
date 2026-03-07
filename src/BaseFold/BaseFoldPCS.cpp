#include "BaseFold/BaseFoldPCS.hpp"

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/ZZ_pX.h>
#include <NTL/mat_ZZ_pE.h>

#include <openssl/evp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "BaseFold/Multilinear.hpp"
#include "BaseFold/Profile.hpp"
#include "GaloisRing/Inverse.hpp"

#if defined(BASEFOLD_USE_OPENMP)
#include <omp.h>
#endif

using NTL::BytesFromZZ;
using NTL::coeff;
using NTL::LogicError;
using NTL::mat_ZZ_pE;
using NTL::NumBits;
using NTL::NumBytes;
using NTL::rep;
using NTL::vec_ZZ_pE;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pEX;
using NTL::ZZ_pX;
using NTL::ZZFromBytes;

namespace basefold {
namespace {

long ParsePositiveEnvLong(const char *name, long fallback) {
  const char *raw = std::getenv(name);
  if (raw == nullptr || *raw == '\0')
    return fallback;
  char *end = nullptr;
  const long v = std::strtol(raw, &end, 10);
  if (end == raw || *end != '\0' || v <= 0)
    return fallback;
  return v;
}

int ParsePositiveEnvInt(const char *name, int fallback) {
  const char *raw = std::getenv(name);
  if (raw == nullptr || *raw == '\0')
    return fallback;
  char *end = nullptr;
  const long v = std::strtol(raw, &end, 10);
  if (end == raw || *end != '\0' || v <= 0 ||
      v > static_cast<long>(std::numeric_limits<int>::max())) {
    return fallback;
  }
  return static_cast<int>(v);
}

VerifierQueryParallelConfig ParseVerifierQueryParallelConfigFromEnv() {
  VerifierQueryParallelConfig cfg;
  cfg.queries_per_thread = ParsePositiveEnvLong(
      "BASEFOLD_VERIFY_QUERY_QUERIES_PER_THREAD", cfg.queries_per_thread);
  cfg.parallel_query_threshold = ParsePositiveEnvLong(
      "BASEFOLD_VERIFY_QUERY_PARALLEL_THRESHOLD", cfg.parallel_query_threshold);
  cfg.max_threads =
      ParsePositiveEnvInt("BASEFOLD_VERIFY_QUERY_MAX_THREADS", cfg.max_threads);
  return cfg;
}

VerifierQueryParallelConfig &MutableVerifierQueryParallelConfig() {
  static VerifierQueryParallelConfig cfg =
      ParseVerifierQueryParallelConfigFromEnv();
  return cfg;
}

VerifierQueryParallelConfig LoadVerifierQueryParallelConfig() {
  VerifierQueryParallelConfig cfg = MutableVerifierQueryParallelConfig();
  if (cfg.queries_per_thread <= 0)
    cfg.queries_per_thread = 1;
  if (cfg.parallel_query_threshold <= 0)
    cfg.parallel_query_threshold = 2;
  if (cfg.max_threads <= 0)
    cfg.max_threads = 8;
  return cfg;
}

bool TryInvertBaseUnit(FieldElement &inv_out, const FieldElement &a);

bool BatchInvertBaseUnits(std::vector<FieldElement> &inverses,
                          const std::vector<FieldElement> &values);

FieldElement EvalLineAtWithInvDenom(const FieldElement &x,
                                    const FieldElement &x1,
                                    const FieldElement &y1,
                                    const FieldElement &y2,
                                    const FieldElement &inv_denom);

void SortAndUniqueIndices(std::vector<long> &indices) {
  std::sort(indices.begin(), indices.end());
  indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
}

std::size_t FindCachedOpeningPosition(const std::vector<long> &indices,
                                      long index, const char *func_name) {
  const auto it = std::lower_bound(indices.begin(), indices.end(), index);
  if (it == indices.end() || *it != index) {
    const std::string msg =
        std::string(func_name) + ": missing cached opening index";
    LogicError(msg.c_str());
  }
  return static_cast<std::size_t>(it - indices.begin());
}

bool SameMerkleOpening(const MerkleOpening &lhs, const MerkleOpening &rhs) {
  return lhs.index == rhs.index && lhs.value == rhs.value &&
         lhs.auth_path.sibling_hashes == rhs.auth_path.sibling_hashes;
}

std::vector<std::vector<long>> CollectBaseQueryIndicesByTree(
    const std::vector<IOPPQueryPlan> &query_plans,
    const FoldableCodeParams &params) {
  std::vector<std::vector<long>> requested;
  requested.resize(static_cast<std::size_t>(params.d + 1));
  for (const IOPPQueryPlan &plan : query_plans) {
    for (long i = 0; i < params.d; ++i) {
      const long mu_i = plan.mu_by_level[static_cast<std::size_t>(i)];
      const long n_i = CodewordLengthAtLevel(params, i);
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

const FieldElement *FindMerkleMultiproofValue(const MerkleMultiproof &proof,
                                              long index) {
  if (static_cast<long>(proof.values.length()) !=
      static_cast<long>(proof.queried_indices.size())) {
    return nullptr;
  }
  const auto it = std::lower_bound(proof.queried_indices.begin(),
                                   proof.queried_indices.end(), index);
  if (it == proof.queried_indices.end() || *it != index) {
    return nullptr;
  }
  const std::size_t pos = static_cast<std::size_t>(it - proof.queried_indices.begin());
  return &proof.values[static_cast<long>(pos)];
}

template <typename Fn>
void ForEachIndexMaybeParallel(long begin, long end, long parallel_threshold,
                               const Fn &fn) {
  if (end <= begin)
    return;
#if defined(BASEFOLD_USE_OPENMP)
  const long work_items = end - begin;
  if (work_items >= parallel_threshold) {
    const int max_threads = omp_get_max_threads();
    int threads_to_use = static_cast<int>(work_items / parallel_threshold);
    if (threads_to_use > max_threads)
      threads_to_use = max_threads;
    if (threads_to_use >= 2) {
      const ZZ base_modulus = NTL::ZZ_p::modulus();
      const ZZ_pX extension_modulus = NTL::ZZ_pE::modulus().val();
#pragma omp parallel num_threads(threads_to_use)
      {
        NTL::ZZ_p::init(base_modulus);
        NTL::ZZ_pE::init(extension_modulus);
#pragma omp for schedule(static)
        for (long i = begin; i < end; ++i) {
          fn(i);
        }
      }
      return;
    }
  }
#endif
  for (long i = begin; i < end; ++i) {
    fn(i);
  }
}

int ChooseQueryVerifyThreads(long num_queries) {
#if defined(BASEFOLD_USE_OPENMP)
  const VerifierQueryParallelConfig cfg = LoadVerifierQueryParallelConfig();
  if (num_queries < cfg.parallel_query_threshold)
    return 1;
  const long blocks =
      (num_queries + cfg.queries_per_thread - 1) / cfg.queries_per_thread;
  int threads_to_use = static_cast<int>(blocks);
  if (threads_to_use > cfg.max_threads)
    threads_to_use = cfg.max_threads;
  const int max_threads = omp_get_max_threads();
  if (threads_to_use > max_threads)
    threads_to_use = max_threads;
  if (threads_to_use < 1)
    threads_to_use = 1;
  return threads_to_use;
#else
  (void)num_queries;
  return 1;
#endif
}

template <typename Fn>
bool VerifyQueriesMaybeParallel(long num_queries, Profile *prof,
                                const Fn &verify_one_query) {
  if (num_queries <= 0)
    return true;

#if defined(BASEFOLD_USE_OPENMP)
  // Keep profile accounting precise: per-thread timer accumulation in parallel
  // would over-count wall-clock time in the breakdown.
  if (prof == nullptr) {
    const int threads_to_use = ChooseQueryVerifyThreads(num_queries);
    if (threads_to_use >= 2) {
      std::vector<unsigned char> query_ok(static_cast<std::size_t>(num_queries),
                                          static_cast<unsigned char>(1));

      const ZZ base_modulus = NTL::ZZ_p::modulus();
      const ZZ_pX extension_modulus = NTL::ZZ_pE::modulus().val();

#pragma omp parallel num_threads(threads_to_use) shared(query_ok)
      {
        NTL::ZZ_p::init(base_modulus);
        NTL::ZZ_pE::init(extension_modulus);

#pragma omp for schedule(static)
        for (long q = 0; q < num_queries; ++q) {
          if (!verify_one_query(q)) {
            query_ok[static_cast<std::size_t>(q)] =
                static_cast<unsigned char>(0);
          }
        }
      }

      for (long q = 0; q < num_queries; ++q) {
        if (query_ok[static_cast<std::size_t>(q)] == 0)
          return false;
      }
      return true;
    }
  }
#endif

  for (long q = 0; q < num_queries; ++q) {
    if (!verify_one_query(q))
      return false;
  }
  return true;
}

long Pow2Checked(long e) {
  if (e < 0)
    LogicError("Pow2Checked: negative exponent");
  if (e >= static_cast<long>(8 * sizeof(long) - 1)) {
    LogicError("Pow2Checked: exponent too large for long");
  }
  return 1L << e;
}

bool IsPowerOfTwoLong(long n) { return n > 0 && (n & (n - 1)) == 0; }

long Log2ExactPowerOfTwoLong(long n) {
  if (!IsPowerOfTwoLong(n))
    LogicError("Log2ExactPowerOfTwoLong: not a power of two");
  long d = 0;
  while (n > 1) {
    n >>= 1;
    ++d;
  }
  return d;
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

void ProverCommitRoundNoValidate(Oracle &pi_i, const Oracle &pi_ip1,
                                 const FieldElement &alpha_i, long level_i,
                                 const FoldableCodeParams &params) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->prover_commit_round_ns : nullptr,
                    prof ? &prof->prover_commit_round_calls : nullptr);

  const long n_i = CodewordLengthAtLevelNoValidate(params, level_i);
  pi_i.SetLength(n_i);

  const Oracle &diag = params.diag_T[static_cast<std::size_t>(level_i)];
  constexpr long kParallelThreshold = 4096;

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
      ForEachIndexMaybeParallel(0, n_i, kParallelThreshold, [&](long j) {
        pi_i[static_cast<std::size_t>(j)] = EvalLineAtWithInvDenom(
            alpha_i, first_x1, pi_ip1[static_cast<std::size_t>(j)],
            pi_ip1[static_cast<std::size_t>(j + n_i)], inv_denom);
      });
      return;
    }
  }

  std::vector<FieldElement> denoms;
  denoms.resize(static_cast<std::size_t>(n_i));
  ForEachIndexMaybeParallel(0, n_i, kParallelThreshold, [&](long j) {
    const FieldElement &x1 = diag[static_cast<std::size_t>(j)];
    denoms[static_cast<std::size_t>(j)] = (params.zeta * x1) - x1;
  });

  std::vector<FieldElement> inv_denoms;
  if (!BatchInvertBaseUnits(inv_denoms, denoms)) {
    inv_denoms.resize(static_cast<std::size_t>(n_i));
    ForEachIndexMaybeParallel(0, n_i, kParallelThreshold, [&](long j) {
      const FieldElement &denom = denoms[static_cast<std::size_t>(j)];
      if (denom == 0) {
        LogicError("ProverCommitRoundNoValidate: x1 must not equal x2");
      }
      if (!TryInvertBaseUnit(inv_denoms[static_cast<std::size_t>(j)], denom)) {
        LogicError("ProverCommitRoundNoValidate: x2-x1 must be a unit");
      }
    });
  }

  ForEachIndexMaybeParallel(0, n_i, kParallelThreshold, [&](long j) {
    const FieldElement &x1 = diag[static_cast<std::size_t>(j)];
    pi_i[static_cast<std::size_t>(j)] = EvalLineAtWithInvDenom(
        alpha_i, x1, pi_ip1[static_cast<std::size_t>(j)],
        pi_ip1[static_cast<std::size_t>(j + n_i)],
        inv_denoms[static_cast<std::size_t>(j)]);
  });
}

IOPPQueryPlan MakeQueryPlanNoValidate(long initial_mu,
                                      const FoldableCodeParams &params) {
  IOPPQueryPlan plan;
  plan.initial_mu = initial_mu;
  plan.mu_by_level.resize(static_cast<std::size_t>(params.d));

  if (params.d == 0)
    return plan;

  const long n_last = CodewordLengthAtLevelNoValidate(params, params.d - 1);
  if (initial_mu < 0 || initial_mu >= n_last) {
    LogicError("MakeQueryPlanNoValidate: initial_mu out of range");
  }

  long mu = initial_mu;
  for (long i = params.d; i-- > 0;) {
    plan.mu_by_level[static_cast<std::size_t>(i)] = mu;
    if (i > 0) {
      const long n_prev = CodewordLengthAtLevelNoValidate(params, i - 1);
      if (mu >= n_prev)
        mu -= n_prev;
    }
  }
  return plan;
}

void AppendU64(Bytes &out, std::uint64_t v) {
  for (int i = 7; i >= 0; --i) {
    out.push_back(static_cast<Byte>((v >> (8 * i)) & 0xff));
  }
}

void StoreU64BigEndian(Byte *dst, std::uint64_t v) {
  for (int i = 0; i < 8; ++i) {
    dst[i] = static_cast<Byte>((v >> (8 * (7 - i))) & 0xff);
  }
}

void AppendSerializedFieldElement(Bytes &out, const FieldElement &x,
                                  const char *func_name) {
  const long r = ZZ_pE::degree();
  if (r <= 0) {
    const std::string msg =
        std::string(func_name) + ": invalid extension degree";
    LogicError(msg.c_str());
  }

  const ZZ_pX &poly = rep(x);
  AppendU64(out, static_cast<std::uint64_t>(r));
  for (long i = 0; i < r; ++i) {
    const ZZ c = rep(coeff(poly, i));
    const long n = NumBytes(c);
    AppendU64(out, static_cast<std::uint64_t>(n));
    if (n > 0) {
      const std::size_t old_size = out.size();
      out.resize(old_size + static_cast<std::size_t>(n));
      BytesFromZZ(reinterpret_cast<unsigned char *>(out.data() + old_size), c,
                  n);
    }
  }
}

Bytes SerializeFieldElement(const FieldElement &x) {
  Bytes out;
  AppendSerializedFieldElement(out, x, "SerializeFieldElement");
  return out;
}

Bytes Sha256(const Byte *data, std::size_t len) {
  Bytes digest;
  digest.resize(32);
  unsigned int out_len = 0;
  const int ok = EVP_Digest(static_cast<const void *>(data), len,
                            reinterpret_cast<unsigned char *>(digest.data()),
                            &out_len, EVP_sha256(), nullptr);
  if (ok != 1) {
    LogicError("Sha256: EVP_Digest failed");
  }
  if (out_len != digest.size()) {
    LogicError("Sha256: unexpected digest size");
  }
  return digest;
}

Bytes Sha256(const Bytes &data) { return Sha256(data.data(), data.size()); }

Digest Sha256Digest(const Byte *data, std::size_t len, const char *func_name) {
  Digest out{};
  unsigned int out_len = 0;
  const int ok = EVP_Digest(static_cast<const void *>(data), len,
                            reinterpret_cast<unsigned char *>(out.data()),
                            &out_len, EVP_sha256(), nullptr);
  if (ok != 1) {
    const std::string msg = std::string(func_name) + ": EVP_Digest failed";
    LogicError(msg.c_str());
  }
  if (out_len != out.size()) {
    const std::string msg = std::string(func_name) + ": unexpected digest size";
    LogicError(msg.c_str());
  }
  return out;
}

Bytes TaggedHash(Byte tag, const Bytes &state, const Bytes &payload) {
  Bytes in;
  in.reserve(1 + state.size() + 8 + payload.size());
  in.push_back(tag);
  in.insert(in.end(), state.begin(), state.end());
  AppendU64(in, static_cast<std::uint64_t>(payload.size()));
  in.insert(in.end(), payload.begin(), payload.end());
  return Sha256(in);
}

Bytes TaggedHash(Byte tag, const Bytes &state, const std::string &payload) {
  Bytes p;
  p.reserve(payload.size());
  p.insert(p.end(), payload.begin(), payload.end());
  return TaggedHash(tag, state, p);
}

class Sha256Transcript {
public:
  Sha256Transcript() {
    Profile *prof = ActiveProfile();
    ScopedTimer timer(prof ? &prof->transcript_absorb_ns : nullptr,
                      prof ? &prof->transcript_absorb_calls : nullptr);

    const std::string domain = "BaseFoldPCS/v1";
    state_ = TaggedHash(static_cast<Byte>(0x42), Bytes{}, domain);
  }

  void AbsorbBytes(const Bytes &data) {
    Profile *prof = ActiveProfile();
    ScopedTimer timer(prof ? &prof->transcript_absorb_ns : nullptr,
                      prof ? &prof->transcript_absorb_calls : nullptr);
    state_ = TaggedHash(static_cast<Byte>(0x01), state_, data);
  }

  void AbsorbDigest(const Digest &digest) {
    const Bytes tmp(digest.begin(), digest.end());
    AbsorbBytes(tmp);
  }

  void AbsorbFieldElement(const FieldElement &x) {
    Profile *prof = ActiveProfile();
    ScopedTimer timer(prof ? &prof->transcript_absorb_ns : nullptr,
                      prof ? &prof->transcript_absorb_calls : nullptr);
    state_ =
        TaggedHash(static_cast<Byte>(0x02), state_, SerializeFieldElement(x));
  }

  void AbsorbQuadraticPoly(const QuadraticPoly &p) {
    AbsorbFieldElement(p.a0);
    AbsorbFieldElement(p.a1);
    AbsorbFieldElement(p.a2);
  }

  FieldElement ChallengeFieldElement(const std::string &label) const {
    Profile *prof = ActiveProfile();
    ScopedTimer timer(prof ? &prof->transcript_challenge_ns : nullptr,
                      prof ? &prof->transcript_challenge_calls : nullptr);

    const long r = ZZ_pE::degree();
    if (r <= 0)
      LogicError("ChallengeFieldElement: invalid extension degree");

    const ZZ modulus = ZZ_p::modulus();
    if (modulus <= 1)
      LogicError("ChallengeFieldElement: invalid base modulus");

    ChallengeStream stream(state_, "fe/" + label);
    ZZ_pX poly;
    NTL::clear(poly);
    for (long i = 0; i < r; ++i) {
      const ZZ c = stream.SampleZZLessThan(modulus);
      ZZ_p c_base;
      NTL::conv(c_base, c);
      NTL::SetCoeff(poly, i, c_base);
    }
    FieldElement out;
    NTL::conv(out, poly);
    return out;
  }

  long ChallengeIndex(const std::string &label, long upper_bound) const {
    Profile *prof = ActiveProfile();
    ScopedTimer timer(prof ? &prof->transcript_challenge_ns : nullptr,
                      prof ? &prof->transcript_challenge_calls : nullptr);

    if (upper_bound <= 0)
      LogicError("ChallengeIndex: upper_bound must be positive");
    if (upper_bound == 1)
      return 0;

    ChallengeStream stream(state_, "idx/" + label);

    std::uint64_t ub = static_cast<std::uint64_t>(upper_bound);
    std::uint64_t t = ub - 1;
    int bits = 0;
    while (t > 0) {
      ++bits;
      t >>= 1;
    }
    if (bits <= 0 || bits > 63)
      LogicError("ChallengeIndex: unsupported range");
    const int byte_len = (bits + 7) / 8;
    const std::uint64_t mask = (bits == 64) ? ~0ULL : ((1ULL << bits) - 1ULL);

    while (true) {
      std::uint8_t buf[8];
      std::memset(buf, 0, sizeof(buf));
      stream.ReadBytes(buf, static_cast<std::size_t>(byte_len));
      std::uint64_t x = 0;
      for (int i = 0; i < byte_len; ++i) {
        x = (x << 8) | static_cast<std::uint64_t>(buf[i]);
      }
      x &= mask;
      if (x < ub)
        return static_cast<long>(x);
    }
  }

private:
  class ChallengeStream {
  public:
    ChallengeStream(const Bytes &state, const std::string &label)
        : state_(state), label_(label) {}

    void ReadBytes(std::uint8_t *out, std::size_t len) {
      std::size_t written = 0;
      while (written < len) {
        if (offset_ == buf_.size()) {
          buf_ = Digest(counter_++);
          offset_ = 0;
        }
        const std::size_t take = std::min(len - written, buf_.size() - offset_);
        std::memcpy(out + written, buf_.data() + offset_, take);
        offset_ += take;
        written += take;
      }
    }

    ZZ SampleZZLessThan(const ZZ &upper_bound) {
      if (upper_bound <= 0)
        LogicError("SampleZZLessThan: upper_bound must be positive");
      if (upper_bound == 1)
        return ZZ(0);

      const ZZ ub_minus_1 = upper_bound - 1;
      const long bits = NumBits(ub_minus_1);
      if (bits <= 0)
        return ZZ(0);
      const long byte_len = (bits + 7) / 8;
      const ZZ two_to_bits = ZZ(1) << bits;

      Bytes tmp;
      tmp.resize(static_cast<std::size_t>(byte_len));
      while (true) {
        ReadBytes(reinterpret_cast<std::uint8_t *>(tmp.data()),
                  static_cast<std::size_t>(byte_len));
        ZZ x = ZZFromBytes(reinterpret_cast<const unsigned char *>(tmp.data()),
                           byte_len);
        x %= two_to_bits;
        if (x < upper_bound)
          return x;
      }
    }

  private:
    Bytes Digest(std::uint64_t ctr) const {
      Bytes payload;
      payload.reserve(8);
      AppendU64(payload, ctr);
      Bytes st = TaggedHash(static_cast<Byte>(0x20), state_, label_);
      return TaggedHash(static_cast<Byte>(0x21), st, payload);
    }

    Bytes state_;
    std::string label_;

    std::uint64_t counter_ = 0;
    Bytes buf_;
    std::size_t offset_ = 0;
  };

  Bytes state_;
};

void ValidateParamsOrThrow(const FoldableCodeParams &params) {
  (void)MessageLength(params);
}

ZZ PositiveMod(const ZZ &a, const ZZ &m) {
  ZZ r = a % m;
  if (r < 0) {
    r += m;
  }
  return r;
}

bool IsPrimePowerOf(const ZZ &n, const ZZ &p) {
  if (n <= 1 || p <= 1) {
    return false;
  }
  ZZ t = n;
  long exp = 0;
  while (true) {
    const ZZ r = t % p;
    if (r != 0) {
      break;
    }
    t /= p;
    ++exp;
  }
  return exp > 0 && t == 1;
}

ZZ_pX ReduceZZpXModPrime(const ZZ_pX &poly_over_pk, const ZZ &p) {
  ZZ_pX out;
  NTL::clear(out);
  const long d = NTL::deg(poly_over_pk);
  for (long i = 0; i <= d; ++i) {
    ZZ_p c_mod_p;
    NTL::conv(c_mod_p, PositiveMod(rep(coeff(poly_over_pk, i)), p));
    if (c_mod_p != 0) {
      NTL::SetCoeff(out, i, c_mod_p);
    }
  }
  out.normalize();
  return out;
}

ZZ_pEX ReduceZZ_pEXToResidueField(const ZZ_pEX &poly_over_pk_ext, const ZZ &p,
                                  long base_degree) {
  ZZ_pEX out;
  NTL::clear(out);
  const long d = NTL::deg(poly_over_pk_ext);
  for (long i = 0; i <= d; ++i) {
    const ZZ_pX coeff_poly_over_pk = rep(coeff(poly_over_pk_ext, i));
    ZZ_pX coeff_poly_mod_p;
    NTL::clear(coeff_poly_mod_p);
    for (long j = 0; j < base_degree; ++j) {
      ZZ_p c_mod_p;
      NTL::conv(c_mod_p, PositiveMod(rep(coeff(coeff_poly_over_pk, j)), p));
      if (c_mod_p != 0) {
        NTL::SetCoeff(coeff_poly_mod_p, j, c_mod_p);
      }
    }
    coeff_poly_mod_p.normalize();
    ZZ_pE coeff_ext;
    NTL::conv(coeff_ext, coeff_poly_mod_p);
    NTL::SetCoeff(out, i, coeff_ext);
  }
  out.normalize();
  return out;
}

void ValidateChallengeConfigOrThrow(
    const BaseFoldPCSChallengeConfig &challenge_cfg,
    const FoldableCodeParams &params) {
  if (!challenge_cfg.use_extension_challenges) {
    return;
  }

  const long ext_degree = NTL::deg(challenge_cfg.challenge_extension_modulus);
  if (ext_degree <= 0) {
    LogicError(
        "ValidateChallengeConfigOrThrow: extension modulus must have positive "
        "degree");
  }

  FieldElement one;
  NTL::set(one);
  if (NTL::LeadCoeff(challenge_cfg.challenge_extension_modulus) != one) {
    LogicError(
        "ValidateChallengeConfigOrThrow: extension modulus must be monic");
  }

  const ZZ modulus = ZZ_p::modulus();
  if (modulus <= 1) {
    LogicError("ValidateChallengeConfigOrThrow: invalid ZZ_p modulus");
  }

  const long base_degree = ZZ_pE::degree();
  if (base_degree <= 0) {
    LogicError("ValidateChallengeConfigOrThrow: invalid ZZ_pE degree");
  }

  ZZ base_prime = params.p;
  if (base_prime <= 1) {
    if (NTL::ProbPrime(modulus)) {
      base_prime = modulus;
    } else {
      LogicError(
          "ValidateChallengeConfigOrThrow: params.p must be set to the base "
          "prime in ring mode");
    }
  }
  if (!NTL::ProbPrime(base_prime)) {
    LogicError("ValidateChallengeConfigOrThrow: params.p must be prime");
  }
  if (!IsPrimePowerOf(modulus, base_prime)) {
    LogicError(
        "ValidateChallengeConfigOrThrow: current ZZ_p modulus must be a power "
        "of params.p");
  }

  if (NTL::ProbPrime(modulus)) {
    if (NTL::IterIrredTest(challenge_cfg.challenge_extension_modulus) != 1) {
      LogicError("ValidateChallengeConfigOrThrow: extension modulus must be "
                 "irreducible over the current base field");
    }
    return;
  }

  const ZZ_pX base_modulus_over_pk = ZZ_pE::modulus().val();

  NTL::ZZ_pBak modulus_bak;
  modulus_bak.save();
  NTL::ZZ_pEBak extension_bak;
  extension_bak.save();

  ZZ_p::init(base_prime);

  const ZZ_pX base_modulus_over_p =
      ReduceZZpXModPrime(base_modulus_over_pk, base_prime);
  if (NTL::deg(base_modulus_over_p) != base_degree) {
    LogicError(
        "ValidateChallengeConfigOrThrow: current ZZ_pE modulus must stay "
        "full-degree after mod-p reduction");
  }
  if (NTL::IterIrredTest(base_modulus_over_p) != 1) {
    LogicError(
        "ValidateChallengeConfigOrThrow: current ZZ_pE modulus must reduce to "
        "an irreducible polynomial modulo p");
  }
  ZZ_pE::init(base_modulus_over_p);

  const ZZ_pEX ext_modulus_over_p = ReduceZZ_pEXToResidueField(
      challenge_cfg.challenge_extension_modulus, base_prime, base_degree);

  if (NTL::deg(ext_modulus_over_p) != ext_degree) {
    LogicError("ValidateChallengeConfigOrThrow: extension modulus must stay "
               "full-degree after mod-p reduction");
  }
  if (NTL::IterIrredTest(ext_modulus_over_p) != 1) {
    LogicError(
        "ValidateChallengeConfigOrThrow: extension modulus must be basic "
        "irreducible (irreducible after mod-p reduction)");
  }
}

void AbsorbPublicInput(Sha256Transcript &transcript,
                       const MerkleRoot &commitment,
                       const std::vector<FieldElement> &z,
                       const FieldElement &y);

Bytes SerializeExtensionPolynomial(const ZZ_pEX &poly) {
  const long d = NTL::deg(poly);
  Bytes out;
  const std::uint64_t coeff_count =
      (d < 0) ? 0ULL : static_cast<std::uint64_t>(d + 1);
  AppendU64(out, coeff_count);
  for (long i = 0; i <= d; ++i) {
    const Bytes coeff_bytes = SerializeFieldElement(coeff(poly, i));
    AppendU64(out, static_cast<std::uint64_t>(coeff_bytes.size()));
    out.insert(out.end(), coeff_bytes.begin(), coeff_bytes.end());
  }
  return out;
}

void AbsorbChallengeConfig(Sha256Transcript &transcript,
                           const BaseFoldPCSChallengeConfig &challenge_cfg) {
  Bytes tag;
  tag.push_back(
      static_cast<Byte>(challenge_cfg.use_extension_challenges ? 1 : 0));
  transcript.AbsorbBytes(tag);
  if (!challenge_cfg.use_extension_challenges) {
    return;
  }
  transcript.AbsorbBytes(
      SerializeExtensionPolynomial(challenge_cfg.challenge_extension_modulus));
}

long ExtensionDegreeOrThrow(const ZZ_pEX &extension_modulus,
                            const char *func_name) {
  const long ext_degree = NTL::deg(extension_modulus);
  if (ext_degree <= 0) {
    const std::string msg =
        std::string(func_name) + ": invalid extension degree";
    LogicError(msg.c_str());
  }
  return ext_degree;
}

const NTL::ZZ_pEXModulus &
ExtensionModulusContextOrThrow(const ZZ_pEX &extension_modulus,
                               const char *func_name) {
  const long ext_degree = ExtensionDegreeOrThrow(extension_modulus, func_name);
  (void)ext_degree;

  struct CachedModulusContext {
    bool initialized = false;
    ZZ base_modulus;
    long base_degree = 0;
    ZZ_pEX modulus_poly;
    NTL::ZZ_pEXModulus modulus_ctx;
  };

  thread_local CachedModulusContext cache;

  const ZZ &cur_base_modulus = ZZ_p::modulus();
  const long cur_base_degree = ZZ_pE::degree();
  if (!cache.initialized || cache.base_modulus != cur_base_modulus ||
      cache.base_degree != cur_base_degree ||
      cache.modulus_poly != extension_modulus) {
    cache.base_modulus = cur_base_modulus;
    cache.base_degree = cur_base_degree;
    cache.modulus_poly = extension_modulus;
    NTL::build(cache.modulus_ctx, cache.modulus_poly);
    cache.initialized = true;
  }
  return cache.modulus_ctx;
}

FieldElement BaseRingOne() {
  FieldElement one;
  NTL::set(one);
  return one;
}

ZZ_pEX LiftBaseToExtension(const FieldElement &x) {
  ZZ_pEX out;
  NTL::clear(out);
  NTL::SetCoeff(out, 0, x);
  return out;
}

FieldElement ProjectExtensionToBaseConstant(const ZZ_pEX &x) {
  return coeff(x, 0);
}

void ReduceExtensionElementInPlace(ZZ_pEX &x, const ZZ_pEX &extension_modulus) {
  const long extension_degree = ExtensionDegreeOrThrow(
      extension_modulus, "ReduceExtensionElementInPlace");
  if (NTL::deg(x) >= extension_degree) {
    const NTL::ZZ_pEXModulus &mod_ctx = ExtensionModulusContextOrThrow(
        extension_modulus, "ReduceExtensionElementInPlace");
    NTL::rem(x, x, mod_ctx);
    x.normalize();
  }
}

ZZ_pEX ExtensionZero() {
  ZZ_pEX out;
  NTL::clear(out);
  return out;
}

ZZ_pEX ExtensionOne() { return LiftBaseToExtension(BaseRingOne()); }

ZZ_pEX MulExtensionByBaseConstant(const ZZ_pEX &a, const FieldElement &scalar) {
  ZZ_pEX out;
  if (NTL::deg(a) < 0 || scalar == 0) {
    NTL::clear(out);
    return out;
  }
  NTL::mul(out, a, scalar);
  out.normalize();
  return out;
}

ZZ_pEX SubBaseConstantFromExtension(const ZZ_pEX &a, const FieldElement &c) {
  ZZ_pEX out = a;
  NTL::SetCoeff(out, 0, coeff(out, 0) - c);
  out.normalize();
  return out;
}

ZZ_pEX AddExtension(const ZZ_pEX &a, const ZZ_pEX &b,
                    const ZZ_pEX &extension_modulus) {
  ZZ_pEX out = a + b;
  const long extension_degree =
      ExtensionDegreeOrThrow(extension_modulus, "AddExtension");
  if (NTL::deg(out) >= extension_degree) {
    const NTL::ZZ_pEXModulus &mod_ctx =
        ExtensionModulusContextOrThrow(extension_modulus, "AddExtension");
    NTL::rem(out, out, mod_ctx);
    out.normalize();
  }
  return out;
}

ZZ_pEX SubExtension(const ZZ_pEX &a, const ZZ_pEX &b,
                    const ZZ_pEX &extension_modulus) {
  ZZ_pEX out = a - b;
  const long extension_degree =
      ExtensionDegreeOrThrow(extension_modulus, "SubExtension");
  if (NTL::deg(out) >= extension_degree) {
    const NTL::ZZ_pEXModulus &mod_ctx =
        ExtensionModulusContextOrThrow(extension_modulus, "SubExtension");
    NTL::rem(out, out, mod_ctx);
    out.normalize();
  }
  return out;
}

ZZ_pEX MulExtension(const ZZ_pEX &a, const ZZ_pEX &b,
                    const ZZ_pEX &extension_modulus) {
  const long extension_degree =
      ExtensionDegreeOrThrow(extension_modulus, "MulExtension");

  const long deg_a = NTL::deg(a);
  const long deg_b = NTL::deg(b);

  ZZ_pEX out;
  if (deg_a <= 0) {
    out = MulExtensionByBaseConstant(b, coeff(a, 0));
  } else if (deg_b <= 0) {
    out = MulExtensionByBaseConstant(a, coeff(b, 0));
  } else {
    const NTL::ZZ_pEXModulus &mod_ctx =
        ExtensionModulusContextOrThrow(extension_modulus, "MulExtension");
    NTL::MulMod(out, a, b, mod_ctx);
  }

  if (NTL::deg(out) >= extension_degree) {
    const NTL::ZZ_pEXModulus &mod_ctx =
        ExtensionModulusContextOrThrow(extension_modulus, "MulExtension");
    NTL::rem(out, out, mod_ctx);
    out.normalize();
  }
  return out;
}

ZZ_pEX EqFactorExtension(const ZZ_pEX &z_i, const ZZ_pEX &x_i,
                         const ZZ_pEX &extension_modulus) {
  // Eq(z, x) = 1 - z - x + 2*z*x, which saves one extension multiplication.
  const ZZ_pEX one = ExtensionOne();
  const ZZ_pEX zx = MulExtension(z_i, x_i, extension_modulus);
  const ZZ_pEX two_zx = AddExtension(zx, zx, extension_modulus);
  const ZZ_pEX linear = SubExtension(SubExtension(one, z_i, extension_modulus),
                                     x_i, extension_modulus);
  return AddExtension(linear, two_zx, extension_modulus);
}

ZZ_pEX EvalExtensionQuadraticPoly(const ExtensionQuadraticPoly &p,
                                  const ZZ_pEX &x,
                                  const ZZ_pEX &extension_modulus) {
  // Horner form: a0 + x*(a1 + a2*x), one fewer multiplication than x^2 path.
  const ZZ_pEX t = AddExtension(p.a1, MulExtension(p.a2, x, extension_modulus),
                                extension_modulus);
  return AddExtension(p.a0, MulExtension(t, x, extension_modulus),
                      extension_modulus);
}

void AppendSerializedExtensionElement(Bytes &out, const ZZ_pEX &x,
                                      const ZZ_pEX &extension_modulus,
                                      const char *func_name) {
  const long ext_degree = ExtensionDegreeOrThrow(extension_modulus, func_name);
  ZZ_pEX reduced = x;
  ReduceExtensionElementInPlace(reduced, extension_modulus);

  AppendU64(out, static_cast<std::uint64_t>(ext_degree));
  for (long i = 0; i < ext_degree; ++i) {
    AppendSerializedFieldElement(out, coeff(reduced, i), func_name);
  }
}

Bytes SerializeExtensionElement(const ZZ_pEX &x,
                                const ZZ_pEX &extension_modulus) {
  Bytes out;
  AppendSerializedExtensionElement(out, x, extension_modulus,
                                   "SerializeExtensionElement");
  return out;
}

void AbsorbExtensionElement(Sha256Transcript &transcript, const ZZ_pEX &x,
                            const ZZ_pEX &extension_modulus) {
  transcript.AbsorbBytes(SerializeExtensionElement(x, extension_modulus));
}

void AbsorbExtensionQuadraticPoly(Sha256Transcript &transcript,
                                  const ExtensionQuadraticPoly &p,
                                  const ZZ_pEX &extension_modulus) {
  AbsorbExtensionElement(transcript, p.a0, extension_modulus);
  AbsorbExtensionElement(transcript, p.a1, extension_modulus);
  AbsorbExtensionElement(transcript, p.a2, extension_modulus);
}

ZZ_pEX SampleExtensionChallenge(const Sha256Transcript &transcript,
                                const std::string &label,
                                const ZZ_pEX &extension_modulus) {
  const long ext_degree =
      ExtensionDegreeOrThrow(extension_modulus, "SampleExtensionChallenge");

  ZZ_pEX sampled;
  NTL::clear(sampled);
  for (long i = 0; i < ext_degree; ++i) {
    const FieldElement c = transcript.ChallengeFieldElement(
        "ext/" + label + "/coeff/" + std::to_string(i));
    NTL::SetCoeff(sampled, i, c);
  }
  ReduceExtensionElementInPlace(sampled, extension_modulus);
  return sampled;
}

bool TryInvertBaseUnit(FieldElement &inv_out, const FieldElement &a) {
  if (a == 0) {
    return false;
  }

  const ZZ modulus = ZZ_p::modulus();
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
    const ZZ a_rep = rep(coeff(rep(a), 0));
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

  const FieldElement one = BaseRingOne();
  return a * inv_out == one;
}

bool BatchInvertBaseUnits(std::vector<FieldElement> &inverses,
                          const std::vector<FieldElement> &values) {
  const long n = static_cast<long>(values.size());
  inverses.resize(static_cast<std::size_t>(n));
  if (n == 0) {
    return true;
  }

  std::vector<FieldElement> prefix;
  prefix.resize(static_cast<std::size_t>(n + 1));
  prefix[0] = BaseRingOne();

  for (long i = 0; i < n; ++i) {
    const FieldElement &v = values[static_cast<std::size_t>(i)];
    if (v == 0) {
      return false;
    }
    prefix[static_cast<std::size_t>(i + 1)] =
        prefix[static_cast<std::size_t>(i)] * v;
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

ZZ_pEX EvalLineAtExtensionWithInvDenom(const ZZ_pEX &x, const FieldElement &x1,
                                       const ZZ_pEX &y1, const ZZ_pEX &y2,
                                       const FieldElement &inv_denom,
                                       const ZZ_pEX &extension_modulus) {
  const ZZ_pEX delta_y = SubExtension(y2, y1, extension_modulus);
  const ZZ_pEX slope = MulExtensionByBaseConstant(delta_y, inv_denom);
  const ZZ_pEX delta_x = SubBaseConstantFromExtension(x, x1);
  const ZZ_pEX correction = MulExtension(slope, delta_x, extension_modulus);
  return AddExtension(y1, correction, extension_modulus);
}

ZZ_pEX EvalLineAtExtension(const ZZ_pEX &x, const FieldElement &x1,
                           const ZZ_pEX &y1, const FieldElement &x2,
                           const ZZ_pEX &y2, const ZZ_pEX &extension_modulus) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->ext_eval_line_at_ns : nullptr,
                    prof ? &prof->ext_eval_line_at_calls : nullptr);

  const FieldElement denom = x2 - x1;
  FieldElement inv_denom;
  if (!TryInvertBaseUnit(inv_denom, denom)) {
    LogicError("EvalLineAtExtension: denominator is not invertible");
  }

  return EvalLineAtExtensionWithInvDenom(x, x1, y1, y2, inv_denom,
                                         extension_modulus);
}

std::vector<ZZ_pEX> LiftOracleToExtension(const Oracle &oracle) {
  std::vector<ZZ_pEX> out;
  out.resize(static_cast<std::size_t>(oracle.length()));
  for (long i = 0; i < oracle.length(); ++i) {
    out[static_cast<std::size_t>(i)] = LiftBaseToExtension(oracle[i]);
  }
  return out;
}

std::vector<ZZ_pEX>
BooleanEvalTableFromMonomialCoeffsExtension(const std::vector<ZZ_pEX> &coeffs,
                                            long k,
                                            const ZZ_pEX &extension_modulus) {
  if (k < 0) {
    LogicError(
        "BooleanEvalTableFromMonomialCoeffsExtension: negative dimension");
  }
  if (static_cast<long>(coeffs.size()) != (1L << k)) {
    LogicError("BooleanEvalTableFromMonomialCoeffsExtension: length mismatch");
  }

  std::vector<ZZ_pEX> eval = coeffs;
  for (long bit = 0; bit < k; ++bit) {
    const long step = 1L << bit;
    for (long mask = 0; mask < static_cast<long>(eval.size()); ++mask) {
      if (mask & step) {
        eval[static_cast<std::size_t>(mask)] = AddExtension(
            eval[static_cast<std::size_t>(mask)],
            eval[static_cast<std::size_t>(mask ^ step)], extension_modulus);
      }
    }
  }
  return eval;
}

std::vector<ZZ_pEX>
Msg0CoeffsAtSuffixChallenges(const vec_ZZ_pE &f_coeffs, long kappa,
                             const std::vector<ZZ_pEX> &r_by_level,
                             const ZZ_pEX &extension_modulus) {
  if (kappa < 0) {
    LogicError("Msg0CoeffsAtSuffixChallenges: negative kappa");
  }
  const long d = static_cast<long>(r_by_level.size());
  const long point_dim = kappa + d;
  if (f_coeffs.length() != (1L << point_dim)) {
    LogicError("Msg0CoeffsAtSuffixChallenges: f_coeffs length mismatch");
  }

  std::vector<ZZ_pEX> cur;
  cur.resize(static_cast<std::size_t>(f_coeffs.length()));
  for (long i = 0; i < f_coeffs.length(); ++i) {
    cur[static_cast<std::size_t>(i)] =
        LiftBaseToExtension(f_coeffs[static_cast<std::size_t>(i)]);
  }

  long cur_len = static_cast<long>(cur.size());
  for (long var = point_dim; var-- > kappa;) {
    const long half = cur_len / 2;
    const ZZ_pEX &r_var = r_by_level[static_cast<std::size_t>(var - kappa)];
    for (long i = 0; i < half; ++i) {
      cur[static_cast<std::size_t>(i)] =
          AddExtension(cur[static_cast<std::size_t>(i)],
                       MulExtension(cur[static_cast<std::size_t>(i + half)],
                                    r_var, extension_modulus),
                       extension_modulus);
    }
    cur.resize(static_cast<std::size_t>(half));
    cur_len = half;
  }
  return cur;
}

std::vector<ZZ_pEX> EncodeC0Extension(const std::vector<ZZ_pEX> &msg0_coeffs,
                                      const FoldableCodeParams &params,
                                      const ZZ_pEX &extension_modulus) {
  if (static_cast<long>(msg0_coeffs.size()) != params.k0) {
    LogicError("EncodeC0Extension: msg0_coeffs has wrong length");
  }

  const long n0 = CodewordLengthAtLevelNoValidate(params, 0);
  std::vector<ZZ_pEX> out;
  out.resize(static_cast<std::size_t>(n0), ExtensionZero());

  for (long j = 0; j < n0; ++j) {
    ZZ_pEX acc = ExtensionZero();
    for (long row = 0; row < params.k0; ++row) {
      const ZZ_pEX g =
          LiftBaseToExtension(params.G0[static_cast<std::size_t>(row)][j]);
      const ZZ_pEX term = MulExtension(
          msg0_coeffs[static_cast<std::size_t>(row)], g, extension_modulus);
      acc = AddExtension(acc, term, extension_modulus);
    }
    out[static_cast<std::size_t>(j)] = acc;
  }

  return out;
}

class ExtensionSumcheckProver {
public:
  ExtensionSumcheckProver(const FieldVec &f_coeffs,
                          const std::vector<FieldElement> &z,
                          const ZZ_pEX &extension_modulus)
      : extension_modulus_(extension_modulus) {
    Profile *prof = ActiveProfile();
    ScopedTimer timer(prof ? &prof->ext_sumcheck_init_ns : nullptr,
                      prof ? &prof->ext_sumcheck_init_calls : nullptr);

    const long n = f_coeffs.length();
    if (!IsPowerOfTwoLong(n)) {
      LogicError("ExtensionSumcheckProver: f_coeffs length must be 2^d");
    }
    d_ = Log2ExactPowerOfTwoLong(n);
    if (static_cast<long>(z.size()) != d_) {
      LogicError("ExtensionSumcheckProver: z dimension mismatch");
    }

    cur_k_ = d_;
    z_.resize(static_cast<std::size_t>(d_));
    for (long i = 0; i < d_; ++i) {
      z_[static_cast<std::size_t>(i)] =
          LiftBaseToExtension(z[static_cast<std::size_t>(i)]);
    }

    std::vector<ZZ_pEX> lifted_coeffs;
    lifted_coeffs.resize(static_cast<std::size_t>(n));
    for (long i = 0; i < n; ++i) {
      lifted_coeffs[static_cast<std::size_t>(i)] =
          LiftBaseToExtension(f_coeffs[static_cast<std::size_t>(i)]);
    }
    f_eval_table_ = BooleanEvalTableFromMonomialCoeffsExtension(
        lifted_coeffs, d_, extension_modulus_);

    prefix_eq_by_vars_.resize(static_cast<std::size_t>(d_));
    if (d_ > 0) {
      prefix_eq_by_vars_[0].resize(1);
      prefix_eq_by_vars_[0][0] = ExtensionOne();

      for (long t = 1; t < d_; ++t) {
        const ZZ_pEX z_var = z_[static_cast<std::size_t>(t - 1)];
        const ZZ_pEX factor0 =
            SubExtension(ExtensionOne(), z_var, extension_modulus_);
        const ZZ_pEX factor1 = z_var;

        const std::vector<ZZ_pEX> &prev =
            prefix_eq_by_vars_[static_cast<std::size_t>(t - 1)];
        const long old = static_cast<long>(prev.size());
        prefix_eq_by_vars_[static_cast<std::size_t>(t)].resize(
            static_cast<std::size_t>(2 * old));
        std::vector<ZZ_pEX> &cur =
            prefix_eq_by_vars_[static_cast<std::size_t>(t)];
        for (long mask = 0; mask < old; ++mask) {
          const ZZ_pEX &base = prev[static_cast<std::size_t>(mask)];
          cur[static_cast<std::size_t>(mask)] =
              MulExtension(base, factor0, extension_modulus_);
          cur[static_cast<std::size_t>(mask + old)] =
              MulExtension(base, factor1, extension_modulus_);
        }
      }
    }

    suffix_eq_prod_ = ExtensionOne();
  }

  ExtensionQuadraticPoly CurrentPolynomial() const {
    Profile *prof = ActiveProfile();
    ScopedTimer timer(prof ? &prof->ext_sumcheck_current_poly_ns : nullptr,
                      prof ? &prof->ext_sumcheck_current_poly_calls : nullptr);

    if (cur_k_ <= 0) {
      LogicError("ExtensionSumcheckProver::CurrentPolynomial: no variables");
    }

    const long k = cur_k_;
    const long n = static_cast<long>(f_eval_table_.size());
    if (n != (1L << k)) {
      LogicError("ExtensionSumcheckProver::CurrentPolynomial: internal length "
                 "mismatch");
    }

    const long half = 1L << (k - 1);
    const std::vector<ZZ_pEX> &prefix =
        prefix_eq_by_vars_[static_cast<std::size_t>(k - 1)];
    if (static_cast<long>(prefix.size()) != half) {
      LogicError(
          "ExtensionSumcheckProver::CurrentPolynomial: prefix size mismatch");
    }

    const FieldElement one_base = BaseRingOne();
    const FieldElement z_k_base =
        ProjectExtensionToBaseConstant(z_[static_cast<std::size_t>(k - 1)]);
    const FieldElement factor0_base = one_base - z_k_base;
    const FieldElement delta_factor_base = z_k_base - factor0_base;

    ExtensionQuadraticPoly out;
    out.a0 = ExtensionZero();
    out.a1 = ExtensionZero();
    out.a2 = ExtensionZero();

    for (long mask = 0; mask < half; ++mask) {
      const ZZ_pEX &prefix_mask = prefix[static_cast<std::size_t>(mask)];
      const FieldElement prefix_mask_base =
          ProjectExtensionToBaseConstant(prefix_mask);
      const ZZ_pEX common =
          MulExtensionByBaseConstant(suffix_eq_prod_, prefix_mask_base);

      const ZZ_pEX eq0 = MulExtensionByBaseConstant(common, factor0_base);
      const ZZ_pEX delta_eq =
          MulExtensionByBaseConstant(common, delta_factor_base);

      const ZZ_pEX &f0 = f_eval_table_[static_cast<std::size_t>(mask)];
      const ZZ_pEX &f1 = f_eval_table_[static_cast<std::size_t>(mask + half)];
      const ZZ_pEX delta_f = SubExtension(f1, f0, extension_modulus_);

      out.a0 = AddExtension(out.a0, MulExtension(f0, eq0, extension_modulus_),
                            extension_modulus_);

      const ZZ_pEX term1 = MulExtension(f0, delta_eq, extension_modulus_);
      const ZZ_pEX term2 = MulExtension(delta_f, eq0, extension_modulus_);
      out.a1 =
          AddExtension(out.a1, AddExtension(term1, term2, extension_modulus_),
                       extension_modulus_);

      out.a2 = AddExtension(out.a2,
                            MulExtension(delta_f, delta_eq, extension_modulus_),
                            extension_modulus_);
    }

    return out;
  }

  void ReceiveChallenge(const ZZ_pEX &r_kminus1) {
    Profile *prof = ActiveProfile();
    ScopedTimer timer(prof ? &prof->ext_sumcheck_receive_challenge_ns : nullptr,
                      prof ? &prof->ext_sumcheck_receive_challenge_calls
                           : nullptr);

    if (cur_k_ <= 0) {
      LogicError("ExtensionSumcheckProver::ReceiveChallenge: no variables");
    }

    const long k = cur_k_;
    const long n = static_cast<long>(f_eval_table_.size());
    if (n != (1L << k)) {
      LogicError("ExtensionSumcheckProver::ReceiveChallenge: internal length "
                 "mismatch");
    }

    const ZZ_pEX eq = EqFactorExtension(z_[static_cast<std::size_t>(k - 1)],
                                        r_kminus1, extension_modulus_);
    suffix_eq_prod_ = MulExtension(suffix_eq_prod_, eq, extension_modulus_);

    const long half = n / 2;
    for (long i = 0; i < half; ++i) {
      const ZZ_pEX &f0 = f_eval_table_[static_cast<std::size_t>(i)];
      const ZZ_pEX &f1 = f_eval_table_[static_cast<std::size_t>(i + half)];
      const ZZ_pEX delta_f = SubExtension(f1, f0, extension_modulus_);
      f_eval_table_[static_cast<std::size_t>(i)] =
          AddExtension(f0, MulExtension(delta_f, r_kminus1, extension_modulus_),
                       extension_modulus_);
    }
    f_eval_table_.resize(static_cast<std::size_t>(half));
    --cur_k_;
  }

private:
  long d_ = 0;
  long cur_k_ = 0;
  ZZ_pEX extension_modulus_;
  std::vector<ZZ_pEX> z_;
  std::vector<ZZ_pEX> f_eval_table_;
  std::vector<std::vector<ZZ_pEX>> prefix_eq_by_vars_;
  ZZ_pEX suffix_eq_prod_;
};

bool CheckExtensionSumcheckRelations(
    const std::vector<ExtensionQuadraticPoly> &h_by_level,
    const std::vector<ZZ_pEX> &r, const FieldElement &claimed_y,
    const ZZ_pEX &extension_modulus) {
  const long d = static_cast<long>(h_by_level.size());
  if (static_cast<long>(r.size()) != d) {
    return false;
  }
  if (d == 0) {
    return true;
  }

  const ZZ_pEX zero = ExtensionZero();
  const ZZ_pEX one = ExtensionOne();
  const ZZ_pEX claimed_y_ext = LiftBaseToExtension(claimed_y);

  const ExtensionQuadraticPoly &h_d =
      h_by_level[static_cast<std::size_t>(d - 1)];
  const ZZ_pEX hd0 = EvalExtensionQuadraticPoly(h_d, zero, extension_modulus);
  const ZZ_pEX hd1 = EvalExtensionQuadraticPoly(h_d, one, extension_modulus);
  if (AddExtension(hd0, hd1, extension_modulus) != claimed_y_ext) {
    return false;
  }

  for (long k = 1; k < d; ++k) {
    const ExtensionQuadraticPoly &h_k =
        h_by_level[static_cast<std::size_t>(k - 1)];
    const ExtensionQuadraticPoly &h_kp1 =
        h_by_level[static_cast<std::size_t>(k)];
    const ZZ_pEX lhs =
        AddExtension(EvalExtensionQuadraticPoly(h_k, zero, extension_modulus),
                     EvalExtensionQuadraticPoly(h_k, one, extension_modulus),
                     extension_modulus);
    const ZZ_pEX rhs = EvalExtensionQuadraticPoly(
        h_kp1, r[static_cast<std::size_t>(k)], extension_modulus);
    if (lhs != rhs) {
      return false;
    }
  }

  return true;
}

void ProverCommitRoundExtensionNoValidate(std::vector<ZZ_pEX> &pi_i,
                                          const std::vector<ZZ_pEX> &pi_ip1,
                                          const ZZ_pEX &alpha_i, long level_i,
                                          const FoldableCodeParams &params,
                                          const ZZ_pEX &extension_modulus) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->ext_prover_commit_round_ns : nullptr,
                    prof ? &prof->ext_prover_commit_round_calls : nullptr);

  const long n_i = CodewordLengthAtLevelNoValidate(params, level_i);
  if (static_cast<long>(pi_ip1.size()) != 2 * n_i) {
    LogicError("ProverCommitRoundExtensionNoValidate: pi_ip1 has wrong length");
  }
  const Oracle &diag = params.diag_T[static_cast<std::size_t>(level_i)];
  if (diag.length() != n_i) {
    LogicError("ProverCommitRoundExtensionNoValidate: diag_T length mismatch");
  }
  pi_i.resize(static_cast<std::size_t>(n_i));

  const FieldElement one = BaseRingOne();
  const FieldElement zeta_minus_one = params.zeta - one;
  if (zeta_minus_one == 0) {
    LogicError("ProverCommitRoundExtensionNoValidate: zeta must not be 1");
  }

  std::vector<FieldElement> denoms;
  denoms.resize(static_cast<std::size_t>(n_i));
  constexpr long kParallelThreshold = 4096;
  ForEachIndexMaybeParallel(0, n_i, kParallelThreshold, [&](long j) {
    const FieldElement &t = diag[static_cast<std::size_t>(j)];
    denoms[static_cast<std::size_t>(j)] = zeta_minus_one * t;
  });

  std::vector<FieldElement> inv_denoms;
  if (!BatchInvertBaseUnits(inv_denoms, denoms)) {
    inv_denoms.resize(static_cast<std::size_t>(n_i));
    ForEachIndexMaybeParallel(0, n_i, kParallelThreshold, [&](long j) {
      if (!TryInvertBaseUnit(inv_denoms[static_cast<std::size_t>(j)],
                             denoms[static_cast<std::size_t>(j)])) {
        LogicError("ProverCommitRoundExtensionNoValidate: denominator is not "
                   "invertible");
      }
    });
  }

  ForEachIndexMaybeParallel(0, n_i, kParallelThreshold, [&](long j) {
    const FieldElement &x1 = diag[static_cast<std::size_t>(j)];
    pi_i[static_cast<std::size_t>(j)] = EvalLineAtExtensionWithInvDenom(
        alpha_i, x1, pi_ip1[static_cast<std::size_t>(j)],
        pi_ip1[static_cast<std::size_t>(j + n_i)],
        inv_denoms[static_cast<std::size_t>(j)], extension_modulus);
  });
}

std::size_t MaxSerializedFieldElementSizeOrThrow(const char *func_name) {
  const long r = ZZ_pE::degree();
  if (r <= 0) {
    const std::string msg =
        std::string(func_name) + ": invalid extension degree";
    LogicError(msg.c_str());
  }
  const long coeff_max_bytes = std::max<long>(1, NumBytes(ZZ_p::modulus()));
  return static_cast<std::size_t>(8) +
         static_cast<std::size_t>(r) *
             (static_cast<std::size_t>(8) +
              static_cast<std::size_t>(coeff_max_bytes));
}

std::size_t
MaxSerializedExtensionElementSizeOrThrow(const ZZ_pEX &extension_modulus,
                                         const char *func_name) {
  const long ext_degree = ExtensionDegreeOrThrow(extension_modulus, func_name);
  return static_cast<std::size_t>(8) +
         static_cast<std::size_t>(ext_degree) *
             MaxSerializedFieldElementSizeOrThrow(func_name);
}

Digest HashExtensionLeaf(long index, const ZZ_pEX &value,
                         const ZZ_pEX &extension_modulus,
                         Bytes *scratch_payload = nullptr) {
  Bytes local_payload;
  Bytes &payload =
      (scratch_payload != nullptr) ? *scratch_payload : local_payload;
  payload.clear();
  payload.push_back(static_cast<Byte>(0x30));
  AppendU64(payload, static_cast<std::uint64_t>(index));

  const std::size_t encoded_len_pos = payload.size();
  AppendU64(payload, 0);
  const std::size_t encoded_start = payload.size();
  AppendSerializedExtensionElement(payload, value, extension_modulus,
                                   "HashExtensionLeaf");
  const std::uint64_t encoded_len =
      static_cast<std::uint64_t>(payload.size() - encoded_start);
  StoreU64BigEndian(payload.data() + encoded_len_pos, encoded_len);

  return Sha256Digest(payload.data(), payload.size(), "HashExtensionLeaf");
}

Digest HashExtensionNode(const Digest &left, const Digest &right) {
  std::array<Byte, 1 + 32 + 32> payload{};
  payload[0] = static_cast<Byte>(0x31);
  std::copy(left.begin(), left.end(), payload.begin() + 1);
  std::copy(right.begin(), right.end(), payload.begin() + 1 + left.size());
  return Sha256Digest(payload.data(), payload.size(), "HashExtensionNode");
}

Digest HashExtensionRawRootEmpty() {
  const Byte payload[1] = {static_cast<Byte>(0x32)};
  return Sha256Digest(payload, sizeof(payload), "HashExtensionRawRootEmpty");
}

Digest HashExtensionRootWithCount(long leaf_count, const Digest &raw_root) {
  if (leaf_count <= 0) {
    LogicError("HashExtensionRootWithCount: invalid leaf_count");
  }
  std::array<Byte, 1 + 8 + 32> payload{};
  payload[0] = static_cast<Byte>(0x33);
  StoreU64BigEndian(payload.data() + 1, static_cast<std::uint64_t>(leaf_count));
  std::copy(raw_root.begin(), raw_root.end(), payload.begin() + 1 + 8);
  return Sha256Digest(payload.data(), payload.size(),
                      "HashExtensionRootWithCount");
}

std::size_t ExpectedExtensionMerkleHeight(long leaf_count) {
  if (leaf_count <= 0) {
    return 0;
  }
  std::size_t height = 0;
  long n = leaf_count;
  while (n > 1) {
    if ((n & 1L) == 1L) {
      ++n;
    }
    n >>= 1;
    ++height;
  }
  return height;
}

class ExtensionMerkleTree {
public:
  ExtensionMerkleTree() = default;

  static ExtensionMerkleTree Build(const std::vector<ZZ_pEX> &oracle,
                                   const ZZ_pEX &extension_modulus) {
    Profile *prof = ActiveProfile();
    ScopedTimer timer(prof ? &prof->ext_merkle_tree_build_ns : nullptr,
                      prof ? &prof->ext_merkle_tree_build_calls : nullptr);

    if (oracle.empty()) {
      LogicError("ExtensionMerkleTree::Build: oracle must be non-empty");
    }

    ExtensionMerkleTree tree;
    tree.leaf_count_ = static_cast<long>(oracle.size());

    Bytes leaf_payload;
    leaf_payload.reserve(static_cast<std::size_t>(1 + 8 + 8) +
                         MaxSerializedExtensionElementSizeOrThrow(
                             extension_modulus, "ExtensionMerkleTree::Build"));

    std::vector<Digest> level;
    level.resize(oracle.size());
    for (long i = 0; i < tree.leaf_count_; ++i) {
      level[static_cast<std::size_t>(i)] =
          HashExtensionLeaf(i, oracle[static_cast<std::size_t>(i)],
                            extension_modulus, &leaf_payload);
    }

    tree.levels_.clear();
    tree.levels_.push_back(level);
    while (level.size() > 1) {
      if ((level.size() % 2U) == 1U) {
        level.push_back(level.back());
      }
      tree.levels_.back() = level;

      std::vector<Digest> next;
      next.reserve(level.size() / 2U);
      for (std::size_t i = 0; i < level.size(); i += 2U) {
        next.push_back(HashExtensionNode(level[i], level[i + 1U]));
      }
      tree.levels_.push_back(next);
      level = std::move(next);
    }

    if (tree.levels_.empty()) {
      tree.raw_root_ = HashExtensionRawRootEmpty();
    } else {
      tree.raw_root_ = tree.levels_.back()[0];
    }
    return tree;
  }

  long LeafCount() const { return leaf_count_; }

  MerkleRoot Root() const {
    return HashExtensionRootWithCount(leaf_count_, raw_root_);
  }

  ExtensionMerkleOpening Open(const std::vector<ZZ_pEX> &oracle,
                              long index) const {
    Profile *prof = ActiveProfile();
    ScopedTimer timer(prof ? &prof->ext_merkle_tree_open_ns : nullptr,
                      prof ? &prof->ext_merkle_tree_open_calls : nullptr);

    if (leaf_count_ <= 0) {
      LogicError("ExtensionMerkleTree::Open: empty tree");
    }
    if (static_cast<long>(oracle.size()) != leaf_count_) {
      LogicError("ExtensionMerkleTree::Open: oracle length mismatch");
    }
    if (index < 0 || index >= leaf_count_) {
      LogicError("ExtensionMerkleTree::Open: index out of range");
    }

    ExtensionMerkleOpening opening;
    opening.index = index;
    opening.value = oracle[static_cast<std::size_t>(index)];

    long cur = index;
    const std::size_t height = levels_.size();
    for (std::size_t h = 0; h + 1U < height; ++h) {
      const std::vector<Digest> &level = levels_[h];
      const long sibling_index = (cur & 1L) ? (cur - 1L) : (cur + 1L);
      opening.auth_path.sibling_hashes.push_back(
          level[static_cast<std::size_t>(sibling_index)]);
      cur >>= 1;
    }
    return opening;
  }

private:
  long leaf_count_ = 0;
  Digest raw_root_{};
  std::vector<std::vector<Digest>> levels_;
};

bool VerifyExtensionMerkleOpening(const MerkleRoot &root, long leaf_count,
                                  const ExtensionMerkleOpening &opening,
                                  const ZZ_pEX &extension_modulus) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->ext_merkle_verify_opening_ns : nullptr,
                    prof ? &prof->ext_merkle_verify_opening_calls : nullptr);

  if (leaf_count <= 0) {
    return false;
  }
  if (opening.index < 0 || opening.index >= leaf_count) {
    return false;
  }
  if (opening.auth_path.sibling_hashes.size() !=
      ExpectedExtensionMerkleHeight(leaf_count)) {
    return false;
  }

  Digest cur =
      HashExtensionLeaf(opening.index, opening.value, extension_modulus);
  long idx = opening.index;
  for (const Digest &sib : opening.auth_path.sibling_hashes) {
    if (idx & 1L) {
      cur = HashExtensionNode(sib, cur);
    } else {
      cur = HashExtensionNode(cur, sib);
    }
    idx >>= 1;
  }
  return HashExtensionRootWithCount(leaf_count, cur) == root;
}

BaseFoldPCSEvalProof ProveEvalWithExtensionChallengesUnchecked(
    const vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const FoldableCodeParams &params,
    const BaseFoldPCSChallengeConfig &challenge_cfg) {
  if (params.d == 0) {
    return BaseFoldPCSProveEvalUnchecked(f_coeffs, z, claimed_y, num_queries,
                                         params);
  }

  const ZZ_pEX &extension_modulus = challenge_cfg.challenge_extension_modulus;
  ExtensionDegreeOrThrow(extension_modulus,
                         "ProveEvalWithExtensionChallengesUnchecked");

  BaseFoldPCSEvalProof proof;
  proof.extension.enabled = true;
  proof.commitments.roots_by_level.clear();
  proof.extension.roots_by_level.resize(static_cast<std::size_t>(params.d));

  Oracle pi_d_base;
  EncodeFoldableUnchecked(pi_d_base, f_coeffs, params);
  const MerkleTree merkle_d = MerkleTree::Build(pi_d_base);
  const MerkleRoot root_d = merkle_d.Root();
  proof.commitments.roots_by_level.push_back(root_d);

  Sha256Transcript transcript;
  AbsorbPublicInput(transcript, root_d, z, claimed_y);
  AbsorbChallengeConfig(transcript, challenge_cfg);

  proof.extension.h_by_level.resize(static_cast<std::size_t>(params.d));
  std::vector<ZZ_pEX> r_by_level;
  r_by_level.resize(static_cast<std::size_t>(params.d));

  std::vector<std::vector<ZZ_pEX>> ext_oracles;
  ext_oracles.resize(static_cast<std::size_t>(params.d + 1));
  ext_oracles[static_cast<std::size_t>(params.d)] =
      LiftOracleToExtension(pi_d_base);

  std::vector<ExtensionMerkleTree> ext_merkle_trees;
  ext_merkle_trees.resize(static_cast<std::size_t>(params.d));

  ExtensionSumcheckProver sumcheck_ext(f_coeffs, z, extension_modulus);
  const ExtensionQuadraticPoly h_d_ext = sumcheck_ext.CurrentPolynomial();
  proof.extension.h_by_level[static_cast<std::size_t>(params.d - 1)] = h_d_ext;
  AbsorbExtensionQuadraticPoly(transcript, h_d_ext, extension_modulus);

  for (long i = params.d; i-- > 0;) {
    const ZZ_pEX r_i_ext = SampleExtensionChallenge(
        transcript, "r/" + std::to_string(i), extension_modulus);
    r_by_level[static_cast<std::size_t>(i)] = r_i_ext;

    ProverCommitRoundExtensionNoValidate(
        ext_oracles[static_cast<std::size_t>(i)],
        ext_oracles[static_cast<std::size_t>(i + 1)], r_i_ext, i, params,
        extension_modulus);
    ext_merkle_trees[static_cast<std::size_t>(i)] = ExtensionMerkleTree::Build(
        ext_oracles[static_cast<std::size_t>(i)], extension_modulus);
    proof.extension.roots_by_level[static_cast<std::size_t>(i)] =
        ext_merkle_trees[static_cast<std::size_t>(i)].Root();
    transcript.AbsorbDigest(
        proof.extension.roots_by_level[static_cast<std::size_t>(i)]);

    sumcheck_ext.ReceiveChallenge(r_i_ext);
    if (i > 0) {
      const ExtensionQuadraticPoly h_i_ext = sumcheck_ext.CurrentPolynomial();
      proof.extension.h_by_level[static_cast<std::size_t>(i - 1)] = h_i_ext;
      AbsorbExtensionQuadraticPoly(transcript, h_i_ext, extension_modulus);
    }
  }

  proof.extension.pi0_full = ext_oracles[0];

  const long kappa = Log2ExactPowerOfTwoLong(params.k0);
  proof.extension.msg0_coeffs = Msg0CoeffsAtSuffixChallenges(
      f_coeffs, kappa, r_by_level, extension_modulus);
  const std::vector<ZZ_pEX> expected_pi0 =
      EncodeC0Extension(proof.extension.msg0_coeffs, params, extension_modulus);
  if (expected_pi0 != proof.extension.pi0_full) {
    LogicError(
        "ProveEvalWithExtensionChallengesUnchecked: internal pi0 mismatch");
  }

  proof.query_proofs.resize(static_cast<std::size_t>(num_queries));
  proof.extension.query_proofs.resize(static_cast<std::size_t>(num_queries));

  const long n_last = CodewordLengthAtLevelNoValidate(params, params.d - 1);
  for (long q = 0; q < num_queries; ++q) {
    const long mu =
        transcript.ChallengeIndex("mu/" + std::to_string(q), n_last);
    const IOPPQueryPlan plan = MakeQueryPlanNoValidate(mu, params);

    BaseFoldPCSQueryProof qp;
    qp.left.resize(1);
    qp.right.resize(1);
    qp.folded.clear();

    BaseFoldPCSQueryProofExtension qp_ext;
    qp_ext.left.resize(static_cast<std::size_t>(params.d));
    qp_ext.right.resize(static_cast<std::size_t>(params.d));
    qp_ext.folded.resize(static_cast<std::size_t>(params.d));

    for (long i = 0; i < params.d; ++i) {
      const long mu_i = plan.mu_by_level[static_cast<std::size_t>(i)];
      const long n_i = CodewordLengthAtLevelNoValidate(params, i);
      qp_ext.folded[static_cast<std::size_t>(i)] =
          ext_merkle_trees[static_cast<std::size_t>(i)].Open(
              ext_oracles[static_cast<std::size_t>(i)], mu_i);

      if (i == params.d - 1) {
        qp_ext.left[static_cast<std::size_t>(i)].index = mu_i;
        qp_ext.left[static_cast<std::size_t>(i)].value =
            ext_oracles[static_cast<std::size_t>(i + 1)]
                       [static_cast<std::size_t>(mu_i)];
        qp_ext.right[static_cast<std::size_t>(i)].index = mu_i + n_i;
        qp_ext.right[static_cast<std::size_t>(i)].value =
            ext_oracles[static_cast<std::size_t>(i + 1)]
                       [static_cast<std::size_t>(mu_i + n_i)];
      } else {
        qp_ext.left[static_cast<std::size_t>(i)] =
            ext_merkle_trees[static_cast<std::size_t>(i + 1)].Open(
                ext_oracles[static_cast<std::size_t>(i + 1)], mu_i);
        qp_ext.right[static_cast<std::size_t>(i)] =
            ext_merkle_trees[static_cast<std::size_t>(i + 1)].Open(
                ext_oracles[static_cast<std::size_t>(i + 1)], mu_i + n_i);
      }
    }

    const long top_i = params.d - 1;
    const long mu_top = plan.mu_by_level[static_cast<std::size_t>(top_i)];
    const long n_top = CodewordLengthAtLevelNoValidate(params, top_i);
    qp.left[0] = merkle_d.Open(pi_d_base, mu_top);
    qp.right[0] = merkle_d.Open(pi_d_base, mu_top + n_top);

    proof.query_proofs[static_cast<std::size_t>(q)] = std::move(qp);
    proof.extension.query_proofs[static_cast<std::size_t>(q)] =
        std::move(qp_ext);
  }

  return proof;
}

bool VerifyEvalWithExtensionChallenges(
    const MerkleRoot &commitment_C, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const BaseFoldPCSEvalProof &proof, const FoldableCodeParams &params,
    const BaseFoldPCSChallengeConfig &challenge_cfg) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->pcs_verify_ns : nullptr,
                    prof ? &prof->pcs_verify_calls : nullptr);

  if (params.d == 0) {
    return BaseFoldPCSVerifyEval(commitment_C, z, claimed_y, num_queries, proof,
                                 params);
  }

  const ZZ_pEX &extension_modulus = challenge_cfg.challenge_extension_modulus;
  ExtensionDegreeOrThrow(extension_modulus,
                         "VerifyEvalWithExtensionChallenges");

  ValidateParamsOrThrow(params);
  if (!IsPowerOfTwoLong(params.k0))
    return false;
  const long kappa = Log2ExactPowerOfTwoLong(params.k0);
  const long point_dim = params.d + kappa;
  if (static_cast<long>(z.size()) != point_dim)
    return false;
  if (num_queries < 0)
    return false;
  if (!proof.extension.enabled)
    return false;
  if (static_cast<long>(proof.extension.roots_by_level.size()) != params.d)
    return false;
  if (static_cast<long>(proof.extension.h_by_level.size()) != params.d)
    return false;
  if (static_cast<long>(proof.extension.r_by_level.size()) != 0 &&
      static_cast<long>(proof.extension.r_by_level.size()) != params.d)
    return false;
  if (static_cast<long>(proof.extension.msg0_coeffs.size()) != params.k0)
    return false;
  const long n0 = CodewordLengthAtLevel(params, 0);
  if (static_cast<long>(proof.extension.pi0_full.size()) != n0)
    return false;
  if (static_cast<long>(proof.extension.query_proofs.size()) != num_queries)
    return false;
  if (static_cast<long>(proof.query_proofs.size()) != num_queries)
    return false;

  const std::size_t base_root_count = proof.commitments.roots_by_level.size();
  if (base_root_count == 1U) {
    if (proof.commitments.roots_by_level[0] != commitment_C) {
      return false;
    }
  } else if (base_root_count == static_cast<std::size_t>(params.d + 1)) {
    if (proof.commitments.roots_by_level[static_cast<std::size_t>(params.d)] !=
        commitment_C) {
      return false;
    }
  } else if (base_root_count != 0U) {
    return false;
  }

  Sha256Transcript transcript;
  AbsorbPublicInput(transcript, commitment_C, z, claimed_y);
  AbsorbChallengeConfig(transcript, challenge_cfg);

  AbsorbExtensionQuadraticPoly(
      transcript,
      proof.extension.h_by_level[static_cast<std::size_t>(params.d - 1)],
      extension_modulus);

  std::vector<ZZ_pEX> r_by_level;
  r_by_level.resize(static_cast<std::size_t>(params.d));
  const bool has_explicit_r =
      (static_cast<long>(proof.extension.r_by_level.size()) == params.d);

  for (long i = params.d; i-- > 0;) {
    const ZZ_pEX r_i = SampleExtensionChallenge(
        transcript, "r/" + std::to_string(i), extension_modulus);
    r_by_level[static_cast<std::size_t>(i)] = r_i;
    if (has_explicit_r &&
        proof.extension.r_by_level[static_cast<std::size_t>(i)] != r_i) {
      return false;
    }
    transcript.AbsorbDigest(
        proof.extension.roots_by_level[static_cast<std::size_t>(i)]);
    if (i > 0) {
      AbsorbExtensionQuadraticPoly(
          transcript,
          proof.extension.h_by_level[static_cast<std::size_t>(i - 1)],
          extension_modulus);
    }
  }

  if (!CheckExtensionSumcheckRelations(proof.extension.h_by_level, r_by_level,
                                       claimed_y, extension_modulus))
    return false;

  const ZZ_pEX r0 = r_by_level[0];
  const ZZ_pEX h1_r0 = EvalExtensionQuadraticPoly(proof.extension.h_by_level[0],
                                                  r0, extension_modulus);

  ZZ_pEX suffix_eq = ExtensionOne();
  for (long i = 0; i < params.d; ++i) {
    const ZZ_pEX zi_ext =
        LiftBaseToExtension(z[static_cast<std::size_t>(kappa + i)]);
    suffix_eq = MulExtension(
        suffix_eq,
        EqFactorExtension(zi_ext, r_by_level[static_cast<std::size_t>(i)],
                          extension_modulus),
        extension_modulus);
  }

  const std::vector<ZZ_pEX> f_eval =
      BooleanEvalTableFromMonomialCoeffsExtension(proof.extension.msg0_coeffs,
                                                  kappa, extension_modulus);

  std::vector<ZZ_pEX> prefix_eq;
  prefix_eq.resize(f_eval.size());
  if (kappa == 0) {
    prefix_eq[0] = ExtensionOne();
  } else {
    prefix_eq.resize(static_cast<std::size_t>(1L << kappa));
    prefix_eq[0] = ExtensionOne();
    const ZZ_pEX one = ExtensionOne();
    for (long var = 0; var < kappa; ++var) {
      const long old = 1L << var;
      const ZZ_pEX zi = LiftBaseToExtension(z[static_cast<std::size_t>(var)]);
      const ZZ_pEX f0 = SubExtension(one, zi, extension_modulus);
      const ZZ_pEX f1 = zi;
      for (long mask = 0; mask < old; ++mask) {
        const ZZ_pEX base = prefix_eq[static_cast<std::size_t>(mask)];
        prefix_eq[static_cast<std::size_t>(mask)] =
            MulExtension(base, f0, extension_modulus);
        prefix_eq[static_cast<std::size_t>(mask + old)] =
            MulExtension(base, f1, extension_modulus);
      }
    }
  }

  ZZ_pEX sum = ExtensionZero();
  for (long mask = 0; mask < static_cast<long>(f_eval.size()); ++mask) {
    sum = AddExtension(sum,
                       MulExtension(f_eval[static_cast<std::size_t>(mask)],
                                    prefix_eq[static_cast<std::size_t>(mask)],
                                    extension_modulus),
                       extension_modulus);
  }

  if (MulExtension(suffix_eq, sum, extension_modulus) != h1_r0)
    return false;

  const std::vector<ZZ_pEX> expected_pi0 =
      EncodeC0Extension(proof.extension.msg0_coeffs, params, extension_modulus);
  if (expected_pi0 != proof.extension.pi0_full) {
    return false;
  }
  if (proof.extension.roots_by_level[0] !=
      ExtensionMerkleTree::Build(proof.extension.pi0_full, extension_modulus)
          .Root()) {
    return false;
  }

  const long n_last = CodewordLengthAtLevel(params, params.d - 1);
  const long n_d = CodewordLengthAtLevel(params, params.d);
  std::vector<IOPPQueryPlan> query_plans;
  query_plans.resize(static_cast<std::size_t>(num_queries));
  for (long q = 0; q < num_queries; ++q) {
    const long mu =
        transcript.ChallengeIndex("mu/" + std::to_string(q), n_last);
    query_plans[static_cast<std::size_t>(q)] = MakeQueryPlan(mu, params);
  }

  auto verify_one_query = [&](long q) -> bool {
    Profile *query_prof = ActiveProfile();
    ScopedTimer query_timer(
        query_prof ? &query_prof->ext_verify_query_merkle_ns : nullptr,
        query_prof ? &query_prof->ext_verify_query_merkle_calls : nullptr);

    const IOPPQueryPlan &plan = query_plans[static_cast<std::size_t>(q)];
    const BaseFoldPCSQueryProof &qp =
        proof.query_proofs[static_cast<std::size_t>(q)];
    const BaseFoldPCSQueryProofExtension &qp_ext =
        proof.extension.query_proofs[static_cast<std::size_t>(q)];
    if (qp.left.size() != qp.right.size())
      return false;
    const bool qp_full = (static_cast<long>(qp.left.size()) == params.d);
    const bool qp_compact = (qp.left.size() == 1U);
    if (!qp_full && !qp_compact)
      return false;
    if (qp_full) {
      if (static_cast<long>(qp.folded.size()) != params.d)
        return false;
    } else {
      if (!qp.folded.empty())
        return false;
    }
    if (static_cast<long>(qp_ext.left.size()) != params.d)
      return false;
    if (static_cast<long>(qp_ext.right.size()) != params.d)
      return false;
    if (static_cast<long>(qp_ext.folded.size()) != params.d)
      return false;

    const long top_i = params.d - 1;
    const long mu_top = plan.mu_by_level[static_cast<std::size_t>(top_i)];
    const long n_top = CodewordLengthAtLevel(params, top_i);
    const std::size_t top_slot = qp_full ? static_cast<std::size_t>(top_i) : 0U;

    const MerkleOpening &left_top = qp.left[top_slot];
    const MerkleOpening &right_top = qp.right[top_slot];
    if (left_top.index != mu_top || right_top.index != mu_top + n_top) {
      return false;
    }
    if (!MerkleVerifyOpening(commitment_C, n_d, left_top) ||
        !MerkleVerifyOpening(commitment_C, n_d, right_top)) {
      return false;
    }

    if (qp_ext.left[static_cast<std::size_t>(top_i)].index != mu_top ||
        qp_ext.right[static_cast<std::size_t>(top_i)].index != mu_top + n_top) {
      return false;
    }
    if (!qp_ext.left[static_cast<std::size_t>(top_i)]
             .auth_path.sibling_hashes.empty() ||
        !qp_ext.right[static_cast<std::size_t>(top_i)]
             .auth_path.sibling_hashes.empty()) {
      return false;
    }
    if (qp_ext.left[static_cast<std::size_t>(top_i)].value !=
            LiftBaseToExtension(left_top.value) ||
        qp_ext.right[static_cast<std::size_t>(top_i)].value !=
            LiftBaseToExtension(right_top.value)) {
      return false;
    }

    for (long i = params.d; i-- > 0;) {
      const long mu_i = plan.mu_by_level[static_cast<std::size_t>(i)];
      const long n_i = CodewordLengthAtLevel(params, i);
      if (mu_i < 0 || mu_i >= n_i)
        return false;

      const ExtensionMerkleOpening &left_ext =
          qp_ext.left[static_cast<std::size_t>(i)];
      const ExtensionMerkleOpening &right_ext =
          qp_ext.right[static_cast<std::size_t>(i)];
      const ExtensionMerkleOpening &folded_ext =
          qp_ext.folded[static_cast<std::size_t>(i)];

      if (folded_ext.index != mu_i)
        return false;
      if (!VerifyExtensionMerkleOpening(
              proof.extension.roots_by_level[static_cast<std::size_t>(i)], n_i,
              folded_ext, extension_modulus)) {
        return false;
      }

      if (left_ext.index != mu_i || right_ext.index != mu_i + n_i)
        return false;
      if (i < params.d - 1) {
        const long n_ip1 = CodewordLengthAtLevel(params, i + 1);
        if (!VerifyExtensionMerkleOpening(
                proof.extension.roots_by_level[static_cast<std::size_t>(i + 1)],
                n_ip1, left_ext, extension_modulus) ||
            !VerifyExtensionMerkleOpening(
                proof.extension.roots_by_level[static_cast<std::size_t>(i + 1)],
                n_ip1, right_ext, extension_modulus)) {
          return false;
        }
      }

      const FieldElement &t = params.diag_T[static_cast<std::size_t>(i)][mu_i];
      const FieldElement x1 = t;
      const FieldElement x2 = params.zeta * t;
      const ZZ_pEX expected_folded = EvalLineAtExtension(
          r_by_level[static_cast<std::size_t>(i)], x1, left_ext.value, x2,
          right_ext.value, extension_modulus);

      if (expected_folded != folded_ext.value)
        return false;

      if (i > 0) {
        const long n_prev = CodewordLengthAtLevel(params, i - 1);
        if (mu_i < n_prev) {
          if (folded_ext.value !=
              qp_ext.left[static_cast<std::size_t>(i - 1)].value) {
            return false;
          }
        } else {
          if (folded_ext.value !=
              qp_ext.right[static_cast<std::size_t>(i - 1)].value) {
            return false;
          }
        }
      } else {
        if (folded_ext.value !=
            proof.extension.pi0_full[static_cast<std::size_t>(mu_i)]) {
          return false;
        }
      }
    }

    return true;
  };

  if (!VerifyQueriesMaybeParallel(num_queries, prof, verify_one_query)) {
    return false;
  }

  return true;
}

void AbsorbPublicInput(Sha256Transcript &transcript,
                       const MerkleRoot &commitment,
                       const std::vector<FieldElement> &z,
                       const FieldElement &y) {
  transcript.AbsorbDigest(commitment);
  for (const FieldElement &zi : z) {
    transcript.AbsorbFieldElement(zi);
  }
  transcript.AbsorbFieldElement(y);
}

} // namespace

void ResetVerifierQueryParallelConfigFromEnv() {
  MutableVerifierQueryParallelConfig() =
      ParseVerifierQueryParallelConfigFromEnv();
}

void SetVerifierQueryParallelConfig(const VerifierQueryParallelConfig &cfg) {
  MutableVerifierQueryParallelConfig() = cfg;
}

VerifierQueryParallelConfig GetVerifierQueryParallelConfig() {
  return LoadVerifierQueryParallelConfig();
}

MerkleRoot BaseFoldPCSCommit(const vec_ZZ_pE &f_coeffs,
                             const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  if (f_coeffs.length() != MessageLength(params))
    LogicError("BaseFoldPCSCommit: f_coeffs has wrong length");

  Oracle pi_d;
  EncodeFoldable(pi_d, f_coeffs, params);
  return MerkleCommitOracle(pi_d);
}

BaseFoldPCSEvalProof BaseFoldPCSProveEval(const vec_ZZ_pE &f_coeffs,
                                          const std::vector<FieldElement> &z,
                                          const FieldElement &claimed_y,
                                          long num_queries,
                                          const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  const long kappa = Log2ExactPowerOfTwoLong(params.k0);
  const long point_dim = params.d + kappa;
  if (static_cast<long>(z.size()) != point_dim)
    LogicError("BaseFoldPCSProveEval: z has wrong dimension");
  if (f_coeffs.length() != MessageLength(params))
    LogicError("BaseFoldPCSProveEval: f_coeffs has wrong length");
  if (num_queries < 0)
    LogicError("BaseFoldPCSProveEval: num_queries must be non-negative");

  if (EvalMultilinearMonomialCoeffs(f_coeffs, z) != claimed_y) {
    LogicError("BaseFoldPCSProveEval: claimed_y != f(z)");
  }

  return BaseFoldPCSProveEvalUnchecked(f_coeffs, z, claimed_y, num_queries,
                                       params);
}

BaseFoldPCSEvalProof
BaseFoldPCSProveEvalUnchecked(const vec_ZZ_pE &f_coeffs,
                              const std::vector<FieldElement> &z,
                              const FieldElement &claimed_y, long num_queries,
                              const FoldableCodeParams &params) {
  BaseFoldPCSEvalProof proof;
  proof.commitments.roots_by_level.resize(
      static_cast<std::size_t>(params.d + 1));
  proof.h_by_level.resize(static_cast<std::size_t>(params.d));

  IOPPOracles oracles;
  oracles.pi.resize(static_cast<std::size_t>(params.d + 1));

  std::vector<MerkleTree> merkle;
  merkle.resize(static_cast<std::size_t>(params.d + 1));

  EncodeFoldableUnchecked(oracles.pi[static_cast<std::size_t>(params.d)],
                          f_coeffs, params);
  merkle[static_cast<std::size_t>(params.d)] =
      MerkleTree::Build(oracles.pi[static_cast<std::size_t>(params.d)]);
  const MerkleRoot root_d = merkle[static_cast<std::size_t>(params.d)].Root();
  proof.commitments.roots_by_level[static_cast<std::size_t>(params.d)] = root_d;

  Sha256Transcript transcript;
  AbsorbPublicInput(transcript, root_d, z, claimed_y);

  if (params.d == 0) {
    proof.pi0_full = oracles.pi[0];
    proof.query_proofs.resize(static_cast<std::size_t>(num_queries));
    return proof;
  }

  SumcheckProver sumcheck(f_coeffs, z);

  // h_d
  const QuadraticPoly h_d = sumcheck.CurrentPolynomial();
  proof.h_by_level[static_cast<std::size_t>(params.d - 1)] = h_d;
  transcript.AbsorbQuadraticPoly(h_d);

  std::vector<FieldElement> r_by_level;
  r_by_level.resize(static_cast<std::size_t>(params.d));

  for (long i = params.d; i-- > 0;) {
    const FieldElement r_i =
        transcript.ChallengeFieldElement("r/" + std::to_string(i));
    r_by_level[static_cast<std::size_t>(i)] = r_i;

    ProverCommitRoundNoValidate(oracles.pi[static_cast<std::size_t>(i)],
                                oracles.pi[static_cast<std::size_t>(i + 1)],
                                r_i, i, params);

    merkle[static_cast<std::size_t>(i)] =
        MerkleTree::Build(oracles.pi[static_cast<std::size_t>(i)]);
    const MerkleRoot root_i = merkle[static_cast<std::size_t>(i)].Root();
    proof.commitments.roots_by_level[static_cast<std::size_t>(i)] = root_i;
    transcript.AbsorbDigest(root_i);

    sumcheck.ReceiveChallenge(r_i);
    if (i > 0) {
      const QuadraticPoly h_i = sumcheck.CurrentPolynomial();
      proof.h_by_level[static_cast<std::size_t>(i - 1)] = h_i;
      transcript.AbsorbQuadraticPoly(h_i);
    }
  }

  proof.pi0_full = oracles.pi[0];

  const long n_last = CodewordLengthAtLevelNoValidate(params, params.d - 1);
  std::vector<IOPPQueryPlan> query_plans;
  query_plans.resize(static_cast<std::size_t>(num_queries));
  for (long q = 0; q < num_queries; ++q) {
    const long mu =
        transcript.ChallengeIndex("mu/" + std::to_string(q), n_last);
    query_plans[static_cast<std::size_t>(q)] = MakeQueryPlanNoValidate(mu, params);
  }

  proof.query_multiproofs.resize(static_cast<std::size_t>(params.d + 1));
  const std::vector<std::vector<long>> requested_indices_by_tree =
      CollectBaseQueryIndicesByTree(query_plans, params);
  for (long tree_level = 0; tree_level <= params.d; ++tree_level) {
    proof.query_multiproofs[static_cast<std::size_t>(tree_level)] =
        merkle[static_cast<std::size_t>(tree_level)].OpenMany(
            oracles.pi[static_cast<std::size_t>(tree_level)],
            requested_indices_by_tree[static_cast<std::size_t>(tree_level)]);
  }

  return proof;
}

bool BaseFoldPCSVerifyEval(const MerkleRoot &commitment_C,
                           const std::vector<FieldElement> &z,
                           const FieldElement &claimed_y, long num_queries,
                           const BaseFoldPCSEvalProof &proof,
                           const FoldableCodeParams &params) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->pcs_verify_ns : nullptr,
                    prof ? &prof->pcs_verify_calls : nullptr);

  ValidateParamsOrThrow(params);
  if (!IsPowerOfTwoLong(params.k0))
    return false;
  const long kappa = Log2ExactPowerOfTwoLong(params.k0);
  const long point_dim = params.d + kappa;
  if (static_cast<long>(z.size()) != point_dim)
    return false;
  if (num_queries < 0)
    return false;
  if (static_cast<long>(proof.commitments.roots_by_level.size()) !=
      params.d + 1)
    return false;
  if (static_cast<long>(proof.h_by_level.size()) != params.d)
    return false;
  const bool has_query_multiproofs = !proof.query_multiproofs.empty();
  if (has_query_multiproofs) {
    if (static_cast<long>(proof.query_multiproofs.size()) != params.d + 1)
      return false;
  } else if (static_cast<long>(proof.query_proofs.size()) != num_queries) {
    return false;
  }

  if (proof.commitments.roots_by_level[static_cast<std::size_t>(params.d)] !=
      commitment_C) {
    return false;
  }

  const long n0 = CodewordLengthAtLevel(params, 0);
  if (proof.pi0_full.length() != n0)
    return false;

  if (params.d == 0) {
    if (!proof.h_by_level.empty())
      return false;
    if (MerkleCommitOracle(proof.pi0_full) != commitment_C)
      return false;

    vec_ZZ_pE msg0;
    if (!DecodeC0(msg0, proof.pi0_full, params))
      return false;
    if (msg0.length() != params.k0)
      return false;

    return EvalMultilinearMonomialCoeffs(msg0, z) == claimed_y;
  }

  if (MerkleCommitOracle(proof.pi0_full) !=
      proof.commitments.roots_by_level[0]) {
    return false;
  }

  Sha256Transcript transcript;
  AbsorbPublicInput(transcript, commitment_C, z, claimed_y);

  // h_d
  transcript.AbsorbQuadraticPoly(
      proof.h_by_level[static_cast<std::size_t>(params.d - 1)]);

  std::vector<FieldElement> r_by_level;
  r_by_level.resize(static_cast<std::size_t>(params.d));

  for (long i = params.d; i-- > 0;) {
    const FieldElement r_i =
        transcript.ChallengeFieldElement("r/" + std::to_string(i));
    r_by_level[static_cast<std::size_t>(i)] = r_i;

    transcript.AbsorbDigest(
        proof.commitments.roots_by_level[static_cast<std::size_t>(i)]);
    if (i > 0) {
      transcript.AbsorbQuadraticPoly(
          proof.h_by_level[static_cast<std::size_t>(i - 1)]);
    }
  }

  if (!CheckSumcheckRelations(proof.h_by_level, r_by_level, claimed_y))
    return false;

  const FieldElement r0 = r_by_level[0];
  const FieldElement h1_r0 = proof.h_by_level[0].Eval(r0);

  vec_ZZ_pE msg0;
  if (!DecodeC0(msg0, proof.pi0_full, params))
    return false;
  if (msg0.length() != params.k0)
    return false;

  // Check the reduced sumcheck claim (BaseFold paper Remark 3):
  //   h1(r0) == Σ_{b∈{0,1}^κ} f(b, r_suffix) * eq_z(b, r_suffix)
  // where msg0 is the monomial-basis coefficient vector of f(·, r_suffix)
  // on the first κ variables, and r_suffix are the d folding challenges.
  FieldElement suffix_eq;
  NTL::set(suffix_eq);
  for (long i = 0; i < params.d; ++i) {
    suffix_eq *= EqFactor(z[static_cast<std::size_t>(kappa + i)],
                          r_by_level[static_cast<std::size_t>(i)]);
  }

  vec_ZZ_pE f_eval = msg0; // in-place subset-sum transform over κ vars
  for (long bit = 0; bit < kappa; ++bit) {
    const long step = 1L << bit;
    for (long mask = 0; mask < f_eval.length(); ++mask) {
      if (mask & step) {
        f_eval[mask] += f_eval[mask ^ step];
      }
    }
  }

  vec_ZZ_pE prefix_eq;
  prefix_eq.SetLength(f_eval.length());
  if (kappa == 0) {
    prefix_eq[0] = FieldElement(1);
  } else {
    prefix_eq.SetLength(1L << kappa);
    prefix_eq[0] = FieldElement(1);
    for (long var = 0; var < kappa; ++var) {
      const long old = 1L << var;
      const FieldElement zi = z[static_cast<std::size_t>(var)];
      const FieldElement f0 = FieldElement(1) - zi; // EqFactor(zi, 0)
      const FieldElement f1 = zi;                   // EqFactor(zi, 1)
      for (long mask = 0; mask < old; ++mask) {
        const FieldElement base = prefix_eq[mask];
        prefix_eq[mask] = base * f0;
        prefix_eq[mask + old] = base * f1;
      }
    }
  }

  FieldElement sum = FieldElement(0);
  for (long mask = 0; mask < f_eval.length(); ++mask) {
    sum += f_eval[mask] * prefix_eq[mask];
  }

  if (suffix_eq * sum != h1_r0)
    return false;

  IOPPChallenges challenges;
  challenges.alphas = r_by_level;

  if (params.d == 0) {
    return true;
  }

  const long n_last = CodewordLengthAtLevel(params, params.d - 1);
  std::vector<IOPPQueryPlan> query_plans;
  query_plans.resize(static_cast<std::size_t>(num_queries));
  for (long q = 0; q < num_queries; ++q) {
    const long mu =
        transcript.ChallengeIndex("mu/" + std::to_string(q), n_last);
    query_plans[static_cast<std::size_t>(q)] = MakeQueryPlan(mu, params);
  }

  const std::vector<std::vector<long>> requested_indices_by_tree =
      CollectBaseQueryIndicesByTree(query_plans, params);

  if (has_query_multiproofs) {
    ScopedTimer query_timer(prof ? &prof->verify_query_merkle_ns : nullptr,
                            prof ? &prof->verify_query_merkle_calls : nullptr);

    for (long tree_level = 0; tree_level <= params.d; ++tree_level) {
      const MerkleMultiproof &multiproof =
          proof.query_multiproofs[static_cast<std::size_t>(tree_level)];
      const std::vector<long> &expected_indices =
          requested_indices_by_tree[static_cast<std::size_t>(tree_level)];
      if (multiproof.queried_indices != expected_indices) {
        return false;
      }
      if (static_cast<long>(multiproof.values.length()) !=
          static_cast<long>(expected_indices.size())) {
        return false;
      }
      const long leaf_count = CodewordLengthAtLevel(params, tree_level);
      if (!MerkleVerifyMultiproof(
              proof.commitments.roots_by_level[static_cast<std::size_t>(tree_level)],
              leaf_count, multiproof)) {
        return false;
      }
    }

    auto verify_one_query = [&](long q) -> bool {
      const IOPPQueryPlan &plan = query_plans[static_cast<std::size_t>(q)];
      for (long i = params.d; i-- > 0;) {
        const long mu = plan.mu_by_level[static_cast<std::size_t>(i)];
        const long n_i = CodewordLengthAtLevel(params, i);
        if (mu < 0 || mu >= n_i) {
          return false;
        }

        const FieldElement *left = FindMerkleMultiproofValue(
            proof.query_multiproofs[static_cast<std::size_t>(i + 1)], mu);
        const FieldElement *right = FindMerkleMultiproofValue(
            proof.query_multiproofs[static_cast<std::size_t>(i + 1)], mu + n_i);
        const FieldElement *folded = FindMerkleMultiproofValue(
            proof.query_multiproofs[static_cast<std::size_t>(i)], mu);
        if (left == nullptr || right == nullptr || folded == nullptr) {
          return false;
        }

        FieldElement x1, x2;
        FoldingPoints(x1, x2, params, i, mu);
        const FieldElement expected = EvalLineAt(
            challenges.alphas[static_cast<std::size_t>(i)], x1, *left, x2, *right);
        if (expected != *folded) {
          return false;
        }
      }
      return true;
    };

    return VerifyQueriesMaybeParallel(num_queries, prof, verify_one_query);
  }

  const bool use_opening_cache =
      (prof != nullptr) || (ChooseQueryVerifyThreads(num_queries) <= 1);
  std::vector<std::vector<MerkleOpening>> cached_openings_by_tree;
  std::vector<std::vector<unsigned char>> cached_openings_seen_by_tree;
  std::vector<std::vector<unsigned char>> cached_openings_valid_by_tree;
  if (use_opening_cache) {
    cached_openings_by_tree.resize(static_cast<std::size_t>(params.d + 1));
    cached_openings_seen_by_tree.resize(static_cast<std::size_t>(params.d + 1));
    cached_openings_valid_by_tree.resize(static_cast<std::size_t>(params.d + 1));
    for (long tree_level = 0; tree_level <= params.d; ++tree_level) {
      cached_openings_by_tree[static_cast<std::size_t>(tree_level)].resize(
          requested_indices_by_tree[static_cast<std::size_t>(tree_level)].size());
      cached_openings_seen_by_tree[static_cast<std::size_t>(tree_level)]
          .assign(requested_indices_by_tree[static_cast<std::size_t>(tree_level)].size(),
                  static_cast<unsigned char>(0));
      cached_openings_valid_by_tree[static_cast<std::size_t>(tree_level)]
          .assign(requested_indices_by_tree[static_cast<std::size_t>(tree_level)].size(),
                  static_cast<unsigned char>(0));
    }
  }

  auto verify_opening_cached = [&](long tree_level, long leaf_count,
                                   const MerkleRoot &root,
                                   const MerkleOpening &opening) -> bool {
    if (!use_opening_cache) {
      return MerkleVerifyOpening(root, leaf_count, opening);
    }
    const std::size_t pos = FindCachedOpeningPosition(
        requested_indices_by_tree[static_cast<std::size_t>(tree_level)],
        opening.index, "BaseFoldPCSVerifyEval/opening-cache");
    unsigned char &seen =
        cached_openings_seen_by_tree[static_cast<std::size_t>(tree_level)][pos];
    unsigned char &valid =
        cached_openings_valid_by_tree[static_cast<std::size_t>(tree_level)][pos];
    MerkleOpening &cached =
        cached_openings_by_tree[static_cast<std::size_t>(tree_level)][pos];
    if (seen != 0) {
      return valid != 0 && SameMerkleOpening(cached, opening);
    }
    cached = opening;
    seen = static_cast<unsigned char>(1);
    valid = MerkleVerifyOpening(root, leaf_count, opening)
                ? static_cast<unsigned char>(1)
                : static_cast<unsigned char>(0);
    return valid != 0;
  };

  auto verify_one_query = [&](long q) -> bool {
    const IOPPQueryPlan &plan = query_plans[static_cast<std::size_t>(q)];
    const BaseFoldPCSQueryProof &qp =
        proof.query_proofs[static_cast<std::size_t>(q)];
    if (static_cast<long>(qp.left.size()) != params.d)
      return false;
    if (static_cast<long>(qp.right.size()) != params.d)
      return false;
    if (static_cast<long>(qp.folded.size()) != params.d)
      return false;

    for (long i = params.d; i-- > 0;) {
      const long mu = plan.mu_by_level[static_cast<std::size_t>(i)];
      const long n_i = CodewordLengthAtLevel(params, i);
      if (mu < 0 || mu >= n_i)
        return false;

      const MerkleOpening &left = qp.left[static_cast<std::size_t>(i)];
      const MerkleOpening &right = qp.right[static_cast<std::size_t>(i)];
      const MerkleOpening &folded = qp.folded[static_cast<std::size_t>(i)];
      if (left.index != mu || right.index != mu + n_i || folded.index != mu) {
        return false;
      }

      const long n_ip1 = 2 * n_i;
      if (!verify_opening_cached(
              i + 1, n_ip1,
              proof.commitments.roots_by_level[static_cast<std::size_t>(i + 1)],
              left)) {
        return false;
      }
      if (!verify_opening_cached(
              i + 1, n_ip1,
              proof.commitments.roots_by_level[static_cast<std::size_t>(i + 1)],
              right)) {
        return false;
      }
      if (!verify_opening_cached(
              i, n_i,
              proof.commitments.roots_by_level[static_cast<std::size_t>(i)],
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

    return true;
  };

  if (!VerifyQueriesMaybeParallel(num_queries, prof, verify_one_query)) {
    return false;
  }

  return true;
}

BaseFoldPCSEvalProof BaseFoldPCSProveEvalWithChallengeConfig(
    const vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const FoldableCodeParams &params,
    const BaseFoldPCSChallengeConfig &challenge_cfg) {
  ValidateChallengeConfigOrThrow(challenge_cfg, params);
  if (!challenge_cfg.use_extension_challenges) {
    return BaseFoldPCSProveEval(f_coeffs, z, claimed_y, num_queries, params);
  }

  ValidateParamsOrThrow(params);
  const long kappa = Log2ExactPowerOfTwoLong(params.k0);
  const long point_dim = params.d + kappa;
  if (static_cast<long>(z.size()) != point_dim)
    LogicError(
        "BaseFoldPCSProveEvalWithChallengeConfig: z has wrong dimension");
  if (f_coeffs.length() != MessageLength(params))
    LogicError(
        "BaseFoldPCSProveEvalWithChallengeConfig: f_coeffs has wrong length");
  if (num_queries < 0)
    LogicError("BaseFoldPCSProveEvalWithChallengeConfig: num_queries must be "
               "non-negative");
  if (EvalMultilinearMonomialCoeffs(f_coeffs, z) != claimed_y) {
    LogicError("BaseFoldPCSProveEvalWithChallengeConfig: claimed_y != f(z)");
  }

  return ProveEvalWithExtensionChallengesUnchecked(
      f_coeffs, z, claimed_y, num_queries, params, challenge_cfg);
}

BaseFoldPCSEvalProof BaseFoldPCSProveEvalWithChallengeConfigUnchecked(
    const vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const FoldableCodeParams &params,
    const BaseFoldPCSChallengeConfig &challenge_cfg) {
  ValidateChallengeConfigOrThrow(challenge_cfg, params);
  if (!challenge_cfg.use_extension_challenges) {
    return BaseFoldPCSProveEvalUnchecked(f_coeffs, z, claimed_y, num_queries,
                                         params);
  }
  return ProveEvalWithExtensionChallengesUnchecked(
      f_coeffs, z, claimed_y, num_queries, params, challenge_cfg);
}

bool BaseFoldPCSVerifyEvalWithChallengeConfig(
    const MerkleRoot &commitment_C, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const BaseFoldPCSEvalProof &proof, const FoldableCodeParams &params,
    const BaseFoldPCSChallengeConfig &challenge_cfg) {
  ValidateChallengeConfigOrThrow(challenge_cfg, params);
  if (!challenge_cfg.use_extension_challenges) {
    return BaseFoldPCSVerifyEval(commitment_C, z, claimed_y, num_queries, proof,
                                 params);
  }
  return VerifyEvalWithExtensionChallenges(
      commitment_C, z, claimed_y, num_queries, proof, params, challenge_cfg);
}

} // namespace basefold
