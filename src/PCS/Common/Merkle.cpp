#include "PCS/Common/Merkle.hpp"

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

#include "PCS/Common/Hash.hpp"
#include "PCS/Common/MerkleMultiproofPlanner.hpp"
#include "PCS/Common/MerkleMultiproofReplay.hpp"
#include "PCS/Common/Profile.hpp"

#if defined(BASEFOLD_USE_OPENMP)
#include <omp.h>
#endif

using NTL::BytesFromZZ;
using NTL::coeff;
using NTL::LogicError;
using NTL::NumBytes;
using NTL::rep;
using NTL::ZZ;
using NTL::ZZ_pE;
using NTL::ZZ_pX;

namespace basefold {
namespace {

constexpr std::size_t kDigestBytes = Digest{}.size();

long ParsePositiveEnvLong(const char *name, long fallback) {
  const char *raw = std::getenv(name);
  if (raw == nullptr || *raw == '\0') return fallback;
  char *end = nullptr;
  const long value = std::strtol(raw, &end, 10);
  if (end == raw || *end != '\0' || value <= 0) return fallback;
  return value;
}

int ParsePositiveEnvInt(const char *name, int fallback) {
  const char *raw = std::getenv(name);
  if (raw == nullptr || *raw == '\0') return fallback;
  char *end = nullptr;
  const long value = std::strtol(raw, &end, 10);
  if (end == raw || *end != '\0' || value <= 0 ||
      value > static_cast<long>(std::numeric_limits<int>::max())) {
    return fallback;
  }
  return static_cast<int>(value);
}

MerkleBuildParallelConfig ParseMerkleBuildParallelConfigFromEnv() {
  MerkleBuildParallelConfig cfg;
  cfg.leaves_per_thread = ParsePositiveEnvLong(
      "BASEFOLD_MERKLE_LEAVES_PER_THREAD", cfg.leaves_per_thread);
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
  if (cfg.leaves_per_thread <= 0) cfg.leaves_per_thread = 32768;
  if (cfg.parallel_level_threshold <= 0) cfg.parallel_level_threshold = 8192;
  if (cfg.max_threads <= 0) cfg.max_threads = 8;
  return cfg;
}

int ChooseMerkleBuildThreads(long leaf_count,
                             const MerkleBuildParallelConfig &cfg) {
#if defined(BASEFOLD_USE_OPENMP)
  if (leaf_count < cfg.leaves_per_thread) return 1;
  const int max_threads = omp_get_max_threads();
  int threads_to_use = static_cast<int>(leaf_count / cfg.leaves_per_thread);
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

void WriteU64BE(Byte *out, std::uint64_t value) {
  for (int i = 7; i >= 0; --i) {
    out[7 - i] = static_cast<Byte>((value >> (8 * i)) & 0xff);
  }
}

void WriteUIntBE(Byte *out, std::uint64_t value, std::size_t bytes) {
  for (std::size_t i = 0; i < bytes; ++i) {
    out[bytes - 1 - i] = static_cast<Byte>((value >> (8 * i)) & 0xff);
  }
}

struct FieldEncodingInfo {
  long degree = 0;
  long coeff_bytes = 0;
};

FieldEncodingInfo CurrentFieldEncodingInfo() {
  const long degree = ZZ_pE::degree();
  if (degree <= 0) {
    LogicError("CurrentFieldEncodingInfo: invalid extension degree");
  }

  static thread_local ZZ cached_modulus;
  static thread_local long cached_degree = 0;
  static thread_local long cached_coeff_bytes = 0;
  static thread_local bool has_cache = false;

  const ZZ &modulus = NTL::ZZ_p::modulus();
  if (!has_cache || degree != cached_degree || modulus != cached_modulus) {
    cached_modulus = modulus;
    cached_degree = degree;
    const ZZ max_coeff = modulus - 1;
    const long coeff_bytes = NumBytes(max_coeff);
    cached_coeff_bytes = (coeff_bytes > 0) ? coeff_bytes : 1;
    has_cache = true;
  }

  FieldEncodingInfo info;
  info.degree = cached_degree;
  info.coeff_bytes = cached_coeff_bytes;
  return info;
}

Digest HashWithDomainTagByte(Byte domain_tag_byte) {
  return HashDigest(&domain_tag_byte, 1, "HashWithDomainTagByte");
}

Digest HashNode(const Digest &left, const Digest &right) {
  std::array<Byte, 1 + 2 * kDigestBytes> input;
  input[0] = static_cast<Byte>(0x01);
  std::memcpy(input.data() + 1, left.data(), left.size());
  std::memcpy(input.data() + 1 + left.size(), right.data(), right.size());
  return HashDigest(input.data(), input.size(), "HashNode");
}

Digest HashLeaf(long index, const FieldElement &value) {
  const FieldEncodingInfo enc = CurrentFieldEncodingInfo();
  const std::size_t coeff_bytes = static_cast<std::size_t>(enc.coeff_bytes);
  const std::size_t payload_bytes =
      1 + 8 + static_cast<std::size_t>(enc.degree) * coeff_bytes;

  static thread_local std::vector<Byte> input;
  if (input.size() != payload_bytes) {
    input.resize(payload_bytes);
  }
  input[0] = static_cast<Byte>(0x00);
  WriteU64BE(input.data() + 1, static_cast<std::uint64_t>(index));

  const ZZ_pX &poly = rep(value);
  std::size_t offset = 1 + 8;
  for (long i = 0; i < enc.degree; ++i) {
    if (enc.coeff_bytes <= 8) {
      unsigned long coeff_value = 0;
      NTL::conv(coeff_value, rep(coeff(poly, i)));
      WriteUIntBE(input.data() + offset,
                  static_cast<std::uint64_t>(coeff_value), coeff_bytes);
    } else {
      const ZZ &coeff_rep = rep(coeff(poly, i));
      BytesFromZZ(reinterpret_cast<unsigned char *>(input.data() + offset),
                  coeff_rep, enc.coeff_bytes);
    }
    offset += coeff_bytes;
  }

  return HashDigest(input.data(), input.size(), "HashLeaf");
}

Digest HashRootWithCount(long leaf_count, const Digest &raw_root) {
  std::array<Byte, 1 + 8 + kDigestBytes> input;
  input[0] = static_cast<Byte>(0x03);
  WriteU64BE(input.data() + 1, static_cast<std::uint64_t>(leaf_count));
  std::memcpy(input.data() + 1 + 8, raw_root.data(), raw_root.size());
  return HashDigest(input.data(), input.size(), "HashRootWithCount");
}

std::size_t ExpectedMerkleHeight(long leaf_count) {
  if (leaf_count <= 0) return 0;
  std::size_t height = 0;
  long width = leaf_count;
  while (width > 1) {
    if ((width & 1L) != 0) width += 1;
    width /= 2;
    ++height;
  }
  return height;
}

using MerkleMultiproofPlan = multiproof_planner::Plan;

Digest MerkleRootRaw(std::vector<Digest> level) {
  if (level.empty()) {
    return HashWithDomainTagByte(static_cast<Byte>(0x04));
  }

  while (level.size() > 1) {
    if ((level.size() & 1U) != 0U) level.push_back(level.back());
    std::vector<Digest> next;
    next.reserve(level.size() / 2U);
    for (std::size_t i = 0; i < level.size(); i += 2U) {
      next.push_back(HashNode(level[i], level[i + 1]));
    }
    level = std::move(next);
  }
  return level[0];
}

bool MerkleVerifyMultiproofNoProfile(const MerkleRoot &root, long leaf_count,
                                     const std::vector<long> &queried_indices,
                                     const MerkleMultiproof &proof) {
  if (leaf_count < 0) return false;
  if (static_cast<long>(proof.values.length()) !=
      static_cast<long>(queried_indices.size())) {
    return false;
  }

  if (leaf_count == 0) {
    if (!queried_indices.empty() || proof.values.length() != 0 ||
        !proof.sibling_hashes.empty()) {
      return false;
    }
    return root ==
           HashRootWithCount(0, HashWithDomainTagByte(static_cast<Byte>(0x04)));
  }

  if (!multiproof_planner::IsSortedUniqueIndicesInRange(
          leaf_count, queried_indices)) {
    return false;
  }

  const MerkleMultiproofPlan plan =
      multiproof_planner::BuildPlanFromSortedUnique(leaf_count,
                                                    queried_indices);
  if (proof.sibling_hashes.size() !=
      static_cast<std::size_t>(plan.stats.unique_sibling_count)) {
    return false;
  }

  if (queried_indices.empty()) {
    return proof.sibling_hashes.empty();
  }

  std::vector<std::pair<long, Digest>> current(queried_indices.size());
  for (std::size_t i = 0; i < queried_indices.size(); ++i) {
    current[i] = {queried_indices[i],
                  HashLeaf(queried_indices[i],
                           proof.values[static_cast<long>(i)])};
  }

  Digest raw_root;
  if (!multiproof_replay::ReplayPlanToRawRoot(
          plan, proof.sibling_hashes, std::move(current),
          [](const Digest &lhs, const Digest &rhs) {
            return HashNode(lhs, rhs);
          },
          &raw_root)) {
    return false;
  }
  return HashRootWithCount(leaf_count, raw_root) == root;
}

}  // namespace

void ResetMerkleBuildParallelConfigFromEnv() {
  MutableMerkleBuildParallelConfig() = ParseMerkleBuildParallelConfigFromEnv();
}

void SetMerkleBuildParallelConfig(const MerkleBuildParallelConfig &cfg) {
  MutableMerkleBuildParallelConfig() = cfg;
}

MerkleBuildParallelConfig GetMerkleBuildParallelConfig() {
  return LoadMerkleBuildParallelConfig();
}

MerkleRoot MerkleCommitOracle(const Oracle &oracle) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->merkle_commit_oracle_ns : nullptr,
                    prof ? &prof->merkle_commit_oracle_calls : nullptr);

  const long leaf_count = oracle.length();
  if (leaf_count < 0) {
    LogicError("MerkleCommitOracle: invalid leaf count");
  }
  if (leaf_count == 0) {
    return HashRootWithCount(0,
                             HashWithDomainTagByte(static_cast<Byte>(0x04)));
  }

  std::vector<Digest> leaf_hashes;
  leaf_hashes.reserve(static_cast<std::size_t>(leaf_count));
  for (long i = 0; i < leaf_count; ++i) {
    leaf_hashes.push_back(HashLeaf(i, oracle[i]));
  }
  const Digest raw_root = MerkleRootRaw(std::move(leaf_hashes));
  return HashRootWithCount(leaf_count, raw_root);
}

MerkleMultiproof MerkleOpenOracleMany(const Oracle &oracle,
                                      const std::vector<long> &queried_indices) {
  const MerkleTree tree = MerkleTree::Build(oracle);
  return tree.OpenMany(oracle, queried_indices);
}

MerkleMultiproofStats PlanMerkleMultiproof(
    long leaf_count, const std::vector<long> &queried_indices) {
  const std::vector<long> unique =
      multiproof_planner::SortAndValidateIndicesOrThrow(
          leaf_count, queried_indices, "PlanMerkleMultiproof");
  return multiproof_planner::BuildPlanFromSortedUnique(leaf_count, unique)
      .stats;
}

bool MerkleVerifyMultiproof(const MerkleRoot &root, long leaf_count,
                            const MerkleMultiproof &proof) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->merkle_verify_opening_ns : nullptr,
                    prof ? &prof->merkle_verify_opening_calls : nullptr);
  return MerkleVerifyMultiproofNoProfile(root, leaf_count,
                                         proof.queried_indices, proof);
}

bool MerkleVerifyMultiproof(const MerkleRoot &root, long leaf_count,
                            const std::vector<long> &queried_indices,
                            const MerkleMultiproof &proof) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->merkle_verify_opening_ns : nullptr,
                    prof ? &prof->merkle_verify_opening_calls : nullptr);
  return MerkleVerifyMultiproofNoProfile(root, leaf_count, queried_indices,
                                         proof);
}

MerkleTree MerkleTree::Build(const Oracle &oracle) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->merkle_tree_build_ns : nullptr,
                    prof ? &prof->merkle_tree_build_calls : nullptr);

  const long leaf_count = oracle.length();
  if (leaf_count < 0) {
    LogicError("MerkleTree::Build: invalid leaf count");
  }

  MerkleTree tree;
  tree.leaf_count_ = leaf_count;
  tree.levels_.clear();
  tree.levels_.reserve(ExpectedMerkleHeight(leaf_count));
  tree.raw_root_ = Digest{};

  if (leaf_count == 0) {
    tree.raw_root_ = HashWithDomainTagByte(static_cast<Byte>(0x04));
    return tree;
  }

  std::vector<Digest> level(static_cast<std::size_t>(leaf_count));
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

#pragma omp parallel num_threads(threads_to_use) shared(level, next, next_size, done, parallel_level, tree)
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
            if ((level.size() & 1U) != 0U) {
              level.push_back(level.back());
            }
            next_size = static_cast<long>(level.size() / 2U);
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
              next[static_cast<std::size_t>(i)] =
                  HashNode(level[j], level[j + 1]);
            }
          }
        }

#pragma omp single
        {
          tree.levels_.push_back(std::move(level));
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
      if ((level.size() & 1U) != 0U) level.push_back(level.back());

      const long next_size = static_cast<long>(level.size() / 2U);
      std::vector<Digest> next(static_cast<std::size_t>(next_size));
      for (long i = 0; i < next_size; ++i) {
        const std::size_t j = static_cast<std::size_t>(2 * i);
        next[static_cast<std::size_t>(i)] = HashNode(level[j], level[j + 1]);
      }
      tree.levels_.push_back(std::move(level));
      level = std::move(next);
    }
  }

  tree.raw_root_ = level[0];
  return tree;
}

MerkleRoot MerkleTree::Root() const {
  return HashRootWithCount(leaf_count_, raw_root_);
}

MerkleMultiproof MerkleTree::OpenMany(
    const Oracle &oracle, const std::vector<long> &queried_indices) const {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->merkle_tree_open_ns : nullptr,
                    prof ? &prof->merkle_tree_open_calls : nullptr);

  if (oracle.length() != leaf_count_) {
    LogicError("MerkleTree::OpenMany: oracle length mismatch");
  }
  const std::vector<long> unique =
      multiproof_planner::SortAndValidateIndicesOrThrow(
          leaf_count_, queried_indices, "MerkleTree::OpenMany");

  MerkleMultiproof proof;
  proof.queried_indices = unique;
  proof.values.SetLength(static_cast<long>(unique.size()));
  if (unique.empty()) {
    return proof;
  }

  const MerkleMultiproofPlan plan =
      multiproof_planner::BuildPlanFromSortedUnique(leaf_count_, unique);
  for (std::size_t i = 0; i < unique.size(); ++i) {
    proof.values[static_cast<long>(i)] = oracle[unique[i]];
  }
  proof.sibling_hashes = multiproof_replay::CollectSiblingHashesForPlan(
      plan, [&](std::size_t level_index, long sibling) {
        return levels_[level_index][static_cast<std::size_t>(sibling)];
      });
  return proof;
}

}  // namespace basefold
