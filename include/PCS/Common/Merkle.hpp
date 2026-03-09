#ifndef BASEFOLD_MERKLE_HPP_
#define BASEFOLD_MERKLE_HPP_

#include <NTL/ZZ_pE.h>
#include <NTL/vec_ZZ_pE.h>

#include <cstdint>
#include <vector>

#include "PCS/Common/Hash.hpp"

namespace basefold {

using FieldElement = NTL::ZZ_pE;
using FieldVec = NTL::vec_ZZ_pE;
using Oracle = NTL::vec_ZZ_pE;

struct MerkleMultiproofStats {
  std::uint64_t opened_leaf_count = 0;
  std::uint64_t unique_sibling_count = 0;
  std::uint64_t verifier_hashes = 0;
};

struct MerkleMultiproof {
  std::vector<long> queried_indices;
  Oracle values;
  std::vector<Digest> sibling_hashes;
};

using MerkleRoot = Digest;

struct MerkleBuildParallelConfig {
  long leaves_per_thread = 32768;
  long parallel_level_threshold = 8192;
  int max_threads = 8;
};

void ResetMerkleBuildParallelConfigFromEnv();
void SetMerkleBuildParallelConfig(const MerkleBuildParallelConfig &cfg);
MerkleBuildParallelConfig GetMerkleBuildParallelConfig();

MerkleRoot MerkleCommitOracle(const Oracle &oracle);

MerkleMultiproof MerkleOpenOracleMany(const Oracle &oracle,
                                      const std::vector<long> &queried_indices);

MerkleMultiproofStats PlanMerkleMultiproof(
    long leaf_count, const std::vector<long> &queried_indices);

bool MerkleVerifyMultiproof(const MerkleRoot &root, long leaf_count,
                            const MerkleMultiproof &proof);

bool MerkleVerifyMultiproof(const MerkleRoot &root, long leaf_count,
                            const std::vector<long> &queried_indices,
                            const MerkleMultiproof &proof);

class MerkleTree {
 public:
  MerkleTree() = default;

  static MerkleTree Build(const Oracle &oracle);

  long LeafCount() const { return leaf_count_; }

  MerkleRoot Root() const;

  MerkleMultiproof OpenMany(const Oracle &oracle,
                            const std::vector<long> &queried_indices) const;

 private:
  long leaf_count_ = 0;
  Digest raw_root_{};
  std::vector<std::vector<Digest>> levels_;
};

}  // namespace basefold

#endif  // BASEFOLD_MERKLE_HPP_
