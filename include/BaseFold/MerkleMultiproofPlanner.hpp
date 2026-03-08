#ifndef BASEFOLD_MERKLEMULTIPROOFPLANNER_HPP_
#define BASEFOLD_MERKLEMULTIPROOFPLANNER_HPP_

#include <algorithm>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include <NTL/ZZ_pE.h>

#include "BaseFold/IOPP.hpp"

namespace basefold {
namespace multiproof_planner {

struct PlanLevel {
  std::vector<long> node_indices;
  std::vector<long> sibling_indices;
  std::vector<long> parent_indices;
};

struct Plan {
  long tree_leaf_count = 0;
  long padded_leaf_count = 0;
  std::vector<long> queried_indices;
  std::vector<PlanLevel> levels;
  MerkleMultiproofStats stats;
};

inline long NextPowerOfTwoLong(long x) {
  if (x <= 1)
    return 1;
  long value = 1;
  while (value < x) {
    if (value > std::numeric_limits<long>::max() / 2) {
      NTL::LogicError("NextPowerOfTwoLong: overflow");
    }
    value *= 2;
  }
  return value;
}

inline bool IsSortedUniqueIndicesInRange(long leaf_count,
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

inline std::vector<long> SortAndValidateIndicesOrThrow(
    long leaf_count, const std::vector<long> &queried_indices,
    const char *func_name) {
  if (leaf_count < 0) {
    const std::string msg = std::string(func_name) + ": invalid leaf count";
    NTL::LogicError(msg.c_str());
  }

  std::vector<long> unique = queried_indices;
  std::sort(unique.begin(), unique.end());
  unique.erase(std::unique(unique.begin(), unique.end()), unique.end());
  for (long index : unique) {
    if (index < 0 || index >= leaf_count) {
      const std::string msg = std::string(func_name) + ": index out of range";
      NTL::LogicError(msg.c_str());
    }
  }
  return unique;
}

inline Plan BuildPlanFromSortedUnique(
    long leaf_count, const std::vector<long> &queried_indices) {
  Plan plan;
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
    PlanLevel level;
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

}  // namespace multiproof_planner
}  // namespace basefold

#endif  // BASEFOLD_MERKLEMULTIPROOFPLANNER_HPP_
