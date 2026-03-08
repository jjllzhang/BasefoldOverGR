#ifndef BASEFOLD_MERKLEMULTIPROOFREPLAY_HPP_
#define BASEFOLD_MERKLEMULTIPROOFREPLAY_HPP_

#include <cstddef>
#include <utility>
#include <vector>

#include "BaseFold/MerkleMultiproofPlanner.hpp"

namespace basefold {
namespace multiproof_replay {

template <typename ValueAtFn>
inline auto FindMultiproofValue(
    const std::vector<long> &queried_indices,
    long value_count, long index,
    const ValueAtFn &value_at) -> decltype(value_at(std::size_t{})) {
  if (value_count != static_cast<long>(queried_indices.size())) {
    return nullptr;
  }
  const auto it =
      std::lower_bound(queried_indices.begin(), queried_indices.end(), index);
  if (it == queried_indices.end() || *it != index) {
    return nullptr;
  }
  const std::size_t pos =
      static_cast<std::size_t>(it - queried_indices.begin());
  return value_at(pos);
}

template <typename GetSiblingHashFn>
inline std::vector<Digest> CollectSiblingHashesForPlan(
    const multiproof_planner::Plan &plan,
    const GetSiblingHashFn &get_sibling_hash) {
  std::vector<Digest> sibling_hashes;
  sibling_hashes.reserve(
      static_cast<std::size_t>(plan.stats.unique_sibling_count));
  for (std::size_t level_index = 0; level_index < plan.levels.size();
       ++level_index) {
    for (long sibling : plan.levels[level_index].sibling_indices) {
      sibling_hashes.push_back(get_sibling_hash(level_index, sibling));
    }
  }
  return sibling_hashes;
}

template <typename NodeHashFn>
inline bool ReplayPlanToRawRoot(
    const multiproof_planner::Plan &plan,
    const std::vector<Digest> &sibling_hashes,
    std::vector<std::pair<long, Digest>> current,
    const NodeHashFn &hash_node,
    Digest *raw_root_out) {
  if (current.size() != plan.queried_indices.size()) {
    return false;
  }
  for (std::size_t i = 0; i < current.size(); ++i) {
    if (current[i].first != plan.queried_indices[i]) {
      return false;
    }
  }
  if (sibling_hashes.size() !=
      static_cast<std::size_t>(plan.stats.unique_sibling_count)) {
    return false;
  }
  if (plan.queried_indices.empty()) {
    return sibling_hashes.empty();
  }

  std::size_t sibling_cursor = 0;
  for (const multiproof_planner::PlanLevel &level : plan.levels) {
    if (current.size() != level.node_indices.size()) {
      return false;
    }
    for (std::size_t i = 0; i < current.size(); ++i) {
      if (current[i].first != level.node_indices[i]) {
        return false;
      }
    }

    std::vector<std::pair<long, Digest>> parents(level.parent_indices.size());
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
        parent_hash = hash_node(hash, current[i + 1U].second);
        i += 2U;
      } else {
        if (sibling_cursor >= sibling_hashes.size()) {
          return false;
        }
        const Digest &sibling = sibling_hashes[sibling_cursor++];
        parent_hash = ((index & 1L) == 0L) ? hash_node(hash, sibling)
                                           : hash_node(sibling, hash);
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

  if (sibling_cursor != sibling_hashes.size() || current.size() != 1U) {
    return false;
  }
  if (raw_root_out != nullptr) {
    *raw_root_out = current.front().second;
  }
  return true;
}

}  // namespace multiproof_replay
}  // namespace basefold

#endif  // BASEFOLD_MERKLEMULTIPROOFREPLAY_HPP_
