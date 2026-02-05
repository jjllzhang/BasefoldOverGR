#ifndef BASEFOLD_PROFILE_HPP_
#define BASEFOLD_PROFILE_HPP_

#include <chrono>
#include <cstdint>
#include <iomanip>
#include <ostream>

namespace basefold {

// Lightweight, opt-in profiler for breaking down verifier time.
//
// Usage (typical in benches):
//   basefold::Profile prof;
//   basefold::ResetProfile(prof);
//   {
//     basefold::ProfileGuard guard(&prof);
//     ... call verifier ...
//   }
//   basefold::PrintProfile(std::cout, prof);
struct Profile {
  std::uint64_t pcs_verify_ns = 0;
  std::uint64_t pcs_verify_calls = 0;

  std::uint64_t verify_query_merkle_ns = 0;
  std::uint64_t verify_query_merkle_calls = 0;

  std::uint64_t merkle_verify_opening_ns = 0;
  std::uint64_t merkle_verify_opening_calls = 0;

  std::uint64_t merkle_commit_oracle_ns = 0;
  std::uint64_t merkle_commit_oracle_calls = 0;

  std::uint64_t eval_line_at_ns = 0;
  std::uint64_t eval_line_at_calls = 0;

  std::uint64_t try_invert_unit_ns = 0;
  std::uint64_t try_invert_unit_calls = 0;
  std::uint64_t try_invert_unit_cache_hits = 0;

  std::uint64_t is_unit_ns = 0;
  std::uint64_t is_unit_calls = 0;

  std::uint64_t inv_fallback_ns = 0;
  std::uint64_t inv_fallback_calls = 0;
};

extern thread_local Profile *g_active_profile;

inline Profile *ActiveProfile() { return g_active_profile; }
inline void SetActiveProfile(Profile *p) { g_active_profile = p; }

inline void ResetProfile(Profile &p) { p = Profile{}; }

inline std::uint64_t NowNs() {
  return static_cast<std::uint64_t>(
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          std::chrono::steady_clock::now().time_since_epoch())
          .count());
}

class ScopedTimer {
 public:
  ScopedTimer(std::uint64_t *ns_acc, std::uint64_t *calls_acc = nullptr)
      : ns_acc_(ns_acc), calls_acc_(calls_acc) {
    if (ns_acc_ != nullptr) start_ns_ = NowNs();
  }

  ~ScopedTimer() {
    if (ns_acc_ != nullptr) *ns_acc_ += (NowNs() - start_ns_);
    if (calls_acc_ != nullptr) ++(*calls_acc_);
  }

  ScopedTimer(const ScopedTimer &) = delete;
  ScopedTimer &operator=(const ScopedTimer &) = delete;

 private:
  std::uint64_t start_ns_ = 0;
  std::uint64_t *ns_acc_ = nullptr;
  std::uint64_t *calls_acc_ = nullptr;
};

class ProfileGuard {
 public:
  explicit ProfileGuard(Profile *p) : prev_(g_active_profile) {
    g_active_profile = p;
  }

  ~ProfileGuard() { g_active_profile = prev_; }

  ProfileGuard(const ProfileGuard &) = delete;
  ProfileGuard &operator=(const ProfileGuard &) = delete;

 private:
  Profile *prev_ = nullptr;
};

inline double NsToMs(std::uint64_t ns) {
  return static_cast<double>(ns) / 1e6;
}

inline void PrintProfile(std::ostream &os, const Profile &p) {
  const double total_ms = NsToMs(p.pcs_verify_ns);
  const double query_ms = NsToMs(p.verify_query_merkle_ns);
  const double merkle_open_ms = NsToMs(p.merkle_verify_opening_ns);
  const double merkle_commit_ms = NsToMs(p.merkle_commit_oracle_ns);
  const double eval_line_ms = NsToMs(p.eval_line_at_ns);
  const double inv_ms = NsToMs(p.try_invert_unit_ns);
  const double is_unit_ms = NsToMs(p.is_unit_ns);
  const double inv_fallback_ms = NsToMs(p.inv_fallback_ns);

  const std::uint64_t query_accounted_ns =
      p.merkle_verify_opening_ns + p.merkle_commit_oracle_ns + p.eval_line_at_ns;
  const double query_other_ms =
      (p.verify_query_merkle_ns > query_accounted_ns)
          ? NsToMs(p.verify_query_merkle_ns - query_accounted_ns)
          : 0.0;
  const double outside_query_ms =
      (p.pcs_verify_ns > p.verify_query_merkle_ns)
          ? NsToMs(p.pcs_verify_ns - p.verify_query_merkle_ns)
          : 0.0;

  os << std::fixed << std::setprecision(3);
  os << "  [profile]\n";
  os << "    BaseFoldPCSVerifyEval total: " << total_ms << " ms"
     << "  (calls " << p.pcs_verify_calls << ")\n";
  os << "    VerifyQueryFromMerkleOpenings: " << query_ms << " ms"
     << "  (calls " << p.verify_query_merkle_calls << ")\n";
  os << "    MerkleVerifyOpening:         " << merkle_open_ms << " ms"
     << "  (calls " << p.merkle_verify_opening_calls << ")\n";
  os << "    MerkleCommitOracle:          " << merkle_commit_ms << " ms"
     << "  (calls " << p.merkle_commit_oracle_calls << ")\n";
  os << "    EvalLineAt:                  " << eval_line_ms << " ms"
     << "  (calls " << p.eval_line_at_calls << ")\n";
  os << "    TryInvertUnit (subset):      " << inv_ms << " ms"
     << "  (calls " << p.try_invert_unit_calls << ", cache hits "
     << p.try_invert_unit_cache_hits << ")\n";
  os << "      IsUnit (subset):           " << is_unit_ms << " ms"
     << "  (calls " << p.is_unit_calls << ")\n";
  os << "      Inv fallback (subset):     " << inv_fallback_ms << " ms"
     << "  (calls " << p.inv_fallback_calls << ")\n";
  os << "    Inside queries other:        " << query_other_ms << " ms\n";
  os << "    Outside queries:             " << outside_query_ms << " ms\n";
}

}  // namespace basefold

#endif  // BASEFOLD_PROFILE_HPP_
