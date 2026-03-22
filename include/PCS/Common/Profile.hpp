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
  std::uint64_t pcs_prove_ns = 0;
  std::uint64_t pcs_prove_calls = 0;

  std::uint64_t pcs_verify_ns = 0;
  std::uint64_t pcs_verify_calls = 0;

  std::uint64_t encode_foldable_unchecked_ns = 0;
  std::uint64_t encode_foldable_unchecked_calls = 0;

  std::uint64_t sumcheck_init_ns = 0;
  std::uint64_t sumcheck_init_calls = 0;

  std::uint64_t sumcheck_current_poly_ns = 0;
  std::uint64_t sumcheck_current_poly_calls = 0;

  std::uint64_t sumcheck_receive_challenge_ns = 0;
  std::uint64_t sumcheck_receive_challenge_calls = 0;

  std::uint64_t prover_commit_round_ns = 0;
  std::uint64_t prover_commit_round_calls = 0;

  std::uint64_t merkle_tree_build_ns = 0;
  std::uint64_t merkle_tree_build_calls = 0;

  std::uint64_t merkle_tree_open_ns = 0;
  std::uint64_t merkle_tree_open_calls = 0;

  std::uint64_t transcript_absorb_ns = 0;
  std::uint64_t transcript_absorb_calls = 0;

  std::uint64_t transcript_challenge_ns = 0;
  std::uint64_t transcript_challenge_calls = 0;

  std::uint64_t z2k_backend_prove_ns = 0;
  std::uint64_t z2k_backend_prove_calls = 0;

  std::uint64_t z2k_backend_verify_ns = 0;
  std::uint64_t z2k_backend_verify_calls = 0;

  std::uint64_t ring_switch_outer_prove_total_ns = 0;
  std::uint64_t ring_switch_outer_prove_total_calls = 0;

  std::uint64_t ring_switch_outer_prove_tensor_build_ns = 0;
  std::uint64_t ring_switch_outer_prove_tensor_build_calls = 0;

  std::uint64_t ring_switch_outer_prove_compute_s_ns = 0;
  std::uint64_t ring_switch_outer_prove_compute_s_calls = 0;

  std::uint64_t ring_switch_outer_prove_recover_partials_ns = 0;
  std::uint64_t ring_switch_outer_prove_recover_partials_calls = 0;

  std::uint64_t ring_switch_outer_prove_batch_prep_ns = 0;
  std::uint64_t ring_switch_outer_prove_batch_prep_calls = 0;

  std::uint64_t ring_switch_outer_prove_sumcheck_ns = 0;
  std::uint64_t ring_switch_outer_prove_sumcheck_calls = 0;

  std::uint64_t ring_switch_outer_prove_final_check_ns = 0;
  std::uint64_t ring_switch_outer_prove_final_check_calls = 0;

  std::uint64_t ring_switch_outer_verify_total_ns = 0;
  std::uint64_t ring_switch_outer_verify_total_calls = 0;

  std::uint64_t ring_switch_outer_verify_tensor_build_ns = 0;
  std::uint64_t ring_switch_outer_verify_tensor_build_calls = 0;

  std::uint64_t ring_switch_outer_verify_recover_partials_ns = 0;
  std::uint64_t ring_switch_outer_verify_recover_partials_calls = 0;

  std::uint64_t ring_switch_outer_verify_prefix_replay_ns = 0;
  std::uint64_t ring_switch_outer_verify_prefix_replay_calls = 0;

  std::uint64_t ring_switch_outer_verify_sumcheck_replay_ns = 0;
  std::uint64_t ring_switch_outer_verify_sumcheck_replay_calls = 0;

  std::uint64_t ring_switch_outer_verify_final_check_ns = 0;
  std::uint64_t ring_switch_outer_verify_final_check_calls = 0;

  std::uint64_t frobenius_outer_prove_total_ns = 0;
  std::uint64_t frobenius_outer_prove_total_calls = 0;

  std::uint64_t frobenius_outer_prove_orbit_build_ns = 0;
  std::uint64_t frobenius_outer_prove_orbit_build_calls = 0;

  std::uint64_t frobenius_outer_prove_compute_s_ns = 0;
  std::uint64_t frobenius_outer_prove_compute_s_calls = 0;

  std::uint64_t frobenius_outer_prove_recover_partials_ns = 0;
  std::uint64_t frobenius_outer_prove_recover_partials_calls = 0;

  std::uint64_t frobenius_outer_prove_batch_prep_ns = 0;
  std::uint64_t frobenius_outer_prove_batch_prep_calls = 0;

  std::uint64_t frobenius_outer_prove_sumcheck_ns = 0;
  std::uint64_t frobenius_outer_prove_sumcheck_calls = 0;

  std::uint64_t frobenius_outer_prove_final_check_ns = 0;
  std::uint64_t frobenius_outer_prove_final_check_calls = 0;

  std::uint64_t frobenius_outer_verify_total_ns = 0;
  std::uint64_t frobenius_outer_verify_total_calls = 0;

  std::uint64_t frobenius_outer_verify_orbit_build_ns = 0;
  std::uint64_t frobenius_outer_verify_orbit_build_calls = 0;

  std::uint64_t frobenius_outer_verify_recover_partials_ns = 0;
  std::uint64_t frobenius_outer_verify_recover_partials_calls = 0;

  std::uint64_t frobenius_outer_verify_prefix_replay_ns = 0;
  std::uint64_t frobenius_outer_verify_prefix_replay_calls = 0;

  std::uint64_t frobenius_outer_verify_sumcheck_replay_ns = 0;
  std::uint64_t frobenius_outer_verify_sumcheck_replay_calls = 0;

  std::uint64_t frobenius_outer_verify_final_check_ns = 0;
  std::uint64_t frobenius_outer_verify_final_check_calls = 0;

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

  std::uint64_t ext_sumcheck_init_ns = 0;
  std::uint64_t ext_sumcheck_init_calls = 0;

  std::uint64_t ext_sumcheck_current_poly_ns = 0;
  std::uint64_t ext_sumcheck_current_poly_calls = 0;

  std::uint64_t ext_sumcheck_receive_challenge_ns = 0;
  std::uint64_t ext_sumcheck_receive_challenge_calls = 0;

  std::uint64_t ext_prover_commit_round_ns = 0;
  std::uint64_t ext_prover_commit_round_calls = 0;

  std::uint64_t ext_merkle_tree_build_ns = 0;
  std::uint64_t ext_merkle_tree_build_calls = 0;

  std::uint64_t ext_merkle_tree_open_ns = 0;
  std::uint64_t ext_merkle_tree_open_calls = 0;

  std::uint64_t ext_msg0_coeffs_suffix_fold_ns = 0;
  std::uint64_t ext_msg0_coeffs_suffix_fold_calls = 0;

  std::uint64_t ext_verify_query_merkle_ns = 0;
  std::uint64_t ext_verify_query_merkle_calls = 0;

  std::uint64_t ext_merkle_verify_opening_ns = 0;
  std::uint64_t ext_merkle_verify_opening_calls = 0;

  std::uint64_t ext_eval_line_at_ns = 0;
  std::uint64_t ext_eval_line_at_calls = 0;
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
  os << std::fixed << std::setprecision(3);

  if (p.pcs_prove_calls > 0) {
    const double total_ms = NsToMs(p.pcs_prove_ns);
    const double encode_ms = NsToMs(p.encode_foldable_unchecked_ns);
    const double sumcheck_init_ms = NsToMs(p.sumcheck_init_ns);
    const double sumcheck_poly_ms = NsToMs(p.sumcheck_current_poly_ns);
    const double sumcheck_recv_ms = NsToMs(p.sumcheck_receive_challenge_ns);
    const double sumcheck_total_ms =
        sumcheck_init_ms + sumcheck_poly_ms + sumcheck_recv_ms;
    const double commit_ms = NsToMs(p.prover_commit_round_ns);
    const double merkle_build_ms = NsToMs(p.merkle_tree_build_ns);
    const double merkle_open_ms = NsToMs(p.merkle_tree_open_ns);
    const double absorb_ms = NsToMs(p.transcript_absorb_ns);
    const double challenge_ms = NsToMs(p.transcript_challenge_ns);
    const double z2k_backend_prove_ms = NsToMs(p.z2k_backend_prove_ns);
    const double ext_sumcheck_init_ms = NsToMs(p.ext_sumcheck_init_ns);
    const double ext_sumcheck_poly_ms = NsToMs(p.ext_sumcheck_current_poly_ns);
    const double ext_sumcheck_recv_ms =
        NsToMs(p.ext_sumcheck_receive_challenge_ns);
    const double ext_sumcheck_total_ms =
        ext_sumcheck_init_ms + ext_sumcheck_poly_ms + ext_sumcheck_recv_ms;
    const double ext_commit_ms = NsToMs(p.ext_prover_commit_round_ns);
    const double ext_merkle_build_ms = NsToMs(p.ext_merkle_tree_build_ns);
    const double ext_merkle_open_ms = NsToMs(p.ext_merkle_tree_open_ns);
    const double ext_msg0_suffix_fold_ms =
        NsToMs(p.ext_msg0_coeffs_suffix_fold_ns);
    const bool has_ext_prover_breakdown =
        (p.ext_sumcheck_init_calls + p.ext_sumcheck_current_poly_calls +
         p.ext_sumcheck_receive_challenge_calls +
         p.ext_prover_commit_round_calls + p.ext_merkle_tree_build_calls +
         p.ext_merkle_tree_open_calls +
         p.ext_msg0_coeffs_suffix_fold_calls) > 0;

    const std::uint64_t accounted_ns = p.encode_foldable_unchecked_ns +
                                      p.sumcheck_init_ns +
                                      p.sumcheck_current_poly_ns +
                                      p.sumcheck_receive_challenge_ns +
                                      p.prover_commit_round_ns +
                                      p.merkle_tree_build_ns +
                                      p.merkle_tree_open_ns +
                                      p.ext_sumcheck_init_ns +
                                      p.ext_sumcheck_current_poly_ns +
                                      p.ext_sumcheck_receive_challenge_ns +
                                      p.ext_prover_commit_round_ns +
                                      p.ext_merkle_tree_build_ns +
                                      p.ext_merkle_tree_open_ns +
                                      p.ext_msg0_coeffs_suffix_fold_ns +
                                      p.z2k_backend_prove_ns +
                                      p.transcript_absorb_ns +
                                      p.transcript_challenge_ns;
    const double other_ms =
        (p.pcs_prove_ns > accounted_ns) ? NsToMs(p.pcs_prove_ns - accounted_ns)
                                        : 0.0;

    os << "  [profile-prover]\n";
    os << "    BaseFoldPCSProveEval total:  " << total_ms << " ms"
       << "  (calls " << p.pcs_prove_calls << ")\n";
    os << "    EncodeFoldableUnchecked:     " << encode_ms << " ms"
       << "  (calls " << p.encode_foldable_unchecked_calls << ")\n";
    os << "    SumcheckProver total:        " << sumcheck_total_ms << " ms\n";
    os << "      init:                      " << sumcheck_init_ms << " ms"
       << "  (calls " << p.sumcheck_init_calls << ")\n";
    os << "      CurrentPolynomial:         " << sumcheck_poly_ms << " ms"
       << "  (calls " << p.sumcheck_current_poly_calls << ")\n";
    os << "      ReceiveChallenge:          " << sumcheck_recv_ms << " ms"
       << "  (calls " << p.sumcheck_receive_challenge_calls << ")\n";
    os << "    ProverCommitRoundNoValidate: " << commit_ms << " ms"
       << "  (calls " << p.prover_commit_round_calls << ")\n";
    os << "    MerkleTree::Build:           " << merkle_build_ms << " ms"
       << "  (calls " << p.merkle_tree_build_calls << ")\n";
    os << "    MerkleTree::Open:            " << merkle_open_ms << " ms"
       << "  (calls " << p.merkle_tree_open_calls << ")\n";
    os << "    Transcript absorb:           " << absorb_ms << " ms"
       << "  (calls " << p.transcript_absorb_calls << ")\n";
    os << "    Transcript challenge:        " << challenge_ms << " ms"
       << "  (calls " << p.transcript_challenge_calls << ")\n";
    if (p.z2k_backend_prove_calls > 0) {
      os << "    Z2kBackendProveEval:        " << z2k_backend_prove_ms << " ms"
         << "  (calls " << p.z2k_backend_prove_calls << ")\n";
    }
    if (has_ext_prover_breakdown) {
      os << "    ExtensionSumcheck total:     " << ext_sumcheck_total_ms << " ms\n";
      os << "      init:                      " << ext_sumcheck_init_ms << " ms"
         << "  (calls " << p.ext_sumcheck_init_calls << ")\n";
      os << "      CurrentPolynomial:         " << ext_sumcheck_poly_ms << " ms"
         << "  (calls " << p.ext_sumcheck_current_poly_calls << ")\n";
      os << "      ReceiveChallenge:          " << ext_sumcheck_recv_ms << " ms"
         << "  (calls " << p.ext_sumcheck_receive_challenge_calls << ")\n";
      os << "    ExtensionCommitRound:        " << ext_commit_ms << " ms"
         << "  (calls " << p.ext_prover_commit_round_calls << ")\n";
      os << "    ExtensionMerkleTree::Build:  " << ext_merkle_build_ms << " ms"
         << "  (calls " << p.ext_merkle_tree_build_calls << ")\n";
      os << "    ExtensionMerkleTree::Open:   " << ext_merkle_open_ms << " ms"
         << "  (calls " << p.ext_merkle_tree_open_calls << ")\n";
      os << "    Msg0CoeffsAtSuffixChallenges:"
         << ((ext_msg0_suffix_fold_ms < 10.0) ? "  " : " ")
         << ext_msg0_suffix_fold_ms << " ms"
         << "  (calls " << p.ext_msg0_coeffs_suffix_fold_calls << ")\n";
    }
    os << "    Other (total - above):       " << other_ms << " ms\n";
  }

  const bool has_ring_switch_outer_prover_breakdown =
      (p.ring_switch_outer_prove_total_calls +
       p.ring_switch_outer_prove_tensor_build_calls +
       p.ring_switch_outer_prove_compute_s_calls +
       p.ring_switch_outer_prove_recover_partials_calls +
       p.ring_switch_outer_prove_batch_prep_calls +
       p.ring_switch_outer_prove_sumcheck_calls +
       p.ring_switch_outer_prove_final_check_calls) > 0;
  if (has_ring_switch_outer_prover_breakdown) {
    const double tensor_ms = NsToMs(p.ring_switch_outer_prove_tensor_build_ns);
    const double compute_s_ms = NsToMs(p.ring_switch_outer_prove_compute_s_ns);
    const double recover_ms =
        NsToMs(p.ring_switch_outer_prove_recover_partials_ns);
    const double batch_prep_ms =
        NsToMs(p.ring_switch_outer_prove_batch_prep_ns);
    const double sumcheck_ms = NsToMs(p.ring_switch_outer_prove_sumcheck_ns);
    const double final_check_ms =
        NsToMs(p.ring_switch_outer_prove_final_check_ns);
    const double profiled_ms = tensor_ms + compute_s_ms + recover_ms +
                               batch_prep_ms + sumcheck_ms + final_check_ms;
    const double total_ms =
        (p.ring_switch_outer_prove_total_calls > 0)
            ? NsToMs(p.ring_switch_outer_prove_total_ns)
            : profiled_ms;
    const double other_ms =
        (total_ms > profiled_ms) ? (total_ms - profiled_ms) : 0.0;

    os << "  [profile-ring-switch-outer-prover]\n";
    os << "    total:                        " << total_ms << " ms\n";
    os << "    BuildComponentTensor:         " << tensor_ms << " ms"
       << "  (calls " << p.ring_switch_outer_prove_tensor_build_calls
       << ")\n";
    os << "    ComputeSByU:                  " << compute_s_ms << " ms"
       << "  (calls " << p.ring_switch_outer_prove_compute_s_calls
       << ")\n";
    os << "    RecoverPartialsEq1:           " << recover_ms << " ms"
       << "  (calls " << p.ring_switch_outer_prove_recover_partials_calls
       << ")\n";
    os << "    BatchPrepPrefix:              " << batch_prep_ms << " ms"
       << "  (calls " << p.ring_switch_outer_prove_batch_prep_calls
       << ")\n";
    os << "    SumcheckSuffixRounds:         " << sumcheck_ms << " ms"
       << "  (calls " << p.ring_switch_outer_prove_sumcheck_calls << ")\n";
    os << "    FinalTstarGstarEq3:           " << final_check_ms << " ms"
       << "  (calls " << p.ring_switch_outer_prove_final_check_calls
       << ")\n";
    os << "    Other/unaccounted:            " << other_ms << " ms\n";
  }

  if (p.pcs_verify_calls > 0) {
    const double total_ms = NsToMs(p.pcs_verify_ns);
    const double query_ms = NsToMs(p.verify_query_merkle_ns);
    const double merkle_open_ms = NsToMs(p.merkle_verify_opening_ns);
    const double merkle_commit_ms = NsToMs(p.merkle_commit_oracle_ns);
    const double eval_line_ms = NsToMs(p.eval_line_at_ns);
    const double z2k_backend_verify_ms = NsToMs(p.z2k_backend_verify_ns);
    const double inv_ms = NsToMs(p.try_invert_unit_ns);
    const double is_unit_ms = NsToMs(p.is_unit_ns);
    const double inv_fallback_ms = NsToMs(p.inv_fallback_ns);
    const double ext_query_ms = NsToMs(p.ext_verify_query_merkle_ns);
    const double ext_merkle_open_ms = NsToMs(p.ext_merkle_verify_opening_ns);
    const double ext_eval_line_ms = NsToMs(p.ext_eval_line_at_ns);
    const bool has_ext_verifier_breakdown =
        (p.ext_verify_query_merkle_calls + p.ext_merkle_verify_opening_calls +
         p.ext_eval_line_at_calls) > 0;

    const std::uint64_t query_accounted_ns = p.merkle_verify_opening_ns +
                                            p.merkle_commit_oracle_ns +
                                            p.eval_line_at_ns;
    const std::uint64_t ext_query_accounted_ns =
        p.ext_merkle_verify_opening_ns + p.ext_eval_line_at_ns;
    const std::uint64_t total_query_ns =
        p.verify_query_merkle_ns + p.ext_verify_query_merkle_ns;
    const std::uint64_t total_query_accounted_ns =
        query_accounted_ns + ext_query_accounted_ns;
    const double query_other_ms =
        (total_query_ns > total_query_accounted_ns)
            ? NsToMs(total_query_ns - total_query_accounted_ns)
            : 0.0;
    const double outside_query_ms =
        (p.pcs_verify_ns > total_query_ns + p.z2k_backend_verify_ns)
            ? NsToMs(p.pcs_verify_ns - total_query_ns -
                     p.z2k_backend_verify_ns)
            : 0.0;

    os << "  [profile-verifier]\n";
    os << "    BaseFoldPCSVerifyEval total: " << total_ms << " ms"
       << "  (calls " << p.pcs_verify_calls << ")\n";
    os << "    VerifyQueryMultiproofs:       " << query_ms << " ms"
       << "  (calls " << p.verify_query_merkle_calls << ")\n";
    os << "    MerkleVerifyMultiproof:      " << merkle_open_ms << " ms"
       << "  (calls " << p.merkle_verify_opening_calls << ")\n";
    os << "    MerkleCommitOracle:          " << merkle_commit_ms << " ms"
       << "  (calls " << p.merkle_commit_oracle_calls << ")\n";
    os << "    EvalLineAt:                  " << eval_line_ms << " ms"
       << "  (calls " << p.eval_line_at_calls << ")\n";
    if (p.z2k_backend_verify_calls > 0) {
      os << "    Z2kBackendVerifyEval:       " << z2k_backend_verify_ms
         << " ms"
         << "  (calls " << p.z2k_backend_verify_calls << ")\n";
    }
    os << "    TryInvertUnit (subset):      " << inv_ms << " ms"
       << "  (calls " << p.try_invert_unit_calls << ", cache hits "
       << p.try_invert_unit_cache_hits << ")\n";
    os << "      IsUnit (subset):           " << is_unit_ms << " ms"
       << "  (calls " << p.is_unit_calls << ")\n";
    os << "      Inv fallback (subset):     " << inv_fallback_ms << " ms"
       << "  (calls " << p.inv_fallback_calls << ")\n";
    if (has_ext_verifier_breakdown) {
      os << "    ExtensionVerifyQuery:        " << ext_query_ms << " ms"
         << "  (calls " << p.ext_verify_query_merkle_calls << ")\n";
      os << "    ExtensionMerkleMultiproof:   " << ext_merkle_open_ms << " ms"
         << "  (calls " << p.ext_merkle_verify_opening_calls << ")\n";
      os << "    EvalLineAtExtension:         " << ext_eval_line_ms << " ms"
         << "  (calls " << p.ext_eval_line_at_calls << ")\n";
    }
    os << "    Inside queries other:        " << query_other_ms << " ms\n";
    os << "    Outside queries:             " << outside_query_ms << " ms\n";
  }

  const bool has_ring_switch_outer_verifier_breakdown =
      (p.ring_switch_outer_verify_total_calls +
       p.ring_switch_outer_verify_tensor_build_calls +
       p.ring_switch_outer_verify_recover_partials_calls +
       p.ring_switch_outer_verify_prefix_replay_calls +
       p.ring_switch_outer_verify_sumcheck_replay_calls +
       p.ring_switch_outer_verify_final_check_calls) > 0;
  if (has_ring_switch_outer_verifier_breakdown) {
    const double tensor_ms = NsToMs(p.ring_switch_outer_verify_tensor_build_ns);
    const double recover_ms =
        NsToMs(p.ring_switch_outer_verify_recover_partials_ns);
    const double prefix_ms =
        NsToMs(p.ring_switch_outer_verify_prefix_replay_ns);
    const double sumcheck_ms =
        NsToMs(p.ring_switch_outer_verify_sumcheck_replay_ns);
    const double final_check_ms =
        NsToMs(p.ring_switch_outer_verify_final_check_ns);
    const double profiled_ms =
        tensor_ms + recover_ms + prefix_ms + sumcheck_ms + final_check_ms;
    const double total_ms =
        (p.ring_switch_outer_verify_total_calls > 0)
            ? NsToMs(p.ring_switch_outer_verify_total_ns)
            : profiled_ms;
    const double other_ms =
        (total_ms > profiled_ms) ? (total_ms - profiled_ms) : 0.0;

    os << "  [profile-ring-switch-outer-verifier]\n";
    os << "    total:                        " << total_ms << " ms\n";
    os << "    BuildComponentTensor:         " << tensor_ms << " ms"
       << "  (calls " << p.ring_switch_outer_verify_tensor_build_calls
       << ")\n";
    os << "    RecoverPartialsEq1:           " << recover_ms << " ms"
       << "  (calls " << p.ring_switch_outer_verify_recover_partials_calls
       << ")\n";
    os << "    ReplayPrefixBatching:         " << prefix_ms << " ms"
       << "  (calls " << p.ring_switch_outer_verify_prefix_replay_calls
       << ")\n";
    os << "    ReplaySumcheckChain:          " << sumcheck_ms << " ms"
       << "  (calls " << p.ring_switch_outer_verify_sumcheck_replay_calls
       << ")\n";
    os << "    FinalGstarEq3:                " << final_check_ms << " ms"
       << "  (calls " << p.ring_switch_outer_verify_final_check_calls
       << ")\n";
    os << "    Other/unaccounted:            " << other_ms << " ms\n";
  }

  const bool has_frobenius_outer_prover_breakdown =
      (p.frobenius_outer_prove_total_calls +
       p.frobenius_outer_prove_orbit_build_calls +
       p.frobenius_outer_prove_compute_s_calls +
       p.frobenius_outer_prove_recover_partials_calls +
       p.frobenius_outer_prove_batch_prep_calls +
       p.frobenius_outer_prove_sumcheck_calls +
       p.frobenius_outer_prove_final_check_calls) > 0;
  if (has_frobenius_outer_prover_breakdown) {
    const double orbit_ms = NsToMs(p.frobenius_outer_prove_orbit_build_ns);
    const double compute_s_ms = NsToMs(p.frobenius_outer_prove_compute_s_ns);
    const double recover_ms =
        NsToMs(p.frobenius_outer_prove_recover_partials_ns);
    const double batch_prep_ms =
        NsToMs(p.frobenius_outer_prove_batch_prep_ns);
    const double sumcheck_ms = NsToMs(p.frobenius_outer_prove_sumcheck_ns);
    const double final_check_ms =
        NsToMs(p.frobenius_outer_prove_final_check_ns);
    const double profiled_ms = orbit_ms + compute_s_ms + recover_ms +
                               batch_prep_ms + sumcheck_ms + final_check_ms;
    const double total_ms =
        (p.frobenius_outer_prove_total_calls > 0)
            ? NsToMs(p.frobenius_outer_prove_total_ns)
            : profiled_ms;
    const double other_ms =
        (total_ms > profiled_ms) ? (total_ms - profiled_ms) : 0.0;

    os << "  [profile-frobenius-outer-prover]\n";
    os << "    total:                        " << total_ms << " ms\n";
    os << "    BuildSuffixOrbitCache:        " << orbit_ms << " ms"
       << "  (calls " << p.frobenius_outer_prove_orbit_build_calls << ")\n";
    os << "    ComputeSByI:                  " << compute_s_ms << " ms"
       << "  (calls " << p.frobenius_outer_prove_compute_s_calls << ")\n";
    os << "    RecoverPartialsEq1:           " << recover_ms << " ms"
       << "  (calls " << p.frobenius_outer_prove_recover_partials_calls
       << ")\n";
    os << "    BatchPrepPrefix:              " << batch_prep_ms << " ms"
       << "  (calls " << p.frobenius_outer_prove_batch_prep_calls << ")\n";
    os << "    SumcheckSuffixRounds:         " << sumcheck_ms << " ms"
       << "  (calls " << p.frobenius_outer_prove_sumcheck_calls << ")\n";
    os << "    FinalTstarGstarEq3:           " << final_check_ms << " ms"
       << "  (calls " << p.frobenius_outer_prove_final_check_calls << ")\n";
    os << "    Other/unaccounted:            " << other_ms << " ms\n";
  }

  const bool has_frobenius_outer_verifier_breakdown =
      (p.frobenius_outer_verify_total_calls +
       p.frobenius_outer_verify_orbit_build_calls +
       p.frobenius_outer_verify_recover_partials_calls +
       p.frobenius_outer_verify_prefix_replay_calls +
       p.frobenius_outer_verify_sumcheck_replay_calls +
       p.frobenius_outer_verify_final_check_calls) > 0;
  if (has_frobenius_outer_verifier_breakdown) {
    const double orbit_ms = NsToMs(p.frobenius_outer_verify_orbit_build_ns);
    const double recover_ms =
        NsToMs(p.frobenius_outer_verify_recover_partials_ns);
    const double prefix_ms =
        NsToMs(p.frobenius_outer_verify_prefix_replay_ns);
    const double sumcheck_ms =
        NsToMs(p.frobenius_outer_verify_sumcheck_replay_ns);
    const double final_check_ms =
        NsToMs(p.frobenius_outer_verify_final_check_ns);
    const double profiled_ms =
        orbit_ms + recover_ms + prefix_ms + sumcheck_ms + final_check_ms;
    const double total_ms =
        (p.frobenius_outer_verify_total_calls > 0)
            ? NsToMs(p.frobenius_outer_verify_total_ns)
            : profiled_ms;
    const double other_ms =
        (total_ms > profiled_ms) ? (total_ms - profiled_ms) : 0.0;

    os << "  [profile-frobenius-outer-verifier]\n";
    os << "    total:                        " << total_ms << " ms\n";
    os << "    BuildSuffixOrbitCache:        " << orbit_ms << " ms"
       << "  (calls " << p.frobenius_outer_verify_orbit_build_calls << ")\n";
    os << "    RecoverPartialsEq1:           " << recover_ms << " ms"
       << "  (calls " << p.frobenius_outer_verify_recover_partials_calls
       << ")\n";
    os << "    ReplayPrefixBatching:         " << prefix_ms << " ms"
       << "  (calls " << p.frobenius_outer_verify_prefix_replay_calls
       << ")\n";
    os << "    ReplaySumcheckChain:          " << sumcheck_ms << " ms"
       << "  (calls " << p.frobenius_outer_verify_sumcheck_replay_calls
       << ")\n";
    os << "    FinalGstarEq3:                " << final_check_ms << " ms"
       << "  (calls " << p.frobenius_outer_verify_final_check_calls
       << ")\n";
    os << "    Other/unaccounted:            " << other_ms << " ms\n";
  }
}

}  // namespace basefold

#endif  // BASEFOLD_PROFILE_HPP_
