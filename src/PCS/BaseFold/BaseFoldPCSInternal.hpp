#ifndef BASEFOLD_SRC_PCS_BASEFOLD_BASEFOLDPCSINTERNAL_HPP_
#define BASEFOLD_SRC_PCS_BASEFOLD_BASEFOLDPCSINTERNAL_HPP_

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

#include <vector>

#include "PCS/BaseFold/BaseFoldPCS.hpp"
#include "PCS/Common/NtlParallel.hpp"
#include "PCS/Common/Profile.hpp"
#include "PCS/Common/Transcript.hpp"

namespace basefold {
namespace basefold_pcs_internal {

long ParsePositiveEnvLong(const char *name, long fallback);
int ParsePositiveEnvInt(const char *name, int fallback);

VerifierQueryParallelConfig ParseVerifierQueryParallelConfigFromEnv();
VerifierQueryParallelConfig &MutableVerifierQueryParallelConfig();
VerifierQueryParallelConfig LoadVerifierQueryParallelConfig();
ProverCommitParallelConfig ParseProverCommitParallelConfigFromEnv();
ProverCommitParallelConfig &MutableProverCommitParallelConfig();
ProverCommitParallelConfig LoadProverCommitParallelConfig();
long LoadEffectiveExtElementsPerThread();

using pcs_common_internal::CaptureNtlThreadContextSnapshot;
using pcs_common_internal::ChooseElementParallelThreads;
using pcs_common_internal::ForEachIndexMaybeParallel;
using pcs_common_internal::InitNtlThreadContext;
using pcs_common_internal::NtlThreadContextSnapshot;

int ChooseQueryVerifyThreads(long num_queries);

template <typename Fn>
bool VerifyQueriesMaybeParallel(long num_queries, Profile *prof,
                                const Fn &verify_one_query) {
  if (num_queries <= 0) {
    return true;
  }

#if defined(BASEFOLD_USE_OPENMP)
  // Keep profile accounting precise: per-thread timer accumulation in parallel
  // would over-count wall-clock time in the breakdown.
  if (prof == nullptr) {
    const int threads_to_use = ChooseQueryVerifyThreads(num_queries);
    if (threads_to_use >= 2) {
      std::vector<unsigned char> query_ok(static_cast<std::size_t>(num_queries),
                                          static_cast<unsigned char>(1));

      const NTL::ZZ base_modulus = NTL::ZZ_p::modulus();
      const NTL::ZZ_pX extension_modulus = NTL::ZZ_pE::modulus().val();

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
        if (query_ok[static_cast<std::size_t>(q)] == 0) {
          return false;
        }
      }
      return true;
    }
  }
#endif

  for (long q = 0; q < num_queries; ++q) {
    if (!verify_one_query(q)) {
      return false;
    }
  }
  return true;
}

void ValidateParamsOrThrow(const FoldableCodeParams &params);

long Pow2Checked(long e);
bool IsPowerOfTwoLong(long n);
long Log2ExactPowerOfTwoLong(long n);
long CodewordLengthAtLevelNoValidate(const FoldableCodeParams &params,
                                     long level);
IOPPQueryPlan MakeQueryPlanNoValidate(long initial_mu,
                                      const FoldableCodeParams &params);

void ValidateCommittedTopOracleArtifactsOrThrow(
    const BaseFoldPCSCommitArtifacts &commit_artifacts, long expected_length,
    const char *func_name);

std::vector<std::vector<long>> CollectBaseQueryIndicesByTree(
    const std::vector<IOPPQueryPlan> &query_plans,
    const FoldableCodeParams &params);

bool HasMerkleMultiproofPayload(const MerkleMultiproof &proof);

HashTranscript MakeBaseFoldTranscript();
void AbsorbPublicInput(HashTranscript &transcript,
                       const MerkleRoot &commitment,
                       const std::vector<FieldElement> &z,
                       const FieldElement &y);

bool TryInvertBaseUnit(FieldElement &inv_out, const FieldElement &a);
bool BatchInvertBaseUnits(std::vector<FieldElement> &inverses,
                          const std::vector<FieldElement> &values);
FieldElement EvalLineAtWithInvDenom(const FieldElement &x,
                                    const FieldElement &x1,
                                    const FieldElement &y1,
                                    const FieldElement &y2,
                                    const FieldElement &inv_denom);

void ProverCommitRoundNoValidate(Oracle &pi_i, const Oracle &pi_ip1,
                                 const FieldElement &alpha_i, long level_i,
                                 const FoldableCodeParams &params);

}  // namespace basefold_pcs_internal
}  // namespace basefold

#endif  // BASEFOLD_SRC_PCS_BASEFOLD_BASEFOLDPCSINTERNAL_HPP_
