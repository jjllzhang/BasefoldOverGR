#ifndef BASEFOLD_PCS_COMMON_NTLPARALLEL_HPP_
#define BASEFOLD_PCS_COMMON_NTLPARALLEL_HPP_

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

#if defined(BASEFOLD_USE_OPENMP)
#include <omp.h>
#endif

namespace basefold {
namespace pcs_common_internal {

struct NtlThreadContextSnapshot {
  NTL::ZZ base_modulus;
  NTL::ZZ_pX extension_modulus;
};

inline NtlThreadContextSnapshot CaptureNtlThreadContextSnapshot() {
  return {NTL::ZZ_p::modulus(), NTL::ZZ_pE::modulus().val()};
}

inline void InitNtlThreadContext(const NtlThreadContextSnapshot &ctx) {
  NTL::ZZ_p::init(ctx.base_modulus);
  NTL::ZZ_pE::init(ctx.extension_modulus);
}

inline int ChooseElementParallelThreads(long work_items,
                                        long parallel_threshold) {
#if defined(BASEFOLD_USE_OPENMP)
  if (parallel_threshold <= 0) {
    parallel_threshold = 1;
  }
  if (work_items >= parallel_threshold) {
    const int max_threads = omp_get_max_threads();
    int threads_to_use = static_cast<int>(work_items / parallel_threshold);
    if (threads_to_use > max_threads) {
      threads_to_use = max_threads;
    }
    if (threads_to_use >= 2) {
      return threads_to_use;
    }
  }
#else
  (void)work_items;
  (void)parallel_threshold;
#endif
  return 1;
}

template <typename Fn>
void ForEachIndexMaybeParallel(long begin, long end, long parallel_threshold,
                               const Fn &fn) {
  if (end <= begin) {
    return;
  }

#if defined(BASEFOLD_USE_OPENMP)
  const long work_items = end - begin;
  const int threads_to_use =
      ChooseElementParallelThreads(work_items, parallel_threshold);
  if (threads_to_use >= 2) {
    const NtlThreadContextSnapshot ntl_ctx = CaptureNtlThreadContextSnapshot();
#pragma omp parallel num_threads(threads_to_use)
    {
      InitNtlThreadContext(ntl_ctx);
#pragma omp for schedule(static)
      for (long i = begin; i < end; ++i) {
        fn(i);
      }
    }
    return;
  }
#endif

  for (long i = begin; i < end; ++i) {
    fn(i);
  }
}

}  // namespace pcs_common_internal
}  // namespace basefold

#endif  // BASEFOLD_PCS_COMMON_NTLPARALLEL_HPP_
