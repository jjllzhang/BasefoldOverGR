#ifndef BASEFOLD_PROOFSIZE_HPP_
#define BASEFOLD_PROOFSIZE_HPP_

#include <cstdint>

#include "PCS/BaseFold/BaseFoldPCS.hpp"

namespace basefold {

struct BaseFoldProofSizeOptions {
  bool include_version_byte = true;
  long challenge_ext_degree = 0;
};

std::uint64_t BaseFoldPCSEvalProofSizeBytes(
    const BaseFoldPCSEvalProof &proof,
    const BaseFoldProofSizeOptions &options = {});

double BaseFoldPCSEvalProofSizeKB(
    const BaseFoldPCSEvalProof &proof,
    const BaseFoldProofSizeOptions &options = {});

}  // namespace basefold

#endif  // BASEFOLD_PROOFSIZE_HPP_
