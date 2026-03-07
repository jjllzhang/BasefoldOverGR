#include "BaseFold/ProofSize.hpp"

#include "BaseFold/ProofSerialize.hpp"

namespace basefold {

std::uint64_t BaseFoldPCSEvalProofSizeBytes(
    const BaseFoldPCSEvalProof &proof,
    const BaseFoldProofSizeOptions &options) {
  FixedProofEncodingOptions encoding_options;
  encoding_options.include_version_byte = options.include_version_byte;
  encoding_options.challenge_ext_degree = options.challenge_ext_degree;
  return CountSerializedBaseFoldPCSEvalProofFixedBytes(proof, encoding_options);
}

double BaseFoldPCSEvalProofSizeKB(
    const BaseFoldPCSEvalProof &proof,
    const BaseFoldProofSizeOptions &options) {
  const std::uint64_t bytes = BaseFoldPCSEvalProofSizeBytes(proof, options);
  return static_cast<double>(bytes) / 1024.0;
}

}  // namespace basefold
