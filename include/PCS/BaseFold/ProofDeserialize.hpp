#ifndef BASEFOLD_PROOFDESERIALIZE_HPP_
#define BASEFOLD_PROOFDESERIALIZE_HPP_

#include "PCS/BaseFold/ProofSerialize.hpp"

namespace basefold {

BaseFoldPCSEvalProof DeserializeBaseFoldPCSEvalProofFixedBytes(
    const Bytes &bytes, const FixedProofEncodingOptions &options = {});

BaseFoldPCSEvalProof DeserializeBaseFoldPCSEvalProofFixedBytes(
    const Byte *data, std::size_t size,
    const FixedProofEncodingOptions &options = {});

}  // namespace basefold

#endif  // BASEFOLD_PROOFDESERIALIZE_HPP_
