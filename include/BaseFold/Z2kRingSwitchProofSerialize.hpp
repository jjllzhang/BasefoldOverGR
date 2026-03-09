#ifndef BASEFOLD_Z2KRINGSWITCHPROOFSERIALIZE_HPP_
#define BASEFOLD_Z2KRINGSWITCHPROOFSERIALIZE_HPP_

#include <cstdint>

#include "BaseFold/Z2kRingSwitchPCS.hpp"

namespace basefold {

struct RingSwitchProofEncodingOptions {
  bool include_version_byte = true;
  Z2kPCSBackendProofEncodingOptions backend_proof_options;
};

Bytes SerializeRingSwitchPCSOuterProofFixedBytes(
    const RingSwitchPCSParams &params, const RingSwitchPCSEvalProof &proof,
    const RingSwitchProofEncodingOptions &options = {});

Bytes SerializeRingSwitchPCSEvalProofFixedBytes(
    const RingSwitchPCSParams &params, const RingSwitchPCSEvalProof &proof,
    const RingSwitchProofEncodingOptions &options = {});

std::uint64_t RingSwitchPCSOuterProofSizeBytes(
    const RingSwitchPCSParams &params, const RingSwitchPCSEvalProof &proof,
    const RingSwitchProofEncodingOptions &options = {});

double RingSwitchPCSOuterProofSizeKB(
    const RingSwitchPCSParams &params, const RingSwitchPCSEvalProof &proof,
    const RingSwitchProofEncodingOptions &options = {});

std::uint64_t RingSwitchPCSEvalProofSizeBytes(
    const RingSwitchPCSParams &params, const RingSwitchPCSEvalProof &proof,
    const RingSwitchProofEncodingOptions &options = {});

double RingSwitchPCSEvalProofSizeKB(
    const RingSwitchPCSParams &params, const RingSwitchPCSEvalProof &proof,
    const RingSwitchProofEncodingOptions &options = {});

}  // namespace basefold

#endif  // BASEFOLD_Z2KRINGSWITCHPROOFSERIALIZE_HPP_
