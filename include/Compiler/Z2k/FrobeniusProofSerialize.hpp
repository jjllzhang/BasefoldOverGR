#ifndef BASEFOLD_Z2K_FROBENIUSPROOFSERIALIZE_HPP_
#define BASEFOLD_Z2K_FROBENIUSPROOFSERIALIZE_HPP_

#include <cstdint>

#include "Compiler/Z2k/FrobeniusPCS.hpp"

namespace basefold {

struct FrobeniusProofEncodingOptions {
  bool include_version_byte = true;
  Z2kPCSBackendProofEncodingOptions backend_proof_options;
};

Bytes SerializeFrobeniusPCSOuterProofFixedBytes(
    const FrobeniusPCSParams &params, const FrobeniusPCSOuterEvalProof &proof,
    const FrobeniusProofEncodingOptions &options = {});

Bytes SerializeFrobeniusPCSOuterProofFixedBytes(
    const FrobeniusPCSParams &params, const FrobeniusPCSEvalProof &proof,
    const FrobeniusProofEncodingOptions &options = {});

Bytes SerializeFrobeniusPCSEvalProofFixedBytes(
    const FrobeniusPCSParams &params, const FrobeniusPCSEvalProof &proof,
    const FrobeniusProofEncodingOptions &options = {});

std::uint64_t FrobeniusPCSOuterProofSizeBytes(
    const FrobeniusPCSParams &params, const FrobeniusPCSOuterEvalProof &proof,
    const FrobeniusProofEncodingOptions &options = {});

std::uint64_t FrobeniusPCSOuterProofSizeBytes(
    const FrobeniusPCSParams &params, const FrobeniusPCSEvalProof &proof,
    const FrobeniusProofEncodingOptions &options = {});

double FrobeniusPCSOuterProofSizeKB(
    const FrobeniusPCSParams &params, const FrobeniusPCSOuterEvalProof &proof,
    const FrobeniusProofEncodingOptions &options = {});

double FrobeniusPCSOuterProofSizeKB(
    const FrobeniusPCSParams &params, const FrobeniusPCSEvalProof &proof,
    const FrobeniusProofEncodingOptions &options = {});

std::uint64_t FrobeniusPCSEvalProofSizeBytes(
    const FrobeniusPCSParams &params, const FrobeniusPCSEvalProof &proof,
    const FrobeniusProofEncodingOptions &options = {});

double FrobeniusPCSEvalProofSizeKB(
    const FrobeniusPCSParams &params, const FrobeniusPCSEvalProof &proof,
    const FrobeniusProofEncodingOptions &options = {});

}  // namespace basefold

#endif  // BASEFOLD_Z2K_FROBENIUSPROOFSERIALIZE_HPP_
